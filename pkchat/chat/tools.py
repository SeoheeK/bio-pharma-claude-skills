"""The PK tool surface exposed to the chatbot.

Each tool has an Anthropic tool-use compatible JSON schema and a Python handler
that maps validated input to a human-readable result string. The same registry
drives both the Claude-backed agent (as function-calling tools) and the local
rule-based router, so behaviour is identical regardless of backend.
"""
from __future__ import annotations

from typing import Callable

from ..models import PKParameters
from ..regimen import DosingRegimen
from ..simulate import simulate, simulate_population
from ..nca import nca
from ..estimation import fit_individual
from ..ml.allometry import allometry_ensemble, SpeciesPK
from ..casedb import recommend_cases


def _params_from(args: dict) -> PKParameters:
    return PKParameters(
        CL=float(args["CL"]),
        V1=float(args["V1"]),
        Q=_opt(args, "Q"),
        V2=_opt(args, "V2"),
        ka=_opt(args, "ka"),
        F=float(args.get("F", 1.0)),
        Vmax=_opt(args, "Vmax"),
        Km=_opt(args, "Km"),
    )


def _opt(args, key):
    v = args.get(key)
    return float(v) if v is not None else None


def _regimen_from(args: dict) -> DosingRegimen:
    dose = float(args["dose"])
    route = args.get("route", "iv_bolus")
    tinf = float(args.get("tinf", 0.0))
    n = int(args.get("n_doses", 1))
    if n > 1:
        interval = float(args["interval"])
        return DosingRegimen.multiple(dose, interval, n, route=route, tinf=tinf)
    return DosingRegimen.single(dose, route=route, tinf=tinf)


# --- Handlers ---------------------------------------------------------------

def _h_simulate(args: dict) -> str:
    p = _params_from(args)
    reg = _regimen_from(args)
    t_end = _opt(args, "t_end")
    res = simulate(p, reg, t_end=t_end, dt=0.02)
    # Report a compact profile plus key metrics.
    d = p.derived()
    lines = ["Simulation result:"]
    lines.append(f"  Cmax = {res.cmax():.4g} mg/L at t = {res.tmax():.3g} h")
    thalf = d.get("t_half_beta", d.get("t_half_1cmt"))
    if thalf is not None:
        lines.append(f"  terminal half-life ~ {thalf:.3g} h")
    # AUC over the simulated window via NCA.
    try:
        n = nca(res.times, res.concentrations,
                dose=reg.events[0].amount, route=reg.events[0].route)
        lines.append(f"  AUC(0-inf) ~ {n.auc_inf:.4g} mg*h/L")
    except Exception:
        pass
    sample_t = _sample_times(res.times, 8)
    lines.append("  profile (t[h]: C[mg/L]): " +
                 ", ".join(f"{t:.2g}:{_interp(res, t):.3g}" for t in sample_t))
    return "\n".join(lines)


def _h_nca(args: dict) -> str:
    times = [float(x) for x in args["times"]]
    concs = [float(x) for x in args["concentrations"]]
    res = nca(times, concs, dose=_opt(args, "dose"),
              route=args.get("route", "iv_bolus"))
    lines = ["Non-compartmental analysis:"]
    lines.append(f"  Cmax = {res.cmax:.4g}, Tmax = {res.tmax:.3g} h")
    lines.append(f"  AUC(0-last) = {res.auc_last:.4g}, AUC(0-inf) = {res.auc_inf:.4g} mg*h/L")
    lines.append(f"  terminal t1/2 = {res.t_half:.3g} h (lambda_z from {res.n_lambda_points} pts, "
                 f"{res.auc_extrap_pct:.1f}% extrapolated)")
    if res.clearance is not None:
        lines.append(f"  CL = {res.clearance:.4g} L/h" +
                     (f", Vz = {res.vz:.4g} L" if res.vz else ""))
    return "\n".join(lines)


def _h_fit(args: dict) -> str:
    times = [float(x) for x in args["times"]]
    concs = [float(x) for x in args["concentrations"]]
    reg = _regimen_from(args)
    estimate = tuple(args.get("estimate", ["CL", "V1"]))
    init = None
    if args.get("ka") is not None:
        init = PKParameters(CL=1.0, V1=10.0, ka=float(args["ka"]))
    fit = fit_individual(times, concs, reg, estimate=estimate, init=init,
                         weighting=args.get("weighting", "proportional"))
    return "Model fit:\n" + fit.summary()


def _h_allometry(args: dict) -> str:
    species = [SpeciesPK(s["species"], float(s["body_weight_kg"]),
                         float(s["clearance"]), _opt(s, "volume"))
               for s in args["species"]]
    res = allometry_ensemble(species, human_bw=float(args.get("human_bw", 70.0)))
    return "Allometric scaling:\n" + res.summary()


def _h_recommend(args: dict) -> str:
    rec = recommend_cases(
        drug_class=args.get("drug_class"),
        modality=args.get("modality"),
        route=args.get("route"),
    )
    return rec.summary()


def _h_population(args: dict) -> str:
    p = _params_from(args)
    reg = _regimen_from(args)
    iiv = args.get("iiv_cv", {"CL": 0.3, "V1": 0.2})
    iiv = {k: float(v) for k, v in iiv.items()}
    res = simulate_population(
        p, reg, n_subjects=int(args.get("n_subjects", 200)),
        iiv_cv=iiv, residual_cv=float(args.get("residual_cv", 0.0)),
        t_end=_opt(args, "t_end"), dt=0.05,
    )
    med = res.median()
    p5 = res.percentile(5)
    p95 = res.percentile(95)
    # Peak of the median profile.
    imax = max(range(len(med)), key=lambda i: med[i]) if med else 0
    lines = [f"Population simulation ({len(res.profiles)} virtual subjects, "
             f"IIV {iiv}):"]
    lines.append(f"  median Cmax = {med[imax]:.4g} mg/L at t = {res.times[imax]:.3g} h")
    lines.append(f"  90% prediction interval at peak: "
                 f"[{p5[imax]:.4g}, {p95[imax]:.4g}] mg/L")
    return "\n".join(lines)


# --- Registry ---------------------------------------------------------------

class Tool:
    def __init__(self, name, description, schema, handler):
        self.name = name
        self.description = description
        self.input_schema = schema
        self.handler = handler


TOOLS: dict[str, Tool] = {}


def _register(name, description, schema, handler):
    TOOLS[name] = Tool(name, description, schema, handler)


_num = {"type": "number"}
_num_arr = {"type": "array", "items": {"type": "number"}}

_register(
    "simulate_pk",
    "Simulate a concentration-time profile for a compartmental PK model given "
    "structural parameters and a dosing regimen. Returns Cmax, Tmax, half-life, "
    "AUC and a sampled profile.",
    {
        "type": "object",
        "properties": {
            "CL": {**_num, "description": "clearance [L/h]"},
            "V1": {**_num, "description": "central volume [L]"},
            "Q": {**_num, "description": "inter-compartmental clearance [L/h] (2-cmt)"},
            "V2": {**_num, "description": "peripheral volume [L] (2-cmt)"},
            "ka": {**_num, "description": "absorption rate [1/h] (oral)"},
            "F": {**_num, "description": "bioavailability 0-1 (oral)"},
            "Vmax": {**_num, "description": "max elim rate [mg/h] (Michaelis-Menten)"},
            "Km": {**_num, "description": "Km [mg/L] (Michaelis-Menten)"},
            "dose": {**_num, "description": "dose amount [mg]"},
            "route": {"type": "string", "enum": ["iv_bolus", "infusion", "oral"]},
            "tinf": {**_num, "description": "infusion duration [h]"},
            "interval": {**_num, "description": "dosing interval [h] for multiple dosing"},
            "n_doses": {"type": "integer", "description": "number of doses"},
            "t_end": {**_num, "description": "simulation end time [h]"},
        },
        "required": ["CL", "V1", "dose"],
    },
    _h_simulate,
)

_register(
    "run_nca",
    "Run non-compartmental analysis on observed concentration-time data. "
    "Returns Cmax, Tmax, AUC (linear-up/log-down), terminal half-life, and "
    "CL/Vz for IV data.",
    {
        "type": "object",
        "properties": {
            "times": {**_num_arr, "description": "sample times [h]"},
            "concentrations": {**_num_arr, "description": "concentrations [mg/L]"},
            "dose": {**_num, "description": "dose [mg] (enables CL/Vz)"},
            "route": {"type": "string", "enum": ["iv_bolus", "infusion", "oral"]},
        },
        "required": ["times", "concentrations"],
    },
    _h_nca,
)

_register(
    "fit_pk_model",
    "Estimate structural PK parameters (CL, V1, and optionally Q, V2, ka) by "
    "fitting a compartmental model to observed concentration data. Returns "
    "point estimates with relative standard errors.",
    {
        "type": "object",
        "properties": {
            "times": {**_num_arr, "description": "sample times [h]"},
            "concentrations": {**_num_arr, "description": "concentrations [mg/L]"},
            "dose": {**_num, "description": "dose [mg]"},
            "route": {"type": "string", "enum": ["iv_bolus", "infusion", "oral"]},
            "ka": {**_num, "description": "fixed absorption rate for oral models"},
            "estimate": {"type": "array", "items": {"type": "string"},
                         "description": "parameters to estimate, e.g. ['CL','V1']"},
            "weighting": {"type": "string",
                          "enum": ["proportional", "uniform", "poisson"]},
        },
        "required": ["times", "concentrations", "dose"],
    },
    _h_fit,
)

_register(
    "allometric_scaling",
    "Predict human clearance (and volume) from animal PK data using an ensemble "
    "of allometric methods (simple allometry, fixed exponent, rule of exponents, "
    "MLP/brain-weight corrections, and a bagged regressor).",
    {
        "type": "object",
        "properties": {
            "species": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "species": {"type": "string"},
                        "body_weight_kg": _num,
                        "clearance": _num,
                        "volume": _num,
                    },
                    "required": ["species", "body_weight_kg", "clearance"],
                },
            },
            "human_bw": {**_num, "description": "human body weight [kg], default 70"},
        },
        "required": ["species"],
    },
    _h_allometry,
)

_register(
    "recommend_pk_case",
    "Find published PopPK case studies similar to a new drug (by class, "
    "modality, route) and propose initial parameter estimates and a model "
    "structure.",
    {
        "type": "object",
        "properties": {
            "drug_class": {"type": "string",
                           "description": "e.g. immune_checkpoint, tki, calcineurin_inhibitor"},
            "modality": {"type": "string", "enum": ["mAb", "small_molecule"]},
            "route": {"type": "string", "enum": ["iv", "oral", "sc"]},
        },
        "required": [],
    },
    _h_recommend,
)

_register(
    "simulate_population",
    "Monte-Carlo population simulation with between-subject variability (IIV) "
    "on parameters and optional residual error. Returns median and 90% "
    "prediction interval at the concentration peak.",
    {
        "type": "object",
        "properties": {
            "CL": _num, "V1": _num, "Q": _num, "V2": _num, "ka": _num, "F": _num,
            "dose": _num,
            "route": {"type": "string", "enum": ["iv_bolus", "infusion", "oral"]},
            "tinf": _num, "interval": _num, "n_doses": {"type": "integer"},
            "n_subjects": {"type": "integer"},
            "iiv_cv": {"type": "object", "description": "param->CV, e.g. {'CL':0.3}"},
            "residual_cv": _num,
            "t_end": _num,
        },
        "required": ["CL", "V1", "dose"],
    },
    _h_population,
)


def tool_schemas() -> list[dict]:
    """Anthropic tool-use tool definitions for all registered tools."""
    return [
        {"name": t.name, "description": t.description, "input_schema": t.input_schema}
        for t in TOOLS.values()
    ]


def dispatch(name: str, args: dict) -> str:
    """Execute a tool by name with a validated argument dict."""
    if name not in TOOLS:
        return f"Unknown tool: {name}"
    try:
        return TOOLS[name].handler(args)
    except Exception as e:  # surface tool errors as readable text
        return f"Tool '{name}' failed: {e}"


# --- helpers ----------------------------------------------------------------

def _sample_times(times, n):
    if len(times) <= n:
        return times
    step = (len(times) - 1) / (n - 1)
    return [times[round(i * step)] for i in range(n)]


def _interp(res, t):
    ts, cs = res.times, res.concentrations
    if t <= ts[0]:
        return cs[0]
    if t >= ts[-1]:
        return cs[-1]
    for i in range(1, len(ts)):
        if ts[i] >= t:
            frac = (t - ts[i - 1]) / (ts[i] - ts[i - 1])
            return cs[i - 1] * (1 - frac) + cs[i] * frac
    return cs[-1]
