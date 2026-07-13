"""Structural pharmacokinetic (PK) models.

Compartmental models expressed both as ODE systems (for arbitrary dosing
regimens, integrated numerically) and as closed-form analytical solutions
(for single-dose IV/PO cases, used for speed and to validate the integrator).

Pure-Python: depends only on the standard library so it runs anywhere. Units
are kept consistent by the caller; the conventional set is CL [L/h], V [L],
dose [mg], time [h] giving concentration [mg/L].
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

Route = str  # "iv_bolus" | "infusion" | "oral"


@dataclass
class PKParameters:
    """Structural PK parameters. Unused fields stay at their defaults.

    CL, V1 are always required. Q/V2 add a peripheral compartment. ka/F add
    first-order oral absorption. Vmax/Km switch elimination to Michaelis-Menten.
    """

    CL: float = 1.0          # clearance [L/h]
    V1: float = 10.0         # central volume [L]
    Q: float | None = None   # inter-compartmental clearance [L/h]
    V2: float | None = None  # peripheral volume [L]
    ka: float | None = None  # absorption rate constant [1/h]
    F: float = 1.0           # oral bioavailability (0-1)
    Vmax: float | None = None  # max elimination rate [mg/h] (Michaelis-Menten)
    Km: float | None = None    # concentration at half Vmax [mg/L]

    @property
    def two_compartment(self) -> bool:
        return self.Q is not None and self.V2 is not None

    @property
    def oral(self) -> bool:
        return self.ka is not None

    @property
    def nonlinear(self) -> bool:
        return self.Vmax is not None and self.Km is not None

    def derived(self) -> dict:
        """Secondary parameters (rate constants, half-lives)."""
        out: dict[str, float] = {}
        k = self.CL / self.V1
        out["k10"] = k
        out["t_half_1cmt"] = math.log(2) / k if k > 0 else float("inf")
        if self.two_compartment:
            k12 = self.Q / self.V1
            k21 = self.Q / self.V2
            out["k12"], out["k21"] = k12, k21
            s = k + k12 + k21
            disc = max(s * s - 4 * k * k21, 0.0)
            alpha = 0.5 * (s + math.sqrt(disc))
            beta = 0.5 * (s - math.sqrt(disc))
            out["alpha"], out["beta"] = alpha, beta
            out["t_half_alpha"] = math.log(2) / alpha if alpha > 0 else float("inf")
            out["t_half_beta"] = math.log(2) / beta if beta > 0 else float("inf")
            out["Vss"] = self.V1 + self.V2
        if self.oral:
            out["ka"] = self.ka
        return out


class CompartmentModel:
    """Compartment layout + ODE right-hand side, driven by dosing events.

    State vector ``amounts`` holds drug *amount* (not concentration) in each
    compartment. Compartment 0 is the depot for oral dosing; otherwise the
    central compartment is index 0.
    """

    def __init__(self, params: PKParameters):
        self.params = params
        # Depot compartment prepended only for oral absorption.
        self.has_depot = params.oral
        self.central = 1 if self.has_depot else 0
        self.peripheral = self.central + 1 if params.two_compartment else None
        self.n = (1 if self.has_depot else 0) + 1 + (1 if params.two_compartment else 0)

    def central_concentration(self, amounts: Sequence[float]) -> float:
        return amounts[self.central] / self.params.V1

    def rhs(self, amounts: Sequence[float], infusion_rate: float) -> list[float]:
        """Derivatives d(amount)/dt. ``infusion_rate`` [mg/h] enters the central
        compartment (zero-order input)."""
        p = self.params
        d = [0.0] * self.n
        c = self.central
        conc = amounts[c] / p.V1

        # Elimination from central compartment.
        if p.nonlinear:
            elim = p.Vmax * conc / (p.Km + conc)
        else:
            elim = p.CL * conc

        d[c] += infusion_rate - elim

        # Oral absorption: depot -> central.
        if self.has_depot:
            ka_flux = p.ka * amounts[0]
            d[0] -= ka_flux
            d[c] += ka_flux

        # Peripheral distribution.
        if self.peripheral is not None:
            k12 = p.Q / p.V1
            k21 = p.Q / p.V2
            central_to_periph = k12 * amounts[c]
            periph_to_central = k21 * amounts[self.peripheral]
            d[c] += periph_to_central - central_to_periph
            d[self.peripheral] += central_to_periph - periph_to_central

        return d


# ---------------------------------------------------------------------------
# Closed-form single-dose solutions (used for validation and fast paths).
# ---------------------------------------------------------------------------

def analytical_concentration(
    params: PKParameters, times: Sequence[float], dose: float,
    route: Route = "iv_bolus", tinf: float = 0.0,
) -> list[float]:
    """Single-dose analytical concentration-time profile [mg/L].

    Supports linear 1- and 2-compartment models with IV bolus, IV infusion, or
    first-order oral absorption. Raises for nonlinear (Michaelis-Menten) models,
    which have no closed form -- integrate those instead.
    """
    if params.nonlinear:
        raise ValueError("nonlinear elimination has no closed form; use simulate()")
    if params.two_compartment:
        return _two_cmt(params, times, dose, route, tinf)
    return _one_cmt(params, times, dose, route, tinf)


def _one_cmt(params, times, dose, route, tinf):
    p = params
    k = p.CL / p.V1
    out = []
    if route == "oral":
        ka = p.ka
        for t in times:
            if t <= 0:
                out.append(0.0)
                continue
            if abs(ka - k) < 1e-9:  # flip-flop degenerate limit
                c = (p.F * dose / p.V1) * k * t * math.exp(-k * t)
            else:
                c = (p.F * dose * ka) / (p.V1 * (ka - k)) * (
                    math.exp(-k * t) - math.exp(-ka * t))
            out.append(max(c, 0.0))
    elif route == "infusion":
        R0 = dose / tinf
        for t in times:
            if t <= 0:
                out.append(0.0)
            elif t < tinf:
                out.append(R0 / p.CL * (1 - math.exp(-k * t)))
            else:
                cmax = R0 / p.CL * (1 - math.exp(-k * tinf))
                out.append(cmax * math.exp(-k * (t - tinf)))
    else:  # iv_bolus
        c0 = dose / p.V1
        out = [c0 * math.exp(-k * t) if t >= 0 else 0.0 for t in times]
    return out


def _two_cmt(params, times, dose, route, tinf):
    p = params
    d = p.derived()
    alpha, beta = d["alpha"], d["beta"]
    k21 = d["k21"]
    out = []
    if route == "iv_bolus":
        A = dose / p.V1 * (alpha - k21) / (alpha - beta)
        B = dose / p.V1 * (k21 - beta) / (alpha - beta)
        for t in times:
            if t < 0:
                out.append(0.0)
            else:
                out.append(A * math.exp(-alpha * t) + B * math.exp(-beta * t))
    elif route == "infusion":
        R0 = dose / tinf
        # Coefficients for infusion input into a 2-cmt model.
        A = (k21 - alpha) / (p.V1 * (beta - alpha))
        B = (k21 - beta) / (p.V1 * (alpha - beta))
        for t in times:
            if t <= 0:
                out.append(0.0)
            elif t < tinf:
                c = R0 * (A / alpha * (1 - math.exp(-alpha * t))
                          + B / beta * (1 - math.exp(-beta * t)))
                out.append(c)
            else:
                ca = R0 * A / alpha * (1 - math.exp(-alpha * tinf)) * math.exp(-alpha * (t - tinf))
                cb = R0 * B / beta * (1 - math.exp(-beta * tinf)) * math.exp(-beta * (t - tinf))
                out.append(ca + cb)
    else:  # oral into a 2-cmt model
        ka = p.ka
        pref = p.F * dose * ka / p.V1
        for t in times:
            if t <= 0:
                out.append(0.0)
                continue
            c = pref * (
                (k21 - alpha) / ((ka - alpha) * (beta - alpha)) * math.exp(-alpha * t)
                + (k21 - beta) / ((ka - beta) * (alpha - beta)) * math.exp(-beta * t)
                + (k21 - ka) / ((alpha - ka) * (beta - ka)) * math.exp(-ka * t))
            out.append(max(c, 0.0))
    return out
