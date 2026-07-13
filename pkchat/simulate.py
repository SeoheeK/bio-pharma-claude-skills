"""Concentration-time simulation.

A fixed-step RK4 integrator advances the compartment ODE system, applying dose
events as they occur: IV boluses add instantaneously to the central amount,
infusions contribute a zero-order rate over their window, and oral doses load
the depot compartment. Population simulation layers log-normal between-subject
variability (IIV) and a residual error model on top.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Sequence

from .models import CompartmentModel, PKParameters
from .regimen import DosingRegimen


@dataclass
class SimulationResult:
    times: list[float]
    concentrations: list[float]

    def cmax(self) -> float:
        return max(self.concentrations) if self.concentrations else 0.0

    def tmax(self) -> float:
        if not self.concentrations:
            return 0.0
        i = max(range(len(self.concentrations)), key=lambda j: self.concentrations[j])
        return self.times[i]

    def as_rows(self) -> list[tuple[float, float]]:
        return list(zip(self.times, self.concentrations))


def _active_infusion_rate(regimen: DosingRegimen, t: float) -> float:
    """Total zero-order input rate [mg/h] active at time ``t``."""
    rate = 0.0
    for e in regimen.events:
        if e.route == "infusion" and e.tinf > 0 and e.time <= t < e.time + e.tinf:
            rate += e.amount / e.tinf
    return rate


def simulate(
    params: PKParameters,
    regimen: DosingRegimen,
    times: Sequence[float] | None = None,
    t_end: float | None = None,
    dt: float = 0.01,
) -> SimulationResult:
    """Integrate the model over a dosing regimen and report central-compartment
    concentration at the requested output ``times``.

    ``dt`` is the internal integration step (small for accuracy); output is
    sampled at ``times`` (defaults to an even grid to ``t_end``).
    """
    model = CompartmentModel(params)
    if t_end is None:
        # Default horizon: a few elimination half-lives past the last dose.
        d = params.derived()
        thalf = d.get("t_half_beta", d.get("t_half_1cmt", 4.0))
        t_end = regimen.duration() + 5 * min(thalf, 100.0) + 1e-9
    if times is None:
        n = max(int(t_end / max(dt, 1e-6)) + 1, 2)
        step = t_end / (n - 1)
        times = [i * step for i in range(n)]
    times = list(times)

    # Instantaneous dose events (bolus, oral) sorted by time.
    bolus_events = sorted(
        [e for e in regimen.events if e.route in ("iv_bolus", "oral")],
        key=lambda e: e.time,
    )

    amounts = [0.0] * model.n
    out_conc: list[float] = []
    out_idx = 0
    sorted_times = sorted(range(len(times)), key=lambda i: times[i])
    order = iter(sorted_times)
    pending = [next(order, None)]
    results: dict[int, float] = {}

    t = 0.0
    horizon = max(t_end, max(times) if times else 0.0)
    be_idx = 0

    def apply_bolus(now: float):
        nonlocal be_idx
        while be_idx < len(bolus_events) and abs(bolus_events[be_idx].time - now) < dt / 2:
            e = bolus_events[be_idx]
            if e.route == "oral":
                amounts[0] += e.amount * params.F  # into depot (index 0)
            else:  # iv_bolus into central
                amounts[model.central] += e.amount
            be_idx += 1

    # Record output at t=0 after any t=0 boluses.
    def record_due(now: float):
        for i, tt in enumerate(times):
            if i not in results and tt <= now + dt / 2:
                results[i] = model.central_concentration(amounts)

    apply_bolus(0.0)
    record_due(0.0)

    steps = int(math.ceil(horizon / dt)) + 1
    for _ in range(steps):
        if t >= horizon:
            break
        h = min(dt, horizon - t)
        if h <= 0:
            break
        # RK4 step with piecewise-constant infusion over the step.
        rate = _active_infusion_rate(regimen, t + h / 2)

        def f(state):
            return model.rhs(state, rate)

        k1 = f(amounts)
        k2 = f([a + h / 2 * d for a, d in zip(amounts, k1)])
        k3 = f([a + h / 2 * d for a, d in zip(amounts, k2)])
        k4 = f([a + h * d for a, d in zip(amounts, k3)])
        for i in range(model.n):
            amounts[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
            if amounts[i] < 0:
                amounts[i] = 0.0
        t += h
        apply_bolus(t)
        record_due(t)

    record_due(horizon)
    for i in range(len(times)):
        out_conc.append(results.get(i, model.central_concentration(amounts)))
    return SimulationResult(times=times, concentrations=out_conc)


@dataclass
class PopulationResult:
    times: list[float]
    profiles: list[list[float]]  # one concentration vector per virtual subject

    def percentile(self, q: float) -> list[float]:
        """Point-wise percentile band across subjects (q in [0,100])."""
        out = []
        for j in range(len(self.times)):
            col = sorted(p[j] for p in self.profiles)
            out.append(_percentile(col, q))
        return out

    def median(self) -> list[float]:
        return self.percentile(50)


def simulate_population(
    params: PKParameters,
    regimen: DosingRegimen,
    n_subjects: int = 100,
    iiv_cv: dict[str, float] | None = None,
    residual_cv: float = 0.0,
    times: Sequence[float] | None = None,
    t_end: float | None = None,
    dt: float = 0.02,
    seed: int | None = 0,
) -> PopulationResult:
    """Monte-Carlo population simulation with log-normal IIV on parameters and
    an optional proportional residual error.

    ``iiv_cv`` maps parameter names (e.g. "CL", "V1") to coefficient of
    variation (e.g. 0.3 for 30%). ``residual_cv`` adds proportional noise to
    each observation.
    """
    rng = random.Random(seed)
    iiv_cv = iiv_cv or {}
    profiles: list[list[float]] = []
    ref_times: list[float] | None = list(times) if times is not None else None

    for _ in range(n_subjects):
        kwargs = dict(
            CL=params.CL, V1=params.V1, Q=params.Q, V2=params.V2,
            ka=params.ka, F=params.F, Vmax=params.Vmax, Km=params.Km,
        )
        for name, cv in iiv_cv.items():
            base = kwargs.get(name)
            if base is None or cv <= 0:
                continue
            # Log-normal: exp(N(0, sigma)) with sigma from CV.
            sigma = math.sqrt(math.log(1 + cv * cv))
            kwargs[name] = base * math.exp(rng.gauss(0, sigma))
        subj = PKParameters(**kwargs)
        res = simulate(subj, regimen, times=ref_times, t_end=t_end, dt=dt)
        if ref_times is None:
            ref_times = res.times
        conc = res.concentrations
        if residual_cv > 0:
            conc = [max(c * (1 + rng.gauss(0, residual_cv)), 0.0) for c in conc]
        profiles.append(conc)

    return PopulationResult(times=ref_times or [], profiles=profiles)


def _percentile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = (q / 100) * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac
