"""Non-compartmental analysis (NCA).

Model-independent PK metrics from a concentration-time profile: AUC (linear and
log-linear trapezoidal), Cmax/Tmax, terminal half-life from log-linear
regression of the elimination phase, and derived clearance/volume for IV data.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Sequence


@dataclass
class NCAResult:
    cmax: float
    tmax: float
    auc_last: float          # AUC to last measurable concentration
    auc_inf: float           # AUC extrapolated to infinity
    lambda_z: float          # terminal elimination rate constant [1/h]
    t_half: float            # terminal half-life [h]
    auc_extrap_pct: float    # % of AUC_inf from extrapolation
    clearance: float | None  # CL = Dose/AUC_inf (IV only)
    vz: float | None         # Vz = Dose/(lambda_z * AUC_inf) (IV only)
    n_lambda_points: int

    def to_dict(self) -> dict:
        return asdict(self)


def _log_linear_terminal(times, concs, min_points=3):
    """Fit log(C) = a - lambda_z * t on the terminal phase.

    Picks the tail window (>= ``min_points`` points ending at the last positive
    concentration) that maximises adjusted R^2, mirroring standard NCA practice.
    """
    pts = [(t, c) for t, c in zip(times, concs) if c > 0]
    if len(pts) < min_points:
        return None
    best = None
    n = len(pts)
    for start in range(0, n - min_points + 1):
        window = pts[start:]
        if len(window) < min_points:
            continue
        xs = [t for t, _ in window]
        ys = [math.log(c) for _, c in window]
        slope, intercept, r2 = _linreg(xs, ys)
        if slope >= 0:  # elimination phase must have negative slope
            continue
        # Prefer later windows with strong fit; adjusted R^2 rewards more points.
        k = len(window)
        adj = 1 - (1 - r2) * (k - 1) / (k - 2) if k > 2 else r2
        score = adj
        if best is None or score > best[0]:
            best = (score, -slope, intercept, k)
    if best is None:
        return None
    _, lambda_z, intercept, k = best
    return lambda_z, intercept, k


def _linreg(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx if sxx else 0.0
    intercept = my - slope * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot else 1.0
    return slope, intercept, r2


def _auc_linear(times, concs):
    auc = 0.0
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        auc += 0.5 * (concs[i] + concs[i - 1]) * dt
    return auc


def _auc_linlog(times, concs):
    """Linear-up / log-down trapezoidal AUC (reduces bias on the decline)."""
    auc = 0.0
    for i in range(1, len(times)):
        c0, c1 = concs[i - 1], concs[i]
        dt = times[i] - times[i - 1]
        if c1 > 0 and c0 > 0 and c1 < c0:  # declining, both positive -> log trap
            auc += dt * (c0 - c1) / math.log(c0 / c1)
        else:
            auc += 0.5 * (c0 + c1) * dt
    return auc


def nca(
    times: Sequence[float],
    concentrations: Sequence[float],
    dose: float | None = None,
    route: str = "iv_bolus",
    method: str = "linlog",
) -> NCAResult:
    """Compute NCA metrics.

    ``method`` = "linear" or "linlog" (linear-up/log-down) for the AUC rule.
    ``dose`` enables CL/Vz for IV routes.
    """
    times = list(times)
    concs = list(concentrations)
    if len(times) != len(concs) or len(times) < 2:
        raise ValueError("need >= 2 aligned time/concentration points")

    cmax = max(concs)
    tmax = times[concs.index(cmax)]
    auc_last = _auc_linlog(times, concs) if method == "linlog" else _auc_linear(times, concs)

    fit = _log_linear_terminal(times, concs)
    if fit is None:
        lambda_z = float("nan")
        t_half = float("inf")
        auc_inf = auc_last
        auc_extrap = 0.0
        n_lambda = 0
    else:
        lambda_z, _intercept, n_lambda = fit
        t_half = math.log(2) / lambda_z
        c_last = next(c for c in reversed(concs) if c > 0)
        auc_tail = c_last / lambda_z
        auc_inf = auc_last + auc_tail
        auc_extrap = 100 * auc_tail / auc_inf if auc_inf else 0.0

    clearance = vz = None
    if dose is not None and route in ("iv_bolus", "infusion") and auc_inf > 0:
        clearance = dose / auc_inf
        if fit is not None and lambda_z > 0:
            vz = clearance / lambda_z

    return NCAResult(
        cmax=cmax, tmax=tmax, auc_last=auc_last, auc_inf=auc_inf,
        lambda_z=lambda_z, t_half=t_half, auc_extrap_pct=auc_extrap,
        clearance=clearance, vz=vz, n_lambda_points=n_lambda,
    )
