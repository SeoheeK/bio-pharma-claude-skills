"""Parameter estimation for compartmental PK models.

Fits structural parameters (CL, V, Q, V2, ka, ...) to observed concentration
data by minimising a weighted residual sum of squares. Parameters are estimated
in log-space so they stay positive. A two-stage population routine fits each
subject individually and summarises the population as geometric-mean typical
values with between-subject CV%.

This is a transparent maximum-likelihood-style fit (individual/two-stage NLME),
not a full FOCE mixed-effects engine -- but it produces the THETA/OMEGA-style
summaries the skills reference, entirely in pure Python.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

from .models import PKParameters, analytical_concentration
from .regimen import DosingRegimen
from .simulate import simulate
from .optim import minimize_restarts


@dataclass
class FitResult:
    params: PKParameters
    estimates: dict[str, float]      # fitted parameter values
    rse_pct: dict[str, float]        # relative standard error (%)
    objective: float                 # final weighted SSR
    n_obs: int
    weighting: str
    predicted: list[float] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"Fit ({self.weighting} weighting), {self.n_obs} obs, OFV={self.objective:.4g}"]
        for k, v in self.estimates.items():
            rse = self.rse_pct.get(k)
            rse_s = f"  (RSE {rse:.1f}%)" if rse is not None and math.isfinite(rse) else ""
            lines.append(f"  {k:5s} = {v:.4g}{rse_s}")
        return "\n".join(lines)


def _predict(param_names, values, template: PKParameters, regimen, times):
    kwargs = dict(
        CL=template.CL, V1=template.V1, Q=template.Q, V2=template.V2,
        ka=template.ka, F=template.F, Vmax=template.Vmax, Km=template.Km,
    )
    for name, val in zip(param_names, values):
        kwargs[name] = val
    p = PKParameters(**kwargs)
    # Analytical fast path for single-dose linear models; else integrate.
    single = len(regimen.events) == 1
    try:
        if single and not p.nonlinear:
            e = regimen.events[0]
            return analytical_concentration(p, times, e.amount, e.route, e.tinf), p
    except Exception:
        pass
    return simulate(p, regimen, times=times).concentrations, p


def fit_individual(
    times: Sequence[float],
    concentrations: Sequence[float],
    regimen: DosingRegimen,
    estimate: Sequence[str] = ("CL", "V1"),
    init: PKParameters | None = None,
    weighting: str = "proportional",
) -> FitResult:
    """Fit structural parameters to one subject's data.

    ``estimate`` names the parameters to fit (others held at ``init``).
    ``weighting``: "proportional" (fit in log space, constant CV error),
    "uniform" (ordinary least squares), or "poisson" (1/pred weighting).
    """
    times = list(times)
    obs = list(concentrations)
    names = list(estimate)
    template = init or _auto_init(times, obs, regimen, names)

    x0 = [math.log(max(getattr(template, n) or 1e-6, 1e-9)) for n in names]

    def objective(logx):
        vals = [math.exp(v) for v in logx]
        pred, _ = _predict(names, vals, template, regimen, times)
        ssr = 0.0
        for o, pr in zip(obs, pred):
            if o <= 0:
                continue
            if weighting == "proportional":
                pr = max(pr, 1e-9)
                ssr += (math.log(o) - math.log(pr)) ** 2
            elif weighting == "poisson":
                pr = max(pr, 1e-9)
                ssr += (o - pr) ** 2 / pr
            else:  # uniform
                ssr += (o - pr) ** 2
        return ssr

    best_log, best_obj = minimize_restarts(objective, x0, restarts=3)
    best_vals = [math.exp(v) for v in best_log]

    pred, fitted = _predict(names, best_vals, template, regimen, times)
    estimates = dict(zip(names, best_vals))
    rse = _rse(objective, best_log, names, len(obs))

    return FitResult(
        params=fitted, estimates=estimates, rse_pct=rse, objective=best_obj,
        n_obs=len([o for o in obs if o > 0]), weighting=weighting, predicted=pred,
    )


def _auto_init(times, obs, regimen, names):
    """Cheap initial guess from NCA-style reasoning."""
    from .nca import nca
    dose = regimen.events[0].amount if regimen.events else 1.0
    route = regimen.events[0].route if regimen.events else "iv_bolus"
    try:
        r = nca(times, obs, dose=dose, route=route)
        cl = r.clearance or (dose / max(r.auc_inf, 1e-6))
        v = (r.vz or (cl * (r.t_half / math.log(2)) if math.isfinite(r.t_half) else cl * 4))
    except Exception:
        cl, v = 1.0, 10.0
    p = PKParameters(CL=max(cl, 1e-3), V1=max(v, 1e-3))
    if "Q" in names:
        p.Q = cl
    if "V2" in names:
        p.V2 = v
    if "ka" in names:
        p.ka = 1.0
    return p


def _rse(objective, best_log, names, n_obs):
    """Approximate relative standard errors from a finite-difference Hessian of
    the (log-parameter) objective. Rough but useful for flagging weak params."""
    k = len(names)
    if n_obs <= k:
        return {n: float("nan") for n in names}
    f0 = objective(best_log)
    sigma2 = f0 / max(n_obs - k, 1)
    h = 1e-4
    hess = [[0.0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i, k):
            xpp = list(best_log); xpp[i] += h; xpp[j] += h
            xpm = list(best_log); xpm[i] += h; xpm[j] -= h
            xmp = list(best_log); xmp[i] -= h; xmp[j] += h
            xmm = list(best_log); xmm[i] -= h; xmm[j] -= h
            val = (objective(xpp) - objective(xpm) - objective(xmp) + objective(xmm)) / (4 * h * h)
            hess[i][j] = hess[j][i] = 0.5 * val / max(sigma2, 1e-12)
    cov = _invert(hess)
    out = {}
    for i, name in enumerate(names):
        if cov is None or cov[i][i] < 0:
            out[name] = float("nan")
        else:
            # SE on log-param ~= relative SE on param.
            out[name] = 100 * math.sqrt(cov[i][i])
    return out


def _invert(mat):
    n = len(mat)
    a = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(mat)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-12:
            return None
        a[col], a[piv] = a[piv], a[col]
        pv = a[col][col]
        a[col] = [x / pv for x in a[col]]
        for r in range(n):
            if r != col and a[r][col] != 0:
                factor = a[r][col]
                a[r] = [x - factor * y for x, y in zip(a[r], a[col])]
    return [row[n:] for row in a]


# ---------------------------------------------------------------------------
# Two-stage population PK.
# ---------------------------------------------------------------------------

@dataclass
class PopPKResult:
    typical: dict[str, float]        # geometric-mean typical values (THETA-like)
    iiv_cv_pct: dict[str, float]     # between-subject CV% (OMEGA-like)
    individual: list[dict[str, float]]
    n_subjects: int

    def summary(self) -> str:
        lines = [f"Two-stage PopPK, N={self.n_subjects} subjects"]
        lines.append("  Typical values (geometric mean) and IIV CV%:")
        for k, v in self.typical.items():
            cv = self.iiv_cv_pct.get(k, float("nan"))
            lines.append(f"    {k:5s} = {v:.4g}   IIV CV {cv:.1f}%")
        return "\n".join(lines)


def fit_population(
    subjects: Sequence[dict],
    estimate: Sequence[str] = ("CL", "V1"),
    weighting: str = "proportional",
) -> PopPKResult:
    """Two-stage NLME. Each ``subjects`` entry is a dict with keys
    ``times``, ``concentrations``, and ``regimen`` (a DosingRegimen).
    """
    names = list(estimate)
    per_subject: list[dict[str, float]] = []
    for s in subjects:
        fit = fit_individual(
            s["times"], s["concentrations"], s["regimen"],
            estimate=names, weighting=weighting,
        )
        per_subject.append(fit.estimates)

    typical: dict[str, float] = {}
    iiv: dict[str, float] = {}
    for name in names:
        vals = [d[name] for d in per_subject if d.get(name, 0) > 0]
        if not vals:
            typical[name] = float("nan"); iiv[name] = float("nan"); continue
        logs = [math.log(v) for v in vals]
        m = sum(logs) / len(logs)
        typical[name] = math.exp(m)
        if len(logs) > 1:
            var = sum((x - m) ** 2 for x in logs) / (len(logs) - 1)
            iiv[name] = 100 * math.sqrt(math.exp(var) - 1)  # CV% from log-variance
        else:
            iiv[name] = float("nan")

    return PopPKResult(
        typical=typical, iiv_cv_pct=iiv, individual=per_subject,
        n_subjects=len(subjects),
    )
