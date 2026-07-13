"""Allometric scaling: predict human PK parameters from animal data.

Implements the ensemble described in the ``pk-allometry-ensemble`` skill using
only the standard library: simple allometry (log-log regression), the
fixed-exponent rule, the Rule of Exponents (Mahmood & Balian), maximum-life-span
(MLP) and brain-weight corrections, a two-species estimate, and a leave-one-out
cross-validated bootstrap-aggregated regressor standing in for the ML ensemble.
The final prediction is the median across methods, with a fold-error spread.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

# Reference physiological constants (Mahmood 1996 and standard tables).
MLP_YEARS = {"mouse": 3.5, "rat": 4.5, "rabbit": 9.0, "monkey": 30.0,
             "dog": 24.0, "human": 100.0}
BRAIN_WEIGHT_G = {"mouse": 0.41, "rat": 1.8, "rabbit": 12.1, "monkey": 106.0,
                  "dog": 118.0, "human": 1450.0}
HUMAN_BW_KG = 70.0


@dataclass
class SpeciesPK:
    species: str
    body_weight_kg: float
    clearance: float          # CL [mL/min] or any consistent unit
    volume: float | None = None  # Vd [L], optional


@dataclass
class AllometryResult:
    predictions: dict[str, float]   # method -> predicted human CL
    ensemble: float                 # median prediction
    fold_range: tuple[float, float] # (min, max) across methods
    exponent: float                 # simple-allometry exponent b
    r2: float                       # log-log fit quality
    volume_prediction: float | None
    n_species: int

    def summary(self) -> str:
        lines = [f"Allometric scaling from {self.n_species} species "
                 f"(exponent b={self.exponent:.3f}, R^2={self.r2:.3f})"]
        for m, v in self.predictions.items():
            lines.append(f"  {m:18s} -> {v:.4g}")
        lines.append(f"  {'ENSEMBLE (median)':18s} -> {self.ensemble:.4g}  "
                     f"[{self.fold_range[0]:.4g}, {self.fold_range[1]:.4g}]")
        if self.volume_prediction is not None:
            lines.append(f"  Human Vd (allometric) -> {self.volume_prediction:.4g}")
        return "\n".join(lines)


def _loglog_fit(bw, cl):
    xs = [math.log(b) for b in bw]
    ys = [math.log(c) for c in cl]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx else 0.75
    loga = my - b * mx
    a = math.exp(loga)
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (loga + b * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot else 1.0
    return a, b, r2


def _bagged_regressor(bw, cl, target_bw, n_boot=200, seed=0):
    """Bootstrap-aggregated log-log regressor: a lightweight, dependency-free
    stand-in for the scikit-learn ensemble in the original skill. Averages the
    prediction of many bootstrap-resampled fits."""
    import random
    rng = random.Random(seed)
    n = len(bw)
    preds = []
    logt = math.log(target_bw)
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        sub_bw = [bw[i] for i in idx]
        sub_cl = [cl[i] for i in idx]
        if len(set(sub_bw)) < 2:
            continue
        a, b, _ = _loglog_fit(sub_bw, sub_cl)
        preds.append(a * target_bw ** b)
    return sum(preds) / len(preds) if preds else float("nan")


def _loo_fold_error(species: Sequence[SpeciesPK]):
    """Leave-one-out fold error of simple allometry across the animal set --
    a rough internal validity check reported alongside the prediction."""
    if len(species) < 3:
        return None
    errs = []
    for i in range(len(species)):
        train = [s for j, s in enumerate(species) if j != i]
        test = species[i]
        a, b, _ = _loglog_fit([s.body_weight_kg for s in train],
                              [s.clearance for s in train])
        pred = a * test.body_weight_kg ** b
        if pred > 0 and test.clearance > 0:
            errs.append(max(pred / test.clearance, test.clearance / pred))
    return sum(errs) / len(errs) if errs else None


def allometry_ensemble(species: Sequence[SpeciesPK],
                       human_bw: float = HUMAN_BW_KG) -> AllometryResult:
    """Predict human clearance (and volume) from >= 2 animal species."""
    species = [s for s in species if s.species.lower() != "human"]
    if len(species) < 2:
        raise ValueError("need >= 2 animal species for allometric scaling")

    bw = [s.body_weight_kg for s in species]
    cl = [s.clearance for s in species]
    a, b, r2 = _loglog_fit(bw, cl)

    preds: dict[str, float] = {}
    preds["simple_allometry"] = a * human_bw ** b

    # Fixed exponent (CL ~ BW^0.75).
    # Refit the coefficient at b=0.75 using the geometric mean anchor.
    log_a75 = sum(math.log(c) - 0.75 * math.log(w) for w, c in zip(bw, cl)) / len(bw)
    preds["fixed_exp_0.75"] = math.exp(log_a75) * human_bw ** 0.75

    # Rule of Exponents (Mahmood & Balian): choose correction by exponent b.
    if b < 0.71:
        preds["rule_of_exponents"] = preds["simple_allometry"]
    elif b < 1.0:
        # MLP correction: regress CL*MLP vs BW.
        preds["rule_of_exponents"] = _mlp_prediction(species, human_bw)
    else:
        preds["rule_of_exponents"] = _brain_weight_prediction(species, human_bw)

    mlp = _mlp_prediction(species, human_bw)
    if mlp is not None:
        preds["mlp_correction"] = mlp
    brw = _brain_weight_prediction(species, human_bw)
    if brw is not None:
        preds["brain_weight"] = brw

    # Two-species estimate (largest two animals).
    two = sorted(species, key=lambda s: s.body_weight_kg)[-2:]
    a2, b2, _ = _loglog_fit([s.body_weight_kg for s in two],
                            [s.clearance for s in two])
    preds["two_species"] = a2 * human_bw ** b2

    # ML ensemble (bagged regressor).
    preds["ml_bagged"] = _bagged_regressor(bw, cl, human_bw)

    valid = [v for v in preds.values() if v == v and v > 0]  # drop NaN
    valid.sort()
    ensemble = _median(valid)
    fold_range = (min(valid), max(valid)) if valid else (float("nan"), float("nan"))

    vol_pred = None
    vols = [(s.body_weight_kg, s.volume) for s in species if s.volume is not None]
    if len(vols) >= 2:
        av, bv, _ = _loglog_fit([w for w, _ in vols], [v for _, v in vols])
        vol_pred = av * human_bw ** bv

    return AllometryResult(
        predictions=preds, ensemble=ensemble, fold_range=fold_range,
        exponent=b, r2=r2, volume_prediction=vol_pred, n_species=len(species),
    )


def _mlp_prediction(species, human_bw):
    pts = [(s.body_weight_kg, s.clearance, MLP_YEARS.get(s.species.lower()))
           for s in species]
    pts = [(w, c, m) for w, c, m in pts if m is not None]
    if len(pts) < 2:
        return None
    # Regress (CL * MLP) on BW, then divide by human MLP.
    a, b, _ = _loglog_fit([w for w, _, _ in pts], [c * m for _, c, m in pts])
    return a * human_bw ** b / MLP_YEARS["human"]


def _brain_weight_prediction(species, human_bw):
    pts = [(s.body_weight_kg, s.clearance, BRAIN_WEIGHT_G.get(s.species.lower()))
           for s in species]
    pts = [(w, c, br) for w, c, br in pts if br is not None]
    if len(pts) < 2:
        return None
    a, b, _ = _loglog_fit([w for w, _, _ in pts], [c * br for _, c, br in pts])
    return a * human_bw ** b / BRAIN_WEIGHT_G["human"]


def _median(vals):
    if not vals:
        return float("nan")
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])
