"""Derivative-free optimisation used for parameter estimation.

A pure-Python Nelder-Mead simplex minimiser with multi-restart support. Kept
dependency-free so the estimation layer runs without SciPy; if NumPy/SciPy are
available the estimation module may prefer them, but this is the guaranteed
fallback.
"""
from __future__ import annotations

import math
from typing import Callable, Sequence


def nelder_mead(
    func: Callable[[list[float]], float],
    x0: Sequence[float],
    step: float = 0.2,
    max_iter: int = 2000,
    xatol: float = 1e-7,
    fatol: float = 1e-9,
) -> tuple[list[float], float]:
    """Minimise ``func`` starting from ``x0``. Returns (best_x, best_f)."""
    n = len(x0)
    x0 = [float(v) for v in x0]

    # Build the initial simplex by perturbing each coordinate.
    simplex = [list(x0)]
    for i in range(n):
        pt = list(x0)
        h = step * (abs(pt[i]) if pt[i] != 0 else 1.0)
        pt[i] += h
        simplex.append(pt)

    fvals = [func(p) for p in simplex]

    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5

    for _ in range(max_iter):
        # Order by objective value.
        order = sorted(range(n + 1), key=lambda i: fvals[i])
        simplex = [simplex[i] for i in order]
        fvals = [fvals[i] for i in order]

        # Convergence: simplex spread small in both x and f.
        fspread = abs(fvals[-1] - fvals[0])
        xspread = max(
            max(abs(simplex[i][j] - simplex[0][j]) for i in range(1, n + 1))
            for j in range(n)
        ) if n else 0.0
        if fspread <= fatol and xspread <= xatol:
            break

        # Centroid of all but the worst point.
        centroid = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]

        # Reflection.
        xr = [centroid[j] + alpha * (centroid[j] - simplex[-1][j]) for j in range(n)]
        fr = func(xr)
        if fvals[0] <= fr < fvals[-2]:
            simplex[-1], fvals[-1] = xr, fr
            continue

        # Expansion.
        if fr < fvals[0]:
            xe = [centroid[j] + gamma * (xr[j] - centroid[j]) for j in range(n)]
            fe = func(xe)
            if fe < fr:
                simplex[-1], fvals[-1] = xe, fe
            else:
                simplex[-1], fvals[-1] = xr, fr
            continue

        # Contraction.
        xc = [centroid[j] + rho * (simplex[-1][j] - centroid[j]) for j in range(n)]
        fc = func(xc)
        if fc < fvals[-1]:
            simplex[-1], fvals[-1] = xc, fc
            continue

        # Shrink towards best.
        best = simplex[0]
        for i in range(1, n + 1):
            simplex[i] = [best[j] + sigma * (simplex[i][j] - best[j]) for j in range(n)]
            fvals[i] = func(simplex[i])

    order = sorted(range(n + 1), key=lambda i: fvals[i])
    return simplex[order[0]], fvals[order[0]]


def minimize_restarts(func, x0, restarts=3, **kwargs):
    """Nelder-Mead with several restarts from the incumbent best, guarding
    against premature convergence to a poor local minimum."""
    best_x, best_f = nelder_mead(func, x0, **kwargs)
    for _ in range(restarts - 1):
        cand_x, cand_f = nelder_mead(func, best_x, **kwargs)
        if cand_f < best_f:
            best_x, best_f = cand_x, cand_f
    return best_x, best_f
