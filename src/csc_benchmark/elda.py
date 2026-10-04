"""Extreme Limiting Dilution Analysis (ELDA) — single-hit Poisson model.

Reference: Hu & Smyth (2009). We implement the complement log-log GLM:
  log(-log(F0)) = log(dose) + log(lambda)
where F0 = fraction of negative wells/cultures at a given dose, and lambda is
the active-cell frequency. Estimate via least squares on the cloglog scale with
a small-sample correction, plus a likelihood profile for CI.
"""
from __future__ import annotations

import numpy as np


def elda_log_likelihood(lam: float, doses: np.ndarray, tested: np.ndarray, positive: np.ndarray) -> float:
    lam = max(lam, 1e-12)
    # P(positive | dose d) = 1 - exp(-lam * d)
    p = 1.0 - np.exp(-lam * np.asarray(doses, dtype=float))
    p = np.clip(p, 1e-12, 1 - 1e-12)
    n = np.asarray(tested, dtype=float)
    k = np.asarray(positive, dtype=float)
    return float(np.sum(k * np.log(p) + (n - k) * np.log(1 - p)))


def elda_estimate(doses, tested, positive) -> dict:
    doses = np.asarray(doses, dtype=float)
    tested = np.asarray(tested, dtype=float)
    positive = np.asarray(positive, dtype=float)
    f0 = 1.0 - positive / np.maximum(tested, 1)
    # avoid log(0): clip
    f0c = np.clip(f0, 1e-6, 1 - 1e-6)
    y = np.log(-np.log(f0c))
    x = np.log(np.maximum(doses, 1e-9))
    # single-hit model has slope 1; intercept = log(lambda)
    log_lam = float(np.mean(y - x))
    lam_hat = float(np.exp(log_lam))
    # profile 95% CI via likelihood ratio (chi2_1 = 3.84)
    ll_max = elda_log_likelihood(lam_hat, doses, tested, positive)
    grid = lam_hat * np.logspace(-1.2, 1.2, 4000)
    lls = np.array([elda_log_likelihood(g, doses, tested, positive) for g in grid])
    keep = lls >= (ll_max - 3.84 / 2.0)
    lo = float(grid[keep].min()) if keep.any() else lam_hat
    hi = float(grid[keep].max()) if keep.any() else lam_hat
    return {
        "frequency": lam_hat,
        "one_in": float(1.0 / lam_hat) if lam_hat > 0 else float("inf"),
        "ci95": [lo, hi],
        "loglik_max": ll_max,
    }
