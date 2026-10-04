"""Synthetic-but-biologically-faithful benchmark.

Ground truth encodes the consensus CSC properties distilled from Rajasekhar
(2014) Ch.1-35 + 2026 reviews (CanSig; Signal Transduction & Targeted Therapy
2026; sc-best-practices):
  - rarity (1% CSC),
  - EMT-coupled stemness,
  - quiescence (low cell-cycle),
  - therapy enrichment (CSC 5x more likely to survive),
  - batch effect across 2 cohorts (tests reproducibility, CanSig-style).

Tasks:
  T1 stemness ranking (AUROC/AUPRC vs true CSC label)
  T2 EMT-stem coupling (Pearson r between stem and EMT scores in CSCs)
  T3 therapy-enrichment recovery (observed CSC fraction pre/post therapy)
  T4 ELDA frequency recovery (estimated 1-in-X vs true 1-in-100)
  T5 cross-batch reproducibility (Spearman of signature ranks across batches)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from .elda import elda_estimate
from .stemness import composite_csc_score, emt_score, quiescence_score, stemness_signature_score

RNG_SEED = 20260214
GENES = (
    ["PROM1", "ALDH1A1", "ALDH1A3", "CD44", "BMI1", "SOX2", "NANOG", "POU5F1", "KLF4", "NOTCH1"]
    + ["CDH1", "EPCAM", "KRT19", "CLDN4", "OCLN"]
    + ["VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "FN1"]
    + ["MKI67", "TOP2A", "CCNB1", "PCNA", "MCM2"]
    + ["CDKN1A", "CDKN1B", "CDKN2A", "HES1"]
    + [f"HOUSE{i:02d}" for i in range(20)]
)


def simulate(n_cells: int = 4000, csc_frac: float = 0.01, seed: int = RNG_SEED, batch_shift: float = 0.35):
    rng = np.random.default_rng(seed)
    n_csc = int(n_cells * csc_frac)
    y = np.zeros(n_cells, dtype=int)
    y[rng.choice(n_cells, size=n_csc, replace=False)] = 1
    batch = (rng.random(n_cells) < 0.5).astype(int)
    base = rng.normal(0, 1, size=(n_cells, len(GENES)))
    df = pd.DataFrame(base, columns=GENES)
    is_csc = y == 1
    # stem program up in CSC (moderate effect => realistic overlap, not separable)
    for g in ["PROM1", "ALDH1A1", "CD44", "BMI1", "SOX2", "NANOG", "POU5F1"]:
        df.loc[is_csc, g] += 1.4
    # EMT coupled in CSC
    for g in ["VIM", "ZEB1", "SNAI1", "TWIST1", "FN1", "CDH2"]:
        df.loc[is_csc, g] += 1.0
    for g in ["CDH1", "EPCAM", "CLDN4", "OCLN"]:
        df.loc[is_csc, g] -= 1.2
    # quiescence: low cycling, high CDKN
    for g in ["MKI67", "TOP2A", "CCNB1", "PCNA"]:
        df.loc[is_csc, g] -= 1.0
    for g in ["CDKN1A", "CDKN1B", "HES1"]:
        df.loc[is_csc, g] += 1.0
    # batch effect on housekeeping + some programs
    df.loc[batch == 1, [f"HOUSE{i:02d}" for i in range(20)]] += batch_shift
    df.loc[batch == 1, ["VIM", "FN1"]] += batch_shift * 0.5
    # therapy: CSC 5x survival
    p_surv = np.where(is_csc, 0.55, 0.11)
    survived = rng.random(n_cells) < p_surv
    return df, pd.Series(y, name="is_csc"), pd.Series(batch, name="batch"), pd.Series(survived, name="survived")


def run_synthetic_benchmark(n_cells: int = 4000, seed: int = RNG_SEED) -> dict:
    expr, y, batch, survived = simulate(n_cells=n_cells, seed=seed)
    s_stem = stemness_signature_score(expr)
    s_emt = emt_score(expr)
    s_qui = quiescence_score(expr)
    s_comp = composite_csc_score(expr)

    def _safe_auc(a, b):
        try:
            return float(roc_auc_score(a, b))
        except Exception:
            return float("nan")

    def _safe_ap(a, b):
        try:
            return float(average_precision_score(a, b))
        except Exception:
            return float("nan")

    t1 = {
        "auroc_stem": _safe_auc(y, s_stem),
        "auprc_stem": _safe_ap(y, s_stem),
        "auroc_composite": _safe_auc(y, s_comp),
        "auprc_composite": _safe_ap(y, s_comp),
    }
    csc_mask = y.values == 1
    t2_r = float(np.corrcoef(s_stem.values[csc_mask], s_emt.values[csc_mask])[0, 1]) if csc_mask.sum() > 3 else float("nan")
    pre = float(y.mean())
    post = float(y[survived.values == 1].mean()) if survived.sum() > 0 else float("nan")
    t3 = {"csc_frac_pre": pre, "csc_frac_post": post, "enrichment": float(post / pre) if pre > 0 else float("nan")}
    # T4: ELDA on simulated LDA (doses around true 1/100)
    doses = np.array([500, 100, 20, 5, 1])
    lam_true = y.mean()
    tested = np.array([24, 24, 48, 48, 96])
    rng = np.random.default_rng(seed + 1)
    p_pos = 1 - np.exp(-lam_true * doses)
    positive = rng.binomial(tested, p_pos)
    elda = elda_estimate(doses, tested, positive)
    t4 = {"true_one_in": float(1 / lam_true), "elda_one_in": elda["one_in"], "elda_ci95_one_in": [float(1 / elda["ci95"][1]), float(1 / elda["ci95"][0])]}
    # T5: rank reproducibility across batches using stem signature gene ranking
    ranks = {}
    for b in (0, 1):
        m = batch.values == b
        means = expr.loc[m].mean(axis=0)
        ranks[b] = means.rank()
    rho, _ = spearmanr(ranks[0], ranks[1])
    t5 = {"spearman_batch_rank": float(rho)}
    # hidden-pattern probes (reported, not tuned)
    qui_csc = float(s_qui.values[csc_mask].mean())
    qui_non = float(s_qui.values[~csc_mask].mean())
    emt_csc = float(s_emt.values[csc_mask].mean())
    emt_non = float(s_emt.values[~csc_mask].mean())
    patterns = {
        "quiescence_gap_csc_minus_bulk": qui_csc - qui_non,
        "emt_gap_csc_minus_bulk": emt_csc - emt_non,
        "therapy_enrichment_fold": t3["enrichment"],
        "batch_rank_rho": float(rho),
    }
    return {"T1_ranking": t1, "T2_emt_coupling_r": t2_r, "T3_therapy": t3, "T4_elda": t4, "T5_reproducibility": t5, "patterns": patterns,
            "n_cells": n_cells, "n_csc": int(y.sum()), "seed": seed}
