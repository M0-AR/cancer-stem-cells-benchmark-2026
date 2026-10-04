"""Stemness / EMT / quiescence scoring.

Implements lightweight, dependency-minimal scores inspired by 2026 best
practices (CanSig signature rediscovery, sc-best-practices QC, ALDH/CD44
literature) without requiring scanpy.

All functions operate on a cells x genes DataFrame of log-normalised counts.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

STEM_GENES = ["PROM1", "ALDH1A1", "ALDH1A3", "CD44", "BMI1", "SOX2", "NANOG", "POU5F1", "KLF4", "NOTCH1"]
EPITHELIAL_GENES = ["CDH1", "EPCAM", "KRT19", "CLDN4", "OCLN"]
MESENCHYMAL_GENES = ["VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "FN1"]
CELL_CYCLE_GENES = ["MKI67", "TOP2A", "CCNB1", "PCNA", "MCM2"]
QUIESCENCE_GENES = ["CDKN1A", "CDKN1B", "CDKN2A", "HES1"]


def _zscore(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=0)
    sd = df.std(axis=0, ddof=0).replace(0, 1.0)
    return (df - mu) / sd


def _mean_available(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    avail = [g for g in genes if g in df.columns]
    if not avail:
        return pd.Series(np.zeros(len(df)), index=df.index)
    z = _zscore(df[avail])
    return z.mean(axis=1)


def stemness_signature_score(expr: pd.DataFrame, stem_genes: list[str] | None = None) -> pd.Series:
    """Mean z-scored stem signature. Higher = more stem-like."""
    return _mean_available(expr, stem_genes or STEM_GENES)


def emt_score(expr: pd.DataFrame) -> pd.Series:
    """Mesenchymal minus epithelial program (z-scored means)."""
    mes = _mean_available(expr, MESENCHYMAL_GENES)
    epi = _mean_available(expr, EPITHELIAL_GENES)
    return mes - epi


def quiescence_score(expr: pd.DataFrame) -> pd.Series:
    """Quiescence program minus cell-cycle program."""
    q = _mean_available(expr, QUIESCENCE_GENES)
    c = _mean_available(expr, CELL_CYCLE_GENES)
    return q - c


def composite_csc_score(expr: pd.DataFrame, w_stem: float = 1.0, w_emt: float = 0.6, w_qui: float = 0.3) -> pd.Series:
    """Weighted composite used as the primary CSC ranker in the benchmark."""
    s = stemness_signature_score(expr)
    e = emt_score(expr)
    q = quiescence_score(expr)

    def _norm(x: pd.Series) -> pd.Series:
        sd = x.std(ddof=0)
        if sd == 0 or np.isnan(sd):
            return x * 0.0
        return (x - x.mean()) / sd

    return w_stem * _norm(s) + w_emt * _norm(e) + w_qui * _norm(q)
