"""End-to-end benchmark entry point.

Runs:
  1. synthetic CSC benchmark (deterministic, seed-fixed),
  2. live public-data verification (GDC + OpenAlex, cached/offline-safe),
Writes:
  results/benchmark_results.json, results/summary.md, data/live_verification.json
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from csc_benchmark.benchmark import run_synthetic_benchmark
from csc_benchmark.live_data import live_verification

ROOT = os.path.join(os.path.dirname(__file__), "..")
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")


def main() -> int:
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)
    bench = run_synthetic_benchmark()
    live = live_verification()
    with open(os.path.join(DATA, "live_verification.json"), "w") as f:
        json.dump(live, f, indent=2)
    payload = {"synthetic_benchmark": bench, "live_verification": live}
    with open(os.path.join(RESULTS, "benchmark_results.json"), "w") as f:
        json.dump(payload, f, indent=2)
    t1 = bench["T1_ranking"]
    lines = [
        "# Benchmark summary (auto-generated, do not hand-edit)",
        "",
        f"- cells={bench['n_cells']} csc={bench['n_csc']} seed={bench['seed']}",
        f"- T1 AUROC stem={t1['auroc_stem']:.3f} composite={t1['auroc_composite']:.3f}",
        f"- T1 AUPRC stem={t1['auprc_stem']:.4f} composite={t1['auprc_composite']:.4f}",
        f"- T2 EMT-stem r in CSC = {bench['T2_emt_coupling_r']:.3f}",
        f"- T3 pre={bench['T3_therapy']['csc_frac_pre']:.4f} post={bench['T3_therapy']['csc_frac_post']:.4f} enrich={bench['T3_therapy']['enrichment']:.2f}x",
        f"- T4 true 1-in-{bench['T4_elda']['true_one_in']:.1f} ELDA 1-in-{bench['T4_elda']['elda_one_in']:.1f}",
        f"- T5 batch rank rho={bench['T5_reproducibility']['spearman_batch_rank']:.3f}",
        f"- live_ok={live['live_ok']} gdc={live['gdc']} openalex={live['openalex']}",
        "",
    ]
    with open(os.path.join(RESULTS, "summary.md"), "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))
    # fail loudly if the synthetic signal is broken (guards against silent edits)
    assert t1["auroc_composite"] > 0.95, "composite ranker broken"
    assert bench["T3_therapy"]["enrichment"] > 2.0, "therapy enrichment broken"
    assert bench["T4_elda"]["elda_one_in"] > 20, "ELDA recovery broken"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
