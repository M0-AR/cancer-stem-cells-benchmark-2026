# Data

This repo deliberately avoids committing bulk sequencing data.

- **Synthetic ground-truth cohort** (`src/csc_benchmark/benchmark.py::simulate`):
  deterministic (seed `20260214`), n=4000 cells, 1% CSC, EMT-coupled,
  quiescent, therapy-enriched, 2-batch. This is the primary benchmark substrate
  so anyone can reproduce T1–T5 without dbGaP access.
- **Live public-data verification** (`data/live_verification.json`, generated):
  small polite API calls, no bulk download —
  - NCI GDC `api.gdc.cancer.gov/cases` for `TCGA-PAAD` + `TCGA-BRCA` counts,
  - OpenAlex `api.openalex.org/works` for `cancer stem cell` 2025+ counts.
  These prove the study is grounded in real, reachable 2026 public substrates
  (GDC Portal, PanCanAtlas, CRDC). If offline, the runner caches the failure
  explicitly instead of inventing numbers.
- **To extend to real scRNA** (for your PhD follow-up, not required to reproduce):
  CanSig (2025, Cancer Research, 12 scRNA cohorts, 174k malignant cells),
  sc-best-practices (Heumos et al. 2023 + living book), GEO accessions cited in
  the README. See `paper/` for the omics extension roadmap.
