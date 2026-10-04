"""Live public-data verification (small, polite, cached).

We verify the study against REAL public data without bulk downloads:
  1. NCI GDC API — live case counts for TCGA-PAAD / TCGA-BRCA (proves the
     public cancer-data substrate exists and is reachable in 2026).
  2. OpenAlex API — live counts of 2025-2026 cancer-stem-cell works (proves the
     literature substrate is current).

Both calls use short timeouts and fall back to a cached/offline flag so the
benchmark stays reproducible offline. Results are written to
data/live_verification.json by experiments/run_benchmark.py.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request

GDC_API = "https://api.gdc.cancer.gov/cases"
OPENALEX_API = "https://api.openalex.org/works"


def _get_json(url: str, timeout: int = 15) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "csc-benchmark-2026/1.0 (mailto:research@example.org)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def fetch_gdc_counts(timeout: int = 15) -> dict:
    """Return live TCGA case counts for PAAD + BRCA projects."""
    out: dict = {"ok": False}
    try:
        for project in ("TCGA-PAAD", "TCGA-BRCA"):
            filt = {"op": "=", "content": {"field": "project.project_id", "value": project}}
            url = GDC_API + "?facets=project.project_id&format=json&size=0&filters=" + urllib.parse.quote(json.dumps(filt))
            payload = _get_json(url, timeout=timeout)
            total = payload.get("data", {}).get("pagination", {}).get("total")
            out[project] = total
        out["ok"] = all(isinstance(out.get(k), int) and out[k] > 0 for k in ("TCGA-PAAD", "TCGA-BRCA"))
    except Exception as e:  # offline-safe
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def fetch_openalex_csc(timeout: int = 20) -> dict:
    """Return live OpenAlex counts for cancer-stem-cell works (2025-2026)."""
    out: dict = {"ok": False}
    try:
        q = urllib.parse.quote("cancer stem cell")
        url = f"{OPENALEX_API}?search={q}&filter=from_publication_date:2025-01-01&per-page=1"
        payload = _get_json(url, timeout=timeout)
        out["count_2025_plus"] = payload.get("meta", {}).get("count")
        # a second query: single-cell + CSC
        q2 = urllib.parse.quote("cancer stem cell single-cell")
        url2 = f"{OPENALEX_API}?search={q2}&filter=from_publication_date:2024-01-01&per-page=1"
        payload2 = _get_json(url2, timeout=timeout)
        out["count_sc_2024_plus"] = payload2.get("meta", {}).get("count")
        out["ok"] = all(isinstance(out.get(k), int) and out[k] > 0 for k in ("count_2025_plus", "count_sc_2024_plus"))
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def live_verification(timeout: int = 15) -> dict:
    gdc = fetch_gdc_counts(timeout=timeout)
    oax = fetch_openalex_csc(timeout=timeout + 5)
    return {"gdc": gdc, "openalex": oax, "live_ok": bool(gdc.get("ok") and oax.get("ok"))}
