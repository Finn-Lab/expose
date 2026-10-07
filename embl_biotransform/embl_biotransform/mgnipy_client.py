"""
mgnipy_client.py
=================

Adapter around `mgnipy` (https://github.com/EBI-Metagenomics/mgnipy), the
official Python client for the **MGnify API v2**
(https://www.ebi.ac.uk/metagenomics/api/v2/), so the rest of the pipeline
(`integrate.py`, `predict.py`) can consume it through the same normalised
record shape as the plain-`requests` `MGnifyClient` in `fetchers.py`.

mgnipy is a young, actively-evolving package (several releases in its first
few months) -- install with `pip install mgnipy`. This adapter targets the
documented Quick Start pattern (search a resource -> `.get()`/`.get_all()`
-> `.search_results`/`.metadata`, and file-backed results such as taxonomy
and functional-annotation tables via `.datasets` / `.stream(alias=...)`).
If a method has moved by the time you run this, check
https://mgnipy.mgnify.org/ for the current call signatures -- the shape of
what comes back into `integrate.py` (a list of dicts with `organism`/
`description`, `source`, `evidence`) is what actually matters, so you can
freely adjust the mgnipy calls in this one file without touching anything
downstream.

Only imported when you actually use it, so `mgnipy` is an optional
dependency -- the rest of the package works without it (e.g. against files
downloaded via the browsable API, or via the plain-requests `MGnifyClient`).
"""

from __future__ import annotations

from typing import Any


class MGnipySource:
    """Pick studies / analyses interactively (e.g. from a notebook) and
    pull their taxonomy + functional annotation, normalised the same way
    as `fetchers.MGnifyClient`."""

    def __init__(self):
        try:
            from mgnipy import MGnipy
        except ImportError as exc:  # pragma: no cover - environment-dependent
            raise ImportError(
                "mgnipy is not installed. Run `pip install mgnipy` to use "
                "MGnipySource, or use fetchers.MGnifyClient (plain requests "
                "against the MGnify v2 REST API) instead."
            ) from exc
        self._mg = MGnipy()

    # -- discovery, for building a notebook picker -------------------------- #

    def search_studies(self, keyword: str, page_size: int = 20) -> list[dict]:
        query = self._mg.studies(search=keyword)
        with self._mg:
            query.get()
        return [dict(r) for r in query.search_results[:page_size]]

    def list_analyses(self, study_accession: str, page_size: int = 50) -> list[dict]:
        query = self._mg.analyses(study_accession=study_accession)
        with self._mg:
            query.get_all()
        return [dict(r) for r in query.search_results[:page_size]]

    # -- pulling the two things the predictor needs -------------------------- #

    def fetch_taxonomy(self, analysis_accession: str) -> list[dict]:
        query = self._mg.analyses(accession=analysis_accession)
        with self._mg:
            query.get()
            query.enrich_details()
        mgazine = query.datasets
        rows = self._stream_first_matching(mgazine, ["taxonomy", "taxa", "ssu"])
        return [{
            "source": "MGnify",
            "organism": row.get("name") or row.get("lineage") or row.get("taxon"),
            "abundance": row.get("count") or row.get("abundance"),
            "evidence": "Taxon detected in the metagenome's taxonomic profile (MGnify API v2 via mgnipy).",
        } for row in rows]

    def fetch_functional_annotations(self, analysis_accession: str) -> list[dict]:
        query = self._mg.analyses(accession=analysis_accession)
        with self._mg:
            query.get()
            query.enrich_details()
        mgazine = query.datasets
        rows = self._stream_first_matching(mgazine, ["interpro", "ko", "go", "function", "pfam"])
        return [{
            "source": "MGnify",
            "annotation_id": row.get("id") or row.get("accession"),
            "description": row.get("description") or row.get("name"),
            "abundance": row.get("count"),
            "evidence": "Functional/enzyme family directly detected in the metagenome's gene annotation "
                        "(MGnify API v2 via mgnipy).",
        } for row in rows]

    @staticmethod
    def _stream_first_matching(mgazine: Any, keywords: list[str]) -> list[dict]:
        """mgnipy exposes an analysis's downloadable files as a `MGazine`;
        pick the first file whose alias looks like the kind of table we
        want and stream it as records. Adjust `keywords` if MGnify renames
        its output files."""
        aliases = getattr(mgazine, "aliases", None) or []
        for alias in aliases:
            if any(k in str(alias).lower() for k in keywords):
                df = mgazine.stream(alias=alias, df_engine="pandas")
                return df.to_dict(orient="records")
        return []
