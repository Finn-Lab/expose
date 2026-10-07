"""
fetchers.py
===========

Thin clients for the EMBL-EBI / SIB resources used to build biotransformation
evidence:

- ChEBI      (compound identity, structure, formula)          https://www.ebi.ac.uk/chebi/backend/api/public
- Rhea       (curated biochemical reactions, ChEBI <-> EC)     https://www.rhea-db.org/help/rest-api
- UniProt    (enzymes / proteins catalysing those EC numbers)  https://rest.uniprot.org
- ChEMBL     (bioactivity / metabolism evidence)                https://www.ebi.ac.uk/chembl/api/data

Every client supports two interchangeable modes, chosen per-call:

1. Live API      -- `fetch(...)` hits the network (requires network access).
2. Local file     -- `load_file(path)` reads a previously-downloaded JSON/TSV
                     response, so the pipeline runs fully offline and is
                     reproducible / testable without hammering the APIs.

Each client normalises its raw response into a small list of plain dicts
("records") with a common shape so `integrate.py` doesn't need to know the
source-specific formats. Every record carries a `source` and `evidence`
field, which is what the graph and the plots use to justify a prediction.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import requests


class EMBLAPIError(RuntimeError):
    """Raised when a live API call fails after retries."""


#: Set by `cache.install_http_cache()`; None means no caching.
HTTP_CACHE = None

#: Minimum seconds between *network* requests to a host, matched on substring.
#: ChEMBL is rate-sensitive and has had long outages, and a big library build
#: is hundreds of calls, so it is paced by default. Cached responses never
#: reach here -- `install_http_cache` returns before the fetch -- so a re-run
#: over cached assays is not slowed at all.
RATE_LIMITS: dict[str, float] = {"ebi.ac.uk/chembl": 1.0}

_LAST_REQUEST: dict[str, float] = {}
_RATE_LOCK = threading.Lock()


def set_rate_limit(host_fragment: str, min_interval: float | None) -> None:
    """Pace (or stop pacing) requests whose URL contains `host_fragment`."""
    if min_interval:
        RATE_LIMITS[host_fragment] = float(min_interval)
    else:
        RATE_LIMITS.pop(host_fragment, None)


def _throttle(url: str) -> None:
    """Block until this host's minimum interval has passed."""
    for fragment, interval in RATE_LIMITS.items():
        if fragment not in url:
            continue
        while True:
            with _RATE_LOCK:
                now = time.monotonic()
                previous = _LAST_REQUEST.get(fragment, 0.0)
                wait = previous + interval - now
                if wait <= 0:
                    _LAST_REQUEST[fragment] = now
                    return
            time.sleep(min(wait, interval))


def _get_json(url: str, params: dict | None = None, timeout: int = 20,
               retries: int = 3, backoff: float = 1.5) -> Any:
    headers = {"User-Agent": "embl-biotransform/0.1 (research tool)"}
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            _throttle(url)
            resp = requests.get(url, params=params, headers=headers, timeout=timeout)
            if 400 <= resp.status_code < 500 and resp.status_code != 429:
                # client error (unknown accession, bad parameter): retrying won't help
                try:
                    detail = resp.json().get("detail")
                except Exception:  # noqa: BLE001
                    detail = None
                raise EMBLAPIError(f"GET {resp.url} -> {resp.status_code}: {detail or resp.reason}")
            resp.raise_for_status()
            return resp.json()
        except EMBLAPIError:
            raise
        except Exception as exc:  # noqa: BLE001 - we want to retry broadly, then surface
            last_exc = exc
            time.sleep(backoff ** attempt)
    raise EMBLAPIError(f"GET {url} failed after {retries} attempts: {last_exc}")


def _get_text(url: str, params: dict | None = None, timeout: int = 20,
              retries: int = 3, backoff: float = 1.5) -> str:
    headers = {"User-Agent": "embl-biotransform/0.1 (research tool)"}
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=timeout)
            resp.raise_for_status()
            return resp.text
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
            time.sleep(backoff ** attempt)
    raise EMBLAPIError(f"GET {url} failed after {retries} attempts: {last_exc}")


def _parse_tsv(text: str) -> list[dict[str, str]]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


# --------------------------------------------------------------------------- #
# ChEBI
# --------------------------------------------------------------------------- #

class ChEBIClient:
    """Compound identity / structure lookups against ChEBI."""

    BASE_URL = "https://www.ebi.ac.uk/chebi/backend/api/public"

    def fetch_compound(self, chebi_id: str) -> dict:
        """Fetch one compound record by ChEBI id (e.g. 'CHEBI:17234' or '17234')."""
        cid = chebi_id.split(":")[-1]
        data = _get_json(f"{self.BASE_URL}/compounds/{cid}/")
        return self._normalise(data)

    def load_file(self, path: str | Path) -> dict:
        data = json.loads(Path(path).read_text())
        return self._normalise(data)

    @staticmethod
    def _normalise(data: dict) -> dict:
        struct = data.get("default_structure") or {}
        chem = data.get("chemical_data") or {}
        return {
            "source": "ChEBI",
            "chebi_id": f"CHEBI:{data.get('id') or data.get('chebi_accession', '')}".replace("CHEBI:CHEBI", "CHEBI"),
            "name": data.get("name") or data.get("ascii_name"),
            "smiles": struct.get("smiles"),
            "formula": chem.get("formula"),
            "mass": chem.get("mass"),
            "evidence": "ChEBI compound record (curated chemical entity).",
        }


# --------------------------------------------------------------------------- #
# Rhea
# --------------------------------------------------------------------------- #

class RheaClient:
    """Reactions a compound participates in, linked to EC numbers."""

    BASE_URL = "https://www.rhea-db.org/rhea"

    def fetch_reactions_for_chebi(self, chebi_id: str, limit: int = 25) -> list[dict]:
        cid = chebi_id if chebi_id.upper().startswith("CHEBI:") else f"CHEBI:{chebi_id}"
        params = {
            "query": cid,
            "columns": "rhea-id,equation,chebi-id,ec",
            "format": "tsv",
            "limit": limit,
        }
        text = _get_text(self.BASE_URL, params=params)
        return self._normalise(_parse_tsv(text))

    def load_file(self, path: str | Path) -> list[dict]:
        text = Path(path).read_text()
        # accept either a raw TSV dump or a JSON list of the same rows
        try:
            rows = json.loads(text)
        except json.JSONDecodeError:
            rows = _parse_tsv(text)
        return self._normalise(rows)

    @staticmethod
    def _normalise(rows: Iterable[dict]) -> list[dict]:
        out = []
        for row in rows:
            ec_field = row.get("EC") or row.get("ec") or ""
            ec_numbers = [e.strip() for e in ec_field.replace("EC:", "").split(";") if e.strip()]
            out.append({
                "source": "Rhea",
                "rhea_id": row.get("RHEA ID") or row.get("rhea-id") or row.get("Reaction identifier"),
                "equation": row.get("Equation") or row.get("equation"),
                "chebi_participants": row.get("ChEBI identifier") or row.get("chebi-id"),
                "ec_numbers": ec_numbers,
                "evidence": "Curated reaction in Rhea (SIB/EMBL-EBI), links compound to EC number(s).",
            })
        return out


# --------------------------------------------------------------------------- #
# UniProt
# --------------------------------------------------------------------------- #

class UniProtClient:
    """Enzymes (reviewed Swiss-Prot entries) catalysing a given EC number."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

    def fetch_enzymes_for_ec(self, ec_number: str, limit: int = 10,
                              reviewed_only: bool = True) -> list[dict]:
        query = f"ec:{ec_number}"
        if reviewed_only:
            query += " AND reviewed:true"
        params = {
            "query": query,
            "fields": "accession,protein_name,gene_names,organism_name,ec,cc_function",
            "format": "json",
            "size": limit,
        }
        data = _get_json(self.BASE_URL, params=params)
        return self._normalise(data.get("results", []))

    def load_file(self, path: str | Path) -> list[dict]:
        data = json.loads(Path(path).read_text())
        results = data.get("results", data) if isinstance(data, dict) else data
        return self._normalise(results)

    @staticmethod
    def _normalise(results: Iterable[dict]) -> list[dict]:
        out = []
        for r in results:
            protein_desc = (r.get("proteinDescription", {}) or {})
            rec_name = (protein_desc.get("recommendedName", {}) or {}).get("fullName", {}) or {}
            genes = r.get("genes", []) or []
            gene_name = (genes[0].get("geneName", {}) or {}).get("value") if genes else None
            organism = (r.get("organism", {}) or {}).get("scientificName")
            out.append({
                "source": "UniProt",
                "accession": r.get("primaryAccession"),
                "protein_name": rec_name.get("value"),
                "gene_name": gene_name,
                "organism": organism,
                "evidence": "Reviewed (Swiss-Prot) enzyme annotated with the matching EC number.",
            })
        return out


# --------------------------------------------------------------------------- #
# ChEMBL
# --------------------------------------------------------------------------- #

class ChEMBLClient:
    """Bioactivity / metabolism-related evidence for a compound."""

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def fetch_metabolism_activities(self, chembl_molecule_id: str, limit: int = 25) -> list[dict]:
        params = {
            "molecule_chembl_id": chembl_molecule_id,
            "format": "json",
            "limit": limit,
        }
        data = _get_json(f"{self.BASE_URL}/activity.json", params=params)
        return self._normalise(data.get("activities", []))

    def load_file(self, path: str | Path) -> list[dict]:
        data = json.loads(Path(path).read_text())
        activities = data.get("activities", data) if isinstance(data, dict) else data
        return self._normalise(activities)

    @staticmethod
    def _normalise(activities: Iterable[dict]) -> list[dict]:
        out = []
        for a in activities:
            out.append({
                "source": "ChEMBL",
                "assay_description": a.get("assay_description"),
                "target_organism": a.get("target_organism"),
                "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
                "evidence": "Bioactivity assay result that may reflect metabolic turnover or clearance.",
            })
        return out

# --------------------------------------------------------------------------- #
# MGnify (metagenome taxonomy + functional annotation)
# --------------------------------------------------------------------------- #

class MGnifyClient:
    """Metagenome analysis results from the MGnify (EMBL-EBI) API v2
    (https://www.ebi.ac.uk/metagenomics/api/v2/): which microbes (taxonomy)
    and which enzyme/functional families (InterPro, Pfam, KEGG Orthologs,
    GO, Rhea) were detected in a sample analysis, plus study search /
    analysis listing for building a dataset picker.

    An "analysis accession" looks like 'MGYA01022557'; a study accession like
    'MGYS00010462'. By default only analyses from MGnify pipeline v6 are used
    (see `pipeline`): the file layouts and annotation sets below are v6's.

    Where the data come from in v2:

    - taxonomy: the JSON endpoint `analyses/{acc}/annotations/taxonomies__<db>`
      (SSU, LSU, PR2, UNITE, ITSoneDB, DADA2). Pipeline-v6 assemblies have no
      rRNA profile there; for those the contig taxonomy (Krona table) listed
      in the analysis' `downloads` is used instead.
    - functions: the per-analysis summary tables listed in `downloads`
      (InterPro / Pfam / KO / GO / Rhea counts). The
      v2 JSON annotation endpoints give descriptions and counts but no
      accessions, and accessions are what `integrate.py` matches on, so the
      JSON is only a fallback (Pfam), matched by name.

    The normalisers accept v2 items (`{"organism"|"description", "count"}`),
    raw v1 JSON:API items (`{"id", "attributes"}`) and records already in this
    package's normalised shape, so a file saved from any of them loads back.
    """

    BASE_URL = "https://www.ebi.ac.uk/metagenomics/api/v2"
    PAGE_SIZE = 100  # the v2 API's maximum

    # UI value -> (label, v2 annotation type or None for the contig table)
    TAXONOMY_SOURCES = {
        "ssu": ("SSU rRNA", "taxonomies__ssu"),
        "lsu": ("LSU rRNA", "taxonomies__lsu"),
        "contigs": ("contig taxonomy", None),
        "pr2": ("PR2", "taxonomies__pr2"),
        "unite": ("UNITE", "taxonomies__unite"),
        "its_one_db": ("ITSoneDB", "taxonomies__its_one_db"),
        "dada2_silva": ("DADA2 SILVA", "taxonomies__dada2_silva"),
        "dada2_pr2": ("DADA2 PR2", "taxonomies__dada2_pr2"),
    }
    AUTO_TAXONOMY_ORDER = ("ssu", "contigs", "lsu", "dada2_silva", "pr2", "unite", "its_one_db", "dada2_pr2")

    # kind -> (label, download alias suffixes of the pipeline-v6 summary tables)
    FUNCTIONAL_KINDS = {
        "interpro": ("InterPro", ("_interpro_summary.tsv.gz",)),
        "kegg": ("KEGG orthologs", ("_ko_summary.tsv.gz",)),
        "pfam": ("Pfam", ("_pfam_summary.tsv.gz",)),
        "go-slim": ("GO slim", ("_goslim_summary.tsv.gz",)),
        "go": ("GO", ("_go_summary.tsv.gz",)),
        "rhea": ("Rhea reactions", ("_proteins2rhea.tsv.gz",)),
    }
    # v1 endpoint names, still accepted
    _KIND_ALIASES = {"interpro-identifiers": "interpro", "kegg-orthologs": "kegg",
                     "pfam-entries": "pfam", "go-terms": "go"}

    def __init__(self, base_url: str | None = None, max_pages: int = 100, pipeline: str | None = "V6"):
        """`pipeline`: only use analyses from this MGnify pipeline version
        (prefix match, so 'V6' covers V6, V6.1, ...); None allows any."""
        self.base_url = (base_url or self.BASE_URL).rstrip("/")
        self.max_pages = max_pages
        self.pipeline = pipeline
        self._details: dict[str, dict] = {}

    def _pipeline_ok(self, version: str | None) -> bool:
        return not self.pipeline or str(version or "").upper().startswith(self.pipeline.upper())

    # -- paging -------------------------------------------------------------- #

    def _get_all(self, url: str, params: dict | None = None, max_items: int | None = None) -> list[dict]:
        """v2 list endpoints return `{"count": N, "items": [...]}`, paged by
        `page` (1-based) and `page_size`."""
        items: list[dict] = []
        page_size = min(max_items or self.PAGE_SIZE, self.PAGE_SIZE)
        for page in range(1, self.max_pages + 1):
            data = _get_json(url, params={**(params or {}), "page": page, "page_size": page_size})
            batch = data.get("items") or []
            items.extend(batch)
            if max_items and len(items) >= max_items:
                return items[:max_items]
            if not batch or len(items) >= (data.get("count") or 0):
                break
        return items

    # -- discovery ----------------------------------------------------------- #

    def search_studies(self, keyword: str, page_size: int = 25, max_candidates: int = 100) -> list[dict]:
        """Studies matching `keyword` (title or accession). With `pipeline`
        set, only studies whose analyses come from that pipeline are kept.
        This is checked per study, because the API's own
        `has_analyses_from_pipeline` filter misses some v6 studies."""
        items = self._get_all(f"{self.base_url}/studies/",
                              params={"search": keyword, "order": "-accession"},
                              max_items=max_candidates if self.pipeline else page_size)
        if self.pipeline:
            with ThreadPoolExecutor(max_workers=8) as pool:
                keep = list(pool.map(lambda it: self._study_has_pipeline(it.get("accession")), items))
            items = [it for it, k in zip(items, keep) if k][:page_size]
        out = []
        for it in items:
            biome = it.get("biome") or {}
            meta = it.get("metadata") or {}
            out.append({
                "accession": it.get("accession"),
                "name": it.get("title") or meta.get("study_title"),
                "biome": biome.get("lineage") or biome.get("biome_name"),
                "ena_accessions": it.get("ena_accessions") or [],
                "abstract": (meta.get("study_description") or meta.get("study_abstract") or "")[:400],
            })
        return out

    def _study_has_pipeline(self, study_accession: str) -> bool:
        """Judged by the study's first listed analysis (MGnify reanalyses with a
        new pipeline get a new study accession, so studies aren't mixed)."""
        try:
            first = self._get_all(f"{self.base_url}/studies/{study_accession}/analyses/", max_items=1)
        except EMBLAPIError:
            return False
        return bool(first) and self._pipeline_ok(first[0].get("pipeline_version"))

    def list_analyses(self, study_accession: str, max_items: int = 500) -> list[dict]:
        items = self._get_all(f"{self.base_url}/studies/{study_accession}/analyses/",
                              max_items=max_items)
        out = []
        for it in items:
            if not self._pipeline_ok(it.get("pipeline_version")):
                continue
            sample = it.get("sample") or {}
            out.append({
                "accession": it.get("accession"),
                "experiment_type": it.get("experiment_type"),
                "pipeline_version": it.get("pipeline_version"),
                "sample": sample.get("accession"),
                "sample_title": sample.get("sample_title"),
                "run": (it.get("run") or {}).get("accession"),
                "assembly": (it.get("assembly") or {}).get("accession"),
            })
        return out

    def find_analysis_for_sample(self, sample_accession: str) -> dict | None:
        """Resolve an ENA / BioSample accession (ERS..., SAMEA...) to the
        MGnify analysis of its assembly, or None if there isn't one.

        This goes via the sample's assemblies rather than via
        `studies/insdc/{project}`, because a project's *reads* study and the
        *assembly* study MGnify derives from it are different accessions with
        different ENA ids: the reads study is often the one INSDC lookup
        returns, and often the one with no analyses.
        """
        try:
            assemblies = self._get_all(f"{self.base_url}/samples/{sample_accession}/assemblies/",
                                       max_items=25)
        except EMBLAPIError:
            return None
        for asm in assemblies:
            accession = asm.get("accession")
            if not accession:
                continue
            try:
                analyses = self._get_all(f"{self.base_url}/assemblies/{accession}/analyses", max_items=25)
            except EMBLAPIError:
                continue
            for an in analyses:
                if self._pipeline_ok(an.get("pipeline_version")):
                    return {"analysis": an.get("accession"),
                            "study": an.get("study_accession") or asm.get("assembly_study_accession"),
                            "assembly": accession,
                            "sample": asm.get("sample_accession") or sample_accession,
                            "pipeline_version": an.get("pipeline_version")}
        return None

    def fetch_annotation_taxonomy(self, analysis_id: str, progress=None) -> dict:
        """Which organism carries each detected enzyme, for one analysis.

        Streams the analysis' annotation-summary GFF and contig-taxonomy table
        and joins them on the contig. Assembly analyses only -- an amplicon
        analysis has neither file.
        """
        from . import contig_taxonomy as ct

        detail = self.get_analysis(analysis_id)
        gff = ct.find_download(detail, "gff")
        taxonomy = ct.find_download(detail, "taxonomy")
        if not gff or not taxonomy:
            raise EMBLAPIError(
                f"{analysis_id} has no annotation-summary GFF and contig taxonomy "
                f"({detail.get('experiment_type') or 'this'} analysis); enzyme-to-organism "
                "attribution needs an assembly analysis.")

        if progress:
            progress(0.1, f"{analysis_id}: contig taxonomy…")
        taxonomy_lines = ct.stream_gzip_lines(taxonomy["url"])
        if progress:
            progress(0.25, f"{analysis_id}: annotations ({gff.get('file_size_bytes') or '~12'} bytes)…")
        gff_lines = ct.stream_gzip_lines(
            gff["url"],
            progress=(lambda n: progress(min(0.95, 0.25 + n / 800_000), f"{analysis_id}: {n:,} CDS rows…"))
            if progress else None)
        links = ct.build_links(taxonomy_lines, gff_lines)
        links["analysis"] = analysis_id
        links["assembly"] = (detail.get("assembly") or {}).get("accession")
        if progress:
            progress(1.0, f"{analysis_id}: {links['pairs']:,} enzyme-taxon pairs")
        return links

    def get_analysis(self, analysis_id: str) -> dict:
        """Analysis detail, including its list of downloadable files. Raises
        EMBLAPIError for an analysis from a pipeline other than `pipeline`."""
        detail = self._details.get(analysis_id)
        if detail is None:
            detail = _get_json(f"{self.base_url}/analyses/{analysis_id}")
            if len(self._details) > 256:
                self._details.clear()
            self._details[analysis_id] = detail
        if not self._pipeline_ok(detail.get("pipeline_version")):
            raise EMBLAPIError(f"{analysis_id} is from MGnify pipeline {detail.get('pipeline_version')}; "
                               f"only {self.pipeline} analyses are supported.")
        return detail

    # -- detections ---------------------------------------------------------- #

    def fetch_taxonomy(self, analysis_id: str, source: str = "auto") -> list[dict]:
        """`source` is a key of TAXONOMY_SOURCES, or 'auto' for the first one
        with data (SSU, then contig taxonomy, then LSU, ...)."""
        if source != "auto" and source not in self.TAXONOMY_SOURCES:
            raise ValueError(f"Unknown taxonomy source {source!r}; "
                             f"use 'auto' or one of {sorted(self.TAXONOMY_SOURCES)}")
        detail = self.get_analysis(analysis_id)
        for src in (self.AUTO_TAXONOMY_ORDER if source == "auto" else (source,)):
            label, annotation_type = self.TAXONOMY_SOURCES[src]
            if annotation_type:
                rows = self._get_all(f"{self.base_url}/analyses/{analysis_id}/annotations/{annotation_type}")
            else:
                f = _find_download(detail, (".krona.txt.gz",))
                rows = _parse_krona(_get_text_file(f["url"])) if f else []
            evidence = f"Taxon detected in the metagenome ({label}, MGnify {analysis_id}, API v2)."
            out = self._normalise_taxonomy({**r, "evidence": evidence} for r in rows)
            if out:
                return out
        tried = "any taxonomy source" if source == "auto" else self.TAXONOMY_SOURCES[source][0]
        raise EMBLAPIError(f"No {tried} data for {analysis_id}.")

    def fetch_functional_annotations(self, analysis_id: str, kind: str = "interpro") -> list[dict]:
        """`kind` is a key of FUNCTIONAL_KINDS: 'interpro', 'kegg', 'pfam',
        'go-slim', 'go' or 'rhea' (Rhea reactions: pipeline-v6 assemblies)."""
        kind = self._KIND_ALIASES.get(kind, kind)
        if kind not in self.FUNCTIONAL_KINDS:
            raise ValueError(f"Unknown functional annotation kind {kind!r}; "
                             f"use one of {sorted(self.FUNCTIONAL_KINDS)}")
        label, suffixes = self.FUNCTIONAL_KINDS[kind]
        detail = self.get_analysis(analysis_id)
        evidence = (f"{label} entry directly detected in the metagenome's gene annotation "
                    f"(MGnify {analysis_id}, API v2).")
        f = _find_download(detail, suffixes)
        if f:
            text = _get_text_file(f["url"])
            rows = _parse_rhea_summary(text) if kind == "rhea" else _parse_summary_table(text)
        elif kind == "pfam":
            # JSON fallback: counts + descriptions only, so matching is by name
            rows = self._get_all(f"{self.base_url}/analyses/{analysis_id}/annotations/pfams")
        else:
            rows = []
        out = self._normalise_functional({**r, "evidence": evidence} for r in rows)
        if not out:
            what = detail.get("experiment_type") or "this"
            raise EMBLAPIError(f"No {label} annotation for {analysis_id} "
                               f"({what} analysis, pipeline {detail.get('pipeline_version') or '?'}).")
        return out

    def load_taxonomy_file(self, path: str | Path) -> list[dict]:
        data = json.loads(Path(path).read_text())
        return self._normalise_taxonomy(_rows_of(data))

    def load_functional_file(self, path: str | Path) -> list[dict]:
        data = json.loads(Path(path).read_text())
        return self._normalise_functional(_rows_of(data))

    @staticmethod
    def _normalise_taxonomy(rows: Iterable[dict]) -> list[dict]:
        out = []
        for item in rows:
            attrs = item.get("attributes", item) or {}
            organism = (item.get("organism") or attrs.get("name") or attrs.get("lineage")
                        or item.get("name") or item.get("lineage") or item.get("taxon"))
            if not organism:
                continue
            organism = _clean_lineage(organism)
            lineage = attrs.get("lineage") or item.get("lineage")
            if not lineage and ";" in organism:  # v2 gives the full lineage as `organism`
                lineage = organism
            abundance = (item.get("abundance") if "abundance" in item else
                         attrs.get("count") if attrs.get("count") is not None else attrs.get("relative-abundance"))
            record = {
                "source": item.get("source", "MGnify"),
                "taxon_id": item.get("taxon_id") or (item.get("id") if "attributes" in item else None),
                "organism": organism,
                "lineage": lineage,
                "rank": attrs.get("rank") or item.get("rank"),
                "abundance": _num(abundance),
                "evidence": item.get("evidence") or
                            "Taxon detected in the metagenome's rRNA-based taxonomic profile (MGnify).",
            }
            # Carry through anything a profiler or the NCBI resolver added. The
            # NCBI `tax_id` especially: dropping it here would silently disable
            # taxon-id matching and fall back to comparing names across
            # taxonomies, which is exactly what it exists to avoid.
            for extra in _TAXON_PASSTHROUGH:
                if item.get(extra) is not None:
                    record[extra] = item[extra]
            out.append(record)
        return out

    @staticmethod
    def _normalise_functional(rows: Iterable[dict]) -> list[dict]:
        out = []
        for item in rows:
            attrs = item.get("attributes", item) or {}
            annotation_id = (item.get("annotation_id") or attrs.get("accession")
                             or item.get("accession") or item.get("id"))
            description = attrs.get("description") or attrs.get("name") or item.get("description")
            if not (annotation_id or description):
                continue
            abundance = (item.get("abundance") if "abundance" in item else
                         attrs.get("count") if "count" in attrs else item.get("count"))
            record = {
                "source": item.get("source", "MGnify"),
                "annotation_id": annotation_id,
                "description": description,
                "abundance": _num(abundance),
                "evidence": item.get("evidence") or
                            "Functional/enzyme family directly detected in the metagenome's gene annotation (MGnify).",
            }
            # same reason as the taxon passthrough: enzyme-to-organism
            # attribution is attached upstream and must not be rebuilt away
            for extra in _FUNCTION_PASSTHROUGH:
                if item.get(extra) is not None:
                    record[extra] = item[extra]
            out.append(record)
        return out


# -- MGnify v2 helpers --------------------------------------------------------- #

_ANNOTATION_ACCESSION = re.compile(r"(IPR\d{6}|PF\d{5}|K\d{5}|GO:\d{7}|RHEA:\d+)", re.I)


def _rows_of(data: Any) -> list:
    """The record list inside a saved API response: v2 `items`, v1 `data`,
    or a bare list."""
    if isinstance(data, dict):
        return data.get("items") or data.get("data") or []
    return data


def _clean_lineage(organism: str) -> str:
    """Legacy-pipeline taxonomies come back from v2 as
    'Bacteria::Bacteroidetes:Flavobacteriia:...:Feifantangia|5.0';
    turn those into the ';'-separated form used everywhere else."""
    s = str(organism).strip()
    if ";" in s or ":" not in s:
        return s
    s = s.split("|", 1)[0]
    return ";".join(p for p in s.split(":") if p.strip())


def _find_download(detail: dict, suffixes: tuple[str, ...]) -> dict | None:
    downloads = detail.get("downloads") or []
    for suffix in suffixes:  # in order of preference
        for f in downloads:
            if str(f.get("alias", "")).lower().endswith(suffix.lower()) and f.get("url"):
                return f
    return None


def _get_text_file(url: str, timeout: int = 120, retries: int = 3, backoff: float = 1.5) -> str:
    """Download a (possibly gzipped) results file from the MGnify FTP site."""
    headers = {"User-Agent": "embl-biotransform/0.1 (research tool)"}
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            _throttle(url)
            resp = requests.get(url, headers=headers, timeout=timeout)
            resp.raise_for_status()
            data = resp.content
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            return data.decode("utf-8", errors="replace")
        except requests.HTTPError as exc:
            if exc.response is not None and 400 <= exc.response.status_code < 500:
                raise EMBLAPIError(f"GET {url} failed: {exc}") from exc
            last_exc = exc
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
        time.sleep(backoff ** attempt)
    raise EMBLAPIError(f"GET {url} failed after {retries} attempts: {last_exc}")


def _parse_summary_table(text: str) -> list[dict]:
    """InterPro / Pfam / KO / GO count tables, in either layout MGnify uses:
    pipeline v6 (TSV with a header, e.g. `interpro_accession description count`)
    or v4/v5 (header-less quoted CSV, e.g. `"588","IPR036388","Winged helix..."`,
    or `"GO:...","term","category","count"`). The accession, the count and the
    description are recognised by content, so column order doesn't matter."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return []
    delimiter = "\t" if "\t" in lines[0] else ","
    rows = []
    for cells in csv.reader(lines, delimiter=delimiter):
        cells = [c.strip() for c in cells]
        acc = next((c for c in cells if _ANNOTATION_ACCESSION.fullmatch(c)), None)
        if not acc:
            continue  # header or malformed
        rest = [c for c in cells if c != acc]
        count = next((c for c in reversed(rest) if _num(c) is not None), None)
        description = next((c for c in rest if c and _num(c) is None), None)
        rows.append({"annotation_id": acc.upper(), "description": description, "abundance": _num(count)})
    return rows


def _parse_rhea_summary(text: str) -> list[dict]:
    """`proteins2rhea` has one row per (protein, reaction); count distinct
    proteins per Rhea reaction."""
    proteins: dict[str, set] = {}
    reaction: dict[str, str] = {}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        rid = (row.get("rhea_id") or "").strip()
        if not rid:
            continue
        proteins.setdefault(rid, set()).add(row.get("protein_id"))
        reaction.setdefault(rid, (row.get("reaction") or "").strip())
    return [{"annotation_id": rid.upper(), "description": reaction[rid], "abundance": len(p)}
            for rid, p in sorted(proteins.items(), key=lambda kv: -len(kv[1]))]


def _parse_krona(text: str) -> list[dict]:
    """Krona text: `count<TAB>d__Bacteria<TAB>k__...<TAB>p__...` per line."""
    rows = []
    for ln in text.splitlines():
        parts = [p.strip() for p in ln.split("\t") if p.strip()]
        if len(parts) < 2 or _num(parts[0]) is None or parts[1].lower() == "unclassified":
            continue
        rows.append({"organism": ";".join(parts[1:]), "count": _num(parts[0])})
    return rows


#: Fields a caller may have attached that must survive normalisation.
_TAXON_PASSTHROUGH = ("tax_id", "gtdb_organism", "gtdb_lineage", "taxonomy_framework",
                      "matched_rank", "genome")

#: Enzyme-to-organism attribution, attached by `pipeline.attach_enzyme_taxonomy`.
_FUNCTION_PASSTHROUGH = ("taxa", "attributed", "unattributed", "attribution_total")


def _num(v):
    try:
        f = float(v)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return None
