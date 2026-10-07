"""
EXPOSE biotransformation viewer -- Flask backend.

Run from the project root (the folder containing `embl_biotransform/`):

    python -m webapp            # http://127.0.0.1:8050
    python -m webapp --port 9000 --debug

The front end (webapp/static/) is a single page with no external
dependencies; everything it needs is served from here. Only the MGnify
endpoints below reach out to the internet (MGnify API v2, plus the
per-analysis summary files it links to on the MGnify FTP site).
"""

from __future__ import annotations

import json
import re
import sys
import threading
import time
from pathlib import Path

from flask import Flask, jsonify, request, send_from_directory

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:  # allow `python webapp/app.py` as well as `python -m webapp`
    sys.path.insert(0, str(ROOT))

from embl_biotransform import (  # noqa: E402
    BioSIFTRFormatError, ChEMBLBiotransformationSource, DetectionFormatError, MGnifyClient,
    ReferenceFormatError, load_biosiftr_run, parse_biosiftr_table, parse_detections_text,
    parse_reference_text, run_analysis, run_multi_analysis, summarise_biosiftr_run,
)
from embl_biotransform.pipeline import DEFAULT_TAXON_RANK  # noqa: E402
from embl_biotransform.cache import Cache, SEARCH_TTL  # noqa: E402
from embl_biotransform.fetchers import EMBLAPIError  # noqa: E402
from embl_biotransform.jobs import JobQueue  # noqa: E402
from embl_biotransform.taxonomy import TaxonomyResolver  # noqa: E402

STATIC_DIR = Path(__file__).resolve().parent / "static"
EXAMPLE_DIR = ROOT / "example_data"
#: No upload cap. This is a local desktop application reading files off the
#: same machine, so the usual reason for a limit -- an untrusted client
#: exhausting a shared server -- does not apply, and a fully built ChEMBL
#: reference set is far larger than any round number worth guessing at.
MAX_UPLOAD_MB = None
#: No cap on a ChEMBL build: the point is to assemble a large reference set
#: over time. Requests are paced (see `fetchers.RATE_LIMITS`) and every assay
#: is cached, so a long run is slow but resumable -- re-running it costs
#: nothing for the assays already held.
#: Analyses per batch download. Each is a real MGnify fetch, so this is a
#: guard against a mistyped request rather than a policy limit; raise it
#: freely.
MAX_ANALYSES = 500


class _MemoryCache:
    """Fallback when the app is created without a persistent cache (tests, or
    `expose --no-cache`). Same two-argument shape as `Cache.get_or_set`."""

    def __init__(self, ttl_seconds: int = 3600):
        self.ttl = ttl_seconds
        self.path = None
        self._data: dict = {}
        self._lock = threading.Lock()

    def get_or_set(self, namespace, key, fn, ttl=None):
        now = time.time()
        composite = (namespace, key)
        with self._lock:
            hit = self._data.get(composite)
            if hit and now - hit[0] < (ttl or self.ttl):
                return hit[1]
        value = fn()
        with self._lock:
            self._data[composite] = (now, value)
        return value

    def set(self, namespace, key, value, ttl=None):
        with self._lock:
            self._data[(namespace, key)] = (time.time(), value)

    def get(self, namespace, key):
        hit = self._data.get((namespace, key))
        return hit[1] if hit else None

    def keys(self, namespace):
        return [k for (ns, k) in self._data if ns == namespace]

    def items(self, namespace, limit=None):
        rows = [(k, v[1], v[0]) for (ns, k), v in self._data.items() if ns == namespace]
        rows.sort(key=lambda r: -r[2])
        return rows[:limit] if limit else rows

    def delete(self, namespace, key):
        with self._lock:
            self._data.pop((namespace, key), None)

    def stats(self):
        return {"path": None, "entries": len(self._data), "file_bytes": 0, "namespaces": {}}


class _NullContext:
    """A JobContext stand-in for the synchronous paths."""

    def progress(self, *_args, **_kwargs): pass
    def step(self, *_args, **_kwargs): pass
    def check_cancelled(self): pass
    cancelled = False


#: Cache namespaces that are keyed by something readable, and so listable.
LIBRARY_NS = "chembl.library"
MGNIFY_LIBRARY_NS = "mgnify.analysis"
ENZYME_TAXONOMY_NS = "mgnify.enzyme-taxonomy"

#: ChEMBL accessions, however they arrive (one per line, a CSV column, prose).
_CHEMBL_ID_RE = re.compile(r"\bCHEMBL\d+\b", re.I)
SETTINGS_NS = "chembl.settings"


def _library_labels(built: dict) -> dict:
    """A human label for a library entry, taken from its first record."""
    first = (built.get("records") or [None])[0] or {}
    community = (built.get("observations") or [None])[0] or {}
    return {"organism": first.get("organism") or community.get("organism"),
            "tier": ("gene" if first.get("entity_type") == "gene"
                     else "microbiome" if community else "microbe")}


def _remember_selection(store, assay_ids, replace: bool = False) -> None:
    current = set() if replace else set(store.get(SETTINGS_NS, "selection") or [])
    store.set(SETTINGS_NS, "selection", sorted(current | set(assay_ids)))


def _merge_counts(dicts) -> dict:
    merged: dict = {}
    for d in dicts:
        for key, value in (d or {}).items():
            merged[key] = merged.get(key, 0) + value
    return merged


def _summarise_records(out: dict) -> dict:
    """The same shape `ChEMBLExtraction.summary()` returns, recomputed after
    merging several per-assay builds."""
    records, observations = out["records"], out["observations"]
    kinds = [r["entity_type"] for r in records]
    return {
        "records": len(records),
        "genes": kinds.count("gene"),
        "microbes": kinds.count("microbe"),
        "positive": sum(1 for r in records if r.get("observed")),
        "negative": sum(1 for r in records if r.get("observed") is False),
        "reaction_classes": len({r["reaction_class"] for r in records}),
        "communities": len(observations),
        "community_observations": sum(len(o.get("observations") or []) for o in observations),
        "warnings": len(out["warnings"]),
    }


def _file_format(filename: str) -> str:
    return Path(filename or "").suffix.lower().lstrip(".") or "json"


def create_app(mgnify_client: MGnifyClient | None = None,
               chembl_source: ChEMBLBiotransformationSource | None = None,
               cache: Cache | None = None, max_workers: int = 4) -> Flask:
    app = Flask(__name__, static_folder=None)
    if MAX_UPLOAD_MB:
        app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_MB * 1024 * 1024
    mgnify = mgnify_client or MGnifyClient()
    chembl = chembl_source or ChEMBLBiotransformationSource(mgnify_client=mgnify)
    store = cache if cache is not None else _MemoryCache()
    jobs = JobQueue(max_workers=max_workers)
    app.config["EXPOSE_CACHE"] = cache
    app.config["EXPOSE_JOBS"] = jobs
    app.config["EXPOSE_STORE"] = store

    def error(message: str, status: int = 400):
        return jsonify({"error": message}), status

    # -- front end ------------------------------------------------------------ #

    @app.get("/")
    def index():
        return send_from_directory(STATIC_DIR, "index.html")

    @app.get("/static/<path:filename>")
    def static_files(filename):
        return send_from_directory(STATIC_DIR, filename)

    # -- API ------------------------------------------------------------------- #

    @app.get("/api/health")
    def health():
        import embl_biotransform
        return jsonify({"status": "ok", "version": embl_biotransform.__version__})

    @app.get("/api/examples")
    def examples():
        try:
            return jsonify({
                "sample_id": "example-sample",
                "taxa": json.loads((EXAMPLE_DIR / "sample_taxonomy.json").read_text()),
                "functions": json.loads((EXAMPLE_DIR / "sample_functions.json").read_text()),
                "references": parse_reference_text(
                    (EXAMPLE_DIR / "reference_associations.json").read_text(), "json"),
            })
        except FileNotFoundError as exc:
            return error(f"Example data not found: {exc}", 404)

    @app.post("/api/parse/reference")
    def parse_reference():
        f = request.files.get("file")
        if f is None:
            return error("No file uploaded (form field 'file').")
        try:
            text = f.read().decode("utf-8-sig")
            records = parse_reference_text(text, _file_format(f.filename))
        except (ReferenceFormatError, UnicodeDecodeError) as exc:
            return error(str(exc))
        return jsonify({"filename": f.filename, "records": records})

    @app.post("/api/parse/detections")
    def parse_detections():
        f = request.files.get("file")
        kind = request.form.get("kind", "")
        if f is None:
            return error("No file uploaded (form field 'file').")
        try:
            text = f.read().decode("utf-8-sig")
            records = parse_detections_text(text, _file_format(f.filename), kind)
        except (DetectionFormatError, UnicodeDecodeError) as exc:
            return error(str(exc))
        return jsonify({"filename": f.filename, "kind": kind, "records": records})

    @app.post("/api/analyze")
    def analyze():
        body = request.get_json(silent=True)
        if not isinstance(body, dict):
            return error("Expected a JSON body.")
        taxa = body.get("taxa") or []
        functions = body.get("functions") or []
        references = body.get("references") or []
        if not references:
            return error("No reference associations supplied.")
        if not taxa and not functions:
            return error("No detections supplied: give a taxonomy and/or a functional annotation table.")
        try:
            result = run_analysis(
                taxa, functions, references,
                sample_id=str(body.get("sample_id") or "sample"),
                allow_name_matching=bool(body.get("allow_name_matching", True)),
                sample_meta=body.get("sample_meta") or {},
            )
        except (ReferenceFormatError, DetectionFormatError) as exc:
            return error(str(exc))
        return jsonify(result)

    @app.post("/api/analyze/multi")
    def analyze_multi():
        """Several samples at once, plus the reaction-class x sample grid."""
        body = request.get_json(silent=True)
        if not isinstance(body, dict):
            return error("Expected a JSON body.")
        samples = body.get("samples") or []
        references = body.get("references") or []
        if not isinstance(samples, list) or not samples:
            return error("No samples supplied (expected a 'samples' list).")
        if not references:
            return error("No reference associations supplied.")
        if not any((s.get("taxa") or s.get("functions")) for s in samples):
            return error("No detections supplied in any sample.")
        try:
            result = run_multi_analysis(
                samples, references,
                allow_name_matching=bool(body.get("allow_name_matching", True)),
                observations=body.get("observations") or [],
                taxon_rank_floor=body.get("taxon_rank_floor") or DEFAULT_TAXON_RANK,
                enzyme_taxonomy=body.get("enzyme_taxonomy") or None,
                require_enzyme_attribution=bool(body.get("require_enzyme_attribution")),
                linked_only=bool(body.get("linked_only", True)))
        except (ReferenceFormatError, DetectionFormatError) as exc:
            return error(str(exc))
        return jsonify(result)

    # -- background jobs -------------------------------------------------------- #

    def submit(name: str, fn, meta: dict | None = None):
        """202 + a job id. The front end polls /api/jobs/<id>."""
        job = jobs.submit(name, fn, meta=meta or {})
        return jsonify(job.to_dict()), 202

    @app.get("/api/jobs/<job_id>")
    def job_status(job_id):
        job = jobs.get(job_id)
        if job is None:
            return error("No such job (it may have expired).", 404)
        return jsonify(job.to_dict())

    @app.get("/api/jobs")
    def job_list():
        return jsonify({"jobs": jobs.list()})

    @app.post("/api/jobs/<job_id>/cancel")
    def job_cancel(job_id):
        if not jobs.cancel(job_id):
            return error("Job is not running.", 409)
        return jsonify({"cancelled": job_id})

    @app.get("/api/cache")
    def cache_stats():
        return jsonify(store.stats())

    @app.delete("/api/cache")
    def cache_clear():
        """Drop cached responses. The front end confirms first: this discards
        hours of downloads, and references that cannot be rebuilt while ChEMBL
        is unavailable."""
        namespace = request.args.get("namespace")
        if cache is None:
            removed = len(store._data) if namespace is None else 0
            store._data.clear() if namespace is None else None
            return jsonify({"removed": removed, "persistent": False})
        return jsonify({"removed": cache.clear(namespace), "persistent": True})

    @app.post("/api/mgnify/enzyme-taxonomy")
    def mgnify_enzyme_taxonomy():
        """Which organism carries each enzyme, per analysis.

        Joins the annotation-summary GFF to the contig taxonomy. ~13 MB per
        analysis to download and a few seconds to stream, so it runs as a job
        and the derived table is cached -- once per analysis, ever.
        """
        body = request.get_json(silent=True) or {}
        accessions = body.get("accessions") or []
        floor = tuple(body.get("floor_ranks") or ("species", "genus"))
        if not isinstance(accessions, list) or not accessions:
            return error("No accessions supplied (expected an 'accessions' list).")
        if len(accessions) > MAX_ANALYSES:
            return error(f"Too many analyses at once (limit {MAX_ANALYSES}).")

        def link(ctx):
            from embl_biotransform import contig_taxonomy as ct
            resolver = TaxonomyResolver()
            out, failed = {}, []
            for i, accession in enumerate(accessions):
                base = i / len(accessions)
                try:
                    links = store.get_or_set(
                        ENZYME_TAXONOMY_NS, accession,
                        lambda a=accession: mgnify.fetch_annotation_taxonomy(
                            a, progress=lambda f, m: ctx.progress(
                                base + f / len(accessions), m)))
                except (EMBLAPIError, ValueError) as exc:
                    failed.append(f"{accession}: {exc}")
                    continue
                ctx.progress(base + 0.9 / len(accessions), f"{accession}: naming taxa…")
                table = links.get("annotations") or {}
                info = store.get_or_set(
                    "mgnify.enzyme-taxa", accession,
                    lambda: resolver.resolve_taxids(ct.observed_taxids(links, table)))
                out[accession] = {
                    "attribution": {a: ct.attribution(links, a, info, floor) for a in table},
                    "stats": {k: links[k] for k in
                              ("contigs", "cds_with_annotations", "cds_on_unclassified_contigs", "pairs")},
                }
            if not out:
                raise EMBLAPIError(" ".join(failed) or "nothing linked")
            return {"analyses": out, "failed": failed}

        return submit(f"Enzyme taxonomy: {len(accessions)} analysis(es)", link,
                      {"kind": "enzyme-taxonomy", "accessions": accessions})

    @app.get("/api/mgnify/library")
    def mgnify_library():
        """Analyses already downloaded, from the same on-disk cache.

        Each is keyed by accession plus the taxonomy/function choice it was
        fetched with, so re-running a comparison costs nothing and works with
        MGnify unreachable.
        """
        entries = []
        for key, value, created in store.items(MGNIFY_LIBRARY_NS):
            # entries written before analyses were keyed by accession are
            # hashes: still valid cache, but nothing readable to show
            if key.count("|") != 2:
                continue
            accession, _, rest = key.partition("|")
            taxonomy, _, functional = rest.partition("|")
            entries.append({
                "accession": accession, "taxonomy": taxonomy, "functional": functional,
                "taxa": len(value.get("taxa") or []),
                "functions": len(value.get("functions") or []),
                "warnings": value.get("warnings") or [],
                "cached_at": created,
            })
        return jsonify({"entries": entries, "persistent": cache is not None})

    @app.delete("/api/mgnify/library/<accession>")
    def mgnify_library_delete(accession):
        removed = [k for k in store.keys(MGNIFY_LIBRARY_NS) if k.split("|")[0] == accession]
        for key in removed:
            store.delete(MGNIFY_LIBRARY_NS, key)
        if not removed:
            return error("Not in the library.", 404)
        return jsonify({"removed": accession, "entries": len(removed)})

    @app.post("/api/chembl/library/match")
    def chembl_library_match():
        """Cross-match a list of ChEMBL assay ids against the library.

        Accepts a JSON `assay_ids` list or an uploaded file (one id per line,
        or any text with CHEMBL ids in it -- a CSV column works). Reports what
        is already held and what would need fetching, so only the missing set
        is downloaded.
        """
        supplied: list[str] = []
        uploaded = request.files.get("file")
        if uploaded is not None:
            try:
                text = uploaded.read().decode("utf-8-sig")
            except UnicodeDecodeError:
                return error("Could not read that file as text.")
            supplied = _CHEMBL_ID_RE.findall(text)
        else:
            body = request.get_json(silent=True) or {}
            raw = body.get("assay_ids")
            if isinstance(raw, str):
                supplied = _CHEMBL_ID_RE.findall(raw)
            elif isinstance(raw, list):
                supplied = [str(a).strip().upper() for a in raw if str(a).strip()]
        supplied = list(dict.fromkeys(a.upper() for a in supplied))
        if not supplied:
            return error("No ChEMBL ids found (expected accessions like CHEMBL5303650).")

        known = set(store.keys(LIBRARY_NS))
        present = [a for a in supplied if a in known]
        missing = [a for a in supplied if a not in known]
        return jsonify({"supplied": supplied, "present": present, "missing": missing,
                        "library_size": len(known)})

    # -- bioSIFTR (shallow-shotgun profiles from reads) ------------------------- #

    @app.post("/api/biosiftr/scan")
    def biosiftr_scan():
        """Read a bioSIFTR output directory on this machine.

        A local path is the natural input for the desktop app: bioSIFTR runs on
        a cluster and writes gigabytes, so the app reads its results in place
        rather than having them uploaded. `/api/parse/biosiftr` takes uploaded
        tables instead, for when the server isn't the machine holding them.
        """
        body = request.get_json(silent=True) or {}
        path = (body.get("path") or "").strip()
        mapper = body.get("mapper") or "sm"
        resolve = bool(body.get("resolve_ncbi", True))
        if not path:
            return error("Give the path to a bioSIFTR output directory.")
        if mapper not in ("sm", "bwa"):
            return error("mapper must be 'sm' (sourmash) or 'bwa'.")
        # fail fast on a bad path, before a job is created
        try:
            run = load_biosiftr_run(path, prefer=mapper)
        except BioSIFTRFormatError as exc:
            return error(str(exc))
        except OSError as exc:
            return error(f"Could not read {path}: {exc}")

        if not resolve:
            return jsonify({"summary": summarise_biosiftr_run(run), "samples": run["samples"],
                            "taxonomy": None})

        def resolve_job(ctx):
            """bioSIFTR lineages are GTDB; the references are NCBI. Resolve the
            profile into NCBI terms so the two can be matched on taxon id."""
            resolver = TaxonomyResolver()
            stats, total = {}, len(run["samples"])
            for i, (sample_id, entry) in enumerate(sorted(run["samples"].items())):
                ctx.step(i, total, f"{sample_id}: resolving taxonomy against NCBI…")
                stats[sample_id] = resolver.annotate(entry["taxa"])
            summary = summarise_biosiftr_run(run)
            summary["taxonomy"] = {
                "framework": "NCBI (resolved from GTDB)",
                "resolved": sum(s["counts"]["resolved"] for s in stats.values()),
                "unresolved": sum(s["counts"]["unresolved"] for s in stats.values()),
                "by_rank": _merge_counts(s["by_rank"] for s in stats.values()),
                "unresolved_examples": sorted({e for s in stats.values()
                                               for e in s["unresolved_examples"] if e})[:8],
            }
            return {"summary": summary, "samples": run["samples"]}

        return submit(f"bioSIFTR: {len(run['samples'])} sample(s)", resolve_job,
                      {"kind": "biosiftr", "path": path, "mapper": mapper})

    @app.post("/api/parse/biosiftr")
    def parse_biosiftr():
        """Upload bioSIFTR tables directly (taxonomy and/or function TSVs)."""
        files = request.files.getlist("file")
        if not files:
            return error("No files uploaded (form field 'file').")
        samples: dict = {}
        warnings: list[str] = []
        for f in files:
            try:
                parsed = parse_biosiftr_table(f.read().decode("utf-8-sig"),
                                              source=f"bioSIFTR ({f.filename})")
            except (BioSIFTRFormatError, UnicodeDecodeError) as exc:
                warnings.append(f"{f.filename}: {exc}")
                continue
            bucket = "taxa" if parsed["kind"] == "taxa" else "functions"
            for sample_id, records in parsed["samples"].items():
                samples.setdefault(sample_id, {"taxa": [], "functions": []})[bucket].extend(records)
        if not samples:
            return error("; ".join(warnings) or "No usable bioSIFTR tables.")
        summary = {"outdir": None, "mapper": None, "sources": {}, "warnings": warnings,
                   "samples": [{"sample_id": k, "taxa": len(v["taxa"]), "functions": len(v["functions"])}
                               for k, v in sorted(samples.items())]}
        return jsonify({"summary": summary, "samples": samples})

    # -- ChEMBL reference associations ------------------------------------------ #

    @app.get("/api/chembl/assays")
    def chembl_assays():
        """Searching ChEMBL can take two filter queries with retries, and the
        API is frequently down, so this runs as a job rather than blocking."""
        q = (request.args.get("search") or "").strip()
        if not q:
            return error("Give a search keyword (?search=...).")

        def search(ctx):
            ctx.progress(-1.0, f"searching ChEMBL for {q!r}…")
            found = store.get_or_set("chembl.search", Cache.make_key(q),
                                     lambda: {"assays": chembl.search_assays(q)}, SEARCH_TTL)
            # flag the hits already built, so a big result set shows at a
            # glance what a build would actually have to fetch
            held = set(store.keys(LIBRARY_NS))
            assays = [dict(a, in_library=a.get("assay_chembl_id") in held)
                      for a in found.get("assays") or []]
            return {"assays": assays, "in_library": sum(a["in_library"] for a in assays)}

        return submit(f"ChEMBL search: {q}", search, {"kind": "chembl-search", "query": q})

    @app.post("/api/chembl/references")
    def chembl_references():
        body = request.get_json(silent=True) or {}
        assay_ids = body.get("assay_ids") or []
        if not isinstance(assay_ids, list) or not assay_ids:
            return error("No assay ids supplied (expected an 'assay_ids' list).")
        assay_ids = list(dict.fromkeys(str(a).strip().upper() for a in assay_ids if str(a).strip()))

        def build(ctx):
            # one assay at a time, keyed by its accession so the library can be
            # listed back; a partial run still leaves every assay it finished
            merged = {"records": [], "observations": [], "warnings": []}
            cached = fetched = 0
            for i, assay_id in enumerate(assay_ids):
                entry = store.get(LIBRARY_NS, assay_id)
                held = entry is not None
                cached += held
                fetched += not held
                ctx.step(i, len(assay_ids),
                         f"{'cached' if held else 'fetching'} {assay_id} "
                         f"({i + 1} of {len(assay_ids)}; {cached} cached, {fetched} fetched)")
                if entry is None:
                    built = chembl.build([assay_id]).to_dict()
                    entry = {"assay_chembl_id": assay_id, "built_at": time.time(),
                             **built, "summary": built.get("summary") or {}}
                    entry.update(_library_labels(built))
                    store.set(LIBRARY_NS, assay_id, entry)
                for field in merged:
                    merged[field].extend(entry.get(field) or [])
            if not merged["records"] and not merged["observations"]:
                raise EMBLAPIError("; ".join(merged["warnings"]) or "No usable assays.")
            merged["summary"] = _summarise_records(merged)
            _remember_selection(store, assay_ids)
            ctx.progress(1.0, "done")
            return merged

        return submit(f"ChEMBL references: {len(assay_ids)} assay(s)", build,
                      {"kind": "chembl-refs", "assay_ids": assay_ids})

    # -- the ChEMBL library: what has been built, kept between runs ------------- #

    @app.get("/api/chembl/library")
    def chembl_library():
        """Assays built previously, with their selection state.

        These survive restarts, so a set of references assembled once can be
        toggled on and off across sessions without rebuilding or re-fetching.
        """
        selected = set(store.get(SETTINGS_NS, "selection") or [])
        entries = []
        for assay_id, entry, created in store.items(LIBRARY_NS):
            summary = entry.get("summary") or {}
            entries.append({
                "assay_chembl_id": assay_id,
                "organism": entry.get("organism"),
                "tier": entry.get("tier"),
                "reaction_classes": summary.get("reaction_classes", 0),
                "positive": summary.get("positive", 0),
                "negative": summary.get("negative", 0),
                "communities": summary.get("communities", 0),
                "warnings": entry.get("warnings") or [],
                "built_at": entry.get("built_at") or created,
                "selected": assay_id in selected,
            })
        return jsonify({"entries": entries, "selected": sorted(selected),
                        "persistent": cache is not None})

    @app.post("/api/chembl/library/selection")
    def chembl_library_selection():
        body = request.get_json(silent=True) or {}
        assay_ids = body.get("assay_ids")
        if not isinstance(assay_ids, list):
            return error("Expected an 'assay_ids' list.")
        known = set(store.keys(LIBRARY_NS))
        unknown = [a for a in assay_ids if a not in known]
        if unknown:
            return error(f"Not in the library: {', '.join(unknown)}")
        _remember_selection(store, assay_ids, replace=True)
        return jsonify({"selected": sorted(set(assay_ids))})

    @app.get("/api/chembl/library/references")
    def chembl_library_references():
        """The merged references for the current selection (or `?ids=`)."""
        requested = request.args.get("ids")
        assay_ids = ([a.strip() for a in requested.split(",") if a.strip()] if requested
                     else list(store.get(SETTINGS_NS, "selection") or []))
        merged = {"records": [], "observations": [], "warnings": [], "assay_ids": []}
        for assay_id in assay_ids:
            entry = store.get(LIBRARY_NS, assay_id)
            if entry is None:
                merged["warnings"].append(f"{assay_id}: not in the library.")
                continue
            merged["assay_ids"].append(assay_id)
            for field in ("records", "observations", "warnings"):
                merged[field].extend(entry.get(field) or [])
        merged["summary"] = _summarise_records(merged)
        return jsonify(merged)

    @app.delete("/api/chembl/library/<assay_id>")
    def chembl_library_delete(assay_id):
        if store.get(LIBRARY_NS, assay_id) is None:
            return error("Not in the library.", 404)
        store.delete(LIBRARY_NS, assay_id)
        remaining = [a for a in (store.get(SETTINGS_NS, "selection") or []) if a != assay_id]
        _remember_selection(store, remaining, replace=True)
        return jsonify({"removed": assay_id})

    # -- MGnify (live) ---------------------------------------------------------- #

    def _mgnify_call(namespace, key, fn, ttl=None):
        # readable keys (accession-based) are kept as-is so the library can list
        # them back; free-text ones are hashed
        cache_key = key if isinstance(key, str) else Cache.make_key(key)
        try:
            return jsonify(store.get_or_set(namespace, cache_key, fn, ttl))
        except EMBLAPIError as exc:
            return error(f"MGnify request failed: {exc}", 502)

    @app.get("/api/mgnify/studies")
    def mgnify_studies():
        q = (request.args.get("search") or "").strip()
        if not q:
            return error("Give a search keyword (?search=...).")
        return _mgnify_call("mgnify.studies", q,
                            lambda: {"studies": mgnify.search_studies(q)}, SEARCH_TTL)

    @app.get("/api/mgnify/studies/<accession>/analyses")
    def mgnify_analyses(accession):
        listed = store.get_or_set("mgnify.analyses", Cache.make_key(accession),
                                  lambda: {"analyses": mgnify.list_analyses(accession)}, SEARCH_TTL)
        # which of these are already on disk, so a long study list shows what
        # is free to load and what needs downloading
        downloaded = {k.split("|")[0] for k in store.keys(MGNIFY_LIBRARY_NS) if "|" in k}
        linked = set(store.keys(ENZYME_TAXONOMY_NS))
        analyses = [dict(a, downloaded=a.get("accession") in downloaded,
                         enzyme_taxonomy=a.get("accession") in linked)
                    for a in listed.get("analyses") or []]
        return jsonify({"analyses": analyses,
                        "downloaded": sum(a["downloaded"] for a in analyses)})

    def _load_analysis(ctx, accession: str, tax: str, kind: str) -> dict:
        """One analysis' detections. Taxonomy and functions are fetched (and
        cached) separately, so a sample whose functional annotation is missing
        still contributes its taxonomy."""
        out = {"accession": accession, "taxa": [], "functions": [], "warnings": []}
        ctx.progress(0.1, f"{accession}: taxonomy…")
        try:
            out["taxa"] = store.get_or_set(
                "mgnify.taxonomy", f"{accession}|{tax}",
                lambda: mgnify.fetch_taxonomy(accession, source=tax))
        except (EMBLAPIError, ValueError) as exc:
            out["warnings"].append(f"Taxonomy unavailable: {exc}")
        ctx.progress(0.55, f"{accession}: functional annotation…")
        try:
            out["functions"] = store.get_or_set(
                "mgnify.functions", f"{accession}|{kind}",
                lambda: mgnify.fetch_functional_annotations(accession, kind=kind))
        except (EMBLAPIError, ValueError) as exc:
            out["warnings"].append(f"Functional annotation unavailable: {exc}")
        if not out["taxa"] and not out["functions"]:
            raise EMBLAPIError(" ".join(out["warnings"]) or "no data returned")
        return out

    @app.get("/api/mgnify/analyses/<accession>")
    def mgnify_analysis_data(accession):
        """Kept synchronous for a single analysis (and for scripts); the
        viewer uses the batch job below when comparing several."""
        tax = request.args.get("taxonomy", "auto")
        kind = request.args.get("functional", "interpro")
        return _mgnify_call(MGNIFY_LIBRARY_NS, f"{accession}|{tax}|{kind}",
                            lambda: _load_analysis(_NullContext(), accession, tax, kind))

    @app.post("/api/mgnify/analyses")
    def mgnify_analyses_batch():
        """Download several analyses as one job: this is the slow step when
        comparing communities (5-15 s each), and it must not block the UI."""
        body = request.get_json(silent=True) or {}
        accessions = body.get("accessions") or []
        tax = body.get("taxonomy", "auto")
        kind = body.get("functional", "interpro")
        if not isinstance(accessions, list) or not accessions:
            return error("No accessions supplied (expected an 'accessions' list).")
        if len(accessions) > MAX_ANALYSES:
            return error(f"Too many analyses at once (limit {MAX_ANALYSES}).")

        def load(ctx):
            loaded, failed = [], []
            for i, accession in enumerate(accessions):
                ctx.step(i, len(accessions), f"{accession} ({i + 1} of {len(accessions)})")
                try:
                    loaded.append(store.get_or_set(
                        MGNIFY_LIBRARY_NS, f"{accession}|{tax}|{kind}",
                        lambda a=accession: _load_analysis(ctx, a, tax, kind)))
                except (EMBLAPIError, ValueError) as exc:
                    failed.append(f"{accession}: {exc}")
            if not loaded:
                raise EMBLAPIError(" ".join(failed) or "nothing loaded")
            return {"analyses": loaded, "failed": failed}

        return submit(f"MGnify: {len(accessions)} analysis(es)", load,
                      {"kind": "mgnify-load", "accessions": accessions})

    @app.errorhandler(413)
    def too_large(_):
        # only reachable if a cap is reinstated (e.g. when serving --host 0.0.0.0)
        return error(f"File too large (limit {MAX_UPLOAD_MB} MB).", 413)

    return app


def main(argv: list[str] | None = None) -> None:
    import argparse
    import webbrowser

    parser = argparse.ArgumentParser(description="EXPOSE biotransformation viewer")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8050)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args(argv)

    url = f"http://{args.host}:{args.port}"
    print(f"EXPOSE viewer running at {url}  (Ctrl+C to stop)")
    if not args.no_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    create_app().run(host=args.host, port=args.port, debug=args.debug, use_reloader=False)


if __name__ == "__main__":
    main()
