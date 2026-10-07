"""
pipeline.py
===========

One-call entry points that tie the pieces together, used by the web app
(and handy from a notebook or script):

- `parse_detections_text` -- read a taxonomy or functional-annotation table
  from JSON (raw MGnify API response or this package's normalised records),
  CSV or TSV.
- `run_analysis`          -- detections + reference associations ->
  JSON-serialisable result: ranked predictions with evidence, the evidence
  graph, and per-detection "was this linked to anything?" flags.
"""

from __future__ import annotations

import csv
import io
import json

from .fetchers import MGnifyClient
from .integrate import EvidenceGraph, rank_is_at_least, taxon_rank
from .predict import BiotransformationPredictor
from .reference import normalise_reference_records


class DetectionFormatError(ValueError):
    """Raised when a detections file can't be parsed."""


_TAXON_COLUMNS = ("organism", "name", "taxon", "lineage", "species", "#SampleID", "taxonomy")
_ABUNDANCE_COLUMNS = ("abundance", "count", "reads", "relative-abundance", "relative_abundance")
_FUNC_ID_COLUMNS = ("annotation_id", "accession", "id", "interpro", "ko", "ec")
_FUNC_DESC_COLUMNS = ("description", "name", "function")


def _pick(row: dict, candidates: tuple[str, ...]):
    lower = {k.lower().strip(): v for k, v in row.items() if k}
    for c in candidates:
        v = lower.get(c.lower())
        if v not in (None, ""):
            return v
    return None


def parse_detections_text(text: str, fmt: str, kind: str) -> list[dict]:
    """`kind` is 'taxa' or 'functions'; `fmt` is 'json', 'csv' or 'tsv'."""
    if kind not in ("taxa", "functions"):
        raise DetectionFormatError(f"kind must be 'taxa' or 'functions', got {kind!r}")
    fmt = fmt.lower().lstrip(".")
    if fmt == "json":
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise DetectionFormatError(f"Invalid JSON: {exc}") from exc
        rows = data
        if isinstance(data, dict):  # MGnify v2 {"count", "items"}, v1 {"data"}, or {"taxa"/"functions"}
            rows = next((data[k] for k in (kind, "items", "data") if isinstance(data.get(k), list)), data)
        if not isinstance(rows, list):
            raise DetectionFormatError("Expected a JSON list of records, or an object with an 'items' or 'data' list.")
    elif fmt in ("csv", "tsv", "txt"):
        delimiter = "," if fmt == "csv" else "\t"
        raw_rows = list(csv.DictReader(io.StringIO(text), delimiter=delimiter))
        rows = []
        for r in raw_rows:
            if kind == "taxa":
                rows.append({"organism": _pick(r, _TAXON_COLUMNS),
                             "abundance": _pick(r, _ABUNDANCE_COLUMNS), "source": "user file"})
            else:
                rows.append({"annotation_id": _pick(r, _FUNC_ID_COLUMNS),
                             "description": _pick(r, _FUNC_DESC_COLUMNS),
                             "abundance": _pick(r, _ABUNDANCE_COLUMNS), "source": "user file"})
    else:
        raise DetectionFormatError(f"Unsupported file type: .{fmt} (use .json, .csv or .tsv)")

    if kind == "taxa":
        out = MGnifyClient._normalise_taxonomy(rows)
    else:
        out = MGnifyClient._normalise_functional(rows)
    if rows and not out:
        raise DetectionFormatError(
            "No usable records found. Taxonomy needs an organism/name/lineage column; "
            "functions need an id/accession or description column.")
    return out


#: Coarser than genus, a detection says an organism *group* is present, which
#: is not evidence that the organism a reference names is there.
DEFAULT_TAXON_RANK = "genus"

#: Rows a microbe/enzyme grid will draw. Generous because those grids show
#: only rows that carry evidence -- with a large reference set every one of
#: them is worth seeing. Anything beyond is reported as `truncated`.
GRID_ROW_LIMIT = 500

#: What a detection row carries out of the graph. Kept in one place because
#: every time this was an inline literal, a field added upstream (the NCBI
#: tax_id, then the enzyme attribution) was silently dropped here.
_DETECTION_FIELDS = (
    "organism", "annotation_id", "description", "abundance", "lineage", "rank",
    "tax_id", "gtdb_organism", "matched_rank", "taxonomy_framework",
    "taxa", "attributed", "unattributed", "attribution_total",
)


def filter_taxa_by_rank(taxa: list[dict], floor: str = DEFAULT_TAXON_RANK) -> tuple[list[dict], dict]:
    """Keep only detections at `floor` rank or finer.

    A contig or read profile reports every level of its lineage, so a sample
    yields rows like `d__Bacteria` and `p__Pseudomonadota` alongside real
    species. Those cannot support a Tier-2 claim about a named organism, and
    they inflate abundance shares, so they are dropped before the graph is
    built. `floor` of "any" keeps everything.
    """
    if floor in (None, "any"):
        return list(taxa), {"floor": "any", "kept": len(taxa), "dropped": 0, "by_rank": {}}
    kept, dropped_by_rank = [], {}
    for record in taxa:
        rank = taxon_rank(record)
        if rank_is_at_least(rank, floor):
            kept.append(record)
        else:
            dropped_by_rank[rank or "unknown"] = dropped_by_rank.get(rank or "unknown", 0) + 1
    return kept, {"floor": floor, "kept": len(kept),
                  "dropped": sum(dropped_by_rank.values()), "by_rank": dropped_by_rank}


def attach_enzyme_taxonomy(functions: list[dict], enzyme_taxonomy: dict | None,
                           require_attribution: bool = False) -> tuple[list[dict], dict]:
    """Attach "which organism carries this enzyme" to each detection.

    `enzyme_taxonomy` is `{accession: attribution}` from
    `contig_taxonomy.attribution()`. Without it nothing changes, so the
    summary-table path keeps working.

    With `require_attribution`, detections that no contig placed at
    genus/species rank are dropped. That is the strict reading: an enzyme
    family found only on contigs classified no further than "Bacteria" shows
    the community has the enzyme, but not that any *named* organism carries
    it -- and treating the two alike is what makes predictions too liberal.
    """
    if not enzyme_taxonomy:
        return functions, {"available": False, "annotated": 0, "dropped": 0}

    kept, annotated, dropped = [], 0, 0
    for record in functions:
        accession = str(record.get("annotation_id") or "").upper()
        hit = enzyme_taxonomy.get(accession)
        if hit:
            annotated += 1
            record["taxa"] = hit.get("taxa") or []
            record["attributed"] = hit.get("attributed") or 0
            record["unattributed"] = hit.get("unattributed") or 0
            record["attribution_total"] = hit.get("total") or 0
            if record["taxa"]:
                names = ", ".join(t["scientific_name"] for t in record["taxa"][:3]
                                  if t.get("scientific_name"))
                record["evidence"] = (
                    f"{record.get('evidence', '')} Carried by contigs assigned to "
                    f"{names} ({record['attributed']} of {record['attribution_total']} "
                    "occurrences attributable).").strip()
        if require_attribution and not (record.get("taxa") or []):
            dropped += 1
            continue
        kept.append(record)
    return kept, {"available": True, "annotated": annotated, "dropped": dropped,
                  "required": require_attribution}


def run_analysis(taxa: list[dict], functions: list[dict], references: list[dict],
                 sample_id: str = "sample", allow_name_matching: bool = True,
                 sample_meta: dict | None = None,
                 taxon_rank_floor: str = DEFAULT_TAXON_RANK,
                 enzyme_taxonomy: dict | None = None,
                 require_enzyme_attribution: bool = False) -> dict:
    references = normalise_reference_records(references)
    # accept raw or normalised detections
    taxa = MGnifyClient._normalise_taxonomy(taxa or [])
    functions = MGnifyClient._normalise_functional(functions or [])
    taxa, rank_filter = filter_taxa_by_rank(taxa, taxon_rank_floor)
    functions, attribution_filter = attach_enzyme_taxonomy(
        functions, enzyme_taxonomy, require_enzyme_attribution)

    graph = EvidenceGraph()
    sample = graph.add_sample(sample_id, **(sample_meta or {}))
    graph.add_detected_taxa(sample, taxa)
    graph.add_detected_functions(sample, functions)
    graph.add_reference_associations(references)
    n_links = graph.link_detections_to_references(allow_name_matching=allow_name_matching)
    predictions = BiotransformationPredictor(graph).predict(sample)

    g = graph.graph
    linked = {u for u, v, d in g.edges(data=True) if d.get("source") == "integration"}
    matched_refs = {v for u, v, d in g.edges(data=True) if d.get("source") == "integration"}

    def _detections(node_type: str) -> list[dict]:
        rows = []
        for n in g.successors(sample):
            d = g.nodes[n]
            if d.get("type") != node_type:
                continue
            refs = sorted({v for _, v, e in g.out_edges(n, data=True) if e.get("source") == "integration"})
            rows.append({"node": n, **{k: d.get(k) for k in _DETECTION_FIELDS},
                         "linked_references": refs})
        rows.sort(key=lambda r: (not r["linked_references"], -(r["abundance"] or 0)))
        return rows

    references_out = []
    for n, d in g.nodes(data=True):
        if d.get("type") == "reference":
            references_out.append({"node": n, "matched": n in matched_refs,
                                   **{k: v for k, v in d.items() if k != "type"}})

    return {
        "sample_id": sample_id,
        "sample_node": sample,
        "summary": {
            "taxa_detected": len([1 for n in g.successors(sample) if g.nodes[n].get("type") == "taxon"]),
            "functions_detected": len([1 for n in g.successors(sample)
                                       if g.nodes[n].get("type") == "functional_annotation"]),
            "detections_linked": len(linked),
            "references_total": len(references_out),
            "references_matched": len(matched_refs),
            "links": n_links,
            "predictions": len(predictions),
            "tier1_predictions": sum(1 for p in predictions if p.tier1_hits),
        },
        "predictions": [p.to_dict() for p in predictions],
        "taxa": _detections("taxon"),
        "functions": _detections("functional_annotation"),
        "references": references_out,
        "graph": _trim_graph(graph.to_dict(), linked),
        "settings": {"allow_name_matching": allow_name_matching,
                     "taxon_rank_floor": rank_filter["floor"],
                     "require_enzyme_attribution": require_enzyme_attribution},
        "rank_filter": rank_filter,
        "attribution_filter": attribution_filter,
    }


def _trim_graph(graph_dict: dict, linked_detections: set) -> dict:
    """Keep the graph small enough to draw: drop detections that aren't
    linked to any reference (a real MGnify sample has thousands)."""
    keep = {n["id"] for n in graph_dict["nodes"]
            if n.get("type") not in ("taxon", "functional_annotation") or n["id"] in linked_detections}
    return {
        "nodes": [n for n in graph_dict["nodes"] if n["id"] in keep],
        "edges": [e for e in graph_dict["edges"]
                  if e["source_node"] in keep and e["target_node"] in keep],
    }


def run_multi_analysis(samples: list[dict], references: list[dict],
                       allow_name_matching: bool = True,
                       observations: list[dict] | None = None,
                       taxon_rank_floor: str = DEFAULT_TAXON_RANK,
                       enzyme_taxonomy: dict | None = None,
                       require_enzyme_attribution: bool = False,
                       linked_only: bool = True) -> dict:
    """Run `run_analysis` over several samples and build a comparison grid.

    `samples` is a list of `{"sample_id", "taxa", "functions", "sample_meta"}`.
    `observations` are ChEMBL community records (see `chembl_references`); where
    one names an MGnify analysis that is also among these samples, its measured
    outcomes are attached to that column, so a cell can say both what was
    predicted and what was actually observed.
    """
    if not samples:
        raise DetectionFormatError("No samples supplied.")
    references = normalise_reference_records(references)

    results = []
    for sample in samples:
        result = run_analysis(
            sample.get("taxa") or [], sample.get("functions") or [], list(references),
            sample_id=str(sample.get("sample_id") or f"sample-{len(results) + 1}"),
            allow_name_matching=allow_name_matching,
            sample_meta=sample.get("sample_meta") or {},
            taxon_rank_floor=taxon_rank_floor,
            enzyme_taxonomy=(enzyme_taxonomy or {}).get(sample.get("sample_id"))
                            or sample.get("enzyme_taxonomy"),
            require_enzyme_attribution=require_enzyme_attribution)
        result["sample_meta"] = sample.get("sample_meta") or {}
        results.append(result)

    return {
        "samples": results,
        "grid": build_grid(results, observations),
        "taxon_grid": build_taxon_grid(results, linked_only=linked_only),
        "enzyme_grid": build_enzyme_grid(results, linked_only=linked_only),
        "settings": {"allow_name_matching": allow_name_matching,
                     "taxon_rank_floor": taxon_rank_floor},
        "summary": {
            "samples": len(results),
            "taxa_dropped_by_rank": sum(r["rank_filter"]["dropped"] for r in results),
            "enzymes_dropped_unattributed": sum(r["attribution_filter"]["dropped"] for r in results),
            "enzymes_attributed": sum(r["attribution_filter"]["annotated"] for r in results),
            "reaction_classes": len({p["reaction_class"] for r in results for p in r["predictions"]}),
            "predictions": sum(len(r["predictions"]) for r in results),
            "references_total": results[0]["summary"]["references_total"] if results else 0,
        },
    }


def _observed_index(observations: list[dict] | None) -> dict:
    """{analysis accession: {reaction_class: observed bool}} for the community
    assays that name an MGnify analysis."""
    index: dict[str, dict[str, bool]] = {}
    for community in observations or []:
        analysis = ((community.get("mgnify") or {}).get("analysis"))
        if not analysis:
            continue
        per_class = index.setdefault(analysis, {})
        for obs in community.get("observations") or []:
            if obs.get("reaction_class") is not None:
                per_class[obs["reaction_class"]] = bool(obs.get("observed"))
    return index


def build_grid(results: list[dict], observations: list[dict] | None = None) -> dict:
    """A reaction-class x sample matrix: one row per reaction class predicted in
    any sample, one cell per sample (None where that sample predicts nothing).

    Rows are ordered by how widely the class is predicted, then by its best
    score, so the classes shared across samples come first -- which is what a
    reader comparing communities is looking for.
    """
    observed_by_analysis = _observed_index(observations)
    sample_ids = [r["sample_id"] for r in results]
    by_class: dict[str, dict] = {}

    for result in results:
        for pred in result["predictions"]:
            row = by_class.setdefault(pred["reaction_class"], {
                "reaction_class": pred["reaction_class"],
                "substrate_chebi": pred.get("substrate_chebi"),
                "cells": {}})
            row["substrate_chebi"] = row["substrate_chebi"] or pred.get("substrate_chebi")
            row["cells"][result["sample_id"]] = {
                "score": pred["score"],
                "tier": 1 if pred["tier1_hits"] else 2,
                "confidence_label": pred.get("confidence_label"),
                "tier1_hits": len(pred["tier1_hits"]),
                "tier2_hits": len(pred["tier2_hits"]),
            }

    # attach measured outcomes, including for classes nothing predicted
    for analysis, per_class in observed_by_analysis.items():
        if analysis not in sample_ids:
            continue
        for reaction_class, observed in per_class.items():
            row = by_class.setdefault(reaction_class, {
                "reaction_class": reaction_class, "substrate_chebi": None, "cells": {}})
            cell = row["cells"].get(analysis)
            if cell is None:
                row["cells"][analysis] = {"score": None, "tier": None, "observed": observed}
            else:
                cell["observed"] = observed

    rows = []
    for row in by_class.values():
        cells = [dict(row["cells"].get(sid) or {}, sample_id=sid) if row["cells"].get(sid) is not None
                 else None for sid in sample_ids]
        scores = [c["score"] for c in cells if c and c.get("score") is not None]
        rows.append({**row, "cells": cells,
                     "samples_predicted": len(scores),
                     "max_score": max(scores) if scores else 0.0})
    rows.sort(key=lambda r: (-r["samples_predicted"], -r["max_score"], r["reaction_class"]))

    return {
        "sample_ids": sample_ids,
        "samples": [{"sample_id": r["sample_id"], "meta": r.get("sample_meta") or {},
                     "predictions": len(r["predictions"]),
                     "has_observations": r["sample_id"] in observed_by_analysis} for r in results],
        "rows": rows,
        "max_score": max((r["max_score"] for r in rows), default=0.0),
        "observed_samples": [s for s in sample_ids if s in observed_by_analysis],
    }


def build_taxon_grid(results: list[dict], limit: int = GRID_ROW_LIMIT, linked_only: bool = True) -> dict:
    """A microbe x sample matrix, the taxonomic counterpart to `build_grid`.

    The reaction grid answers "which communities could do this?"; this one
    answers "which organisms are they, and which of them carry the evidence?"
    -- so the two can be read against each other.

    A gut profile has hundreds to thousands of taxa, far more than is readable,
    so rows are ranked with the ones that matter first: taxa linked to a
    reference association, then those present in the most samples, then the
    most abundant. `limit` caps the rows and `truncated` reports what was cut.

    Abundances are not comparable between samples as they stand -- bioSIFTR
    gives fractions, MGnify gives read counts, and sequencing depth differs --
    so each cell also carries `share`, the taxon's fraction of that sample's
    total, which is what the heatmap should colour by.
    """
    sample_ids = [r["sample_id"] for r in results]
    totals = {r["sample_id"]: sum((t.get("abundance") or 0) for t in r["taxa"]) or 0
              for r in results}

    by_taxon: dict[str, dict] = {}
    for result in results:
        sample_id = result["sample_id"]
        for taxon in result["taxa"]:
            name = taxon.get("organism")
            if not name:
                continue
            # a taxon id is a stabler key than a name across samples, and two
            # profilers can spell the same organism differently
            key = f"taxid:{taxon['tax_id']}" if taxon.get("tax_id") else f"name:{name}"
            row = by_taxon.setdefault(key, {
                "taxon": name,
                "tax_id": taxon.get("tax_id"),
                "rank": taxon.get("matched_rank") or taxon.get("rank"),
                "lineage": taxon.get("lineage"),
                "gtdb_organism": taxon.get("gtdb_organism"),
                "linked_references": [],
                "cells": {}})
            abundance = taxon.get("abundance")
            total = totals.get(sample_id) or 0
            row["cells"][sample_id] = {
                "abundance": abundance,
                "share": (abundance / total) if (abundance and total) else None,
                "linked": bool(taxon.get("linked_references")),
            }
            for ref in taxon.get("linked_references") or []:
                if ref not in row["linked_references"]:
                    row["linked_references"].append(ref)

    rows = []
    for row in by_taxon.values():
        cells = [dict(row["cells"][s], sample_id=s) if s in row["cells"] else None for s in sample_ids]
        shares = [c["share"] for c in cells if c and c["share"] is not None]
        rows.append({**row, "cells": cells,
                     "samples_present": sum(1 for c in cells if c),
                     "max_share": max(shares) if shares else 0.0,
                     "linked": bool(row["linked_references"])})
    rows.sort(key=lambda r: (not r["linked"], -r["samples_present"], -r["max_share"], r["taxon"]))

    linked_count = sum(1 for r in rows if r["linked"])
    shown = [r for r in rows if r["linked"]] if linked_only else rows
    kept, truncated = shown[:limit], max(0, len(shown) - limit)
    return {
        "sample_ids": sample_ids,
        "rows": kept,
        "total_taxa": len(rows),
        "truncated": truncated,
        "linked_taxa": linked_count,
        "linked_only": linked_only,
        "hidden_unlinked": len(rows) - linked_count if linked_only else 0,
        # scale the ramp over what is on screen, not over taxa that were hidden
        "max_share": max((r["max_share"] for r in kept), default=0.0),
    }


def build_enzyme_grid(results: list[dict], limit: int = GRID_ROW_LIMIT, linked_only: bool = True) -> dict:
    """An enzyme x sample matrix: the functional counterpart to the other two.

    Rows are the functional annotations detected (InterPro / Pfam / KO / EC /
    Rhea ids). Those linked to a gene reference come first, because they are
    the Tier-1 evidence behind the reaction grid -- this view is where you see
    *which* enzyme produced a prediction, and whether it was found by id or
    only by name.
    """
    sample_ids = [r["sample_id"] for r in results]
    totals = {r["sample_id"]: sum((f.get("abundance") or 0) for f in r["functions"]) or 0
              for r in results}

    by_enzyme: dict[str, dict] = {}
    for result in results:
        sample_id = result["sample_id"]
        # how each reference was matched, so the grid can flag name-only hits
        match_types: dict[str, set] = {}
        for prediction in result["predictions"]:
            for hit in prediction["tier1_hits"]:
                for detection in hit["detections"]:
                    key = detection.get("annotation_id") or detection.get("label")
                    if key:
                        match_types.setdefault(key, set()).add(detection.get("match_type"))

        for function in result["functions"]:
            key = function.get("annotation_id") or function.get("description")
            if not key:
                continue
            row = by_enzyme.setdefault(key, {
                "annotation_id": function.get("annotation_id"),
                "description": function.get("description"),
                "linked_references": [],
                "match_types": [],
                "taxa": [],
                "attributed": 0,
                "unattributed": 0,
                "cells": {}})
            row["description"] = row["description"] or function.get("description")
            row["attributed"] += function.get("attributed") or 0
            row["unattributed"] += function.get("unattributed") or 0
            for taxon in function.get("taxa") or []:
                existing = next((t for t in row["taxa"] if t["tax_id"] == taxon["tax_id"]), None)
                if existing:
                    existing["count"] += taxon.get("count") or 0
                else:
                    row["taxa"].append(dict(taxon))
            abundance = function.get("abundance")
            total = totals.get(sample_id) or 0
            row["cells"][sample_id] = {
                "abundance": abundance,
                "share": (abundance / total) if (abundance and total) else None,
                "linked": bool(function.get("linked_references")),
            }
            for ref in function.get("linked_references") or []:
                if ref not in row["linked_references"]:
                    row["linked_references"].append(ref)
            for kind in match_types.get(key, ()):
                if kind and kind not in row["match_types"]:
                    row["match_types"].append(kind)

    rows = []
    for row in by_enzyme.values():
        cells = [dict(row["cells"][s], sample_id=s) if s in row["cells"] else None for s in sample_ids]
        shares = [c["share"] for c in cells if c and c["share"] is not None]
        row["taxa"].sort(key=lambda t: -(t.get("count") or 0))
        rows.append({**row, "cells": cells,
                     "samples_present": sum(1 for c in cells if c),
                     "max_share": max(shares) if shares else 0.0,
                     "linked": bool(row["linked_references"]),
                     "name_only": bool(row["match_types"]) and set(row["match_types"]) == {"name"},
                     "unattributable": bool(row["attributed"] + row["unattributed"]) and not row["taxa"]})
    rows.sort(key=lambda r: (not r["linked"], -r["samples_present"], -r["max_share"],
                             str(r["annotation_id"] or "")))

    linked_count = sum(1 for r in rows if r["linked"])
    shown = [r for r in rows if r["linked"]] if linked_only else rows
    kept, truncated = shown[:limit], max(0, len(shown) - limit)
    return {
        "sample_ids": sample_ids,
        "rows": kept,
        "total_enzymes": len(rows),
        "truncated": truncated,
        "linked_enzymes": linked_count,
        "linked_only": linked_only,
        "hidden_unlinked": len(rows) - linked_count if linked_only else 0,
        "name_only": sum(1 for r in kept if r["name_only"]),
        "with_taxonomy": sum(1 for r in kept if r["taxa"]),
        "unattributable": sum(1 for r in kept if r["unattributable"]),
        "max_share": max((r["max_share"] for r in kept), default=0.0),
    }
