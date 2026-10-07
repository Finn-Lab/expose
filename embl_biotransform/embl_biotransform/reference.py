"""
reference.py
============

Loads the curated "ground truth" biotransformation associations that anchor
everything else: statements of the form "this microbe" or "this gene/enzyme"
is known to carry out "this biotransformation". No single EMBL-EBI API
packages this directly, so it is supplied as a file (JSON or CSV) -- built,
for example, by combining literature curation with UniProt (gene -> EC ->
organism) and Rhea (EC -> reaction) lookups.

Each association is one of two kinds, matching the two ways biotransformation
evidence was described:

- entity_type == "gene"    : tied to a specific enzyme/gene. Matched against
                              a metagenome's *functional* annotation
                              (EC / KO / InterPro id). This is the strongest
                              (Tier 1) evidence when it's directly detected.
- entity_type == "microbe" : tied to an organism/taxon. Matched against a
                              metagenome's *taxonomic* profile. This is
                              weaker (Tier 2) evidence -- the organism is
                              present, but its specific enzyme wasn't
                              independently confirmed in this sample.

Expected JSON shape (a list of records):
[
  {
    "entity_type": "gene",
    "entity_id": "EC:1.14.14.1",        # or a KO/InterPro id
    "alt_ids": ["IPR000000", "K00000"],  # optional: other ids for the same enzyme,
                                         # matched exactly against detected annotations
    "entity_name": "CYP153 alkane hydroxylase",
    "reaction_class": "Alkane hydroxylation",
    "chebi_substrate": "CHEBI:138366",   # optional, e.g. dodecane
    "source": "Literature / UniProt",
    "evidence": "Characterised alkane omega-hydroxylase (PMID:...)."
  },
  {
    "entity_type": "microbe",
    "entity_id": "Pseudomonas putida",
    "entity_name": "Pseudomonas putida",
    "reaction_class": "Alkane hydroxylation",
    "chebi_substrate": "CHEBI:138366",
    "source": "Literature",
    "evidence": "Reported to degrade medium-chain alkanes (PMID:...)."
  }
]

CSV files use the same column names as the JSON keys above; `alt_ids` is a
semicolon-separated string in CSV.
"""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path


REQUIRED_FIELDS = {"entity_type", "entity_id", "reaction_class"}


class ReferenceFormatError(ValueError):
    """Raised when a reference file or record is malformed."""


def normalise_reference_records(records: list[dict]) -> list[dict]:
    """Validate and fill defaults for a list of reference records (in place
    and returned). Shared by the file loader and the web API."""
    if not isinstance(records, list):
        raise ReferenceFormatError("Reference data must be a list of records.")
    out = []
    for i, raw in enumerate(records, start=1):
        if not isinstance(raw, dict):
            raise ReferenceFormatError(f"Record {i} is not an object: {raw!r}")
        rec = {k.strip(): (v.strip() if isinstance(v, str) else v)
               for k, v in raw.items() if k is not None}
        # CSV rows give "" for empty cells -- treat as missing
        rec = {k: v for k, v in rec.items() if v not in ("", None)}
        missing = REQUIRED_FIELDS - rec.keys()
        if missing:
            raise ReferenceFormatError(
                f"Record {i} is missing required field(s) {sorted(missing)}: {raw}")
        rec["entity_type"] = str(rec["entity_type"]).lower()
        if rec["entity_type"] not in ("gene", "microbe"):
            raise ReferenceFormatError(
                f"Record {i}: entity_type must be 'gene' or 'microbe', got {rec['entity_type']!r}")
        alt = rec.get("alt_ids", [])
        if isinstance(alt, str):
            alt = [a.strip() for a in alt.replace(",", ";").split(";") if a.strip()]
        rec["alt_ids"] = list(alt)
        rec.setdefault("source", "user-supplied reference file")
        rec.setdefault("evidence", f"Curated association: {rec['entity_id']} -> {rec['reaction_class']}.")
        rec.setdefault("entity_name", rec["entity_id"])
        out.append(rec)
    return out


def parse_reference_text(text: str, fmt: str) -> list[dict]:
    """Parse reference data from a string. `fmt` is 'json' or 'csv'/'tsv'."""
    fmt = fmt.lower().lstrip(".")
    if fmt == "json":
        try:
            records = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ReferenceFormatError(f"Invalid JSON: {exc}") from exc
    elif fmt in ("csv", "tsv"):
        delimiter = "\t" if fmt == "tsv" else ","
        records = list(csv.DictReader(io.StringIO(text), delimiter=delimiter))
    else:
        raise ReferenceFormatError(f"Unsupported reference file type: .{fmt} (use .json, .csv or .tsv)")
    return normalise_reference_records(records)


def load_reference_associations(path: str | Path) -> list[dict]:
    path = Path(path)
    return parse_reference_text(path.read_text(), path.suffix)
