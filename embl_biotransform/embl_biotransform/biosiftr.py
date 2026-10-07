"""
biosiftr.py
===========

Read the output of **bioSIFTR** (https://github.com/EBI-Metagenomics/biosiftr)
and turn it into the detections this pipeline predicts from.

bioSIFTR profiles *shallow shotgun reads* (< 10 M) against MGnify's
biome-specific genome catalogues, which is a better fit for this pipeline than
what MGnify's own V6 assembly analyses can offer: those carry no rRNA profile
at all, so the taxonomy falls back to contig assignments and sees only what
assembled. Reads-based profiling keeps the community, and the functional
tables come from the same catalogue genomes, so taxonomy and function are
finally two views of one profile.

Running bioSIFTR needs Nextflow, containers and a sizeable reference database,
so this module deliberately *consumes* its output rather than running it. Run
it wherever suits (cluster, workstation), then point EXPOSE at the results.

What it reads
-------------

Every table bioSIFTR writes has the same shape -- a feature id in the first
column, one column per sample -- so one reader handles them all, whether it is
a single sample or the integrated matrix across many:

| file | first column | values |
|------|--------------|--------|
| `taxonomy_tables/<s>_sm_species.tsv` (or `_bwa_`) | `lineage` | relative abundance |
| `function_tables/<s>_*_community_kegg.tsv` | `ko_id` | counts |
| `function_tables/<s>_*_community_pfams.tsv` | `pfam_id` | counts |
| `integrated_annotation/*_matrix.tsv` | `feature_id` | counts |

The lineage carries the catalogue's representative **genome accession** as a
final element (`...;s__Escherichia_coli;MGYG000000001`), which is not a rank.
It is split off into `genome` so the deepest *taxonomic* name is what the
matcher sees -- otherwise every organism would read as "MGYG000000001".
"""

from __future__ import annotations

import csv
import io
import re
from pathlib import Path

#: The first-column header each kind of table uses.
TAXONOMY_HEADERS = ("lineage",)
KEGG_HEADERS = ("ko_id",)
PFAM_HEADERS = ("pfam_id",)
GENERIC_HEADERS = ("feature_id",)

#: MGnify genome accessions, appended to the lineage by bioSIFTR.
_GENOME_RE = re.compile(r"^MGYG\d+$", re.I)
_KO_RE = re.compile(r"^K\d{5}$")
_PFAM_RE = re.compile(r"^PF\d{5}$", re.I)

#: Where the pipeline publishes each kind of table, relative to its outdir.
TAXONOMY_DIR = "taxonomy_tables"
FUNCTION_DIR = "function_tables"
INTEGRATED_DIR = "integrated_annotation"


class BioSIFTRFormatError(ValueError):
    """Raised when a file doesn't look like a bioSIFTR table."""


def split_lineage(lineage: str) -> tuple[str, str | None]:
    """`('d__Bacteria;...;s__Escherichia_coli', 'MGYG000000001')`.

    bioSIFTR builds the lineage as the catalogue lineage plus the
    representative genome, and replaces spaces with underscores; the
    underscores are left alone here because `integrate.clean_taxon_name`
    already turns them back into spaces.
    """
    parts = [p.strip() for p in str(lineage or "").split(";") if p.strip()]
    if parts and _GENOME_RE.match(parts[-1]):
        return ";".join(parts[:-1]), parts[-1]
    return ";".join(parts), None


def detect_kind(header: str, sample_ids: list[str]) -> str:
    """'taxa', 'kegg' or 'pfam' for a table's first-column header, falling back
    to the shape of the feature ids for the integrated matrices, whose header
    is just `feature_id`."""
    name = (header or "").strip().lower()
    if name in TAXONOMY_HEADERS:
        return "taxa"
    if name in KEGG_HEADERS:
        return "kegg"
    if name in PFAM_HEADERS:
        return "pfam"
    if name in GENERIC_HEADERS or True:
        for feature in sample_ids:
            if _KO_RE.match(feature):
                return "kegg"
            if _PFAM_RE.match(feature):
                return "pfam"
            if ";" in feature:
                return "taxa"
    raise BioSIFTRFormatError(
        f"Could not tell what kind of table this is (first column {header!r}).")


def parse_table(text: str, source: str = "bioSIFTR") -> dict:
    """Parse one bioSIFTR table.

    Returns `{"kind": ..., "samples": {sample_id: [records]}}`, with one entry
    per sample column -- so a single-sample table and an integrated matrix are
    handled identically.
    """
    rows = [r for r in csv.reader(io.StringIO(text), delimiter="\t") if r and any(c.strip() for c in r)]
    if len(rows) < 2:
        raise BioSIFTRFormatError("Table is empty or has no data rows.")
    header, *data = rows
    if len(header) < 2:
        raise BioSIFTRFormatError(
            "Expected a feature column plus at least one sample column (tab-separated).")

    feature_col, sample_names = header[0], [h.strip() for h in header[1:]]
    kind = detect_kind(feature_col, [r[0].strip() for r in data[:25] if r])

    samples: dict[str, list[dict]] = {name: [] for name in sample_names}
    for row in data:
        feature = (row[0] or "").strip()
        if not feature:
            continue
        for i, name in enumerate(sample_names, start=1):
            value = _number(row[i]) if i < len(row) else None
            # the matrices are sparse: a zero means "not detected here"
            if value in (None, 0):
                continue
            samples[name].append(_record(kind, feature, value, source))
    return {"kind": kind, "feature_column": feature_col, "samples": samples}


def _record(kind: str, feature: str, value, source: str) -> dict:
    if kind == "taxa":
        lineage, genome = split_lineage(feature)
        return {
            "source": source,
            "organism": lineage,
            "lineage": lineage,
            "genome": genome,
            "abundance": value,
            "evidence": ("Species detected by read-mapping against the MGnify genome "
                         f"catalogue{f' (representative genome {genome})' if genome else ''} "
                         "(bioSIFTR shallow-shotgun profile)."),
        }
    label = "KEGG ortholog" if kind == "kegg" else "Pfam entry"
    return {
        "source": source,
        "annotation_id": feature.upper(),
        "description": None,       # bioSIFTR emits ids only; matching is by id
        "abundance": value,
        "evidence": (f"{label} inferred from the genomes detected in the sample "
                     "(bioSIFTR shallow-shotgun profile)."),
    }


def _number(value):
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return int(number) if number.is_integer() else number


# --------------------------------------------------------------------------- #
# whole-run loading
# --------------------------------------------------------------------------- #

def _pick_tables(directory: Path, patterns: tuple[str, ...]) -> list[Path]:
    found: list[Path] = []
    for pattern in patterns:
        found.extend(sorted(directory.glob(pattern)))
    return found


def load_run(outdir: str | Path, prefer: str = "sm") -> dict:
    """Read a whole bioSIFTR output directory.

    `prefer` picks the mapper when a run produced both: 'sm' (sourmash, always
    present) or 'bwa' (only with `--run_bwa`).

    Returns `{"samples": {id: {"taxa": [...], "functions": [...]}}, "warnings": [...],
    "sources": {...}}` -- ready to hand to `run_multi_analysis`.
    """
    root = Path(outdir).expanduser()
    if not root.is_dir():
        raise BioSIFTRFormatError(f"Not a directory: {root}")

    # tolerate being pointed at the outdir or at one of its subdirectories
    taxonomy_dir = root / TAXONOMY_DIR if (root / TAXONOMY_DIR).is_dir() else root
    function_dir = root / FUNCTION_DIR if (root / FUNCTION_DIR).is_dir() else root
    integrated_dir = root / INTEGRATED_DIR if (root / INTEGRATED_DIR).is_dir() else root

    other = "bwa" if prefer == "sm" else "sm"
    samples: dict[str, dict] = {}
    warnings: list[str] = []
    sources: dict[str, list[str]] = {"taxonomy": [], "functions": []}

    def absorb(path: Path, bucket: str) -> None:
        try:
            parsed = parse_table(path.read_text(), source=f"bioSIFTR ({path.name})")
        except (BioSIFTRFormatError, OSError, UnicodeDecodeError) as exc:
            warnings.append(f"{path.name}: {exc}")
            return
        for sample_id, records in parsed["samples"].items():
            entry = samples.setdefault(sample_id, {"taxa": [], "functions": []})
            entry[bucket].extend(records)
        sources[bucket if bucket == "taxonomy" else "functions"].append(path.name)

    taxonomy_files = (_pick_tables(taxonomy_dir, (f"*_{prefer}_species.tsv",))
                      or _pick_tables(taxonomy_dir, (f"*_{other}_species.tsv", "*_species.tsv")))
    for path in taxonomy_files:
        absorb(path, "taxa")
    sources["taxonomy"] = [p.name for p in taxonomy_files]

    function_files = (_pick_tables(function_dir, (f"*_{prefer}_community_kegg.tsv",
                                                  f"*_{prefer}_community_pfams.tsv"))
                      or _pick_tables(function_dir, ("*_community_kegg.tsv", "*_community_pfams.tsv")))
    for path in function_files:
        absorb(path, "functions")
    sources["functions"] = [p.name for p in function_files]

    # integrated matrices cover every sample at once; use them only to fill gaps
    if not taxonomy_files or not function_files:
        for path in _pick_tables(integrated_dir, ("*_matrix.tsv",)):
            try:
                parsed = parse_table(path.read_text(), source=f"bioSIFTR ({path.name})")
            except (BioSIFTRFormatError, OSError, UnicodeDecodeError) as exc:
                warnings.append(f"{path.name}: {exc}")
                continue
            bucket = "taxa" if parsed["kind"] == "taxa" else "functions"
            if (bucket == "taxa" and taxonomy_files) or (bucket == "functions" and function_files):
                continue
            for sample_id, records in parsed["samples"].items():
                samples.setdefault(sample_id, {"taxa": [], "functions": []})[bucket].extend(records)
            sources["taxonomy" if bucket == "taxa" else "functions"].append(path.name)

    if not samples:
        raise BioSIFTRFormatError(
            f"No bioSIFTR tables found under {root}. Expected {TAXONOMY_DIR}/ and "
            f"{FUNCTION_DIR}/ (or the tables themselves) -- is this the pipeline's --outdir?")

    for sample_id, entry in samples.items():
        if not entry["taxa"]:
            warnings.append(f"{sample_id}: no taxonomic profile found.")
        if not entry["functions"]:
            warnings.append(f"{sample_id}: no functional profile found.")

    return {"samples": samples, "warnings": warnings, "sources": sources,
            "outdir": str(root), "mapper": prefer}


def summarise_run(run: dict) -> dict:
    """Counts for the UI, without shipping every record."""
    return {
        "outdir": run["outdir"],
        "mapper": run["mapper"],
        "samples": [
            {"sample_id": sample_id,
             "taxa": len(entry["taxa"]),
             "functions": len(entry["functions"])}
            for sample_id, entry in sorted(run["samples"].items())
        ],
        "sources": run["sources"],
        "warnings": run["warnings"],
    }
