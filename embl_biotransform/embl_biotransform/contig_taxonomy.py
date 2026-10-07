"""
contig_taxonomy.py
==================

Attribute each detected enzyme to the organism whose contig carries it.

The problem this solves
-----------------------

The per-analysis InterPro/Pfam/KO summary tables are community-wide counts:
they say an enzyme family occurs *somewhere* in the sample, with no link to
which organism carries it. A prediction built on them claims only "this
community contains the enzyme, and separately contains the microbe" -- which
is far weaker than it looks, and makes predictions too liberal.

Pipeline v6 publishes the pieces to do better. Every CDS in
`<assembly>_annotation_summary.gff.gz` is on a named contig:

    ERZ29588192_1  Pyrodigal  CDS  760  1821  ...  ID=ERZ29588192_1_1;kegg=ko:K03453;pfam=PF01758;interpro=IPR002657,IPR038770

and `<assembly>_contigs_taxonomy.tsv.gz` assigns each contig an **NCBI taxid
lineage**:

    ERZ29588192_1   taxid assigned   based on 543/547 ORFs   1;131567;2;1783272;201174;1760;85004;31953;1678   1.00;0.96;...

Joining on the contig gives enzyme -> taxon, in NCBI taxids, with the
classifier's per-rank confidence.

Why not InterProScan
--------------------

The same join can be made from `<assembly>_interproscan.tsv.gz`, whose protein
ids (`ERZ29588192_2784_12`) contain the contig (`ERZ29588192_2784`). But that
file is **~1 GB compressed per analysis**, against **12 MB** for the
annotation summary -- and the summary additionally carries Pfam, KEGG and GO,
not InterPro alone. Same answer, ~80x less to download, so this module reads
the summary. (The InterProScan file remains the only source for *which part*
of a protein matched, which nothing here needs.)

What comes out is small: on a real gut assembly, 392k CDS rows reduce to
24k annotations over 320k (annotation, taxid) pairs -- about 4 MB, cached, so
the download happens once per analysis.

The catch, and why this makes predictions stricter
--------------------------------------------------

Many contigs cannot be classified below domain. On a real assembly the
enzymes behind several predictions sat on contigs assigned only to taxid 2
(Bacteria) or 131567 (cellular organisms). Those occurrences are real but
unattributable, so `attribution()` separates them: counts at genus/species
rank are evidence that a *particular* organism carries the enzyme; the rest
is not.
"""

from __future__ import annotations

import collections
import gzip
import io
import re
from typing import Callable, Iterable, Iterator

import requests

from .fetchers import EMBLAPIError

#: Attribute keys in the annotation summary GFF that name a functional entry.
GFF_ANNOTATION_KEYS = ("interpro", "pfam", "kegg")

_ATTR_RE = re.compile(r"(?:^|;)(%s)=([^;]+)" % "|".join(GFF_ANNOTATION_KEYS), re.I)

#: Taxids at or above "cellular organisms" carry no information about which
#: organism a contig came from.
UNINFORMATIVE_TAXIDS = {"1", "131567"}

DOWNLOAD_ALIASES = {
    "gff": ("_annotation_summary.gff.gz",),
    "taxonomy": ("_contigs_taxonomy.tsv.gz",),
}


def stream_gzip_lines(url: str, timeout: int = 300,
                      progress: Callable[[int], None] | None = None) -> Iterator[str]:
    """Decompress a remote gzip line by line.

    Streamed rather than downloaded: the annotation summary is ~100 MB once
    decompressed, and nothing here needs more than one line at a time.
    """
    headers = {"User-Agent": "embl-biotransform/0.3 (research tool)"}
    try:
        response = requests.get(url, headers=headers, timeout=timeout, stream=True)
        response.raise_for_status()
    except Exception as exc:  # noqa: BLE001
        raise EMBLAPIError(f"GET {url} failed: {exc}") from exc

    seen = 0
    with response:
        raw = response.raw
        raw.decode_content = True
        with gzip.GzipFile(fileobj=raw) as gz:
            for line in io.TextIOWrapper(gz, encoding="utf-8", errors="replace"):
                if progress is not None:
                    seen += 1
                    if seen % 50_000 == 0:
                        progress(seen)
                yield line


def parse_contig_taxonomy(lines: Iterable[str]) -> dict[str, list[str]]:
    """`{contig: [taxid, ...]}` -- the NCBI lineage, root first."""
    out: dict[str, list[str]] = {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        contig, lineage = parts[0].strip(), parts[3].strip()
        if not contig or not lineage:
            continue
        taxids = [t.strip() for t in lineage.split(";") if t.strip()]
        if taxids:
            out[contig] = taxids
    return out


def iter_gff_annotations(lines: Iterable[str]) -> Iterator[tuple[str, list[str]]]:
    """`(contig, [accession, ...])` per CDS row of an annotation summary GFF."""
    for line in lines:
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        accessions = []
        for _key, value in _ATTR_RE.findall(cols[8]):
            for accession in value.split(","):
                accession = accession.strip()
                if accession.lower().startswith("ko:"):
                    accession = accession[3:]
                accession = accession.upper()
                if accession and accession != "-":
                    accessions.append(accession)
        if accessions:
            yield cols[0].strip(), accessions


def build_links(taxonomy_lines: Iterable[str], gff_lines: Iterable[str]) -> dict:
    """Join the two files into `{accession: {taxid: count}}` plus totals."""
    contigs = parse_contig_taxonomy(taxonomy_lines)
    links: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    cds = unclassified = 0
    for contig, accessions in iter_gff_annotations(gff_lines):
        cds += 1
        lineage = contigs.get(contig)
        if not lineage:
            unclassified += 1
            continue
        deepest = lineage[-1]
        for accession in accessions:
            links[accession][deepest] += 1
    return {
        "annotations": {a: dict(c) for a, c in links.items()},
        "contigs": len(contigs),
        "cds_with_annotations": cds,
        "cds_on_unclassified_contigs": unclassified,
        "pairs": sum(len(c) for c in links.values()),
    }


def find_download(detail: dict, kind: str) -> dict | None:
    """The annotation-summary GFF or contig-taxonomy file of an analysis."""
    for suffix in DOWNLOAD_ALIASES[kind]:
        for f in detail.get("downloads") or []:
            if str(f.get("alias", "")).lower().endswith(suffix.lower()) and f.get("url"):
                return f
    return None


def attribution(links: dict, accession: str, taxon_info: dict,
                floor_ranks: tuple[str, ...] = ("species", "genus")) -> dict:
    """Split one annotation's occurrences into attributable and not.

    `taxon_info` maps taxid -> `{"scientific_name", "rank"}`. A count only
    supports "this organism carries this enzyme" when its contig was
    classified to `floor_ranks`; everything else (domain, phylum, or
    unclassified) is counted separately rather than quietly treated the same.
    """
    counts = (links.get("annotations") or {}).get(accession) or {}
    attributed, coarse = [], 0
    for taxid, count in counts.items():
        info = taxon_info.get(str(taxid)) or {}
        rank = (info.get("rank") or "").lower()
        if taxid in UNINFORMATIVE_TAXIDS or rank not in floor_ranks:
            coarse += count
            continue
        attributed.append({"tax_id": str(taxid), "count": count,
                           "scientific_name": info.get("scientific_name"), "rank": rank})
    attributed.sort(key=lambda t: -t["count"])
    return {
        "accession": accession,
        "taxa": attributed,
        "attributed": sum(t["count"] for t in attributed),
        "unattributed": coarse,
        "total": sum(counts.values()),
    }


def observed_taxids(links: dict, accessions: Iterable[str]) -> list[str]:
    """Every taxid seen for the given accessions, for one batched lookup."""
    seen: set[str] = set()
    table = links.get("annotations") or {}
    for accession in accessions:
        seen.update(str(t) for t in (table.get(accession) or {}))
    return sorted(seen - UNINFORMATIVE_TAXIDS)
