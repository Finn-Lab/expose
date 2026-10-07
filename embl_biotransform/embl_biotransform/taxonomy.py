"""
taxonomy.py
===========

Move a taxonomic profile from the **GTDB** framework into the **NCBI** one, so
detections and reference associations can be compared at all.

Why this is needed
------------------

bioSIFTR profiles against MGnify's genome catalogues, whose lineages are GTDB
(`d__Bacteria;p__Firmicutes_A;...;s__Bacteroides_fragilis_A`). Reference
associations from ChEMBL are NCBI: an organism name plus `assay_tax_id`, a
real NCBI taxon id. Matching GTDB names against NCBI names fails silently and
asymmetrically, because the two frameworks disagree in both directions:

- GTDB renames that NCBI has since adopted -- `Ruminococcus gnavus` is now
  `Mediterraneibacter gnavus` (taxid 33038), so the *old* name (which is what
  a 2019 paper used) no longer resolves.
- GTDB-only decorations -- `Bacteroides fragilis_A`, `Enterococcus_G italicus`
  carry a suffix marking a polyphyletic split; strip it and NCBI knows them.
- GTDB placeholders -- `CAG-272 sp900556615`, `UBA11517 sp900768545`. In the
  human-gut catalogue **about three quarters** of species are of this kind:
  uncultured taxa with no NCBI species at all. They resolve at genus level or
  not at all, which is correct -- an uncultured species will not be the
  subject of a drug-metabolism assay either.

So the fix is not to translate names but to resolve both sides to a **taxon
id** and match on that, falling back to names only where no id exists.

How a name is resolved
----------------------

Against ENA's taxonomy service (EBI, same estate as everything else here),
walking up the lineage until something matches:

1. the name as given, exact scientific name;
2. the same name with GTDB suffixes stripped (`_A`, `_G`);
3. `any-name`, which also matches synonyms and so recovers the old names
   (`Clostridium bolteae` -> 208479 `Enterocloster bolteae`);
4. the next rank up, and so on to phylum.

`matched_rank` records where on that ladder the answer came from, so a
prediction can say whether a match was species- or genus-level. Unknown names
answer `200` with an empty list rather than an error, and that answer is
cached too -- an uncultured species is looked up once, ever.
"""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

from .fetchers import EMBLAPIError, _get_json

ENA_TAXONOMY_URL = "https://www.ebi.ac.uk/ena/taxonomy/rest"

#: GTDB rank prefixes, deepest first.
_RANK_PREFIXES = ("s__", "g__", "f__", "o__", "c__", "p__", "k__", "d__")
_RANK_NAMES = {"s__": "species", "g__": "genus", "f__": "family", "o__": "order",
               "c__": "class", "p__": "phylum", "k__": "kingdom", "d__": "domain"}

#: GTDB marks polyphyletic splits with a trailing capital suffix:
#: `Firmicutes_A`, `Bacteroides fragilis_A`, `Enterococcus_G italicus`.
_GTDB_SUFFIX_RE = re.compile(r"_[A-Z]+(?=\b|_|$)")

#: GTDB placeholder species: `CAG-272 sp900556615`, `Blautia sp900541345`.
_PLACEHOLDER_RE = re.compile(r"\bsp\d{5,}\b|^(?:CAG|UBA|RUG|HGM|QALS|BX|RC)\d*[-\d]", re.I)


def strip_gtdb_suffix(name: str) -> str:
    """`'Bacteroides fragilis_A'` -> `'Bacteroides fragilis'`."""
    return _GTDB_SUFFIX_RE.sub("", str(name or "")).strip()


def is_placeholder(name: str) -> bool:
    """True for GTDB's uncultured placeholders, which have no NCBI species."""
    return bool(_PLACEHOLDER_RE.search(str(name or "")))


def split_ranks(lineage: str) -> list[tuple[str, str]]:
    """`[('species', 'Bacteroides fragilis_A'), ('genus', 'Bacteroides'), ...]`,
    deepest first, from a GTDB lineage string."""
    out: list[tuple[str, str]] = []
    for part in str(lineage or "").split(";"):
        part = part.strip()
        if not part:
            continue
        for prefix in _RANK_PREFIXES:
            if part.lower().startswith(prefix):
                # strip the polyphyly suffix while the underscores are still
                # there: bioSIFTR writes 'Bacteroides_fragilis_A', and once the
                # underscores become spaces the '_A' is no longer recognisable
                label = strip_gtdb_suffix(part[len(prefix):]).replace("_", " ").strip()
                # a rank prefix with nothing after it means "unassigned here"
                if label:
                    out.append((_RANK_NAMES[prefix], label))
                break
        else:
            out.append(("", strip_gtdb_suffix(part).replace("_", " ").strip()))
    out.reverse()
    return out


@dataclass
class TaxonRef:
    """One resolved taxon, in NCBI terms."""

    tax_id: str | None = None
    scientific_name: str | None = None
    rank: str | None = None
    lineage: str | None = None          # NCBI lineage, as ENA reports it
    matched_rank: str | None = None     # which GTDB rank produced the hit
    matched_name: str | None = None     # the string that actually matched
    query: str | None = None

    @property
    def resolved(self) -> bool:
        return self.tax_id is not None

    def to_dict(self) -> dict:
        return {"tax_id": self.tax_id, "scientific_name": self.scientific_name,
                "rank": self.rank, "lineage": self.lineage,
                "matched_rank": self.matched_rank, "matched_name": self.matched_name}


class TaxonomyResolver:
    """GTDB (or any) taxon names -> NCBI taxon ids, cached and batched.

    Every lookup goes through `fetchers._get_json`, so when an on-disk cache is
    installed each distinct name is fetched once across all runs.
    """

    def __init__(self, base_url: str = ENA_TAXONOMY_URL, max_workers: int = 6,
                 deepest_rank: str = "phylum"):
        self.base_url = base_url.rstrip("/")
        self.max_workers = max_workers
        self.deepest_rank = deepest_rank
        self._memo: dict[str, TaxonRef] = {}

    # -- single lookups -------------------------------------------------------- #

    def _lookup(self, name: str) -> dict | None:
        """First hit for a name, trying the exact scientific name before
        `any-name` (which also matches synonyms, and so is looser)."""
        if not name:
            return None
        for endpoint in ("scientific-name", "any-name"):
            try:
                data = _get_json(f"{self.base_url}/{endpoint}/{_quote(name)}", retries=2)
            except EMBLAPIError:
                continue
            if isinstance(data, list) and data:
                return data[0]
        return None

    def resolve_name(self, name: str) -> TaxonRef:
        """Resolve one name, trying it as given and then suffix-stripped."""
        if name in self._memo:
            return self._memo[name]
        candidates = [name]
        stripped = strip_gtdb_suffix(name)
        if stripped and stripped != name:
            candidates.append(stripped)
        ref = TaxonRef(query=name)
        for candidate in candidates:
            hit = self._lookup(candidate)
            if hit:
                ref = TaxonRef(tax_id=str(hit.get("taxId")) if hit.get("taxId") else None,
                               scientific_name=hit.get("scientificName"),
                               rank=hit.get("rank"), lineage=hit.get("lineage"),
                               matched_name=candidate, query=name)
                break
        self._memo[name] = ref
        return ref

    def resolve_lineage(self, lineage: str) -> TaxonRef:
        """Walk a GTDB lineage from species upwards until a rank resolves."""
        ranks = split_ranks(lineage)
        allowed = _ranks_at_or_above(self.deepest_rank)
        for rank, label in ranks:
            if rank and rank not in allowed:
                continue
            if rank == "species" and is_placeholder(label):
                # 'Blautia sp900541345' will never resolve as a species; its
                # genus is the useful answer, so don't spend a lookup on it
                continue
            ref = self.resolve_name(label)
            if ref.resolved:
                ref.matched_rank = rank or None
                return ref
        return TaxonRef(query=lineage)

    def resolve_taxid(self, tax_id: str) -> TaxonRef:
        """A taxon id straight to its NCBI record -- no name guessing needed.

        Contig taxonomy gives taxids, so this is the direct path; the name
        ladder above is only for sources that report names (GTDB lineages).
        """
        key = f"taxid:{tax_id}"
        if key in self._memo:
            return self._memo[key]
        ref = TaxonRef(query=str(tax_id))
        try:
            data = _get_json(f"{self.base_url}/tax-id/{_quote(tax_id)}", retries=2)
        except EMBLAPIError:
            data = None
        if isinstance(data, dict) and data.get("taxId"):
            ref = TaxonRef(tax_id=str(data["taxId"]), scientific_name=data.get("scientificName"),
                           rank=(data.get("rank") or "").lower() or None,
                           lineage=data.get("lineage"), matched_rank=(data.get("rank") or "").lower() or None,
                           matched_name=data.get("scientificName"), query=str(tax_id))
        self._memo[key] = ref
        return ref

    def resolve_taxids(self, tax_ids: list[str]) -> dict[str, dict]:
        """`{taxid: {"scientific_name", "rank", "lineage"}}`, fetched in parallel."""
        unique = [str(t) for t in dict.fromkeys(tax_ids) if str(t).strip()]
        if not unique:
            return {}
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            refs = list(pool.map(self.resolve_taxid, unique))
        return {t: {"scientific_name": r.scientific_name, "rank": r.rank, "lineage": r.lineage}
                for t, r in zip(unique, refs) if r.resolved}

    # -- batch ----------------------------------------------------------------- #

    def resolve_lineages(self, lineages: list[str]) -> dict[str, TaxonRef]:
        """Resolve many lineages concurrently, de-duplicated."""
        unique = list(dict.fromkeys(l for l in lineages if l))
        if not unique:
            return {}
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            refs = list(pool.map(self.resolve_lineage, unique))
        return dict(zip(unique, refs))

    def annotate(self, taxa: list[dict], progress=None) -> dict:
        """Rewrite detection records into the NCBI framework, in place.

        The GTDB name is kept as `gtdb_organism` / `gtdb_lineage` so the
        evidence trail can still show what the profiler actually reported,
        while `organism`, `lineage` and `tax_id` become NCBI's -- which is what
        the matcher and the reference associations speak.
        """
        lineages = [t.get("lineage") or t.get("organism") or "" for t in taxa]
        if progress:
            progress(0.0, f"resolving {len(set(l for l in lineages if l))} taxa against NCBI…")
        resolved = self.resolve_lineages(lineages)

        counts = {"resolved": 0, "unresolved": 0}
        by_rank: dict[str, int] = {}
        for record, key in zip(taxa, lineages):
            ref = resolved.get(key)
            if ref is None or not ref.resolved:
                counts["unresolved"] += 1
                record.setdefault("gtdb_organism", record.get("organism"))
                record.setdefault("gtdb_lineage", record.get("lineage"))
                record["taxonomy_framework"] = "GTDB (unresolved against NCBI)"
                continue
            counts["resolved"] += 1
            by_rank[ref.matched_rank or "?"] = by_rank.get(ref.matched_rank or "?", 0) + 1
            record.setdefault("gtdb_organism", record.get("organism"))
            record.setdefault("gtdb_lineage", record.get("lineage"))
            record["organism"] = ref.scientific_name
            record["lineage"] = ref.lineage or ref.scientific_name
            record["rank"] = ref.rank
            record["tax_id"] = ref.tax_id
            record["matched_rank"] = ref.matched_rank
            record["taxonomy_framework"] = "NCBI"
            if ref.matched_rank and ref.matched_rank != "species":
                record["evidence"] = (
                    f"{record.get('evidence', '')} GTDB "
                    f"'{record.get('gtdb_organism')}' has no NCBI species; matched at "
                    f"{ref.matched_rank} level ({ref.scientific_name}).").strip()
        if progress:
            progress(1.0, f"{counts['resolved']} of {len(taxa)} taxa resolved")
        return {"counts": counts, "by_rank": by_rank,
                "unresolved_examples": [t.get("gtdb_organism") for t in taxa
                                        if not t.get("tax_id")][:10]}


def _ranks_at_or_above(deepest: str) -> set[str]:
    order = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
    if deepest not in order:
        return set(order)
    return set(order[: order.index(deepest) + 1])


def _quote(value: str) -> str:
    from urllib.parse import quote
    return quote(str(value), safe="")
