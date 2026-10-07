"""
integrate.py
============

Builds one evidence graph out of:

- a metagenome **sample** and what MGnify detected in it (taxa + functional
  annotations),
- the curated **reference associations** (per-microbe or per-gene
  biotransformation statements) from `reference.py`,
- optionally, ChEBI/Rhea/UniProt records for extra compound/reaction
  context.

Node types: sample, taxon, functional_annotation, reference (a curated
microbe- or gene-level association), reaction_class, compound.

Edges always carry `source` and `evidence`, so any prediction can be traced
back to exactly which detection(s) and which curated statement(s) produced
it -- this is what `predict.py` reasons over and `visualize.py` / the web
app display. Detection -> reference edges additionally carry `match_type`,
recording *how* the link was made (see `link_detections_to_references`).
"""

from __future__ import annotations

import re

import networkx as nx


# --------------------------------------------------------------------------- #
# Identifier / name helpers
# --------------------------------------------------------------------------- #

_ID_PREFIXES = ("EC:", "EC ", "INTERPRO:", "KEGG:", "KO:", "PFAM:", "IPR:")

# Words that appear in many functional-family names and so say nothing about
# *which* enzyme it is. They are ignored when matching by name.
NAME_STOPWORDS = {
    "family", "superfamily", "subfamily", "subunit", "alpha", "beta", "gamma", "delta",
    "protein", "proteins", "domain", "domains", "type", "like", "putative", "probable",
    "enzyme", "chain", "large", "small", "component", "dependent", "containing",
    "region", "terminal", "n-terminal", "c-terminal", "and", "the", "with", "for",
    "related", "associated", "class", "group", "uncharacterised", "uncharacterized",
    "conserved", "hypothetical", "site", "binding", "motif", "system",
}


def normalise_id(identifier: str) -> str:
    """'EC:1.14.14.1' -> '1.14.14.1', 'ipr001279' -> 'IPR001279'."""
    s = str(identifier).strip()
    up = s.upper()
    for p in _ID_PREFIXES:
        if up.startswith(p):
            s = s[len(p):]
            break
    return s.strip().upper()


def _is_ec(identifier: str) -> bool:
    return bool(re.fullmatch(r"\d+(\.(\d+|-|n\d+)){0,3}", identifier, flags=re.I))


def ids_match(reference_id: str, detected_id: str) -> bool:
    """Exact id match, with EC wildcards: reference '1.14.13.-' matches a
    detected '1.14.13.39'."""
    ref, det = normalise_id(reference_id), normalise_id(detected_id)
    if not ref or not det:
        return False
    if ref == det:
        return True
    if _is_ec(ref) and _is_ec(det) and "-" in ref:
        prefix = [p for p in ref.split(".") if p != "-"]
        return det.split(".")[: len(prefix)] == prefix
    return False


def _tokens(text: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9\-]+", str(text).lower()) if t]


def distinctive_tokens(name: str) -> list[str]:
    out = []
    for t in _tokens(name):
        for part in {t, *t.split("-")}:
            if len(part) >= 4 and part not in NAME_STOPWORDS and not part.isdigit():
                out.append(part)
    return sorted(set(out))


def name_match(reference_name: str, detected_description: str) -> list[str]:
    """Return the distinctive words of `reference_name` found (as whole words)
    in `detected_description`, if enough of them match to be meaningful:
    at least 2, or all of them when the name has only one distinctive word."""
    ref_tokens = distinctive_tokens(reference_name)
    if not ref_tokens:
        return []
    det_tokens = set()
    for t in _tokens(detected_description):
        det_tokens.add(t)
        det_tokens.update(t.split("-"))
    hits = [t for t in ref_tokens if t in det_tokens]
    return hits if len(hits) >= min(2, len(ref_tokens)) else []


#: Ranks from most to least specific. Anything below `genus` is too coarse to
#: be evidence that a *particular* organism is present.
RANK_ORDER = ("species", "genus", "family", "order", "class", "phylum", "kingdom", "domain")

#: Lineage prefixes, as GTDB and MGnify write them.
_RANK_BY_PREFIX = {"s": "species", "g": "genus", "f": "family", "o": "order",
                   "c": "class", "p": "phylum", "k": "kingdom", "sk": "domain", "d": "domain"}


def taxon_rank(record: dict) -> str | None:
    """The rank of a detected taxon, or None when it can't be told.

    Checked in order of authority: the rank NCBI resolution settled on, an
    explicit `rank` field, then the deepest prefix in the lineage
    (`...;g__Blautia` is a genus). A bare binomial with no rank information at
    all is read as a species, which is what a two-word organism name means
    everywhere this pipeline gets its input.
    """
    for key in ("matched_rank", "rank"):
        value = str(record.get(key) or "").strip().lower()
        if value in RANK_ORDER:
            return value

    lineage = str(record.get("lineage") or record.get("organism") or "")
    if ";" in lineage or "__" in lineage:
        for part in reversed([p.strip() for p in lineage.split(";") if p.strip()]):
            prefix, sep, label = part.partition("__")
            if sep and label.strip():
                rank = _RANK_BY_PREFIX.get(prefix.strip().lower())
                if rank:
                    return rank
        return None

    name = clean_taxon_name(lineage)
    return "species" if len(name.split()) >= 2 else None


def rank_is_at_least(rank: str | None, floor: str) -> bool:
    """True when `rank` is `floor` or more specific. An unknown rank passes:
    dropping a detection because its rank wasn't stated would quietly discard
    perfectly good input (an uploaded `organism` column, say)."""
    if rank is None:
        return True
    if floor not in RANK_ORDER or rank not in RANK_ORDER:
        return True
    return RANK_ORDER.index(rank) <= RANK_ORDER.index(floor)


def clean_taxon_name(name: str) -> str:
    """Turn MGnify lineage strings such as
    'sk__Bacteria;p__Proteobacteria;g__Pseudomonas;s__Pseudomonas_putida'
    into 'Pseudomonas putida' (the most specific named rank)."""
    if not name:
        return ""
    parts = [p.strip() for p in str(name).split(";") if p.strip()]
    for part in reversed(parts):
        label = re.sub(r"^[a-z]{1,2}__", "", part).replace("_", " ").strip()
        if label:
            return label
    return str(name).strip()


def taxon_match(reference_taxon: str, detected_taxon: str) -> str | None:
    """'exact' if the same taxon, 'within-genus' if the reference names a genus
    and the detected taxon belongs to it; otherwise None. A species-level
    reference is *not* matched by a genus-only detection."""
    ref = clean_taxon_name(reference_taxon).lower()
    det = clean_taxon_name(detected_taxon).lower()
    if not ref or not det:
        return None
    if ref == det:
        return "exact"
    if " " not in ref and det.split(" ")[0] == ref:
        return "within-genus"
    return None


def taxon_match_record(reference: dict, detection: dict) -> str | None:
    """Match a microbe reference to a detected taxon, preferring taxon ids.

    Names are a poor key across taxonomies: a profile built on GTDB and a
    reference curated from an older paper can describe the same organism as
    'Mediterraneibacter gnavus' and 'Ruminococcus gnavus' and never match.
    When both sides carry an NCBI taxon id (ChEMBL supplies `assay_tax_id`;
    `taxonomy.TaxonomyResolver` puts one on each detection) that id decides it,
    and the name comparison is only the fallback.

    A detection resolved no deeper than genus keeps `matched_rank`, so a
    genus-level id equal to a *species* reference's id can't happen -- the ids
    simply differ, which is the correct answer.
    """
    ref_id = _clean_taxid(reference.get("tax_id"))
    det_id = _clean_taxid(detection.get("tax_id"))
    if ref_id and det_id:
        if ref_id == det_id:
            return "taxid"
        # the reference's organism may be a genus containing the detected
        # species; fall through to names, which are now both NCBI's
        lineage = detection.get("lineage") or ""
        ref_name = clean_taxon_name(reference.get("entity_id") or "").lower()
        if ref_name and ref_name in [p.strip().lower() for p in str(lineage).split(";")]:
            return "within-lineage"
    return taxon_match(reference.get("entity_id") or "", detection.get("organism") or "")


def _clean_taxid(value) -> str | None:
    text = str(value).strip() if value not in (None, "") else ""
    return text or None


# --------------------------------------------------------------------------- #
# The graph
# --------------------------------------------------------------------------- #

class EvidenceGraph:
    def __init__(self) -> None:
        self.graph = nx.MultiDiGraph()

    # -- metagenome side ---------------------------------------------------- #

    def add_sample(self, sample_id: str, **meta) -> str:
        node = f"sample:{sample_id}"
        self.graph.add_node(node, type="sample", sample_id=sample_id, **meta)
        return node

    def add_detected_taxa(self, sample_node: str, taxa: list[dict]) -> list[str]:
        nodes = []
        for t in taxa:
            name = clean_taxon_name(t.get("organism") or "")
            if not name:
                continue
            tnode = f"taxon:{name}"
            attrs = {**t, "organism": name}
            if self.graph.has_node(tnode):  # same taxon reported twice: sum abundance
                prev = self.graph.nodes[tnode].get("abundance") or 0
                attrs["abundance"] = (prev or 0) + (t.get("abundance") or 0)
            self.graph.add_node(tnode, type="taxon", **attrs)
            if not self.graph.has_edge(sample_node, tnode):
                self.graph.add_edge(sample_node, tnode, source=t.get("source", "user"),
                                    evidence=t.get("evidence", "Taxon detected in sample."))
            nodes.append(tnode)
        return nodes

    def add_detected_functions(self, sample_node: str, annotations: list[dict]) -> list[str]:
        nodes = []
        for a in annotations:
            key = a.get("annotation_id") or a.get("description")
            if not key:
                continue
            anode = f"func:{key}"
            self.graph.add_node(anode, type="functional_annotation", **a)
            if not self.graph.has_edge(sample_node, anode):
                self.graph.add_edge(sample_node, anode, source=a.get("source", "user"),
                                    evidence=a.get("evidence", "Function detected in sample."))
            nodes.append(anode)
        return nodes

    # -- reference (curated) side -------------------------------------------- #

    def add_reference_associations(self, references: list[dict]) -> list[str]:
        """Adds curated microbe-/gene-level biotransformation statements and
        links them to a shared `reaction_class` node so multiple lines of
        evidence for the same reaction converge on one place."""
        nodes = []
        for ref in references:
            rnode = f"ref:{ref['entity_type']}:{ref['entity_id']}:{ref['reaction_class']}"
            self.graph.add_node(rnode, type="reference",
                                tier=1 if ref["entity_type"] == "gene" else 2, **ref)

            rc_node = f"reaction_class:{ref['reaction_class']}"
            self.graph.add_node(rc_node, type="reaction_class", name=ref["reaction_class"])
            self.graph.add_edge(rnode, rc_node, source=ref.get("source", ""), evidence=ref.get("evidence", ""))

            if ref.get("chebi_substrate"):
                cnode = f"compound:{ref['chebi_substrate']}"
                self.graph.add_node(cnode, type="compound", chebi_id=ref["chebi_substrate"],
                                    name=ref.get("substrate_name") or ref["chebi_substrate"])
                if not self.graph.has_edge(rc_node, cnode):
                    self.graph.add_edge(rc_node, cnode, source=ref.get("source", ""),
                                        evidence=f"Reaction class acts on {ref['chebi_substrate']}.")
            nodes.append(rnode)
        return nodes

    # -- linking detections to references ------------------------------------ #

    def link_detections_to_references(self, allow_name_matching: bool = True) -> int:
        """Connects sample-level detections (taxa / functional annotations)
        to matching reference nodes. Returns the number of links made.

        Each link records a `match_type`:

        - gene references
            * ``id``   -- the detected annotation id equals the reference's
                          `entity_id` or one of its `alt_ids` (EC wildcards
                          such as 1.14.13.- are honoured). Strongest.
            * ``name`` -- no id match, but >= 2 distinctive words of the
                          reference's enzyme name appear in the detected
                          family description (generic words like "family",
                          "subunit" are ignored). Weaker; can be disabled.
        - microbe references
            * ``taxid``          -- both sides carry the same NCBI taxon id.
                                    Strongest, and immune to naming drift
                                    between GTDB and NCBI.
            * ``within-lineage`` -- the reference taxon appears in the
                                    detection's resolved NCBI lineage.
            * ``exact``          -- same taxon name (no ids to compare).
            * ``within-genus``   -- the reference names a genus and the
                                    detected taxon belongs to it.
        """
        g = self.graph
        taxon_nodes = [(n, d) for n, d in g.nodes(data=True) if d.get("type") == "taxon"]
        func_nodes = [(n, d) for n, d in g.nodes(data=True) if d.get("type") == "functional_annotation"]
        ref_nodes = [(n, d) for n, d in g.nodes(data=True) if d.get("type") == "reference"]
        n_links = 0

        # Index the detections once, rather than comparing every reference
        # against every detection. A fully built ChEMBL set is tens of
        # thousands of references and a metagenome has ~15,000 annotations, so
        # the pairwise form is hundreds of millions of comparisons; the indexes
        # turn the common cases (exact id, exact taxon, taxon id) into lookups.
        index = _DetectionIndex(taxon_nodes, func_nodes)

        for rnode, rdata in ref_nodes:
            if rdata["entity_type"] == "microbe":
                for tnode, tdata in index.candidate_taxa(rdata):
                    kind = taxon_match_record(rdata, tdata)
                    if kind:
                        how = {"taxid": f"NCBI taxon {tdata.get('tax_id')}",
                               "within-lineage": "reference taxon appears in the detected NCBI lineage"}.get(kind, kind)
                        origin = (f" [reported by the profiler as '{tdata['gtdb_organism'].split(';')[-1]}']"
                                  if tdata.get("gtdb_organism") and tdata.get("taxonomy_framework") == "NCBI" else "")
                        g.add_edge(tnode, rnode, source="integration", match_type=kind,
                                   evidence=(f"Detected taxon '{tdata['organism']}'{origin} matches reference "
                                             f"microbe '{rdata['entity_id']}' ({how})."))
                        n_links += 1
            else:  # gene
                ref_ids = [rdata["entity_id"], *rdata.get("alt_ids", [])]
                matched_by_id = set()
                for fnode, fdata, matched_id in index.id_matches(ref_ids):
                    matched_by_id.add(fnode)
                    g.add_edge(fnode, rnode, source="integration", match_type="id",
                               evidence=(f"Detected annotation {fdata.get('annotation_id')} matches "
                                         f"reference id {matched_id}."))
                    n_links += 1
                if allow_name_matching:
                    for fnode, fdata in index.name_candidates(rdata.get("entity_name", "")):
                        if fnode in matched_by_id:
                            continue
                        words = name_match(rdata.get("entity_name", ""), fdata.get("description") or "")
                        if words:
                            g.add_edge(fnode, rnode, source="integration", match_type="name",
                                       matched_words=words,
                                       evidence=(f"Detected family '{fdata.get('description')}' shares the "
                                                 f"words {', '.join(words)} with reference enzyme "
                                                 f"'{rdata.get('entity_name')}' (name match, no id match)."))
                            n_links += 1
        return n_links

    # -- querying ------------------------------------------------------------ #

    def evidence_trail(self, source_node: str, target_node: str, cutoff: int = 4) -> list[list[dict]]:
        trails = []
        for path in nx.all_simple_edge_paths(self.graph, source_node, target_node, cutoff=cutoff):
            trail = []
            for u, v, k in path:
                edge = self.graph.edges[u, v, k]
                trail.append({"from": u, "to": v, "source": edge["source"], "evidence": edge["evidence"]})
            trails.append(trail)
        return trails

    def summary(self) -> dict:
        by_type: dict[str, int] = {}
        for _, data in self.graph.nodes(data=True):
            by_type[data.get("type", "?")] = by_type.get(data.get("type", "?"), 0) + 1
        return {"nodes": self.graph.number_of_nodes(), "edges": self.graph.number_of_edges(), "by_type": by_type}

    def to_dict(self) -> dict:
        """JSON-serialisable node/edge lists (used by the web app)."""
        nodes = [{"id": n, **{k: v for k, v in d.items() if _jsonable(v)}}
                 for n, d in self.graph.nodes(data=True)]
        edges = [{"source_node": u, "target_node": v, **{k: val for k, val in d.items() if _jsonable(val)}}
                 for u, v, d in self.graph.edges(data=True)]
        return {"nodes": nodes, "edges": edges}


class _DetectionIndex:
    """Lookup tables over a sample's detections.

    Only narrows the candidate set; every candidate still goes through the
    same `ids_match` / `name_match` / `taxon_match_record` decision, so what
    counts as a match is unchanged.
    """

    def __init__(self, taxon_nodes, func_nodes):
        self.taxon_nodes = taxon_nodes
        self.func_nodes = func_nodes

        # -- functional annotations ------------------------------------------ #
        self._by_id: dict[str, list] = {}
        self._ec_nodes: list = []
        self._by_token: dict[str, list] = {}
        for fnode, fdata in func_nodes:
            detected = normalise_id(fdata.get("annotation_id") or "")
            if detected:
                self._by_id.setdefault(detected, []).append((fnode, fdata))
                if _is_ec(detected):
                    # EC wildcards ('1.14.13.-') need a prefix scan, but only
                    # over ids that are themselves EC numbers
                    self._ec_nodes.append((fnode, fdata, detected))
            for token in _description_tokens(fdata.get("description") or ""):
                self._by_token.setdefault(token, []).append((fnode, fdata))

        # -- taxa -------------------------------------------------------------- #
        self._by_taxid: dict[str, list] = {}
        self._by_name: dict[str, list] = {}
        self._by_genus: dict[str, list] = {}
        self._lineage_taxa: list = []
        for tnode, tdata in taxon_nodes:
            tax_id = _clean_taxid(tdata.get("tax_id"))
            if tax_id:
                self._by_taxid.setdefault(tax_id, []).append((tnode, tdata))
                # a reference may name a genus that appears in this lineage
                if tdata.get("lineage"):
                    self._lineage_taxa.append((tnode, tdata))
            name = clean_taxon_name(tdata.get("organism") or "").lower()
            if name:
                self._by_name.setdefault(name, []).append((tnode, tdata))
                self._by_genus.setdefault(name.split(" ")[0], []).append((tnode, tdata))

    # -- genes ----------------------------------------------------------------- #

    def id_matches(self, reference_ids):
        """`(node, data, matched_reference_id)` for every exact or EC-wildcard hit."""
        seen = set()
        for reference_id in reference_ids:
            normalised = normalise_id(reference_id)
            if not normalised:
                continue
            for fnode, fdata in self._by_id.get(normalised, ()):
                if fnode not in seen:
                    seen.add(fnode)
                    yield fnode, fdata, reference_id
            if _is_ec(normalised) and "-" in normalised:
                for fnode, fdata, detected in self._ec_nodes:
                    if fnode not in seen and ids_match(reference_id, detected):
                        seen.add(fnode)
                        yield fnode, fdata, reference_id

    def name_candidates(self, reference_name: str):
        """Detections sharing at least one distinctive word with the reference
        name -- a superset of what `name_match` will accept."""
        tokens = distinctive_tokens(reference_name)
        if not tokens:
            return
        seen = set()
        for token in tokens:
            for fnode, fdata in self._by_token.get(token, ()):
                if fnode not in seen:
                    seen.add(fnode)
                    yield fnode, fdata

    # -- microbes ---------------------------------------------------------------- #

    def candidate_taxa(self, reference: dict):
        """Detections that could match this microbe reference: the same taxon
        id, the same name, anything in that genus, and -- when the reference
        carries a taxon id -- taxa with a resolved lineage to search."""
        seen, out = set(), []

        def add(items):
            for tnode, tdata in items:
                if tnode not in seen:
                    seen.add(tnode)
                    out.append((tnode, tdata))

        tax_id = _clean_taxid(reference.get("tax_id"))
        if tax_id:
            add(self._by_taxid.get(tax_id, ()))
            add(self._lineage_taxa)
        name = clean_taxon_name(reference.get("entity_id") or "").lower()
        if name:
            add(self._by_name.get(name, ()))
            if " " not in name:                      # a genus-level reference
                add(self._by_genus.get(name, ()))
        return out


def _description_tokens(description: str) -> set[str]:
    out = set()
    for token in _tokens(description):
        out.add(token)
        out.update(token.split("-"))
    return out


def _jsonable(v) -> bool:
    return isinstance(v, (str, int, float, bool, type(None), list, dict))
