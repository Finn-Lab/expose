"""
predict.py
==========

Scores candidate biotransformations for a metagenome sample using the
two-tier evidence model:

  Tier 1 (strongest) -- the specific ENZYME/gene behind a reference
                         association was directly detected in the sample's
                         functional annotation.
  Tier 2 (weaker)     -- a MICROBE known (from the reference data) to carry
                         out the reaction is present in the sample's
                         taxonomic profile, but its enzyme was not
                         independently confirmed there.

A reaction class can accumulate both kinds of evidence (e.g. three
supporting microbes and one directly-detected enzyme); agreement across
independent references and across tiers raises confidence, matching the
idea that convergent, independent evidence is stronger than any single
line of it.

Each "hit" is one curated reference statement supported by one or more
detections in the sample; counts in the score are of *distinct reference
statements*, so five detected strains of one genus that all match the same
genus-level reference still count once.

Nothing here claims certainty -- every prediction keeps its full evidence
trail (which detections + which curated reference statements produced it)
so a human can audit the basis for the call.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .integrate import EvidenceGraph

# How much a Tier-1 hit counts towards the score, by how it was matched.
TIER1_MATCH_WEIGHT = {"id": 1.0, "name": 0.5}


@dataclass
class BiotransformationPrediction:
    reaction_class: str
    tier1_hits: list[dict] = field(default_factory=list)  # enzyme-level, directly detected
    tier2_hits: list[dict] = field(default_factory=list)  # microbe-level, inferred
    substrate_chebi: str | None = None

    @staticmethod
    def _best_match(hit: dict) -> str:
        types = {d["match_type"] for d in hit.get("detections", [])}
        return "id" if "id" in types else ("name" if "name" in types else next(iter(types), "?"))

    @property
    def tier1_strength(self) -> float:
        return sum(TIER1_MATCH_WEIGHT.get(self._best_match(h), 0.5) for h in self.tier1_hits)

    @property
    def score(self) -> float:
        """0-1 confidence.
        - Any Tier-1 (enzyme-detected) hit dominates the score; multiple
          independent Tier-1 hits push it toward 1.0. A Tier-1 hit made only
          by enzyme-name similarity (no id match) counts half.
        - Tier-2-only (microbe-present) support caps out lower, reflecting
          that the specific catalytic machinery wasn't itself confirmed.
        - Having both tiers agree on the same reaction class adds a small
          convergence bonus on top of the Tier-1 base.
        """
        t1 = min(self.tier1_strength, 3) / 3          # 0..1
        t2 = min(len(self.tier2_hits), 3) / 3          # 0..1
        if self.tier1_hits:
            base = 0.6 + 0.3 * t1                       # 0.6 - 0.9
            bonus = 0.1 if self.tier2_hits else 0.0      # convergence bonus
            return round(min(base + bonus, 1.0), 3)
        # microbe evidence only
        return round(0.15 + 0.35 * t2, 3)                # 0.15 - 0.5

    @property
    def confidence_label(self) -> str:
        if self.tier1_hits:
            by_name_only = all(self._best_match(h) == "name" for h in self.tier1_hits)
            return ("high - enzyme directly detected in metagenome"
                    + (" (matched by name only)" if by_name_only else "")
                    + (", corroborated by microbe presence" if self.tier2_hits else ""))
        return "moderate/low - inferred only from microbe presence (enzyme not confirmed)"

    @property
    def confidence_level(self) -> str:
        return "high" if self.tier1_hits else "low"

    @property
    def evidence_summary(self) -> list[str]:
        lines = []
        for h in self.tier1_hits:
            dets = "; ".join(f"'{d['label']}' ({d['match_type']} match)" for d in h["detections"])
            lines.append(f"[Tier 1 / enzyme] {h['entity_name']} ({h['entity_id']}) detected in sample "
                         f"via functional annotation {dets}.")
        for h in self.tier2_hits:
            dets = "; ".join(f"'{d['label']}' ({d['match_type']})" for d in h["detections"])
            lines.append(f"[Tier 2 / microbe] {h['entity_name']} present in sample's taxonomic profile "
                         f"as {dets} (enzyme itself not confirmed).")
        return lines

    def to_dict(self) -> dict:
        return {
            "reaction_class": self.reaction_class,
            "score": self.score,
            "confidence_level": self.confidence_level,
            "confidence_label": self.confidence_label,
            "substrate_chebi": self.substrate_chebi,
            "tier1_hits": self.tier1_hits,
            "tier2_hits": self.tier2_hits,
            "evidence_summary": self.evidence_summary,
        }


class BiotransformationPredictor:
    """Given an EvidenceGraph that already has a sample's detections and the
    reference associations loaded and linked (`link_detections_to_references`),
    produce ranked, evidence-backed biotransformation predictions."""

    def __init__(self, graph: EvidenceGraph):
        self.graph = graph

    def predict(self, sample_node: str) -> list[BiotransformationPrediction]:
        g = self.graph.graph
        detected = {n for n in g.successors(sample_node)
                    if g.nodes[n].get("type") in ("taxon", "functional_annotation")}

        by_reaction: dict[str, BiotransformationPrediction] = {}

        for ref_node, rdata in g.nodes(data=True):
            if rdata.get("type") != "reference":
                continue
            reaction_class = rdata["reaction_class"]
            pred = by_reaction.setdefault(
                reaction_class,
                BiotransformationPrediction(reaction_class=reaction_class,
                                            substrate_chebi=rdata.get("chebi_substrate")),
            )
            if not pred.substrate_chebi and rdata.get("chebi_substrate"):
                pred.substrate_chebi = rdata["chebi_substrate"]

            detections = []
            for det_node in g.predecessors(ref_node):
                if det_node not in detected:
                    continue
                for e in g.get_edge_data(det_node, ref_node).values():
                    if e.get("source") != "integration":
                        continue
                    ddata = g.nodes[det_node]
                    detections.append({
                        "node": det_node,
                        "label": ddata.get("organism") or ddata.get("description") or ddata.get("annotation_id"),
                        "annotation_id": ddata.get("annotation_id"),
                        "abundance": ddata.get("abundance"),
                        "match_type": e.get("match_type", "?"),
                        "evidence": e.get("evidence"),
                    })
            if not detections:
                continue

            hit = {k: v for k, v in rdata.items() if k != "type"}
            hit["node"] = ref_node
            hit["detections"] = detections
            if rdata["entity_type"] == "gene":
                pred.tier1_hits.append(hit)
            else:
                pred.tier2_hits.append(hit)

        # only keep reaction classes with at least one hit in *this* sample
        predictions = [p for p in by_reaction.values() if p.tier1_hits or p.tier2_hits]
        return sorted(predictions, key=lambda p: (-p.score, p.reaction_class))
