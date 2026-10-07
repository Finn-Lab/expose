"""
chembl_references.py
====================

Turn ChEMBL "Bacterial Biotransformation" assays into the reference
associations this pipeline predicts *from*, plus -- for the rare community
assays -- the measured outcomes to validate those predictions *against*.

ChEMBL records these experiments at three levels of resolution, which line up
with the pipeline's evidence tiers:

| tier         | how the assay looks in ChEMBL                                   | becomes                    |
|--------------|-----------------------------------------------------------------|----------------------------|
| `gene`       | `confidence_score` 9, a real single-protein target resolving to a | a `gene` reference         |
|              | UniProt accession with InterPro / Pfam cross-references           | (Tier 1 evidence)          |
| `microbe`    | `confidence_score` 0, target is the generic ADMET placeholder;    | a `microbe` reference      |
|              | the organism is in `assay_organism` / `assay_tax_id`              | (Tier 2 evidence)          |
| `microbiome` | as `microbe`, but the organism is a metagenome and the assay      | an *observation* tied to   |
|              | parameters carry ENA study / sample accessions                    | one MGnify analysis        |

Community (microbiome) assays are rare -- almost no metagenome has measured
drug-metabolism data. Where one does exist it is worth more than a reference:
it says which biotransformations that specific community was *observed* to
perform, so predictions made from gene and microbe evidence can be checked
against it, and the confirmed ones become the strongest evidence available.

Both positive and negative results are kept. A negative ("no biotransformation
occurred") is genuine evidence of absence and is rare in curated sets, so it is
recorded with `observed: false` rather than dropped. Only `observed: true`
records are associations to predict from, which is what `references()` returns
by default.

ChEMBL never states *what chemistry* a biotransformation is, only that one
happened, so `reaction_class` stays at that level -- "Biotransformation of
<drug>", one class per substrate.

Substrate ChEBI ids come from UniChem: ChEMBL molecule records have no ChEBI
cross-reference of their own.
"""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import requests

from .fetchers import EMBLAPIError, MGnifyClient, _get_json

CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
UNICHEM_URL = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
INTERPRO_BASE_URL = "https://www.ebi.ac.uk/interpro/api/entry/interpro"

#: ChEMBL's placeholder target, used when an assay has no protein target.
GENERIC_TARGET = "CHEMBL612558"

#: The activity `standard_type` that marks these experiments.
ACTIVITY_STANDARD_TYPE = "Bacterial Biotransformation"

#: InterPro entry types specific enough to be evidence. A
#: `homologous_superfamily` such as IPR029058 (Alpha/Beta hydrolase fold)
#: spans ~1.7 million proteins and would match almost any metagenome, so it is
#: excluded by default: matching on it would manufacture Tier-1 hits.
USEFUL_INTERPRO_TYPES = ("family", "domain")

# Each published screen phrases its result differently, so outcomes are read
# from two places: `standard_text_value`, a short structured verdict ("Compound
# metabolized"), and failing that the free-text `activity_comment`.
#
# Negatives must be tested before positives in both, because the negative
# wording contains the positive as a substring ("No Biotransformation Occurred"
# contains "biotransformation occurred"; "not metabolized" contains
# "metabolized"). ChEMBL's comments also carry the typo "Oberved" for
# "Observed", hence `obs?erved`.
_NEGATIVE_RE = re.compile(
    r"\bnot\s+(?:metabol|biotransformed)|"
    r"no biotransformation|no metabolit|could not be mediated|was not detected", re.I)
_POSITIVE_RE = re.compile(
    r"\bis biotransformed\b|\bmetabol(?:ized|ised)\b|"
    r"biotransformation occurred|proven to be mediated|obs?erved biotransformation", re.I)
_METABOLITE_RE = re.compile(r"metabolite with m/z\s*([\d.]+)", re.I)


class ChEMBLReferenceError(RuntimeError):
    """Raised when an assay can't be turned into reference records."""


def classify_assay(assay: dict) -> str:
    """'gene', 'microbe' or 'microbiome' -- see the table above."""
    target = assay.get("target_chembl_id")
    if assay.get("confidence_score") and target and target != GENERIC_TARGET:
        return "gene"
    organism = (assay.get("assay_organism") or "").lower()
    if "metagenome" in organism or "microbiome" in organism or _ena_accessions(assay).get("sample"):
        return "microbiome"
    return "microbe"


def observation_from_comment(comment: str | None) -> bool | None:
    """True / False for a stated outcome, None when the text says neither."""
    text = comment or ""
    if _NEGATIVE_RE.search(text):
        return False
    if _POSITIVE_RE.search(text):
        return True
    return None


def observation_from_activity(activity: dict) -> bool | None:
    """The outcome of one ChEMBL activity row.

    `standard_text_value` is a curated verdict ("Compound metabolized") and is
    trusted first; the free-text comment is the fallback for the screens that
    leave it empty.
    """
    verdict = observation_from_comment(activity.get("standard_text_value"))
    if verdict is not None:
        return verdict
    return observation_from_comment(activity.get("activity_comment"))


def _ena_accessions(assay: dict) -> dict:
    """The ENA study / sample accessions a community assay carries in its
    free-text assay parameters."""
    out = {}
    for p in assay.get("assay_parameters") or []:
        name = (p.get("type") or "").lower()
        value = (p.get("value") or p.get("text_value") or "").strip()
        if not value:
            continue
        if "study accession" in name:
            out["study"] = value
        elif "sample accession" in name:
            out["sample"] = value
    return out


@dataclass
class ChEMBLExtraction:
    """What `ChEMBLBiotransformationSource.build()` produces."""

    records: list[dict] = field(default_factory=list)
    observations: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def references(self, observed_only: bool = True) -> list[dict]:
        """Reference associations ready for `run_analysis`.

        Negative records are excluded by default: `run_analysis` treats every
        record it is given as evidence *for* a reaction, so passing a
        "no biotransformation occurred" record would invert its meaning.
        """
        return [{k: v for k, v in r.items() if k != "observed"}
                for r in self.records if r.get("observed") or not observed_only]

    def summary(self) -> dict:
        kinds = [r["entity_type"] for r in self.records]
        return {
            "records": len(self.records),
            "genes": kinds.count("gene"),
            "microbes": kinds.count("microbe"),
            "positive": sum(1 for r in self.records if r.get("observed")),
            "negative": sum(1 for r in self.records if r.get("observed") is False),
            "reaction_classes": len({r["reaction_class"] for r in self.records}),
            "communities": len(self.observations),
            "community_observations": sum(len(o["observations"]) for o in self.observations),
            "warnings": len(self.warnings),
        }

    def to_dict(self) -> dict:
        return {"summary": self.summary(), "records": self.records,
                "observations": self.observations, "warnings": self.warnings}


class ChEMBLBiotransformationSource:
    """Build reference associations from ChEMBL biotransformation assays.

        source = ChEMBLBiotransformationSource()
        out = source.build(["CHEMBL5303650", "CHEMBL5724809", "CHEMBL5725002"])
        run_analysis(taxa, functions, out.references(), sample_id="...")
    """

    def __init__(self, base_url: str = CHEMBL_BASE_URL, resolve_chebi: bool = True,
                 interpro_types: tuple[str, ...] = USEFUL_INTERPRO_TYPES,
                 mgnify_client: MGnifyClient | None = None, max_workers: int = 8):
        self.base_url = base_url.rstrip("/")
        self.resolve_chebi = resolve_chebi
        self.interpro_types = tuple(t.lower() for t in interpro_types) if interpro_types else ()
        self.mgnify = mgnify_client
        self.max_workers = max_workers
        self._chebi: dict[str, str | None] = {}
        self._interpro_type: dict[str, str | None] = {}

    # -- ChEMBL ------------------------------------------------------------- #

    def fetch_assay(self, assay_id: str) -> dict:
        return _get_json(f"{self.base_url}/assay/{assay_id}.json")

    def fetch_target(self, target_id: str) -> dict:
        return _get_json(f"{self.base_url}/target/{target_id}.json")

    def search_assays(self, keyword: str, limit: int = 50) -> list[dict]:
        """Find biotransformation assays by accession, organism or description.

        There is no assay-level flag for "this is a biotransformation
        experiment" -- it is a property of the *activities* -- so this searches
        activities with `standard_type` "Bacterial Biotransformation" and
        groups them back up into their assays.
        """
        keyword = (keyword or "").strip()
        if not keyword:
            return []
        if re.fullmatch(r"CHEMBL\d+", keyword, re.I):
            try:
                assay = self.fetch_assay(keyword.upper())
            except EMBLAPIError:
                return []
            return [self._assay_summary(assay, None)]

        found: dict[str, dict] = {}
        for field in ("assay_organism__icontains", "assay_description__icontains"):
            params = {"standard_type": ACTIVITY_STANDARD_TYPE, field: keyword, "limit": 1000,
                      "only": "assay_chembl_id,assay_organism,assay_description,target_chembl_id"}
            try:
                data = _get_json(f"{self.base_url}/activity.json", params=params)
            except EMBLAPIError:
                continue
            for act in data.get("activities") or []:
                accession = act.get("assay_chembl_id")
                if not accession:
                    continue
                entry = found.setdefault(accession, {
                    "assay_chembl_id": accession,
                    "assay_organism": act.get("assay_organism"),
                    "description": act.get("assay_description"),
                    "target_chembl_id": act.get("target_chembl_id"),
                    "activities": 0})
                entry["activities"] += 1
            if len(found) >= limit:
                break
        out = sorted(found.values(), key=lambda e: -e["activities"])[:limit]
        for entry in out:
            entry["tier"] = self._tier_from_summary(entry)
        return out

    @staticmethod
    def _tier_from_summary(entry: dict) -> str:
        """Tier from a search hit, which carries no confidence_score: a real
        target means 'gene', a metagenome organism means 'microbiome'."""
        target = entry.get("target_chembl_id")
        organism = (entry.get("assay_organism") or "").lower()
        if target and target != GENERIC_TARGET:
            return "gene"
        if "metagenome" in organism or "microbiome" in organism:
            return "microbiome"
        return "microbe"

    def _assay_summary(self, assay: dict, activities: int | None) -> dict:
        return {"assay_chembl_id": assay.get("assay_chembl_id"),
                "assay_organism": assay.get("assay_organism"),
                "description": assay.get("description"),
                "target_chembl_id": assay.get("target_chembl_id"),
                "activities": activities,
                "tier": classify_assay(assay)}

    def fetch_activities(self, assay_id: str, max_items: int = 5000) -> list[dict]:
        """Every activity measured in one assay (paged)."""
        items: list[dict] = []
        url = f"{self.base_url}/activity.json"
        params: dict | None = {"assay_chembl_id": assay_id, "limit": 1000}
        while url and len(items) < max_items:
            data = _get_json(url, params=params)
            items.extend(data.get("activities") or [])
            nxt = (data.get("page_meta") or {}).get("next")
            url, params = (f"https://www.ebi.ac.uk{nxt}" if nxt else None), None
        return items[:max_items]

    # -- cross-references ---------------------------------------------------- #

    def _target_protein(self, target_id: str) -> dict | None:
        """UniProt accession, name and the annotation ids a metagenome could
        actually carry for the assay's protein target."""
        target = self.fetch_target(target_id)
        components = target.get("target_components") or []
        if not components:
            return None
        component = components[0]
        alt_ids, dropped = [], []
        for x in component.get("target_component_xrefs") or []:
            src, xid = (x.get("xref_src_db") or "").lower(), x.get("xref_id")
            if not xid:
                continue
            if src == "pfam":
                alt_ids.append(xid)
            elif src == "interpro":
                kind = self._interpro_entry_type(xid)
                if not self.interpro_types or kind in self.interpro_types:
                    alt_ids.append(xid)
                else:
                    dropped.append(f"{xid} ({kind})")
        return {
            "accession": component.get("accession"),
            "name": component.get("component_description") or target.get("pref_name"),
            "alt_ids": alt_ids,
            "dropped_ids": dropped,
            "organism": target.get("organism"),
            "tax_id": target.get("tax_id"),
        }

    def _interpro_entry_type(self, entry_id: str) -> str | None:
        """'family', 'domain', 'homologous_superfamily', ... (cached)."""
        if entry_id not in self._interpro_type:
            try:
                data = _get_json(f"{INTERPRO_BASE_URL}/{entry_id}")
                self._interpro_type[entry_id] = ((data.get("metadata") or {}).get("type") or "").lower()
            except EMBLAPIError:
                self._interpro_type[entry_id] = None  # unknown: keep it rather than lose evidence
        return self._interpro_type[entry_id]

    def chebi_for_molecules(self, molecule_ids: list[str]) -> dict[str, str | None]:
        """ChEMBL molecule id -> ChEBI id, via UniChem (ChEMBL's own molecule
        records carry no ChEBI cross-reference)."""
        todo = [m for m in dict.fromkeys(molecule_ids) if m and m not in self._chebi]
        if todo and self.resolve_chebi:
            with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
                for mol, chebi in zip(todo, pool.map(self._unichem_chebi, todo)):
                    self._chebi[mol] = chebi
        return {m: self._chebi.get(m) for m in molecule_ids}

    @staticmethod
    def _unichem_chebi(molecule_id: str) -> str | None:
        try:
            resp = requests.post(UNICHEM_URL, timeout=30,
                                 json={"type": "sourceID", "compound": molecule_id, "sourceID": 1})
            resp.raise_for_status()
            data = resp.json()
        except Exception:  # noqa: BLE001 - a missing mapping must not fail the build
            return None
        for compound in data.get("compounds") or []:
            for src in compound.get("sources") or []:
                if (src.get("shortName") or "").lower() == "chebi":
                    return src.get("compoundId")
        return None

    # -- building ------------------------------------------------------------ #

    def build(self, assay_ids: list[str]) -> ChEMBLExtraction:
        out = ChEMBLExtraction()
        for assay_id in assay_ids:
            try:
                self._build_one(assay_id, out)
            except (EMBLAPIError, ChEMBLReferenceError) as exc:
                out.warnings.append(f"{assay_id}: {exc}")
        return out

    def _build_one(self, assay_id: str, out: ChEMBLExtraction) -> None:
        assay = self.fetch_assay(assay_id)
        if not assay.get("assay_chembl_id"):
            raise ChEMBLReferenceError("not an assay record")
        tier = classify_assay(assay)
        activities = self.fetch_activities(assay_id)
        outcomes = self._collapse_activities(activities)
        if not outcomes:
            raise ChEMBLReferenceError(self._why_no_outcomes(activities))
        chebi = self.chebi_for_molecules([o["molecule_chembl_id"] for o in outcomes])
        citation = self._citation(assay)

        if tier == "microbiome":
            out.observations.append(self._community(assay, outcomes, chebi, citation, out))
            return

        protein = None
        if tier == "gene":
            protein = self._target_protein(assay["target_chembl_id"])
            if not protein or not protein.get("accession"):
                out.warnings.append(f"{assay_id}: protein target could not be resolved; recorded as a microbe")
                tier = "microbe"
            elif protein["dropped_ids"]:
                out.warnings.append(
                    f"{assay_id}: ignored over-general InterPro entries {', '.join(protein['dropped_ids'])}")

        for outcome in outcomes:
            out.records.append(self._record(tier, assay, protein, outcome, chebi, citation))

    def _record(self, tier: str, assay: dict, protein: dict | None,
                outcome: dict, chebi: dict, citation: str) -> dict:
        drug = outcome["drug"]
        organism = assay.get("assay_organism")
        strain = assay.get("assay_strain")
        record = {
            "entity_type": "gene" if tier == "gene" else "microbe",
            "reaction_class": f"Biotransformation of {drug.lower()}",
            "observed": outcome["observed"],
            "molecule_chembl_id": outcome["molecule_chembl_id"],
            "organism": organism,
            "tax_id": assay.get("assay_tax_id"),
            "assay_chembl_id": assay.get("assay_chembl_id"),
            "source": f"ChEMBL {assay.get('assay_chembl_id')}{citation}",
        }
        if chebi.get(outcome["molecule_chembl_id"]):
            record["chebi_substrate"] = chebi[outcome["molecule_chembl_id"]]
        if tier == "gene":
            record["entity_id"] = protein["accession"]
            record["entity_name"] = protein["name"]
            record["alt_ids"] = list(protein["alt_ids"])
            subject = f"{protein['name']} ({protein['accession']})"
        else:
            record["entity_id"] = organism
            record["entity_name"] = f"{organism} {strain}".strip() if strain else organism
            subject = record["entity_name"]
        verb = "was observed to transform" if outcome["observed"] else "did not transform"
        detail = f" Metabolites at m/z {', '.join(outcome['metabolites'])}." if outcome["metabolites"] else ""
        record["evidence"] = (f"{subject} {verb} {drug.lower()} in vitro.{detail} "
                              f"ChEMBL assay {assay.get('assay_chembl_id')}{citation}.")
        return record

    def _community(self, assay: dict, outcomes: list[dict], chebi: dict,
                   citation: str, out: ChEMBLExtraction) -> dict:
        """A measured community: which drugs this metagenome did and didn't
        transform, tied where possible to the MGnify analysis of that sample."""
        ena = _ena_accessions(assay)
        community = {
            "assay_chembl_id": assay.get("assay_chembl_id"),
            "community": (assay.get("aidx") or "").split("_")[1] if "_" in (assay.get("aidx") or "") else None,
            "organism": assay.get("assay_organism"),
            "tax_id": assay.get("assay_tax_id"),
            "ena_study": ena.get("study"),
            "ena_sample": ena.get("sample"),
            "mgnify": None,
            "source": f"ChEMBL {assay.get('assay_chembl_id')}{citation}",
            "observations": [{
                "reaction_class": f"Biotransformation of {o['drug'].lower()}",
                "molecule_chembl_id": o["molecule_chembl_id"],
                "chebi_substrate": chebi.get(o["molecule_chembl_id"]),
                "observed": o["observed"],
            } for o in outcomes],
        }
        if ena.get("sample"):
            client = self.mgnify or MGnifyClient()
            match = client.find_analysis_for_sample(ena["sample"])
            if match:
                community["mgnify"] = match
            else:
                out.warnings.append(
                    f"{assay.get('assay_chembl_id')}: no {client.pipeline or 'MGnify'} analysis "
                    f"for ENA sample {ena['sample']}")
        return community

    @staticmethod
    def _why_no_outcomes(activities: list[dict]) -> str:
        """Say *why* an assay yielded nothing, quoting a real comment: these
        screens are published by different groups and each phrases its result
        differently, so an unrecognised wording is the likeliest cause and the
        message needs to carry enough to fix the patterns."""
        if not activities:
            return "no activities returned for this assay"
        sample = next((a.get("activity_comment") or a.get("standard_text_value")
                       for a in activities if a.get("activity_comment") or a.get("standard_text_value")), None)
        types = sorted({a.get("standard_type") for a in activities if a.get("standard_type")})
        if sample is None:
            return (f"{len(activities)} activities, none with a comment or verdict "
                    f"(standard_type: {', '.join(map(str, types)) or 'none'})")
        return (f"{len(activities)} activities, but no recognisable outcome wording; "
                f"example comment: {sample[:120]!r}")

    @staticmethod
    def _collapse_activities(activities: list[dict]) -> list[dict]:
        """One outcome per drug. ChEMBL lists a headline row plus one row per
        putative metabolite, so the rows for a drug are merged: observed if any
        row says so, with the metabolite m/z values kept as detail."""
        merged: dict[str, dict] = {}
        for act in activities:
            molecule = act.get("molecule_chembl_id")
            observed = observation_from_activity(act)
            if not molecule or observed is None:
                continue
            drug = (act.get("molecule_pref_name") or molecule).strip()
            entry = merged.setdefault(molecule, {
                "molecule_chembl_id": molecule, "drug": drug,
                "observed": False, "metabolites": [], "activities": 0})
            entry["activities"] += 1
            entry["observed"] = entry["observed"] or observed
            found = _METABOLITE_RE.search(act.get("activity_comment") or "")
            if found and found.group(1) not in entry["metabolites"]:
                entry["metabolites"].append(found.group(1))
        return sorted(merged.values(), key=lambda e: (not e["observed"], e["drug"]))

    def _citation(self, assay: dict) -> str:
        doc_id = assay.get("document_chembl_id")
        if not doc_id:
            return ""
        try:
            doc = _get_json(f"{self.base_url}/document/{doc_id}.json")
        except EMBLAPIError:
            return ""
        pmid, year = doc.get("pubmed_id"), doc.get("year")
        if pmid:
            return f" (PMID:{pmid}{f', {year}' if year else ''})"
        return f" ({year})" if year else ""


def build_reference_associations(assay_ids: list[str], **kwargs) -> ChEMBLExtraction:
    """Convenience wrapper: `ChEMBLBiotransformationSource(**kwargs).build(assay_ids)`."""
    return ChEMBLBiotransformationSource(**kwargs).build(assay_ids)
