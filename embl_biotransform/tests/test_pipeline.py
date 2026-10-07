"""Tests for the core pipeline and the web API.

Run from the project root:  python -m pytest   (or: python -m unittest discover tests)
"""

from __future__ import annotations

import io
import json
import sys
import threading
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from embl_biotransform import (  # noqa: E402
    EvidenceGraph, BiotransformationPredictor, ReferenceFormatError,
    load_reference_associations, parse_detections_text, parse_reference_text, run_analysis,
    run_multi_analysis, DetectionFormatError,
)
from embl_biotransform.pipeline import build_taxon_grid  # noqa: E402
from embl_biotransform.reference import normalise_reference_records  # noqa: E402
from embl_biotransform.fetchers import EMBLAPIError  # noqa: E402
from embl_biotransform.integrate import clean_taxon_name, ids_match, name_match, taxon_match  # noqa: E402
from embl_biotransform import (  # noqa: E402
    BioSIFTRFormatError, load_biosiftr_run, parse_biosiftr_table, summarise_biosiftr_run,
)

EX = ROOT / "example_data"


def _example():
    return (json.loads((EX / "sample_taxonomy.json").read_text()),
            json.loads((EX / "sample_functions.json").read_text()),
            load_reference_associations(EX / "reference_associations.json"))


class MatchingTests(unittest.TestCase):
    def test_ids_match_exact_and_prefixes(self):
        self.assertTrue(ids_match("EC:1.14.14.1", "1.14.14.1"))
        self.assertTrue(ids_match("ipr001279", "IPR001279"))
        self.assertFalse(ids_match("IPR001279", "IPR001280"))

    def test_ids_match_ec_wildcard(self):
        self.assertTrue(ids_match("EC:1.14.13.-", "1.14.13.39"))
        self.assertTrue(ids_match("1.14.-.-", "EC:1.14.99.1"))
        self.assertFalse(ids_match("EC:1.14.13.-", "1.14.14.1"))

    def test_name_match_ignores_generic_words(self):
        # only generic words shared -> no match
        self.assertEqual(name_match("Haloalkane dehalogenase family protein", "Glutathione transferase family protein"), [])
        # distinctive words shared -> match
        self.assertTrue(name_match("Alkane hydroxylase (AlkB)", "Alkane hydroxylase family (AlkB)"))

    def test_taxon_names(self):
        self.assertEqual(clean_taxon_name("sk__Bacteria;p__Proteobacteria;g__Pseudomonas;s__Pseudomonas_putida"),
                         "Pseudomonas putida")
        self.assertEqual(clean_taxon_name("sk__Bacteria;g__Pseudomonas;s__"), "Pseudomonas")
        self.assertEqual(taxon_match("Pseudomonas putida", "pseudomonas putida"), "exact")
        self.assertEqual(taxon_match("Pseudomonas", "Pseudomonas putida"), "within-genus")
        # a species-level reference must not be matched by a genus-only detection
        self.assertIsNone(taxon_match("Pseudomonas putida", "Pseudomonas"))


class PipelineTests(unittest.TestCase):
    def test_example_predictions(self):
        taxa, functions, refs = _example()
        g = EvidenceGraph()
        s = g.add_sample("demo")
        g.add_detected_taxa(s, taxa)
        g.add_detected_functions(s, functions)
        g.add_reference_associations(refs)
        g.link_detections_to_references()
        preds = {p.reaction_class: p for p in BiotransformationPredictor(g).predict(s)}
        self.assertEqual(set(preds), {"Alkane hydroxylation", "Aromatic ring hydroxylation", "Dehalogenation"})
        self.assertTrue(preds["Alkane hydroxylation"].tier1_hits and preds["Alkane hydroxylation"].tier2_hits)
        self.assertFalse(preds["Dehalogenation"].tier1_hits)
        self.assertGreater(preds["Alkane hydroxylation"].score, preds["Dehalogenation"].score)

    def test_name_matching_can_be_disabled(self):
        taxa, functions, refs = _example()
        r = run_analysis(taxa, functions, refs, allow_name_matching=False)
        self.assertEqual(r["summary"]["tier1_predictions"], 0)

    def test_id_match_scores_higher_than_name_match(self):
        refs = [{"entity_type": "gene", "entity_id": "EC:3.8.1.5", "alt_ids": ["IPR000639"],
                 "entity_name": "Haloalkane dehalogenase", "reaction_class": "Dehalogenation"}]
        by_id = run_analysis([], [{"annotation_id": "IPR000639", "description": "Epoxide hydrolase-like"}], refs)
        by_name = run_analysis([], [{"annotation_id": "IPR999999", "description": "Haloalkane dehalogenase"}], refs)
        self.assertEqual(by_id["predictions"][0]["tier1_hits"][0]["detections"][0]["match_type"], "id")
        self.assertEqual(by_name["predictions"][0]["tier1_hits"][0]["detections"][0]["match_type"], "name")
        self.assertGreater(by_id["predictions"][0]["score"], by_name["predictions"][0]["score"])

    def test_genus_reference_counts_once(self):
        refs = [{"entity_type": "microbe", "entity_id": "Pseudomonas", "reaction_class": "X"}]
        taxa = [{"organism": "Pseudomonas putida", "abundance": 3}, {"organism": "Pseudomonas stutzeri", "abundance": 2}]
        r = run_analysis(taxa, [], refs)
        p = r["predictions"][0]
        self.assertEqual(len(p["tier2_hits"]), 1)
        self.assertEqual(len(p["tier2_hits"][0]["detections"]), 2)

    def test_result_is_json_and_graph_trimmed(self):
        taxa, functions, refs = _example()
        r = run_analysis(taxa, functions, refs)
        json.dumps(r)
        ids = {n["id"] for n in r["graph"]["nodes"]}
        self.assertNotIn("taxon:Bacillus subtilis", ids)  # detected but not linked -> not drawn
        self.assertIn("taxon:Pseudomonas putida", ids)


class ParsingTests(unittest.TestCase):
    def test_reference_csv(self):
        text = ("entity_type,entity_id,reaction_class,alt_ids\n"
                "gene,EC:1.1.1.1,Oxidation,IPR1;K2\n"
                "microbe,Pseudomonas putida,Oxidation,\n")
        recs = parse_reference_text(text, "csv")
        self.assertEqual(recs[0]["alt_ids"], ["IPR1", "K2"])
        self.assertEqual(recs[1]["alt_ids"], [])
        self.assertEqual(recs[1]["entity_name"], "Pseudomonas putida")

    def test_reference_errors(self):
        with self.assertRaises(ReferenceFormatError):
            parse_reference_text('[{"entity_type": "virus", "entity_id": "x", "reaction_class": "y"}]', "json")
        with self.assertRaises(ReferenceFormatError):
            parse_reference_text('[{"entity_id": "x"}]', "json")

    def test_mgnify_raw_json(self):
        raw = {"data": [{"id": "IPR000001", "type": "interpro-identifiers",
                         "attributes": {"accession": "IPR000001", "description": "Kringle", "count": 12}}]}
        recs = parse_detections_text(json.dumps(raw), "json", "functions")
        self.assertEqual(recs[0]["annotation_id"], "IPR000001")
        self.assertEqual(recs[0]["abundance"], 12)
        raw_tax = {"data": [{"id": "1", "attributes": {"lineage": "sk__Bacteria;g__Pseudomonas;s__Pseudomonas_putida",
                                                      "count": 5, "rank": "species"}}]}
        recs = parse_detections_text(json.dumps(raw_tax), "json", "taxa")
        r = run_analysis(recs, [], [{"entity_type": "microbe", "entity_id": "Pseudomonas putida", "reaction_class": "X"}])
        self.assertEqual(r["summary"]["predictions"], 1)

    def test_mgnify_v2_json(self):
        raw = {"count": 2, "items": [
            {"count": 7, "description": None, "organism": "sk__Bacteria;g__Pseudomonas;s__Pseudomonas_putida"},
            {"count": 1, "description": None,
             "organism": "Bacteria::Proteobacteria:Gammaproteobacteria:Pseudomonadales:Pseudomonadaceae:Pseudomonas|5.0"}]}
        recs = parse_detections_text(json.dumps(raw), "json", "taxa")
        self.assertEqual([r["abundance"] for r in recs], [7, 1])
        self.assertEqual(clean_taxon_name(recs[0]["organism"]), "Pseudomonas putida")
        self.assertEqual(clean_taxon_name(recs[1]["organism"]), "Pseudomonas")
        funcs = parse_detections_text(json.dumps({"count": 1, "items": [
            {"count": 3, "description": "ABC transporter", "organism": None}]}), "json", "functions")
        self.assertEqual((funcs[0]["description"], funcs[0]["abundance"]), ("ABC transporter", 3))

    def test_normalised_roundtrip(self):
        taxa, functions, _ = _example()
        self.assertEqual(parse_detections_text(json.dumps(taxa), "json", "taxa")[0]["organism"], "Pseudomonas putida")
        self.assertEqual(parse_detections_text(json.dumps(functions), "json", "functions")[0]["annotation_id"], "IPR001279")

    def test_tsv_table(self):
        text = "accession\tdescription\tcount\nIPR1\tSomething\t4\n"
        recs = parse_detections_text(text, "tsv", "functions")
        self.assertEqual((recs[0]["annotation_id"], recs[0]["abundance"]), ("IPR1", 4))


class _FakeMGnify:
    def search_studies(self, q):
        return [{"accession": "MGYS1", "name": f"study about {q}"}]

    def list_analyses(self, acc):
        return [{"accession": "MGYA1"}]

    def fetch_taxonomy(self, acc, source="auto"):
        return [{"organism": "Pseudomonas putida", "abundance": 1, "source": "MGnify", "evidence": "x"}]

    def fetch_functional_annotations(self, acc, kind="interpro"):
        raise EMBLAPIError("boom")


class _FakeChEMBL:
    def search_assays(self, keyword, limit=50):
        return [{"assay_chembl_id": "CHEMBL1", "assay_organism": f"bug {keyword}",
                 "description": "d", "target_chembl_id": "CHEMBLT", "activities": 3, "tier": "gene"}]

    def build(self, assay_ids):
        from embl_biotransform.chembl_references import ChEMBLExtraction
        if assay_ids == ["EMPTY"]:
            return ChEMBLExtraction(warnings=["EMPTY: no biotransformation activities"])
        return ChEMBLExtraction(
            records=[{"entity_type": "gene", "entity_id": "Q1", "reaction_class": "Biotransformation of d",
                      "observed": True, "alt_ids": ["IPR001279"]},
                     {"entity_type": "gene", "entity_id": "Q1", "reaction_class": "Biotransformation of e",
                      "observed": False, "alt_ids": ["IPR001279"]}],
            observations=[{"assay_chembl_id": "CHEMBL3", "mgnify": {"analysis": "MGYA1"},
                           "observations": [{"reaction_class": "Biotransformation of d", "observed": True}]}],
            warnings=["ignored over-general InterPro entries IPR029058 (homologous_superfamily)"])


class WebAPITests(unittest.TestCase):
    def setUp(self):
        # a fresh app per test: the ChEMBL and MGnify libraries live in the
        # app's cache, so a shared one would leak entries between tests
        try:
            from webapp.app import create_app
        except ImportError as exc:  # flask not installed
            raise unittest.SkipTest(f"web app deps missing: {exc}")
        self.client = create_app(mgnify_client=_FakeMGnify(),
                                 chembl_source=_FakeChEMBL()).test_client()

    def _await_job(self, response, timeout=10.0):
        """ChEMBL and batch-MGnify endpoints return 202 + a job id; poll it."""
        self.assertEqual(response.status_code, 202, response.json)
        job_id = response.json["id"]
        deadline = time.time() + timeout
        while time.time() < deadline:
            job = self.client.get(f"/api/jobs/{job_id}").json
            if job["state"] in ("done", "failed", "cancelled"):
                return job
            time.sleep(0.02)
        self.fail(f"job {job_id} did not finish within {timeout}s")

    def test_index_and_health(self):
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.client.get("/static/app.js").status_code, 200)
        self.assertEqual(self.client.get("/api/health").json["status"], "ok")

    def test_analyze_example(self):
        ex = self.client.get("/api/examples").json
        r = self.client.post("/api/analyze", json=ex)
        self.assertEqual(r.status_code, 200)
        self.assertEqual(r.json["summary"]["predictions"], 3)

    def test_analyze_errors(self):
        self.assertEqual(self.client.post("/api/analyze", json={"taxa": [{"organism": "x"}]}).status_code, 400)
        r = self.client.post("/api/analyze", json={"taxa": [{"organism": "x"}],
                                                   "references": [{"entity_type": "bad", "entity_id": "x", "reaction_class": "y"}]})
        self.assertEqual(r.status_code, 400)
        self.assertIn("entity_type", r.json["error"])

    def test_large_uploads_are_accepted(self):
        """No upload cap: a fully built ChEMBL reference set is bigger than any
        round number worth guessing at, and this reads local files."""
        self.assertIsNone(self.client.application.config.get("MAX_CONTENT_LENGTH"))
        rows = b"entity_type,entity_id,reaction_class\n" + \
               b"".join(b"gene,IPR%06d,Biotransformation of drug %d\n" % (i, i)
                        for i in range(60_000))
        self.assertGreater(len(rows), 2_000_000)
        r = self.client.post("/api/parse/reference",
                             data={"file": (io.BytesIO(rows), "big_refs.csv")},
                             content_type="multipart/form-data")
        self.assertEqual(r.status_code, 200, str(r.json)[:200])
        self.assertEqual(len(r.json["records"]), 60_000)

    def test_upload_parsing(self):
        data = {"file": (io.BytesIO(b"entity_type,entity_id,reaction_class\ngene,EC:1.1.1.1,Ox\n"), "refs.csv")}
        r = self.client.post("/api/parse/reference", data=data, content_type="multipart/form-data")
        self.assertEqual(r.status_code, 200, r.json)
        self.assertEqual(len(r.json["records"]), 1)
        data = {"file": (io.BytesIO(b"organism,count\nPseudomonas putida,4\n"), "t.csv"), "kind": "taxa"}
        r = self.client.post("/api/parse/detections", data=data, content_type="multipart/form-data")
        self.assertEqual(r.json["records"][0]["abundance"], 4)

    def test_analyze_multi(self):
        ex = self.client.get("/api/examples").json
        r = self.client.post("/api/analyze/multi", json={
            "samples": [{"sample_id": "A", "taxa": ex["taxa"], "functions": ex["functions"]},
                        {"sample_id": "B", "taxa": ex["taxa"], "functions": []}],
            "references": ex["references"]})
        self.assertEqual(r.status_code, 200, r.json)
        self.assertEqual(len(r.json["samples"]), 2)
        self.assertEqual(r.json["grid"]["sample_ids"], ["A", "B"])
        self.assertTrue(r.json["grid"]["rows"])

    def test_analyze_multi_errors(self):
        ex = self.client.get("/api/examples").json
        self.assertEqual(self.client.post("/api/analyze/multi",
                         json={"references": ex["references"]}).status_code, 400)
        self.assertEqual(self.client.post("/api/analyze/multi",
                         json={"samples": [{"sample_id": "A"}], "references": ex["references"]}).status_code, 400)
        self.assertEqual(self.client.post("/api/analyze/multi",
                         json={"samples": [{"sample_id": "A", "taxa": ex["taxa"]}]}).status_code, 400)

    def test_chembl_routes(self):
        job = self._await_job(self.client.get("/api/chembl/assays?search=Bacteroides"))
        self.assertEqual(job["state"], "done", job)
        self.assertEqual(job["result"]["assays"][0]["assay_chembl_id"], "CHEMBL1")
        self.assertEqual(self.client.get("/api/chembl/assays").status_code, 400)

        job = self._await_job(self.client.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}))
        self.assertEqual(job["state"], "done", job)
        out = job["result"]
        self.assertEqual(out["summary"]["positive"], 1)
        self.assertEqual(out["summary"]["negative"], 1)
        self.assertEqual(len(out["observations"]), 1)
        self.assertTrue(out["warnings"])

        # validation still rejects before a job is created
        self.assertEqual(self.client.post("/api/chembl/references", json={}).status_code, 400)
        # an assay with nothing usable fails the job rather than the request
        job = self._await_job(self.client.post("/api/chembl/references", json={"assay_ids": ["EMPTY"]}))
        self.assertEqual(job["state"], "failed")
        self.assertIn("no biotransformation activities", job["error"])

    def test_large_assay_lists_are_accepted(self):
        """No cap: building a large reference set is the point, and it is
        paced and cached rather than refused."""
        r = self.client.post("/api/chembl/references",
                             json={"assay_ids": ["CHEMBL1"] * 3 + ["EMPTY"] * 60})
        self.assertEqual(r.status_code, 202)
        job = self._await_job(r, timeout=20)
        # duplicates collapse, so this is 2 distinct assays, not 63
        self.assertEqual(len(job["meta"]["assay_ids"]), 2)

    def test_search_flags_assays_already_built(self):
        found = self._await_job(self.client.get("/api/chembl/assays?search=x"))["result"]
        self.assertFalse(found["assays"][0]["in_library"])
        self.assertEqual(found["in_library"], 0)
        self._await_job(self.client.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}))
        # the search result is cached, so ask for a different term to re-flag
        found = self._await_job(self.client.get("/api/chembl/assays?search=y"))["result"]
        self.assertTrue(found["assays"][0]["in_library"])
        self.assertEqual(found["in_library"], 1)

    def test_job_listing_and_cancel(self):
        job = self._await_job(self.client.get("/api/chembl/assays?search=x"))
        listed = self.client.get("/api/jobs").json["jobs"]
        self.assertIn(job["id"], [j["id"] for j in listed])
        self.assertNotIn("result", listed[0])          # results only on the detail route
        self.assertEqual(self.client.get("/api/jobs/nosuch").status_code, 404)
        self.assertEqual(self.client.post(f"/api/jobs/{job['id']}/cancel").status_code, 409)

    def test_mgnify_batch_job(self):
        r = self.client.post("/api/mgnify/analyses", json={"accessions": ["MGYA1", "MGYA2"]})
        job = self._await_job(r)
        self.assertEqual(job["state"], "done", job)
        # the fake client serves taxonomy but fails functions, so both load partially
        self.assertEqual(len(job["result"]["analyses"]), 2)
        self.assertEqual(self.client.post("/api/mgnify/analyses", json={}).status_code, 400)
        self.assertEqual(self.client.post("/api/mgnify/analyses",
                         json={"accessions": ["X"] * 501}).status_code, 400)

    def test_cache_endpoint_and_clear(self):
        self._await_job(self.client.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}))
        stats = self.client.get("/api/cache").json
        self.assertIn("entries", stats)
        self.assertTrue(stats["entries"])
        self.assertTrue(self.client.get("/api/chembl/library").json["entries"])

        removed = self.client.delete("/api/cache").json["removed"]
        self.assertTrue(removed)
        self.assertEqual(self.client.get("/api/cache").json["entries"], 0)
        self.assertEqual(self.client.get("/api/chembl/library").json["entries"], [])

    def test_chembl_references_feed_the_grid(self):
        """The positives build a grid, and the community assay marks its column."""
        built = self._await_job(
            self.client.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}))["result"]
        refs = [{k: v for k, v in rec.items() if k != "observed"}
                for rec in built["records"] if rec["observed"]]
        r = self.client.post("/api/analyze/multi", json={
            "samples": [{"sample_id": "MGYA1",
                         "functions": [{"annotation_id": "IPR001279", "description": "x", "abundance": 2}]}],
            "references": refs, "observations": built["observations"]})
        self.assertEqual(r.status_code, 200, r.json)
        grid = r.json["grid"]
        self.assertEqual(grid["observed_samples"], ["MGYA1"])
        row = next(x for x in grid["rows"] if x["reaction_class"] == "Biotransformation of d")
        self.assertIs(row["cells"][0]["observed"], True)
        self.assertEqual(row["cells"][0]["tier"], 1)

    def test_chembl_library_persists_and_toggles(self):
        """Built assays stay in the library, and their selection is remembered
        so a set of references can be toggled on and off between runs."""
        self.assertEqual(self.client.get("/api/chembl/library").json["entries"], [])
        self._await_job(self.client.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}))

        lib = self.client.get("/api/chembl/library").json
        self.assertEqual([e["assay_chembl_id"] for e in lib["entries"]], ["CHEMBL1"])
        entry = lib["entries"][0]
        self.assertEqual((entry["positive"], entry["negative"]), (1, 1))
        self.assertEqual(entry["tier"], "gene")
        self.assertTrue(entry["selected"])          # building selects it

        # references for the current selection, without rebuilding
        refs = self.client.get("/api/chembl/library/references").json
        self.assertEqual(refs["assay_ids"], ["CHEMBL1"])
        self.assertEqual(refs["summary"]["positive"], 1)

        # toggle off
        r = self.client.post("/api/chembl/library/selection", json={"assay_ids": []})
        self.assertEqual(r.json["selected"], [])
        self.assertFalse(self.client.get("/api/chembl/library").json["entries"][0]["selected"])
        self.assertEqual(self.client.get("/api/chembl/library/references").json["assay_ids"], [])

        # and back on
        self.client.post("/api/chembl/library/selection", json={"assay_ids": ["CHEMBL1"]})
        self.assertTrue(self.client.get("/api/chembl/library").json["entries"][0]["selected"])

        self.assertEqual(self.client.post("/api/chembl/library/selection",
                                          json={"assay_ids": ["NOPE"]}).status_code, 400)
        self.assertEqual(self.client.delete("/api/chembl/library/CHEMBL1").status_code, 200)
        self.assertEqual(self.client.get("/api/chembl/library").json["entries"], [])
        self.assertEqual(self.client.delete("/api/chembl/library/CHEMBL1").status_code, 404)

    def test_library_survives_a_restart(self):
        """The point of persisting it: a new app on the same cache dir sees it."""
        import tempfile
        from embl_biotransform.cache import Cache
        from webapp.app import create_app
        path = Path(tempfile.mkdtemp()) / "c.sqlite"
        cache = Cache(path)
        self.addCleanup(cache.close)
        first = create_app(mgnify_client=_FakeMGnify(), chembl_source=_FakeChEMBL(),
                           cache=cache).test_client()
        job_id = first.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}).json["id"]
        deadline = time.time() + 10
        while time.time() < deadline:
            if first.get(f"/api/jobs/{job_id}").json["state"] in ("done", "failed"):
                break
            time.sleep(0.02)

        reopened = Cache(path)
        self.addCleanup(reopened.close)
        second = create_app(mgnify_client=_FakeMGnify(), chembl_source=_FakeChEMBL(),
                            cache=reopened).test_client()
        lib = second.get("/api/chembl/library").json
        self.assertEqual([e["assay_chembl_id"] for e in lib["entries"]], ["CHEMBL1"])
        self.assertTrue(lib["persistent"])
        self.assertEqual(second.get("/api/chembl/library/references").json["summary"]["positive"], 1)

    def test_analysis_list_flags_downloaded(self):
        listed = self.client.get("/api/mgnify/studies/MGYS1/analyses").json
        self.assertFalse(listed["analyses"][0]["downloaded"])
        self.assertEqual(listed["downloaded"], 0)
        self._await_job(self.client.post("/api/mgnify/analyses", json={"accessions": ["MGYA1"]}))
        listed = self.client.get("/api/mgnify/studies/MGYS1/analyses").json
        by_id = {a["accession"]: a for a in listed["analyses"]}
        self.assertTrue(by_id["MGYA1"]["downloaded"])
        self.assertFalse(by_id["MGYA1"]["enzyme_taxonomy"])
        self.assertEqual(listed["downloaded"], 1)

    def test_mgnify_library(self):
        """Downloaded analyses are listed back from the cache, by accession."""
        self.assertEqual(self.client.get("/api/mgnify/library").json["entries"], [])
        self._await_job(self.client.post("/api/mgnify/analyses", json={"accessions": ["MGYA1"]}))
        entries = self.client.get("/api/mgnify/library").json["entries"]
        self.assertEqual([e["accession"] for e in entries], ["MGYA1"])
        self.assertEqual(entries[0]["taxonomy"], "auto")
        self.assertEqual(entries[0]["taxa"], 1)
        self.assertEqual(self.client.delete("/api/mgnify/library/MGYA1").json["removed"], "MGYA1")
        self.assertEqual(self.client.get("/api/mgnify/library").json["entries"], [])
        self.assertEqual(self.client.delete("/api/mgnify/library/MGYA1").status_code, 404)

    def test_mgnify_library_ignores_legacy_hashed_keys(self):
        """Analyses cached before they were keyed by accession are hashes;
        they are still usable cache, but must not show as garbage rows."""
        from webapp.app import MGNIFY_LIBRARY_NS
        self._await_job(self.client.post("/api/mgnify/analyses", json={"accessions": ["MGYA1"]}))
        store = self.client.application.config["EXPOSE_STORE"]
        store.set(MGNIFY_LIBRARY_NS, "a" * 64, {"taxa": [1], "functions": []})
        entries = self.client.get("/api/mgnify/library").json["entries"]
        self.assertEqual([e["accession"] for e in entries], ["MGYA1"])

    def test_chembl_id_list_matching(self):
        """Upload a list of ids: say what is held and what must be fetched."""
        self._await_job(self.client.post("/api/chembl/references", json={"assay_ids": ["CHEMBL1"]}))
        data = {"file": (io.BytesIO(b"CHEMBL1\nCHEMBL2\n# a comment\nCHEMBL1\n"), "ids.txt")}
        r = self.client.post("/api/chembl/library/match", data=data,
                             content_type="multipart/form-data")
        self.assertEqual(r.json["present"], ["CHEMBL1"])
        self.assertEqual(r.json["missing"], ["CHEMBL2"])
        self.assertEqual(r.json["supplied"], ["CHEMBL1", "CHEMBL2"])     # de-duplicated

        # ids embedded in a CSV column are found too
        data = {"file": (io.BytesIO(b"assay,note\nCHEMBL2,x\nCHEMBL3,y\n"), "ids.csv")}
        r = self.client.post("/api/chembl/library/match", data=data,
                             content_type="multipart/form-data")
        self.assertEqual(r.json["missing"], ["CHEMBL2", "CHEMBL3"])

        r = self.client.post("/api/chembl/library/match", json={"assay_ids": ["chembl1"]})
        self.assertEqual(r.json["present"], ["CHEMBL1"])                 # case-insensitive
        self.assertEqual(self.client.post("/api/chembl/library/match",
                                          json={"assay_ids": []}).status_code, 400)
        data = {"file": (io.BytesIO(b"nothing useful here"), "empty.txt")}
        self.assertEqual(self.client.post("/api/chembl/library/match", data=data,
                                          content_type="multipart/form-data").status_code, 400)

    def test_analyze_multi_honours_rank_floor(self):
        r = self.client.post("/api/analyze/multi", json={
            "samples": [{"sample_id": "A", "taxa": [
                {"organism": "d__Bacteria;p__Bacillota", "lineage": "d__Bacteria;p__Bacillota"},
                {"organism": "d__Bacteria;g__Blautia", "lineage": "d__Bacteria;g__Blautia"}]}],
            "references": [{"entity_type": "microbe", "entity_id": "Blautia",
                            "reaction_class": "Biotransformation of x"}],
            "taxon_rank_floor": "genus"})
        self.assertEqual(r.status_code, 200, r.json)
        self.assertEqual(r.json["summary"]["taxa_dropped_by_rank"], 1)
        self.assertEqual(r.json["settings"]["taxon_rank_floor"], "genus")
        self.assertIn("enzyme_grid", r.json)

    def test_mgnify_routes(self):
        self.assertEqual(self.client.get("/api/mgnify/studies?search=soil").json["studies"][0]["accession"], "MGYS1")
        self.assertEqual(self.client.get("/api/mgnify/studies/MGYS1/analyses").json["analyses"][0]["accession"], "MGYA1")
        r = self.client.get("/api/mgnify/analyses/MGYA1").json
        self.assertEqual(len(r["taxa"]), 1)
        self.assertEqual(r["functions"], [])
        self.assertTrue(r["warnings"])


class MGnifyV2ClientTests(unittest.TestCase):
    """The v2 client against canned responses (no network)."""

    DETAIL = {"accession": "MGYA1", "experiment_type": "Assembly", "pipeline_version": "V6", "downloads": [
        {"alias": "X_interpro_summary.tsv.gz", "url": "https://ftp/X_interpro_summary.tsv.gz"},
        {"alias": "X_proteins2rhea.tsv.gz", "url": "https://ftp/X_proteins2rhea.tsv.gz"},
        {"alias": "X.krona.txt.gz", "url": "https://ftp/X.krona.txt.gz"},
        {"alias": "X_ko_summary.tsv.gz", "url": "https://ftp/X_ko_summary.tsv.gz"}]}
    OLD = {"accession": "MGYA2", "experiment_type": "Assembly", "pipeline_version": "V5", "downloads": []}
    FILES = {
        "https://ftp/X_interpro_summary.tsv.gz":
            "interpro_accession\tdescription\tcount\nIPR001128\tCytochrome P450\t42\n",
        "https://ftp/X_proteins2rhea.tsv.gz":
            "contig_id\tprotein_id\tprotein_hash\trhea_id\tchebi_reaction\treaction\ttop_hit\n"
            "c1\tp1\th\tRHEA:10000\tx\tA = B\tTrue\nc1\tp2\th\tRHEA:10000\tx\tA = B\tTrue\n"
            "c1\tp2\th\tRHEA:10000\tx\tA = B\tFalse\n",
        "https://ftp/X.krona.txt.gz":
            "100\tunclassified\n12\td__Bacteria\tp__Pseudomonadota\tg__Pseudomonas\n",
        "https://ftp/X_ko_summary.tsv.gz": "ko\tdescription\tcount\nK03406\tmethyl-accepting chemotaxis protein\t114\n",
    }

    def setUp(self):
        from unittest import mock
        from embl_biotransform import fetchers
        self.client = fetchers.MGnifyClient()
        self.calls = []

        def fake_json(url, params=None, **_):
            self.calls.append((url, params))
            if url.endswith("/analyses/MGYA1"):
                return self.DETAIL
            if url.endswith("/analyses/MGYA2"):
                return self.OLD
            if url.endswith("/studies/"):
                return {"count": 2, "items": [{"accession": "MGYS2", "title": "old"}, {"accession": "MGYS1", "title": "v6"}]}
            if url.endswith("/studies/MGYS2/analyses/"):
                return {"count": 1, "items": [{"accession": "MGYA2", "pipeline_version": "V5"}]}
            if "/annotations/taxonomies__" in url:
                return {"count": 0, "items": []}
            if url.endswith("/studies/MGYS1/analyses/"):
                page, size = params["page"], params["page_size"]
                total = 130
                return {"count": total, "items": [{"accession": f"MGYA{i}", "pipeline_version": "V5" if i == 129 else "V6"}
                                                  for i in range((page - 1) * size, min(page * size, total))]}
            raise fetchers.EMBLAPIError(f"unexpected {url}")

        patches = [mock.patch.object(fetchers, "_get_json", fake_json),
                   mock.patch.object(fetchers, "_get_text_file", lambda url, **_: self.FILES[url])]
        for p in patches:
            p.start()
            self.addCleanup(p.stop)

    def test_paging(self):
        analyses = self.client.list_analyses("MGYS1")
        self.assertEqual(len(analyses), 129)  # the one V5 analysis is dropped
        self.assertEqual([c[1]["page"] for c in self.calls], [1, 2])

    def test_only_pipeline_v6(self):
        self.assertEqual([s["accession"] for s in self.client.search_studies("x")], ["MGYS1"])
        for fetch in (self.client.fetch_taxonomy, self.client.fetch_functional_annotations):
            with self.assertRaisesRegex(EMBLAPIError, "only V6"):
                fetch("MGYA2")

    def test_legacy_summary_layout(self):
        from embl_biotransform.fetchers import _parse_summary_table
        rows = _parse_summary_table('"114","K03406","methyl-accepting chemotaxis protein"\n')
        self.assertEqual(rows, [{"annotation_id": "K03406", "description": "methyl-accepting chemotaxis protein",
                                 "abundance": 114}])

    def test_functional_summary_tables(self):
        recs = self.client.fetch_functional_annotations("MGYA1", "interpro")
        self.assertEqual((recs[0]["annotation_id"], recs[0]["abundance"]), ("IPR001128", 42))
        self.assertIn("MGYA1", recs[0]["evidence"])
        ko = self.client.fetch_functional_annotations("MGYA1", "kegg-orthologs")  # v1 name still accepted
        self.assertEqual((ko[0]["annotation_id"], ko[0]["description"], ko[0]["abundance"]),
                         ("K03406", "methyl-accepting chemotaxis protein", 114))
        rhea = self.client.fetch_functional_annotations("MGYA1", "rhea")
        self.assertEqual((rhea[0]["annotation_id"], rhea[0]["abundance"]), ("RHEA:10000", 2))
        with self.assertRaises(EMBLAPIError):
            self.client.fetch_functional_annotations("MGYA1", "go-slim")
        with self.assertRaises(ValueError):
            self.client.fetch_functional_annotations("MGYA1", "nonsense")

    def test_taxonomy_falls_back_to_contigs(self):
        taxa = self.client.fetch_taxonomy("MGYA1")
        self.assertEqual(len(taxa), 1)
        self.assertEqual(clean_taxon_name(taxa[0]["organism"]), "Pseudomonas")
        self.assertIn("contig taxonomy", taxa[0]["evidence"])
        with self.assertRaises(EMBLAPIError):
            self.client.fetch_taxonomy("MGYA1", "ssu")


class ChEMBLReferenceTests(unittest.TestCase):
    """Building reference associations from canned ChEMBL assay records."""

    GENE_ASSAY = {"assay_chembl_id": "CHEMBL5303650", "confidence_score": 9,
                  "target_chembl_id": "CHEMBL5303331", "assay_organism": "Bacteroides thetaiotaomicron",
                  "assay_strain": "VPI-5482", "assay_tax_id": 818, "document_chembl_id": "DOC1",
                  "aidx": "Zimmermann_VPI-5482_Bacteroides thetaiotaomicron_biotransformation_Q8ABF8"}
    MICROBE_ASSAY = {"assay_chembl_id": "CHEMBL5724809", "confidence_score": 0,
                     "target_chembl_id": "CHEMBL612558", "assay_organism": "Bifidobacterium adolescentis",
                     "assay_strain": "DSM20083", "assay_tax_id": 1680, "document_chembl_id": "DOC1",
                     "aidx": "BaghaiArassi_DSM20083_Bifidobacterium adolescentis_biotransformation_anaerobic"}
    COMMUNITY_ASSAY = {"assay_chembl_id": "CHEMBL5725002", "confidence_score": 0,
                       "target_chembl_id": "CHEMBL612558", "assay_organism": "human gut metagenome",
                       "assay_tax_id": 408170, "document_chembl_id": "DOC1",
                       "aidx": "Mastrorilli_WLS-007_biotransformation",
                       "assay_parameters": [
                           {"type": "ENA study accession number", "text_value": "PRJEB37062"},
                           {"type": "ENA sample accession number", "text_value": "ERS12561742"}]}
    TARGET = {"pref_name": "Acetyl esterase (Acetylxylosidase)", "organism": "Bacteroides thetaiotaomicron",
              "tax_id": 818, "target_components": [{
                  "accession": "Q8ABF8", "component_description": "Acetyl esterase (Acetylxylosidase)",
                  "target_component_xrefs": [
                      {"xref_src_db": "InterPro", "xref_id": "IPR000801"},
                      {"xref_src_db": "InterPro", "xref_id": "IPR029058"},
                      {"xref_src_db": "Pfam", "xref_id": "PF00756"},
                      {"xref_src_db": "AlphaFoldDB", "xref_id": "Q8ABF8"},
                      {"xref_src_db": "GoFunction", "xref_id": "GO:0016747"}]}]}
    INTERPRO = {"IPR000801": "family", "IPR029058": "homologous_superfamily"}
    ACTIVITIES = {
        "CHEMBL5303650": [
            {"molecule_chembl_id": "CHEMBL942", "molecule_pref_name": "BISACODYL",
             "activity_comment": "The Oberved Biotransformation Of Drug Bisacodyl Was Proven To Be Mediated By Bt_0152."},
            {"molecule_chembl_id": "CHEMBL942", "molecule_pref_name": "BISACODYL",
             "activity_comment": "Putative Metabolite Identification: The Observed Biotransformation Of Drug "
                                 "Bisacodyl Into Metabolite With M/Z 183.0685 Was Proven To Be Mediated By Bt_0152."},
            {"molecule_chembl_id": "CHEMBL269671", "molecule_pref_name": "ARTEMISININ",
             "activity_comment": "The Biotransformation Of Drug Artemisinin Could Not Be Mediated By Bt_0152."},
            {"molecule_chembl_id": "CHEMBL9999", "molecule_pref_name": "NOISE", "activity_comment": None}],
        "CHEMBL5724809": [
            {"molecule_chembl_id": "CHEMBL1542", "molecule_pref_name": "AZATHIOPRINE",
             "activity_comment": "Biotransformation Occurred For Drug Azathioprine, For The Bacteria B. Adolescentis."},
            {"molecule_chembl_id": "CHEMBL24", "molecule_pref_name": "ATENOLOL",
             "activity_comment": "No Biotransformation Occurred For Drug Atenolol, For The Bacteria B. Adolescentis."}],
        "CHEMBL5725002": [
            {"molecule_chembl_id": "CHEMBL421", "molecule_pref_name": "SULFASALAZINE",
             "activity_comment": "Biotransformation Occurred For Drug Sulfasalazine, Within Community, Wls-007."},
            {"molecule_chembl_id": "CHEMBL30", "molecule_pref_name": "CIMETIDINE",
             "activity_comment": "No Biotransformation Occurred For Drug Cimetidine, Within Community, Wls-007."}],
    }
    ASSAYS = {"CHEMBL5303650": GENE_ASSAY, "CHEMBL5724809": MICROBE_ASSAY, "CHEMBL5725002": COMMUNITY_ASSAY}

    class _FakeMGnify:
        pipeline = "V6"

        def find_analysis_for_sample(self, accession):
            return ({"analysis": "MGYA01030891", "study": "MGYS00010462", "assembly": "ERZ29587455",
                     "sample": "SAMEA110463714", "pipeline_version": "V6"}
                    if accession == "ERS12561742" else None)

    def setUp(self):
        from unittest import mock
        from embl_biotransform import chembl_references as cr

        def fake_json(url, params=None, **_):
            if "/assay/" in url:
                acc = url.rsplit("/", 1)[-1].replace(".json", "")
                if acc not in self.ASSAYS:
                    raise cr.EMBLAPIError(f"GET {url} -> 404: not found")
                return self.ASSAYS[acc]
            if "/target/" in url:
                return self.TARGET
            if "/document/" in url:
                return {"pubmed_id": 31158845, "year": 2019}
            if "/interpro/" in url:
                return {"metadata": {"type": self.INTERPRO[url.rsplit("/", 1)[-1]]}}
            if "activity.json" in url:
                return {"activities": self.ACTIVITIES[params["assay_chembl_id"]], "page_meta": {"next": None}}
            raise AssertionError(f"unexpected {url}")

        patches = [mock.patch.object(cr, "_get_json", fake_json),
                   mock.patch.object(cr.ChEMBLBiotransformationSource, "_unichem_chebi",
                                     staticmethod(lambda m: {"CHEMBL942": "CHEBI:3125",
                                                             "CHEMBL1542": "CHEBI:2948",
                                                             "CHEMBL421": "CHEBI:9334"}.get(m)))]
        for p in patches:
            p.start()
            self.addCleanup(p.stop)
        self.source = cr.ChEMBLBiotransformationSource(mgnify_client=self._FakeMGnify())

    def test_classification(self):
        from embl_biotransform import classify_assay
        self.assertEqual(classify_assay(self.GENE_ASSAY), "gene")
        self.assertEqual(classify_assay(self.MICROBE_ASSAY), "microbe")
        self.assertEqual(classify_assay(self.COMMUNITY_ASSAY), "microbiome")

    def test_observation_from_comment(self):
        from embl_biotransform.chembl_references import observation_from_comment
        self.assertTrue(observation_from_comment("Biotransformation Occurred For Drug X."))
        self.assertTrue(observation_from_comment("The Oberved Biotransformation ... Was Proven To Be Mediated"))
        # "No Biotransformation Occurred" contains "biotransformation occurred"
        self.assertFalse(observation_from_comment("No Biotransformation Occurred For Drug X."))
        self.assertFalse(observation_from_comment("The Biotransformation Of Drug X Could Not Be Mediated."))
        self.assertIsNone(observation_from_comment(None))

    def test_observation_across_screens(self):
        """Each published screen words its result differently; the curated
        `standard_text_value` wins over the free-text comment."""
        from embl_biotransform.chembl_references import observation_from_activity as read
        # Zimmermann gene screen (CHEMBL5303650)
        self.assertTrue(read({"activity_comment":
            "The Oberved Biotransformation Of Drug Bisacodyl Was Proven To Be Mediated By Bt_0152."}))
        # Baghai Arassi strain screen (CHEMBL5724809)
        self.assertFalse(read({"activity_comment":
            "No Biotransformation Occurred For Drug Atenolol, For The Bacteria B. Adolescentis."}))
        # Bacteroides dorei screen (CHEMBL5303585): different wording again, plus a verdict field
        self.assertTrue(read({"standard_text_value": "Compound metabolized", "activity_comment":
            "Putative Metabolite Identification: The Drug  Artemisinin Is Biotransformed By Bacteroides "
            "Dorei Dsm17855 And A Corresponding Drug Metabolite With M/Z 108.094 Is Detected."}))
        self.assertTrue(read({"activity_comment":
            "The Drug Danazol Is Biotransformed By Bacteroides Dorei Dsm17855."}))
        # "not metabolized" contains "metabolized", so negatives must win
        self.assertFalse(read({"standard_text_value": "Compound not metabolized"}))
        self.assertIsNone(read({"activity_comment": "Culture prepared anaerobically."}))

    def test_real_dorei_activities_build(self):
        """The captured CHEMBL5303585 rows (a screen whose wording the first
        patterns missed) now collapse into outcomes with metabolite detail."""
        from embl_biotransform.chembl_references import ChEMBLBiotransformationSource as S
        rows = [{"molecule_chembl_id": "CHEMBL269671", "molecule_pref_name": "ARTEMISININ",
                 "standard_text_value": "Compound metabolized", "activity_comment":
                 "Putative Metabolite Identification: The Drug  Artemisinin Is Biotransformed By Bacteroides "
                 "Dorei Dsm17855 And A Corresponding Drug Metabolite With M/Z 108.094 At Retention Time "
                 "3.542 Minutes Is Detected."},
                {"molecule_chembl_id": "CHEMBL269671", "molecule_pref_name": "ARTEMISININ",
                 "standard_text_value": "Compound metabolized", "activity_comment":
                 "Putative Metabolite Identification: The Drug  Artemisinin Is Biotransformed By Bacteroides "
                 "Dorei Dsm17855 And A Corresponding Drug Metabolite With M/Z 148.125 Is Detected."},
                {"molecule_chembl_id": "CHEMBL1479", "molecule_pref_name": "DANAZOL",
                 "standard_text_value": "Compound not metabolized", "activity_comment": None}]
        out = S._collapse_activities(rows)
        self.assertEqual([(o["drug"], o["observed"]) for o in out],
                         [("ARTEMISININ", True), ("DANAZOL", False)])
        self.assertEqual(out[0]["metabolites"], ["108.094", "148.125"])
        self.assertEqual(out[0]["activities"], 2)

    def test_gene_records(self):
        out = self.source.build(["CHEMBL5303650"])
        self.assertEqual(out.summary()["genes"], 2)
        pos = [r for r in out.records if r["observed"]]
        self.assertEqual(len(pos), 1)
        r = pos[0]
        self.assertEqual(r["entity_id"], "Q8ABF8")
        self.assertEqual(r["reaction_class"], "Biotransformation of bisacodyl")
        self.assertEqual(r["chebi_substrate"], "CHEBI:3125")
        # the homologous superfamily is too general to be evidence; Pfam/family kept
        self.assertEqual(sorted(r["alt_ids"]), ["IPR000801", "PF00756"])
        self.assertIn("183.0685", r["evidence"])
        self.assertIn("PMID:31158845", r["source"])
        self.assertTrue(any("IPR029058" in w for w in out.warnings))
        neg = [x for x in out.records if not x["observed"]][0]
        self.assertEqual(neg["reaction_class"], "Biotransformation of artemisinin")
        self.assertIn("did not transform", neg["evidence"])

    def test_references_excludes_negatives_by_default(self):
        out = self.source.build(["CHEMBL5303650"])
        refs = out.references()
        self.assertEqual(len(refs), 1)
        self.assertNotIn("observed", refs[0])
        self.assertEqual(len(out.references(observed_only=False)), 2)
        normalise_reference_records(refs)  # must satisfy the reference schema

    def test_microbe_records(self):
        out = self.source.build(["CHEMBL5724809"])
        self.assertEqual(out.summary()["microbes"], 2)
        r = [x for x in out.records if x["observed"]][0]
        self.assertEqual(r["entity_type"], "microbe")
        self.assertEqual(r["entity_id"], "Bifidobacterium adolescentis")
        self.assertEqual(r["entity_name"], "Bifidobacterium adolescentis DSM20083")
        self.assertEqual(r["tax_id"], 1680)

    def test_community_observations(self):
        out = self.source.build(["CHEMBL5725002"])
        self.assertEqual(out.records, [])  # a community is not a reference
        self.assertEqual(len(out.observations), 1)
        c = out.observations[0]
        self.assertEqual(c["community"], "WLS-007")
        self.assertEqual((c["ena_study"], c["ena_sample"]), ("PRJEB37062", "ERS12561742"))
        self.assertEqual(c["mgnify"]["analysis"], "MGYA01030891")
        self.assertEqual([o["observed"] for o in c["observations"]], [True, False])

    def test_build_is_resilient(self):
        out = self.source.build(["CHEMBL5303650", "NOSUCH"])
        self.assertEqual(out.summary()["genes"], 2)
        self.assertTrue(any("NOSUCH" in w for w in out.warnings))


class GridTests(unittest.TestCase):
    """The reaction-class x sample matrix behind the viewer's heatmap."""

    REFS = [{"entity_type": "gene", "entity_id": "IPR001279", "reaction_class": "Biotransformation of a"},
            {"entity_type": "microbe", "entity_id": "Pseudomonas putida", "reaction_class": "Biotransformation of b"}]
    SAMPLES = [
        {"sample_id": "S1", "functions": [{"annotation_id": "IPR001279", "description": "d", "abundance": 5}], "taxa": []},
        {"sample_id": "S2", "taxa": [{"organism": "Pseudomonas putida", "abundance": 3}], "functions": []},
        {"sample_id": "S3", "taxa": [], "functions": []},
    ]

    def test_grid_shape_and_order(self):
        out = run_multi_analysis(self.SAMPLES, self.REFS)
        grid = out["grid"]
        self.assertEqual(grid["sample_ids"], ["S1", "S2", "S3"])
        self.assertEqual(len(out["samples"]), 3)
        # one cell per sample, in sample order, None where nothing was predicted
        for row in grid["rows"]:
            self.assertEqual(len(row["cells"]), 3)
        row_a = next(r for r in grid["rows"] if r["reaction_class"] == "Biotransformation of a")
        self.assertEqual(row_a["cells"][0]["tier"], 1)
        self.assertIsNone(row_a["cells"][1])
        self.assertIsNone(row_a["cells"][2])

    def test_rows_sorted_by_breadth_then_score(self):
        samples = self.SAMPLES[:2] + [dict(self.SAMPLES[0], sample_id="S3")]
        grid = run_multi_analysis(samples, self.REFS)["grid"]
        # 'a' is predicted in two samples, 'b' in one, so 'a' comes first
        self.assertEqual([r["samples_predicted"] for r in grid["rows"]], [2, 1])
        self.assertEqual(grid["rows"][0]["reaction_class"], "Biotransformation of a")

    def test_observations_overlay(self):
        observations = [{"mgnify": {"analysis": "S1"}, "observations": [
            {"reaction_class": "Biotransformation of a", "observed": True},
            {"reaction_class": "Biotransformation of unpredicted", "observed": False}]}]
        grid = run_multi_analysis(self.SAMPLES, self.REFS, observations=observations)["grid"]
        self.assertEqual(grid["observed_samples"], ["S1"])
        self.assertTrue(grid["samples"][0]["has_observations"])
        self.assertFalse(grid["samples"][1]["has_observations"])
        row_a = next(r for r in grid["rows"] if r["reaction_class"] == "Biotransformation of a")
        self.assertIs(row_a["cells"][0]["observed"], True)
        # measured but never predicted: the row still appears, with no score
        row_u = next(r for r in grid["rows"] if r["reaction_class"] == "Biotransformation of unpredicted")
        self.assertIsNone(row_u["cells"][0]["score"])
        self.assertIs(row_u["cells"][0]["observed"], False)
        self.assertEqual(row_u["samples_predicted"], 0)

    def test_observations_for_other_samples_ignored(self):
        observations = [{"mgnify": {"analysis": "SOMETHING_ELSE"},
                         "observations": [{"reaction_class": "Biotransformation of a", "observed": True}]}]
        grid = run_multi_analysis(self.SAMPLES, self.REFS, observations=observations)["grid"]
        self.assertEqual(grid["observed_samples"], [])
        self.assertNotIn("observed", grid["rows"][0]["cells"][0])

    def test_empty_samples_rejected(self):
        with self.assertRaises(DetectionFormatError):
            run_multi_analysis([], self.REFS)


class CacheTests(unittest.TestCase):
    """The on-disk cache: the thing that makes an API outage survivable."""

    def setUp(self):
        import tempfile
        from embl_biotransform.cache import Cache
        self.dir = tempfile.mkdtemp()
        self.cache = Cache(Path(self.dir) / "c.sqlite")
        self.addCleanup(self.cache.close)

    def test_roundtrip_and_namespacing(self):
        self.cache.set("a", "k", {"v": 1})
        self.assertEqual(self.cache.get("a", "k"), {"v": 1})
        self.assertIsNone(self.cache.get("b", "k"))       # same key, other namespace
        self.assertIsNone(self.cache.get("a", "missing"))

    def test_get_or_set_calls_once(self):
        calls = []
        fn = lambda: calls.append(1) or {"n": len(calls)}
        first = self.cache.get_or_set("ns", "k", fn)
        second = self.cache.get_or_set("ns", "k", fn)
        self.assertEqual(first, second)
        self.assertEqual(len(calls), 1)

    def test_failures_are_not_cached(self):
        """A failed fetch must be retried next time, not remembered."""
        def boom():
            raise RuntimeError("api down")
        with self.assertRaises(RuntimeError):
            self.cache.get_or_set("ns", "k", boom)
        self.assertIsNone(self.cache.get("ns", "k"))
        self.assertEqual(self.cache.get_or_set("ns", "k", lambda: "recovered"), "recovered")

    def test_expiry(self):
        self.cache.set("ns", "k", "v", ttl=-1)             # already expired
        self.assertIsNone(self.cache.get("ns", "k"))
        self.cache.set("ns", "keep", "v", ttl=None)
        self.cache.set("ns", "gone", "v", ttl=-1)
        self.assertEqual(self.cache.purge_expired(), 1)
        self.assertEqual(self.cache.get("ns", "keep"), "v")

    def test_survives_reopen(self):
        from embl_biotransform.cache import Cache
        self.cache.set("ns", "k", {"kept": True})
        reopened = Cache(Path(self.dir) / "c.sqlite")
        self.addCleanup(reopened.close)
        self.assertEqual(reopened.get("ns", "k"), {"kept": True})

    def test_make_key_is_stable_and_order_independent(self):
        from embl_biotransform.cache import Cache
        self.assertEqual(Cache.make_key("u", {"a": 1, "b": 2}), Cache.make_key("u", {"b": 2, "a": 1}))
        self.assertNotEqual(Cache.make_key("u", {"a": 1}), Cache.make_key("u", {"a": 2}))

    def test_http_cache_serves_during_an_outage(self):
        """The point of the whole thing: fetch once while the API is up, then
        keep working when it starts failing."""
        from embl_biotransform import fetchers
        from embl_biotransform.cache import install_http_cache
        install_http_cache(self.cache)
        self.addCleanup(setattr, fetchers, "HTTP_CACHE", None)

        state = {"up": True}

        class _Resp:
            status_code = 200
            url = "http://x/api"
            def raise_for_status(self): pass
            def json(self): return {"payload": "real"}

        def fake_get(url, **kwargs):
            if not state["up"]:
                raise RuntimeError("500 Server Error")
            return _Resp()

        original = fetchers.requests.get
        fetchers.requests.get = fake_get
        try:
            self.assertEqual(fetchers._get_json("http://x/api")["payload"], "real")
            state["up"] = False                                   # the API goes down
            self.assertEqual(fetchers._get_json("http://x/api")["payload"], "real")
            with self.assertRaises(fetchers.EMBLAPIError):        # uncached call still fails
                fetchers._get_json("http://x/other", retries=1, backoff=0)
        finally:
            fetchers.requests.get = original


class JobQueueTests(unittest.TestCase):
    def setUp(self):
        from embl_biotransform.jobs import JobQueue
        self.q = JobQueue(max_workers=3)
        self.addCleanup(self.q.shutdown)

    def _drain(self, *jobs, timeout=5.0):
        from embl_biotransform.jobs import TERMINAL
        deadline = time.time() + timeout
        while time.time() < deadline:
            if all(j.state in TERMINAL for j in jobs):
                return
            time.sleep(0.01)
        self.fail("jobs did not finish")

    def test_success_failure_and_progress(self):
        def work(ctx):
            ctx.step(1, 2, "half")
            ctx.step(2, 2, "all")
            return "ok"
        def boom(ctx):
            raise ValueError("nope")
        good, bad = self.q.submit("good", work), self.q.submit("bad", boom)
        self._drain(good, bad)
        self.assertEqual((good.state, good.result, good.progress), ("done", "ok", 1.0))
        self.assertEqual(bad.state, "failed")
        self.assertIn("ValueError: nope", bad.error)

    def test_cancellation(self):
        started = threading.Event()
        def slow(ctx):
            started.set()
            for i in range(500):
                ctx.step(i, 500)
                time.sleep(0.01)
            return "should not finish"
        job = self.q.submit("slow", slow)
        self.assertTrue(started.wait(2))
        self.assertTrue(self.q.cancel(job.id))
        self._drain(job)
        self.assertEqual(job.state, "cancelled")
        self.assertIsNone(job.result)
        self.assertFalse(self.q.cancel(job.id))        # already finished

    def test_rapid_submits_keep_every_job(self):
        """Regression: a job finishing while another was being submitted could
        be pruned as 'stale', because it was terminal before its finish time
        was recorded."""
        jobs = [self.q.submit(f"j{i}", (lambda ctx: 1) if i % 2 else
                              (lambda ctx: (_ for _ in ()).throw(ValueError("x"))))
                for i in range(12)]
        self._drain(*jobs)
        self.assertEqual(len(self.q.list()), 12)
        self.assertTrue(all(self.q.get(j.id) is not None for j in jobs))


class BioSIFTRTests(unittest.TestCase):
    """Reading bioSIFTR output: taxonomy from reads, in the GTDB framework."""

    SPECIES = ("lineage\tsampleA\n"
               "d__Bacteria;p__Bacillota_A;g__Mediterraneibacter;s__Mediterraneibacter_gnavus;MGYG000001234\t0.0421\n"
               "d__Bacteria;p__Bacteroidota;g__Bacteroides;s__Bacteroides_fragilis_A;MGYG000002345\t0.0187\n"
               "d__Bacteria;p__Bacillota;g__Blautia_A;s__Blautia_A_sp900541345;MGYG000003456\t0\n")
    KEGG = "ko_id\tsampleA\nK01580\t42\nK00929\t7\n"
    PFAM = "pfam_id\tsampleA\nPF00756\t12\n"
    INTEGRATED = "feature_id\tsampleA\tsampleB\nPF00756\t12\t0\nPF00005\t3\t8\n"

    def _run_dir(self):
        import tempfile
        root = Path(tempfile.mkdtemp())
        (root / "taxonomy_tables").mkdir()
        (root / "function_tables").mkdir()
        (root / "taxonomy_tables" / "sampleA_sm_species.tsv").write_text(self.SPECIES)
        (root / "function_tables" / "sampleA_sm_community_kegg.tsv").write_text(self.KEGG)
        (root / "function_tables" / "sampleA_sm_community_pfams.tsv").write_text(self.PFAM)
        return root

    def test_species_table(self):
        out = parse_biosiftr_table(self.SPECIES)
        self.assertEqual(out["kind"], "taxa")
        records = out["samples"]["sampleA"]
        self.assertEqual(len(records), 2)                      # the zero row is dropped
        self.assertEqual(records[0]["genome"], "MGYG000001234")
        # the representative genome is split off, so the deepest name is a taxon
        self.assertEqual(clean_taxon_name(records[0]["organism"]), "Mediterraneibacter gnavus")
        self.assertEqual(records[0]["abundance"], 0.0421)

    def test_function_tables(self):
        kegg = parse_biosiftr_table(self.KEGG)
        self.assertEqual(kegg["kind"], "kegg")
        self.assertEqual([r["annotation_id"] for r in kegg["samples"]["sampleA"]], ["K01580", "K00929"])
        pfam = parse_biosiftr_table(self.PFAM)
        self.assertEqual(pfam["kind"], "pfam")
        self.assertIsNone(pfam["samples"]["sampleA"][0]["description"])   # ids only: id matching

    def test_integrated_matrix_fans_out(self):
        out = parse_biosiftr_table(self.INTEGRATED)
        self.assertEqual(out["kind"], "pfam")                  # inferred from the feature ids
        self.assertEqual(sorted(out["samples"]), ["sampleA", "sampleB"])
        self.assertEqual([r["annotation_id"] for r in out["samples"]["sampleB"]], ["PF00005"])

    def test_load_run(self):
        run = load_biosiftr_run(self._run_dir())
        self.assertEqual(sorted(run["samples"]), ["sampleA"])
        entry = run["samples"]["sampleA"]
        self.assertEqual(len(entry["taxa"]), 2)
        self.assertEqual(len(entry["functions"]), 3)           # 2 KO + 1 Pfam
        self.assertEqual(summarise_biosiftr_run(run)["samples"][0]["functions"], 3)

    def test_load_run_rejects_other_directories(self):
        import tempfile
        with self.assertRaises(BioSIFTRFormatError):
            load_biosiftr_run(tempfile.mkdtemp())
        with self.assertRaises(BioSIFTRFormatError):
            load_biosiftr_run(Path(tempfile.mkdtemp()) / "does-not-exist")

    def test_bad_table(self):
        with self.assertRaises(BioSIFTRFormatError):
            parse_biosiftr_table("lineage\n")                  # no sample column
        with self.assertRaises(BioSIFTRFormatError):
            parse_biosiftr_table("")

    def test_feeds_the_pipeline(self):
        """A bioSIFTR profile drops straight into run_analysis."""
        run = load_biosiftr_run(self._run_dir())
        entry = run["samples"]["sampleA"]
        result = run_analysis(entry["taxa"], entry["functions"],
                              [{"entity_type": "microbe", "entity_id": "Mediterraneibacter gnavus",
                                "reaction_class": "Biotransformation of x"},
                               {"entity_type": "gene", "entity_id": "K01580",
                                "reaction_class": "Biotransformation of y"}],
                              sample_id="sampleA")
        self.assertEqual(result["summary"]["predictions"], 2)


class TaxonomyResolverTests(unittest.TestCase):
    """GTDB -> NCBI. The lookups are stubbed; the ladder logic is the point."""

    NCBI = {
        "Mediterraneibacter gnavus": ("33038", "Mediterraneibacter gnavus", "species",
                                      "Bacteria; Bacillota; Clostridia; Lachnospiraceae; Mediterraneibacter;"),
        "Bacteroides fragilis": ("817", "Bacteroides fragilis", "species", "Bacteria; Bacteroidota;"),
        "Blautia": ("572511", "Blautia", "genus", "Bacteria; Bacillota;"),
        "Bacteroides": ("816", "Bacteroides", "genus", "Bacteria; Bacteroidota;"),
    }

    def setUp(self):
        from unittest import mock
        from embl_biotransform import taxonomy as tx
        self.queries = []

        def fake_json(url, **_):
            from urllib.parse import unquote
            name = unquote(url.rsplit("/", 1)[-1])
            self.queries.append(name)
            hit = self.NCBI.get(name)
            if not hit:
                return []                       # ENA answers 200 + [] for unknown names
            tax_id, sci, rank, lineage = hit
            return [{"taxId": tax_id, "scientificName": sci, "rank": rank, "lineage": lineage}]

        patch = mock.patch.object(tx, "_get_json", fake_json)
        patch.start()
        self.addCleanup(patch.stop)
        self.resolver = tx.TaxonomyResolver()

    def test_helpers(self):
        from embl_biotransform.taxonomy import is_placeholder, split_ranks, strip_gtdb_suffix
        self.assertEqual(strip_gtdb_suffix("Bacteroides fragilis_A"), "Bacteroides fragilis")
        self.assertEqual(strip_gtdb_suffix("Enterococcus_G italicus"), "Enterococcus italicus")
        self.assertTrue(is_placeholder("Blautia sp900541345"))
        self.assertFalse(is_placeholder("Coprococcus eutactus"))
        # the polyphyly suffix must go before underscores become spaces
        self.assertEqual(split_ranks("s__Bacteroides_fragilis_A")[0], ("species", "Bacteroides fragilis"))

    def test_resolves_species(self):
        ref = self.resolver.resolve_lineage(
            "d__Bacteria;p__Bacillota_A;g__Mediterraneibacter;s__Mediterraneibacter_gnavus")
        self.assertEqual((ref.tax_id, ref.matched_rank), ("33038", "species"))

    def test_gtdb_suffix_recovered_at_species_level(self):
        ref = self.resolver.resolve_lineage("d__Bacteria;g__Bacteroides;s__Bacteroides_fragilis_A")
        self.assertEqual((ref.tax_id, ref.matched_rank), ("817", "species"))

    def test_placeholder_falls_back_to_genus(self):
        ref = self.resolver.resolve_lineage("d__Bacteria;g__Blautia_A;s__Blautia_A_sp900541345")
        self.assertEqual((ref.tax_id, ref.matched_rank), ("572511", "genus"))
        # a placeholder species is never looked up: it cannot resolve
        self.assertNotIn("Blautia sp900541345", self.queries)

    def test_unresolvable_lineage(self):
        ref = self.resolver.resolve_lineage("d__Nothing;g__CAG-272;s__CAG-272_sp900556615")
        self.assertFalse(ref.resolved)

    def test_annotate_rewrites_into_ncbi(self):
        taxa = [{"organism": "d__Bacteria;g__Bacteroides;s__Bacteroides_fragilis_A",
                 "lineage": "d__Bacteria;g__Bacteroides;s__Bacteroides_fragilis_A",
                 "evidence": "bioSIFTR profile."},
                {"organism": "d__X;g__CAG-272;s__CAG-272_sp900556615",
                 "lineage": "d__X;g__CAG-272;s__CAG-272_sp900556615",
                 "evidence": "bioSIFTR profile."}]
        stats = self.resolver.annotate(taxa)
        self.assertEqual(stats["counts"], {"resolved": 1, "unresolved": 1})
        self.assertEqual(taxa[0]["tax_id"], "817")
        self.assertEqual(taxa[0]["organism"], "Bacteroides fragilis")
        self.assertEqual(taxa[0]["taxonomy_framework"], "NCBI")
        self.assertIn("s__Bacteroides_fragilis_A", taxa[0]["gtdb_organism"])   # provenance kept
        self.assertIsNone(taxa[1].get("tax_id"))
        self.assertIn("unresolved", taxa[1]["taxonomy_framework"])

    def test_taxid_survives_normalisation(self):
        """Regression: run_analysis normalises detections through
        MGnifyClient._normalise_taxonomy, which rebuilt each record from a
        fixed field list and dropped the resolved NCBI tax_id -- silently
        disabling taxon-id matching for every bioSIFTR profile."""
        from embl_biotransform.fetchers import MGnifyClient
        taxa = [{"organism": "d__Bacteria;g__Mediterraneibacter;s__Mediterraneibacter_gnavus",
                 "lineage": "d__Bacteria;g__Mediterraneibacter;s__Mediterraneibacter_gnavus",
                 "abundance": 0.04, "evidence": ""}]
        self.resolver.annotate(taxa)
        normalised = MGnifyClient._normalise_taxonomy(taxa)[0]
        self.assertEqual(normalised["tax_id"], "33038")
        self.assertIn("gtdb_organism", normalised)
        self.assertEqual(normalised["taxonomy_framework"], "NCBI")

    def test_end_to_end_old_name_reference_still_matches(self):
        """The whole point, through run_analysis: a reference curated under the
        old name predicts from a GTDB-profiled sample."""
        taxa = [{"organism": "d__Bacteria;g__Mediterraneibacter;s__Mediterraneibacter_gnavus",
                 "lineage": "d__Bacteria;g__Mediterraneibacter;s__Mediterraneibacter_gnavus",
                 "abundance": 0.04, "evidence": "bioSIFTR profile."}]
        self.resolver.annotate(taxa)
        references = [{"entity_type": "microbe", "entity_id": "Ruminococcus gnavus",
                       "tax_id": "33038", "reaction_class": "Biotransformation of bisacodyl"}]
        result = run_analysis(taxa, [], references, sample_id="gutA")
        self.assertEqual(result["summary"]["predictions"], 1)
        hit = result["predictions"][0]["tier2_hits"][0]
        self.assertEqual(hit["detections"][0]["match_type"], "taxid")
        # without ids it would not match at all
        plain = run_analysis([{"organism": "Mediterraneibacter gnavus", "abundance": 0.04}], [],
                             [{"entity_type": "microbe", "entity_id": "Ruminococcus gnavus",
                               "reaction_class": "Biotransformation of bisacodyl"}], sample_id="gutA")
        self.assertEqual(plain["summary"]["predictions"], 0)

    def test_taxid_matching_beats_naming_drift(self):
        """The reason this module exists: a reference curated under the old
        name still matches a GTDB profile once both carry a taxon id."""
        from embl_biotransform.integrate import taxon_match_record
        taxa = [{"organism": "d__Bacteria;g__Mediterraneibacter;s__Mediterraneibacter_gnavus",
                 "lineage": "d__Bacteria;g__Mediterraneibacter;s__Mediterraneibacter_gnavus",
                 "evidence": ""}]
        self.resolver.annotate(taxa)
        old_name_reference = {"entity_id": "Ruminococcus gnavus", "tax_id": "33038"}
        self.assertIsNone(taxon_match(old_name_reference["entity_id"], taxa[0]["organism"]))
        self.assertEqual(taxon_match_record(old_name_reference, taxa[0]), "taxid")
        # a genus-level reference matches through the resolved NCBI lineage
        self.assertEqual(taxon_match_record({"entity_id": "Mediterraneibacter", "tax_id": "1"}, taxa[0]),
                         "within-lineage")
        # an unrelated organism must not match
        self.assertIsNone(taxon_match_record({"entity_id": "Bacteroides fragilis", "tax_id": "817"}, taxa[0]))


class TaxonGridTests(unittest.TestCase):
    """The microbe x metagenome matrix."""

    REFS = [{"entity_type": "microbe", "entity_id": "Bacteroides fragilis", "tax_id": "817",
             "reaction_class": "Biotransformation of x"}]

    def _run(self, limit=120, linked_only=False):
        samples = [
            {"sample_id": "S1", "functions": [], "taxa": [
                {"organism": "Bacteroides fragilis", "tax_id": "817", "abundance": 30},
                {"organism": "Blautia", "tax_id": "572511", "abundance": 70}]},
            {"sample_id": "S2", "functions": [], "taxa": [
                {"organism": "Blautia", "tax_id": "572511", "abundance": 10},
                {"organism": "Escherichia coli", "tax_id": "562", "abundance": 90}]},
        ]
        return run_multi_analysis(samples, self.REFS, linked_only=linked_only)["taxon_grid"]

    def test_shape_and_share(self):
        grid = self._run()
        self.assertEqual(grid["sample_ids"], ["S1", "S2"])
        self.assertEqual(grid["total_taxa"], 3)
        rows = {r["taxon"]: r for r in grid["rows"]}
        # share is within-sample, so columns are comparable despite depth
        self.assertAlmostEqual(rows["Bacteroides fragilis"]["cells"][0]["share"], 0.3)
        self.assertIsNone(rows["Bacteroides fragilis"]["cells"][1])
        self.assertAlmostEqual(rows["Blautia"]["cells"][1]["share"], 0.1)

    def test_linked_taxa_come_first(self):
        grid = self._run()
        self.assertEqual(grid["rows"][0]["taxon"], "Bacteroides fragilis")
        self.assertTrue(grid["rows"][0]["linked"])
        self.assertTrue(grid["rows"][0]["linked_references"])
        self.assertEqual(grid["linked_taxa"], 1)
        # then by how many samples hold them
        self.assertEqual(grid["rows"][1]["taxon"], "Blautia")
        self.assertEqual(grid["rows"][1]["samples_present"], 2)

    def test_rows_are_capped(self):
        grid = self._run()
        self.assertEqual(grid["truncated"], 0)
        samples = [{"sample_id": "S1", "functions": [],
                    "taxa": [{"organism": f"Taxon {i}", "abundance": i} for i in range(1, 200)]}]
        big = run_multi_analysis(samples, self.REFS, linked_only=False)["taxon_grid"]
        self.assertEqual(len(big["rows"]), 199)          # under the limit
        capped = build_taxon_grid(run_multi_analysis(samples, self.REFS,
                                                     linked_only=False)["samples"],
                                  limit=50, linked_only=False)
        self.assertEqual(len(capped["rows"]), 50)
        self.assertEqual(capped["total_taxa"], 199)
        self.assertEqual(capped["truncated"], 149)

    def test_linked_only_is_the_default(self):
        """A profile has thousands of taxa; only the ones carrying evidence
        belong in a summary, so unlinked rows are hidden unless asked for."""
        grid = self._run(linked_only=True)
        self.assertEqual([r["taxon"] for r in grid["rows"]], ["Bacteroides fragilis"])
        self.assertTrue(all(r["linked"] for r in grid["rows"]))
        self.assertEqual(grid["total_taxa"], 3)
        self.assertEqual(grid["linked_taxa"], 1)
        self.assertEqual(grid["hidden_unlinked"], 2)
        self.assertTrue(grid["linked_only"])
        # the ramp scales over what is shown, not over hidden rows
        self.assertAlmostEqual(grid["max_share"], 0.3)

    def test_taxa_merge_on_taxon_id(self):
        """Two profilers spelling one organism differently still make one row."""
        samples = [
            {"sample_id": "S1", "functions": [],
             "taxa": [{"organism": "Mediterraneibacter gnavus", "tax_id": "33038", "abundance": 5}]},
            {"sample_id": "S2", "functions": [],
             "taxa": [{"organism": "Ruminococcus gnavus", "tax_id": "33038", "abundance": 5}]},
        ]
        grid = run_multi_analysis(samples, self.REFS, linked_only=False)["taxon_grid"]
        self.assertEqual(grid["total_taxa"], 1)
        self.assertEqual(grid["rows"][0]["samples_present"], 2)


class RankFilterTests(unittest.TestCase):
    """Only genus/species detections should support a claim about an organism."""

    TAXA = [
        {"organism": "d__Bacteria", "lineage": "d__Bacteria", "abundance": 500},
        {"organism": "d__Bacteria;p__Bacillota", "lineage": "d__Bacteria;p__Bacillota", "abundance": 300},
        {"organism": "d__Bacteria;g__Blautia", "lineage": "d__Bacteria;g__Blautia", "abundance": 100},
        {"organism": "d__Bacteria;g__Bacteroides;s__Bacteroides_fragilis",
         "lineage": "d__Bacteria;g__Bacteroides;s__Bacteroides_fragilis", "abundance": 50},
    ]

    def test_taxon_rank_detection(self):
        from embl_biotransform.integrate import taxon_rank
        ranks = [taxon_rank(t) for t in self.TAXA]
        self.assertEqual(ranks, ["domain", "phylum", "genus", "species"])
        # explicit fields win over the lineage
        self.assertEqual(taxon_rank({"organism": "x", "matched_rank": "genus"}), "genus")
        self.assertEqual(taxon_rank({"organism": "Pseudomonas putida"}), "species")
        # an empty rank slot (MGnify writes 'k__') must not be read as a rank
        self.assertEqual(taxon_rank({"lineage": "sk__Bacteria;k__;p__Acidobacteriota"}), "phylum")
        self.assertIsNone(taxon_rank({"organism": "Bacteria"}))

    def test_rank_is_at_least(self):
        from embl_biotransform.integrate import rank_is_at_least
        self.assertTrue(rank_is_at_least("species", "genus"))
        self.assertTrue(rank_is_at_least("genus", "genus"))
        self.assertFalse(rank_is_at_least("phylum", "genus"))
        self.assertFalse(rank_is_at_least("genus", "species"))
        # an unknown rank is kept rather than silently discarded
        self.assertTrue(rank_is_at_least(None, "species"))

    def test_filter_levels(self):
        from embl_biotransform.pipeline import filter_taxa_by_rank
        kept, stats = filter_taxa_by_rank(self.TAXA, "genus")
        self.assertEqual(len(kept), 2)
        self.assertEqual(stats["by_rank"], {"domain": 1, "phylum": 1})
        kept, stats = filter_taxa_by_rank(self.TAXA, "species")
        self.assertEqual(len(kept), 1)
        self.assertEqual(stats["dropped"], 3)
        kept, stats = filter_taxa_by_rank(self.TAXA, "any")
        self.assertEqual((len(kept), stats["dropped"]), (4, 0))

    def test_filtering_changes_the_analysis(self):
        references = [{"entity_type": "microbe", "entity_id": "Blautia",
                       "reaction_class": "Biotransformation of x"}]
        loose = run_analysis(self.TAXA, [], references, taxon_rank_floor="any")
        genus = run_analysis(self.TAXA, [], references, taxon_rank_floor="genus")
        strict = run_analysis(self.TAXA, [], references, taxon_rank_floor="species")
        self.assertEqual(loose["summary"]["taxa_detected"], 4)
        self.assertEqual(genus["summary"]["taxa_detected"], 2)
        self.assertEqual(genus["rank_filter"]["dropped"], 2)
        self.assertEqual(genus["settings"]["taxon_rank_floor"], "genus")
        # the genus-level reference loses its support under species-only
        self.assertEqual(genus["summary"]["predictions"], 1)
        self.assertEqual(strict["summary"]["predictions"], 0)


class EnzymeGridTests(unittest.TestCase):
    """The enzyme x sample matrix behind the Enzymes view."""

    def _grid(self, allow_name_matching=True, linked_only=False):
        references = [{"entity_type": "gene", "entity_id": "IPR001279",
                       "entity_name": "Cytochrome P450 monooxygenase",
                       "reaction_class": "Biotransformation of x"}]
        samples = [
            {"sample_id": "S1", "taxa": [], "functions": [
                {"annotation_id": "IPR001279", "description": "Cytochrome P450", "abundance": 40},
                {"annotation_id": "IPR000001", "description": "Kringle", "abundance": 10}]},
            {"sample_id": "S2", "taxa": [], "functions": [
                {"annotation_id": "IPR000001", "description": "Kringle", "abundance": 5}]},
        ]
        return run_multi_analysis(samples, references, linked_only=linked_only,
                                  allow_name_matching=allow_name_matching)["enzyme_grid"]

    def test_shape_and_ordering(self):
        grid = self._grid()
        self.assertEqual(grid["sample_ids"], ["S1", "S2"])
        self.assertEqual(grid["total_enzymes"], 2)
        self.assertEqual(grid["linked_enzymes"], 1)
        # the linked enzyme sorts first
        self.assertEqual(grid["rows"][0]["annotation_id"], "IPR001279")
        self.assertTrue(grid["rows"][0]["linked"])
        self.assertEqual(grid["rows"][0]["match_types"], ["id"])
        self.assertAlmostEqual(grid["rows"][0]["cells"][0]["share"], 0.8)
        self.assertIsNone(grid["rows"][0]["cells"][1])

    def test_name_only_matches_are_flagged(self):
        """A name match is much weaker than an id match and must be visible."""
        references = [{"entity_type": "gene", "entity_id": "EC:1.14.14.1",
                       "entity_name": "Cytochrome P450 monooxygenase",
                       "reaction_class": "Biotransformation of x"}]
        samples = [{"sample_id": "S1", "taxa": [], "functions": [
            {"annotation_id": "IPR999999", "description": "Cytochrome P450 monooxygenase", "abundance": 3}]}]
        grid = run_multi_analysis(samples, references)["enzyme_grid"]
        self.assertEqual(grid["rows"][0]["match_types"], ["name"])
        self.assertTrue(grid["rows"][0]["name_only"])
        self.assertEqual(grid["name_only"], 1)

    def test_unlinked_enzymes_listed_only_when_asked(self):
        shown = self._grid(linked_only=False)
        kringle = next(r for r in shown["rows"] if r["annotation_id"] == "IPR000001")
        self.assertFalse(kringle["linked"])
        self.assertEqual(kringle["samples_present"], 2)

        # by default a sample's ~15,000 annotations collapse to the few that
        # actually support a reference
        default = self._grid(linked_only=True)
        self.assertEqual([r["annotation_id"] for r in default["rows"]], ["IPR001279"])
        self.assertEqual(default["total_enzymes"], 2)
        self.assertEqual(default["hidden_unlinked"], 1)


class ContigTaxonomyTests(unittest.TestCase):
    """Enzyme -> contig -> organism, the join that makes Tier 1 specific."""

    TAXONOMY = [
        "# contig\tclassification\treason\tlineage\tlineage scores",
        "c1\ttaxid assigned\tbased on 543/547 ORFs\t1;131567;2;1783272;1678\t1.00;0.96;0.96;0.69;0.61",
        "c2\ttaxid assigned\tbased on 24/24 ORFs\t1;131567;2\t1.00;0.58;0.58",
        "c3\tno taxid\t-\t\t-",
    ]
    GFF = [
        "##gff-version 3",
        "c1\tPyrodigal\tCDS\t760\t1821\t.\t-\t.\tID=c1_1;eggNOG=1680.BADO_0541;kegg=ko:K03453;pfam=PF01758;interpro=IPR002657,IPR038770",
        "c2\tPyrodigal\tCDS\t1\t500\t.\t+\t.\tID=c2_1;interpro=IPR002657",
        "c4\tPyrodigal\tCDS\t1\t500\t.\t+\t.\tID=c4_1;interpro=IPR999999",
        "c1\tFragGeneScanRS\tCDS\t131\t740\t.\t+\t.\tID=c1_131_740_+",
    ]
    INFO = {"1678": {"scientific_name": "Bifidobacterium", "rank": "genus"},
            "2": {"scientific_name": "Bacteria", "rank": "superkingdom"}}

    def _links(self):
        from embl_biotransform.contig_taxonomy import build_links
        return build_links(self.TAXONOMY, self.GFF)

    def test_parsing_and_join(self):
        from embl_biotransform.contig_taxonomy import parse_contig_taxonomy
        contigs = parse_contig_taxonomy(self.TAXONOMY)
        self.assertEqual(contigs["c1"][-1], "1678")
        self.assertNotIn("c3", contigs)              # no lineage: not a contig we can use

        links = self._links()
        self.assertEqual(links["annotations"]["IPR002657"], {"1678": 1, "2": 1})
        self.assertEqual(links["annotations"]["K03453"], {"1678": 1})   # 'ko:' prefix stripped
        self.assertEqual(links["cds_with_annotations"], 3)   # the row with no annotations is skipped
        self.assertEqual(links["cds_on_unclassified_contigs"], 1)       # c4 has no taxonomy

    def test_attribution_separates_coarse_hits(self):
        """The point: an enzyme seen only on 'Bacteria' contigs is not evidence
        that a named organism carries it."""
        from embl_biotransform.contig_taxonomy import attribution
        links = self._links()
        shared = attribution(links, "IPR002657", self.INFO)
        self.assertEqual(shared["total"], 2)
        self.assertEqual(shared["attributed"], 1)
        self.assertEqual(shared["unattributed"], 1)
        self.assertEqual(shared["taxa"][0]["scientific_name"], "Bifidobacterium")

        # the same accession with only a domain-level contig attributes nothing
        only_coarse = attribution({"annotations": {"X": {"2": 5}}}, "X", self.INFO)
        self.assertEqual((only_coarse["attributed"], only_coarse["unattributed"]), (0, 5))
        self.assertEqual(only_coarse["taxa"], [])

    def test_uninformative_taxids_never_attribute(self):
        from embl_biotransform.contig_taxonomy import attribution, observed_taxids
        links = {"annotations": {"X": {"1": 3, "131567": 2}}}
        self.assertEqual(observed_taxids(links, ["X"]), [])
        self.assertEqual(attribution(links, "X", {})["attributed"], 0)

    def test_attach_to_detections(self):
        from embl_biotransform.contig_taxonomy import attribution
        from embl_biotransform.pipeline import attach_enzyme_taxonomy
        links = self._links()
        table = {a: attribution(links, a, self.INFO) for a in links["annotations"]}
        functions = [{"annotation_id": "IPR002657", "description": "x", "abundance": 2, "evidence": "seen."},
                     {"annotation_id": "IPR038770", "description": "y", "abundance": 1, "evidence": "seen."},
                     {"annotation_id": "NOTLINKED", "description": "z", "abundance": 1, "evidence": "seen."}]

        kept, stats = attach_enzyme_taxonomy([dict(f) for f in functions], table)
        self.assertEqual(len(kept), 3)                    # nothing dropped by default
        self.assertEqual(stats["annotated"], 2)
        self.assertEqual(kept[0]["attributed"], 1)
        self.assertIn("Bifidobacterium", kept[0]["evidence"])

        strict, stats = attach_enzyme_taxonomy([dict(f) for f in functions], table,
                                               require_attribution=True)
        self.assertEqual([f["annotation_id"] for f in strict], ["IPR002657", "IPR038770"])
        self.assertEqual(stats["dropped"], 1)

    def test_strict_mode_changes_predictions(self):
        """End to end: an enzyme with no organism behind it stops predicting."""
        from embl_biotransform.contig_taxonomy import attribution
        links = {"annotations": {"IPR000391": {"2": 4}}}      # only 'Bacteria' contigs
        table = {"IPR000391": attribution(links, "IPR000391", self.INFO)}
        references = [{"entity_type": "gene", "entity_id": "IPR000391",
                       "reaction_class": "Biotransformation of x"}]
        functions = [{"annotation_id": "IPR000391", "description": "dioxygenase", "abundance": 4}]

        loose = run_analysis([], functions, references, enzyme_taxonomy=table)
        strict = run_analysis([], functions, references, enzyme_taxonomy=table,
                              require_enzyme_attribution=True)
        self.assertEqual(loose["summary"]["predictions"], 1)
        self.assertEqual(strict["summary"]["predictions"], 0)
        self.assertEqual(strict["attribution_filter"]["dropped"], 1)

    def test_enzyme_grid_reports_attribution(self):
        from embl_biotransform.contig_taxonomy import attribution
        links = self._links()
        table = {a: attribution(links, a, self.INFO) for a in links["annotations"]}
        samples = [{"sample_id": "S1", "taxa": [], "enzyme_taxonomy": table, "functions": [
            {"annotation_id": "IPR002657", "description": "x", "abundance": 2},
            {"annotation_id": "IPR999999", "description": "unplaced", "abundance": 1}]}]
        grid = run_multi_analysis(samples, [{"entity_type": "gene", "entity_id": "IPR002657",
                                             "reaction_class": "Biotransformation of x"}])["enzyme_grid"]
        placed = next(r for r in grid["rows"] if r["annotation_id"] == "IPR002657")
        self.assertEqual(placed["taxa"][0]["scientific_name"], "Bifidobacterium")
        self.assertEqual(grid["with_taxonomy"], 1)


class MatchIndexEquivalenceTests(unittest.TestCase):
    """The detection index only narrows candidates; it must never change which
    pairs match. Checked against a brute-force pass over every pair."""

    @staticmethod
    def _brute_force(taxa, functions, references, allow_name_matching=True):
        """What the matcher did before indexing: compare everything."""
        from embl_biotransform.integrate import ids_match, name_match, taxon_match_record
        links = set()
        for i, ref in enumerate(references):
            if ref["entity_type"] == "microbe":
                for t in taxa:
                    if taxon_match_record(ref, t):
                        links.add((i, "taxon", t["organism"]))
            else:
                ref_ids = [ref["entity_id"], *ref.get("alt_ids", [])]
                for f in functions:
                    det = str(f.get("annotation_id") or "")
                    if any(ids_match(r, det) for r in ref_ids):
                        links.add((i, "id", det))
                    elif allow_name_matching and name_match(ref.get("entity_name", ""),
                                                            f.get("description") or ""):
                        links.add((i, "name", det))
        return links

    @staticmethod
    def _via_graph(taxa, functions, references, allow_name_matching=True):
        result = run_analysis(taxa, functions, references, taxon_rank_floor="any",
                              allow_name_matching=allow_name_matching)
        order = {f"ref:{r['entity_type']}:{r['entity_id']}:{r['reaction_class']}": i
                 for i, r in enumerate(references)}
        links = set()
        for row in result["taxa"]:
            for ref in row["linked_references"]:
                links.add((order[ref], "taxon", row["organism"]))
        graph = result["graph"]
        by_node = {n["id"]: n for n in graph["nodes"]}
        for edge in graph["edges"]:
            source, target = by_node.get(edge["source_node"]), by_node.get(edge["target_node"])
            if not source or source.get("type") != "functional_annotation":
                continue
            if target and target.get("type") == "reference":
                links.add((order[target["id"]], edge.get("match_type"), source.get("annotation_id")))
        return links

    def _case(self, seed):
        """Inputs with no duplicate reference statements or detections.

        The graph deliberately merges references sharing
        (entity_type, entity_id, reaction_class) -- they are the same
        statement -- and detections sharing an id or organism name. Both are
        correct, but they would confound a comparison against a brute-force
        pass over the raw records, so the generator avoids them.
        """
        import random
        rng = random.Random(seed)
        names = ["Cytochrome P450 monooxygenase", "Aromatic ring hydroxylating dioxygenase",
                 "Acetyl esterase", "Nitroreductase family protein", "Azoreductase"]
        genera = ["Bacteroides", "Blautia", "Eggerthella", "Dorea"]
        references, taxa, functions = [], [], []
        for i in range(30):
            if rng.random() < 0.5:
                references.append({"entity_type": "gene",
                                   "entity_id": rng.choice(["IPR00%04d" % rng.randrange(20),
                                                            "EC:1.14.13.-", "K0%04d" % rng.randrange(20)]),
                                   "alt_ids": ["PF00%03d" % rng.randrange(10)] if rng.random() < 0.4 else [],
                                   "entity_name": rng.choice(names),
                                   "reaction_class": "R%d" % (i % 5)})
            else:
                genus = rng.choice(genera)
                references.append({"entity_type": "microbe",
                                   "entity_id": genus if rng.random() < 0.4 else f"{genus} species{rng.randrange(5)}",
                                   "tax_id": str(rng.randrange(800, 810)) if rng.random() < 0.5 else None,
                                   "reaction_class": "R%d" % (i % 5)})
        for i in range(40):
            functions.append({"annotation_id": rng.choice(["IPR00%04d" % rng.randrange(20),
                                                           "1.14.13.%d" % rng.randrange(5),
                                                           "PF00%03d" % rng.randrange(10)]),
                              "description": rng.choice(names + ["Hypothetical protein", "ABC transporter"]),
                              "abundance": rng.randrange(1, 50)})
        for i in range(25):
            genus = rng.choice(genera)
            taxa.append({"organism": f"{genus} species{rng.randrange(5)}",
                         "lineage": f"Bacteria; {genus};",
                         "tax_id": str(rng.randrange(800, 810)) if rng.random() < 0.5 else None,
                         "abundance": rng.randrange(1, 40)})
        references = [{k: v for k, v in r.items() if v is not None} for r in references]
        taxa = [{k: v for k, v in t.items() if v is not None} for t in taxa]

        def unique(items, key):
            seen, out = set(), []
            for item in items:
                if key(item) not in seen:
                    seen.add(key(item))
                    out.append(item)
            return out

        references = unique(references, lambda r: (r["entity_type"], r["entity_id"], r["reaction_class"]))
        functions = unique(functions, lambda f: f["annotation_id"])
        taxa = unique(taxa, lambda t: t["organism"])
        return taxa, functions, references

    def test_index_matches_brute_force(self):
        for seed in range(12):
            taxa, functions, references = self._case(seed)
            for allow_names in (True, False):
                self.assertEqual(
                    self._via_graph(taxa, functions, references, allow_names),
                    self._brute_force(taxa, functions, references, allow_names),
                    f"index and brute force disagree (seed={seed}, name_matching={allow_names})")

    def test_scales_to_a_large_reference_set(self):
        """A fully built ChEMBL set against a real-sized metagenome: the
        pairwise form took minutes, so guard against it creeping back."""
        import time
        functions = [{"annotation_id": "IPR%06d" % i, "description": "family %d" % i,
                      "abundance": i % 50} for i in range(5000)]
        references = [{"entity_type": "gene", "entity_id": "IPR%06d" % (i % 9000),
                       "entity_name": "Enzyme family %d" % i,
                       "reaction_class": "R%d" % (i % 200)} for i in range(20000)]
        started = time.time()
        result = run_analysis([], functions, references)
        elapsed = time.time() - started
        self.assertEqual(result["summary"]["predictions"], 200)
        self.assertLess(elapsed, 20, f"matching 20k references took {elapsed:.1f}s")


if __name__ == "__main__":
    unittest.main()
