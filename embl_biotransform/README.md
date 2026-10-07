# embl_biotransform

A small, extensible pipeline that integrates evidence from EMBL-EBI /
SIB resources to predict the biotransformation capabilities of a
metagenome sample, and shows the evidence behind each prediction.

## The data model

Three kinds of input, as specified:

1. **Reference biotransformation data**, tied to either an individual
   **microbe** or an individual **gene/enzyme**. There's no single EMBL
   API that packages this directly, so it's supplied as a file
   (`reference.py` loads JSON or CSV) -- see
   `example_data/reference_associations.json`. In practice you'd build
   this by combining literature curation with lookups against
   [UniProt](https://www.uniprot.org) (gene → EC number → organism) and
   [Rhea](https://www.rhea-db.org) (EC number → reaction).
2. **A metagenome analysis**, giving both:
   - which **microbes** are present (taxonomic profile), and
   - which **enzymes/functional families** are present (functional
     annotation) -- from [MGnify](https://www.ebi.ac.uk/metagenomics),
     EMBL-EBI's metagenomics resource.
3. **Evidence tiering**, exactly as described: an **enzyme detected
   directly** in the metagenome is the strongest evidence for
   biotransformation activity (**Tier 1**); a **microbe's presence**
   alone, without its specific enzyme being confirmed, is weaker,
   inferred evidence (**Tier 2**).

## Architecture

```
fetchers.py       -- ChEBI / Rhea / UniProt / ChEMBL / MGnify (API v2, plain requests) clients.
                      Each supports both a live fetch() and a load_file() for offline/reproducible runs.
                      MGnifyClient also does study search / analysis listing (paged).
mgnipy_client.py   -- MGnipySource: adapter over `mgnipy` (https://github.com/EBI-Metagenomics/mgnipy),
                      the official client for MGnify API v2. Optional dependency.
reference.py       -- loads + validates curated per-microbe / per-gene associations (JSON/CSV/TSV).
chembl_references.py -- builds those associations from ChEMBL "Bacterial Biotransformation" assays
                      (gene / microbe / microbiome tiers), with ChEBI substrates via UniChem.
integrate.py       -- EvidenceGraph: a typed, evidence-annotated networkx graph tying together
                      sample -> detected taxa / functional annotations -> reference associations
                      -> reaction classes (-> substrate compound, where known).
predict.py         -- BiotransformationPredictor: scores each reaction class per sample using the
                      Tier-1 (enzyme) / Tier-2 (microbe) evidence model, with the full trail kept.
pipeline.py        -- run_analysis(): one call from detections + references to a JSON-ready result;
                      parse_detections_text(): read MGnify JSON or simple CSV/TSV tables.
visualize.py       -- matplotlib plots of the graph and the ranked scores (for notebooks / reports).
contig_taxonomy.py -- joins an assembly's annotation GFF to its contig taxonomy, so each
                      detected enzyme can be attributed to the organism carrying it.
biosiftr.py        -- reads bioSIFTR output (shallow-shotgun profiles from reads) into the same
                      detection records as everything else.
taxonomy.py        -- resolves GTDB taxon names to NCBI taxon ids, so profiles and references
                      can be matched on id rather than on name.
cache.py           -- persistent on-disk (SQLite) cache; install_http_cache() puts every API call
                      and results file behind it, so an outage is survivable.
jobs.py            -- small thread-backed job queue with progress + cancellation, so slow or
                      flaky work never blocks a request.
webapp/            -- the EXPOSE viewer: Flask API (app.py) + a dependency-free single-page front end.
expose_app/        -- the application shell: the `expose` CLI and the native desktop window.
tests/             -- unit tests for matching, parsing, scoring and the web API.
```

Every edge in the graph carries a `source` and an `evidence` string, so
every prediction is traceable back to exactly which detection(s) and
which curated statement(s) produced it -- nothing is a black box.

### How detections are linked to reference statements

Each link records a `match_type`, shown in the viewer and in exports:

| Reference kind | match_type     | Rule |
|----------------|----------------|------|
| gene           | `id`           | detected annotation id equals `entity_id` or one of `alt_ids` (EC wildcards like `1.14.13.-` honoured) |
| gene           | `name`         | no id match, but >= 2 distinctive words of the enzyme name appear in the detected family description (generic words such as "family", "subunit", "protein" ignored). Weaker: counts half in the score, and can be switched off |
| microbe        | `exact`        | same taxon (MGnify lineage strings like `g__Pseudomonas;s__Pseudomonas_putida` are cleaned first) |
| microbe        | `within-genus` | the reference names a genus and the detected taxon belongs to it |

Add InterPro / KO ids to a gene reference's `alt_ids` to get exact `id`
matches against MGnify's functional annotation.

Scores count *distinct reference statements*: five detected strains of a
genus that all match one genus-level reference count once.

## The EXPOSE application

```bash
# once (Python >= 3.10)
python3.11 -m venv ~/.venvs/expose
source ~/.venvs/expose/bin/activate
pip install -e '.[desktop]'      # drop [desktop] for a browser tab instead

# every time, from anywhere
expose                            # native window
expose --browser                  # browser tab
expose --headless --host 0.0.0.0  # serve to another machine
```

`python -m webapp` still works as an alias.

### Caching

Every API response and results file is cached on disk, keyed by URL, along
with composites such as a built set of ChEMBL references. Accession-keyed
records never expire (an MGnify analysis is immutable once published);
searches expire after a day. **Only successes are cached**, so a failed call
is retried rather than remembered as a failure.

This is what makes the tool usable against unreliable services: anything
fetched before an outage stays available during one, and re-running a
comparison of the same analyses is instant instead of a minute of downloads.

```bash
expose cache stats     # what is cached, and how big
expose cache clear     # drop it all, or one namespace
expose cache purge     # drop only expired entries
expose --no-cache      # ignore it for one run
```

Location: `~/Library/Caches/expose` (macOS), `$XDG_CACHE_HOME/expose`
(Linux), `%LOCALAPPDATA%\expose\cache` (Windows), or `$EXPOSE_CACHE_DIR`.

### Background jobs

Slow or outage-prone work -- ChEMBL searches and builds, bulk MGnify
downloads, bioSIFTR taxonomy resolution -- runs as a background job. The
request returns a job id immediately and the UI polls it, showing progress and
a cancel button, so a ChEMBL outage can no longer hang the page. Because every
step writes through the cache, a job that fails part-way leaves its completed
fetches behind and a re-run continues from there.

| route | purpose |
|-------|---------|
| `GET /api/jobs/<id>` | state, progress, result |
| `GET /api/jobs` | recent jobs |
| `POST /api/jobs/<id>/cancel` | ask a job to stop |
| `GET /api/cache`, `DELETE /api/cache` | cache stats / clear |

In the viewer you:

1. Pick a sample: the bundled example, **bioSIFTR** output (see below), your
   own files (MGnify JSON, or
   CSV/TSV with `organism,abundance` / `annotation_id,description,abundance`
   columns), or live from **MGnify API v2**, restricted to **pipeline v6**
   assembly analyses with contig-based taxonomy (the assignment that ties each
   enzyme to its organism)
   analyses (the MGnify tab opens on study MGYS00010462, a v6 assembly study
   known to work end to end; search studies -> pick an analysis ->
   choose a taxonomy (best available, SSU/LSU rRNA, contig taxonomy, PR2, UNITE,
   ...) and InterPro / KEGG / Pfam / Rhea / GO-slim functions). Functions come
   from the per-analysis summary tables v2 links to, because v2's JSON
   annotation endpoints carry no accessions to match on. Amplicon analyses
   have taxonomy only.
2. Pick reference associations: the bundled example or your own JSON/CSV.
3. Run, then explore. With **several analyses selected** the results open on a
   comparison grid (see below). For a single analysis it goes straight to the
   per-sample view: ranked predictions with the evidence behind each one;
   an evidence graph (click a prediction or node to trace its evidence;
   large sets of matching taxa collapse into one node); tables of every
   detected taxon/function and whether it was linked; export to JSON or CSV.
   Identifiers link out to InterPro, KEGG, ExPASy ENZYME, ChEBI, NCBI
   Taxonomy and MGnify.

### Comparing several analyses

Select several analyses (ctrl/cmd-click or shift-click) and they are all
loaded, scored against the same references, and laid out as a
**reaction class x analysis** grid, in two views:

- **Heatmap** -- one row per reaction class predicted in *any* analysis, cell
  shade = confidence. Rows are ordered by how many analyses share the class,
  so what the communities have in common floats to the top. `▲` means the
  enzyme itself was detected (Tier 1), `○` microbe presence only (Tier 2),
  `·` nothing predicted. Click a cell to jump to that analysis's evidence.
All three grids show **only rows where a link was made** -- the enzymes and
microbes that actually carry evidence, not the whole profile. On a real gut
analysis that is 9 enzymes rather than 14,958, and the grid becomes a summary
you can read. A checkbox in the legend reveals the rest when you want it.

- **Microbes** -- the same layout for organisms: one row per taxon, cell
  shaded by that taxon's **share of its own sample**. Abundances are not
  comparable as reported (bioSIFTR gives fractions, MGnify read counts, and
  depth differs), so the share is what is coloured, scaled against the grid's
  largest -- otherwise every cell lands in the lightest bucket, since one
  organism rarely exceeds a few per cent of a community. Taxa linked to a
  reference sort first and carry a `●`, so the microbes actually carrying the
  evidence are visible against the rest of the community. Rows merge on
  **taxon id**, so two profilers spelling one organism differently still make
  one row. Capped at 120 rows, with the remainder reported.
- **Enzymes** -- one row per linked functional annotation (InterPro / Pfam /
  KO / EC / Rhea), with a **Found in** column naming the organisms whose
  contigs carry it. Cells show the **count**, shaded by share of that sample's
  annotations: a sample has tens of thousands of entries, so one family's
  share rounds to 0.00% and says nothing, while the count does. Anything
  matched only by enzyme *name* rather than by id is badged **name only**, and
  an entry no contig placed below domain says so in red -- between them, those
  two labels show how much of a prediction rests on weak evidence.
- **Small multiples** -- the top predictions per analysis, side by side.

Confidence uses a single-hue sequential ramp in five steps (magnitude, so one
hue, light to dark); tier is carried by the glyph rather than by colour, and
dark mode has its own steps chosen against the dark surface. Every step's
label ink clears 4.5:1.

Where a ChEMBL **community assay** names one of the loaded analyses, its
measured outcomes are overlaid: `✓` observed in vitro, `✗` not observed, and a
`◆` on that column's header. A class that was *measured but never predicted*
still gets a row, with no score -- which is the false negative you want to see.
"Export grid CSV" writes out whichever grid is on screen -- reactions with a
score / tier / observed triple per analysis, microbes and enzymes with
abundance and share -- as `expose_<view>_<n>_analyses.csv`.

Only the MGnify search/download (the v2 API and the result files on the MGnify FTP site) goes to the internet; MGnify responses are
cached in memory for an hour.


## Taxonomic profiles from reads (bioSIFTR)

MGnify's V6 *assembly* analyses carry no rRNA profile, so their taxonomy falls
back to contig assignments and sees only what assembled -- on one gut analysis
the `unclassified` bucket alone was larger than every classified taxon
combined. [bioSIFTR](https://github.com/EBI-Metagenomics/biosiftr) profiles
shallow-shotgun *reads* against MGnify's biome-specific genome catalogues
instead, and derives its functional tables from the same catalogue genomes, so
taxonomy and function become two views of one profile.

Running it needs Nextflow, containers and a large reference database, so
EXPOSE **consumes its output** rather than running it. Run the pipeline
wherever suits, then in the viewer pick **bioSIFTR** and give it the
pipeline's `--outdir`:

```
<outdir>/taxonomy_tables/<sample>_sm_species.tsv           lineage + relative abundance
<outdir>/function_tables/<sample>_sm_community_kegg.tsv    ko_id + counts
<outdir>/function_tables/<sample>_sm_community_pfams.tsv   pfam_id + counts
<outdir>/integrated_annotation/*_matrix.tsv                all samples at once
```

Every one of these is "a feature id, then one column per sample", so a single
reader handles them all and an integrated matrix fans out into several samples
automatically. Uploading the tables works too, for when the results sit on
another machine. `--run_bwa` output is used instead if you pick that mapper.

Two details that matter:

- The species lineage ends with the catalogue's **representative genome**
  (`...;s__Escherichia_coli;MGYG000000001`). That is not a rank, so it is split
  off into `genome`; left in place, every organism would read as
  "MGYG000000001".
- bioSIFTR gives ids without descriptions, so functional matching is by **id
  only** -- which is the precise path anyway.

### Attributing enzymes to organisms

The per-analysis InterPro/Pfam/KO **summary** tables are community-wide
counts: they say a family occurs *somewhere* in the sample, with no link to
which organism carries it. A prediction built on them claims only "this
community contains the enzyme, and separately contains the microbe" -- which
is weaker than it looks, and makes predictions too liberal.

Pipeline v6 publishes the pieces to do better. Each CDS in
`<assembly>_annotation_summary.gff.gz` sits on a named contig, and
`<assembly>_contigs_taxonomy.tsv.gz` gives each contig an **NCBI taxid
lineage** with the classifier's per-rank confidence. Joining on the contig
gives enzyme -> organism. Tick **Attribute enzymes to organisms** under
*Options*; the result is cached per analysis.

**Why the annotation summary rather than InterProScan.** The same join can be
made from `<assembly>_interproscan.tsv.gz`, whose protein ids
(`ERZ29588192_2784_12`) contain the contig (`ERZ29588192_2784`). But that file
is **~1 GB compressed per analysis** against **12 MB** for the annotation
summary, which additionally carries Pfam, KEGG and GO rather than InterPro
alone. Same answer, ~80x less to download. (InterProScan remains the only
source of *which part* of a protein matched, which nothing here needs.)

Measured on one real gut assembly: 136,548 contigs and 392,435 CDS rows stream
and join in about 6 seconds, reducing to 24,336 annotations over 319,850
(annotation, taxid) pairs -- roughly 4 MB, cached, so it happens once.

**What it shows, and why it tightens things.** Many contigs cannot be
classified below domain. On two real analyses, of the 9 InterPro entries
supporting a prediction, **7 sat only on contigs assigned no further than
"Bacteria"**:

```
IPR004183   on Clostridium sp. LIBA-8841 x2, Eggerthella
IPR001663   no organism below genus
IPR000391   no organism below genus
IPR014436   on Acidaminococcus
IPR015879   no organism below genus            ... and four more
```

Those occurrences are real, but they are not evidence that a *named* organism
carries the enzyme. The Enzymes grid labels each row with its carriers, or
flags it in red. **Use only enzymes placed on a genus/species contig** makes
that strict: on the same run it dropped 7,991 enzyme detections and left 2 of
the 9 linked entries standing.

### Working at scale

A fully built ChEMBL reference set is far larger than the worked examples, so
the limits that used to guard the small case are gone:

- **No upload cap.** This reads files off the same machine, so the usual
  reason for one does not apply. A 67 MB / 700,000-record reference CSV parses
  in about 5 seconds. (Reinstate `MAX_UPLOAD_MB` in `webapp/app.py` if you
  ever serve this with `--host 0.0.0.0`.)
- **No cap on assays per build**, and up to 500 MGnify analyses per batch.
- Matching is **indexed**, not pairwise. Detections are put into lookup tables
  by annotation id, EC number, description word, taxon id, taxon name and
  genus, so a reference is compared only against plausible candidates. The
  pairwise form took 103 s for 2,000 references against a 15,000-annotation
  metagenome and would have taken roughly 45 minutes for 50,000; indexed, the
  same 50,000 take **1.7 s**.

The index only narrows the candidate set -- every candidate still goes through
the same `ids_match` / `name_match` / `taxon_match_record` decision -- and a
test compares it against a brute-force pass over every pair, across a dozen
randomised inputs, so the speed-up cannot quietly change what matches.

Two things the graph merges by design, worth knowing when reading counts:
references sharing `(entity_type, entity_id, reaction_class)` are the same
statement and become one node, as do detections sharing an annotation id or
organism name.

### Which taxa count as evidence

A contig or read profile reports **every level** of its lineage, so a sample
yields `d__Bacteria` and `p__Pseudomonadota` rows alongside real species.
Anything coarser than genus says a *group* is present, which cannot support a
claim about the organism a reference names, and it inflates the abundance
shares of everything else. So detections are filtered by rank before the graph
is built:

| setting | keeps | effect on a real gut analysis (2,837 taxa) |
|---------|-------|--------------------------------------------|
| Genus and species *(default)* | genus, species | 156 dropped |
| Species only | species | 436 dropped |
| Any rank | everything | nothing dropped |

The rank is taken from NCBI resolution where it happened, else an explicit
`rank`, else the deepest prefix in the lineage (`...;g__Blautia` is a genus).
A detection whose rank cannot be determined is **kept** -- discarding an
uploaded `organism` column because it carried no rank would lose good data --
and a bare binomial is read as a species. The count dropped is reported in the
result, and the control sits under *Options*.

### GTDB to NCBI

bioSIFTR lineages are **GTDB**; reference associations are **NCBI** (ChEMBL
supplies `assay_tax_id`). Comparing the two by name fails silently and in both
directions, so `taxonomy.py` resolves the profile into NCBI terms and matching
happens on **taxon id** (`match_type: taxid`), with names only as a fallback.

Neither the MGnify API nor the catalogue metadata carries an NCBI lineage --
both are GTDB-only -- so names are resolved against ENA's taxonomy service,
walking up the lineage until something matches:

| case | example | outcome |
|------|---------|---------|
| adopted rename | reference says `Ruminococcus gnavus`, profile says `Mediterraneibacter gnavus` | both are taxid **33038** -- matches |
| polyphyly suffix | `Bacteroides fragilis_A` | suffix stripped -> **817**, species level |
| GTDB placeholder | `Blautia_A sp900541345` | no NCBI species; genus **572511** |
| GTDB-only | `CAG-272 sp900556615` | unresolved, and reported as such |

In the human-gut catalogue about **three quarters** of species are
placeholders of the last two kinds. That is expected rather than a failure: an
uncultured species is not going to be the subject of a drug-metabolism assay
either. What matters is that the *named, cultured* organisms -- the ones
references are actually about -- match despite the naming drift. On a worked
example, 2 of 4 microbe references matched only by taxon id and would have
been missed entirely by name.

The GTDB name is kept as `gtdb_organism` / `gtdb_lineage`, and the evidence
line says when a match was made above species level.

## Building reference associations from ChEMBL

Curating associations by hand is the slow part, so `chembl_references.py`
builds them from ChEMBL's `Bacterial Biotransformation` assays (~47,000
activities, largely from gut-microbiome drug-metabolism screens).

In the viewer this is the **ChEMBL** tab under *Reference associations*:
search assays by organism, drug or accession, tick the ones you want, and
"Build references". It reports how many positive and negative records came
back, how many reaction classes they span, and anything dropped along the way.
From a script:

```python
from embl_biotransform import ChEMBLBiotransformationSource, run_analysis

source = ChEMBLBiotransformationSource()
out = source.build(["CHEMBL5303650", "CHEMBL5724809", "CHEMBL5725002"])

print(out.summary())
result = run_analysis(taxa, functions, out.references(), sample_id="MGYA01030891")
```

ChEMBL records these experiments at three levels, which map onto the
pipeline's evidence tiers. Which one an assay is, is decided by
`confidence_score` and whether its target is the generic `CHEMBL612558`
ADMET placeholder:

| tier | example | yields |
|------|---------|--------|
| **gene** | `CHEMBL5303650` -> UniProt `Q8ABF8` | a `gene` reference; the target's InterPro / Pfam cross-references become `alt_ids`, so it matches MGnify by **id** rather than by name |
| **microbe** | `CHEMBL5724809` -> *Bifidobacterium adolescentis* DSM20083 | a `microbe` reference (Tier 2) |
| **microbiome** | `CHEMBL5725002` -> community WLS-007 | not a reference but an **observation**: what this community actually did |

### The library

Every assay built is kept in the on-disk cache, keyed by its accession, and
listed under the ChEMBL tab with its tier, organism and counts. Tick and
untick them to change which references the analysis uses; the selection is
saved with the library, so it survives closing the app.

Nothing is re-fetched to do this -- toggling an assay back on costs one
lookup in SQLite -- which matters when ChEMBL is unavailable, as it often is.
A `×` drops an assay from the library.

There is **no upload limit** and **no cap** on how many assays a build may cover: assembling a large
reference set is the point. Requests to ChEMBL are paced at one per second
(`fetchers.RATE_LIMITS`), so a few hundred new assays takes a while -- but the
pacing applies only to calls that actually reach the network, and every assay
is cached by accession, so re-running the same list costs nothing. Search
results and the build button both say how many of the selected assays are
already held and how many must be fetched.

You can also **upload a list of ChEMBL ids** (one per line, or any text or CSV
with accessions in it). The list is cross-matched against the library first,
so only the assays not already held are fetched, and the whole list is then
selected.

MGnify analyses are cached the same way and listed under the MGnify tab, with
the taxonomy/function choice they were fetched with, so re-running a
comparison costs nothing and works with MGnify unreachable:
`GET /api/mgnify/library`, `DELETE /api/mgnify/library/<accession>`. The
analysis picker marks each entry **cached** (already downloaded) and
**linked** (enzyme-to-organism links already built), so it is clear what is
free to load.

The header carries the cache size and a **Cache…** button to clear it. That
discards hours of downloads -- and, while ChEMBL is down, references that
cannot be rebuilt at all -- so it asks for a second, explicit confirmation
naming how many entries and how much disk is about to go.

| route | purpose |
|-------|---------|
| `GET /api/chembl/library` | what has been built, and what is selected |
| `DELETE /api/cache` | clear the cache (the UI confirms first) |
| `POST /api/chembl/library/match` | cross-match an id list; say what is held and what is missing |
| `POST /api/chembl/library/selection` | set the selection |
| `GET /api/chembl/library/references` | merged references for the selection |
| `DELETE /api/chembl/library/<assay>` | drop one |

### Positives and negatives

Both are kept, with an `observed` flag. A "no biotransformation occurred"
result is real evidence of absence and rare in curated sets, so it is
recorded rather than dropped. `references()` returns only the positives by
default, because `run_analysis` treats every record it is given as evidence
*for* a reaction -- passing a negative would invert its meaning. Use
`references(observed_only=False)` to get everything.

`reaction_class` stays at the level ChEMBL actually asserts -- ChEMBL says
*that* a biotransformation happened, never what chemistry it was -- so it is
`"Biotransformation of <drug>"`, one class per substrate, with
`chebi_substrate` resolved through UniChem (ChEMBL molecule records carry no
ChEBI cross-reference of their own).

### Reading the outcome

Each published screen words its result differently, so an outcome is read from
two places: `standard_text_value`, a short curated verdict ("Compound
metabolized"), which is trusted first, and failing that the free-text
`activity_comment`. Three wordings seen so far:

| screen | wording |
|--------|---------|
| Zimmermann gene screen | "...Was **Proven To Be Mediated** By Bt_0152" / "**Could Not Be Mediated**" |
| Baghai Arassi strain screen | "**Biotransformation Occurred** For Drug ..." / "**No Biotransformation Occurred**" |
| *Bacteroides dorei* screen | "The Drug ... **Is Biotransformed By** ..." + verdict "Compound metabolized" |

Negatives are matched before positives, because every negative contains its
positive as a substring ("No Biotransformation Occurred" contains
"biotransformation occurred"; "not metabolized" contains "metabolized"). An
unrecognised wording yields no record rather than a guessed one, and the
warning quotes a real comment from that assay so the pattern can be extended.

### Over-general identifiers are dropped

A protein target's InterPro cross-references include homologous superfamilies
such as `IPR029058` (Alpha/Beta hydrolase fold, ~1.7 million proteins), which
would match nearly any metagenome and manufacture Tier-1 hits. Only `family`
and `domain` entries are kept; anything dropped is reported in
`out.warnings`. Change this with `interpro_types=(...)`, or `()` to keep all.

### Community assays: validation, and stronger evidence

Community assays are rare -- almost no metagenome has measured
drug-metabolism data. Where one exists it is worth more than a reference,
because it can be tied to a real MGnify analysis: the assay parameters carry
ENA study and sample accessions, and
`MGnifyClient.find_analysis_for_sample()` resolves the sample to the v6
analysis of its assembly:

```
ERS12561742 -> sample SAMEA110463714 -> assembly ERZ29587455 -> analysis MGYA01030891 (MGYS00010462)
```

Note this goes via the *sample*, not via `studies/insdc/{project}`: a
project's reads study and the assembly study MGnify derives from it have
different accessions, and the INSDC lookup returns the reads study
(`PRJEB37062` -> `MGYS00010430`), which has no analyses.

So for such a sample you can predict from gene and microbe evidence, then
check the prediction against what was measured -- and where a prediction is
confirmed, the gene, the organism and the community observation form three
independent, agreeing lines of evidence.

## Tests

```bash
python -m pytest            # or: python -m unittest discover tests
```

## Interface: Jupyter notebook

`biotransformation_explorer.ipynb` is the interactive front end. It lets
a user:

1. Pick a data source: live MGnify (via `mgnipy`, API v2) with a study
   search → analysis picker, or point at previously-downloaded JSON
   files (`example_data/sample_taxonomy.json`,
   `example_data/sample_functions.json` are provided as a worked
   example).
2. Point at a reference-associations file.
3. Run the integration + prediction step.
4. See the evidence graph and the ranked, annotated confidence chart
   inline.

`mgnipy` is a young, actively-evolving package -- if a call in
`mgnipy_client.py` no longer matches its current API, check
https://mgnipy.mgnify.org/. The shape of what comes back downstream (a
list of `{organism, source, evidence}` / `{annotation_id, description,
source, evidence}` dicts) is what `integrate.py` actually depends on, so
you can adjust the mgnipy calls without touching anything else.

## Quick start (offline, using the bundled example data)

```python
from embl_biotransform import (
    EvidenceGraph, BiotransformationPredictor, load_reference_associations,
    plot_evidence_graph, plot_evidence_scores,
)
import json

taxa = json.load(open("example_data/sample_taxonomy.json"))
functions = json.load(open("example_data/sample_functions.json"))
references = load_reference_associations("example_data/reference_associations.json")

graph = EvidenceGraph()
sample = graph.add_sample("demo-sample")
graph.add_detected_taxa(sample, taxa)
graph.add_detected_functions(sample, functions)
graph.add_reference_associations(references)
graph.link_detections_to_references()

predictions = BiotransformationPredictor(graph).predict(sample)
for p in predictions:
    print(p.reaction_class, p.score, p.confidence_label)

plot_evidence_graph(graph, sample_node=sample, save_path="graph.png")
plot_evidence_scores(predictions, save_path="scores.png")
```

Or open `biotransformation_explorer.ipynb` for the interactive version.

## Extending this

- **Matching logic**: see the table above. The most useful next step is
  filling `alt_ids` (InterPro / KO) for each gene reference so matches are
  by id rather than by name.
- **Scoring**: `BiotransformationPrediction.score` in `predict.py` is a
  simple, documented heuristic (Tier-1 dominates, Tier-2 caps lower,
  convergence gets a bonus) -- swap in a different weighting, or a proper
  statistical model, without touching the graph or the notebook.
- **More reference sources**: `fetchers.py` already includes ChEBI (compound
  structure) and Rhea (EC ↔ reaction) clients so you can auto-populate
  parts of the reference file from UniProt + Rhea instead of hand-curating
  it, for genes with a known EC number.
- **Substrate structure**: `ChEBIClient` returns SMILES for a compound;
  if you install `rdkit` you can add structural-alert matching (e.g. via
  SMARTS patterns for known biotransformation sites) as a third,
  independent line of evidence alongside Tier 1 / Tier 2.

## Requirements

Python >= 3.10; see `requirements.txt`. `rdkit` and `mgnipy` are optional:
the pipeline and the web viewer run without either of them installed.
