/* EXPOSE viewer front end: no dependencies, talks to the Flask API. */
"use strict";

const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];

const state = {
  source: "example",       // example | files | mgnify
  refSource: "example",    // example | file
  example: null,           // {taxa, functions, references}
  files: { taxa: null, functions: null },
  mgnify: { study: null, analyses: [], data: {} },   // data keyed by `${acc}|${tax}|${fun}`
  chembl: { hits: [], selected: new Set(), built: null, library: [] },
  biosiftr: { samples: {}, selected: new Set(), summary: null },
  enzymeTaxonomy: {},      // analysis accession -> {accession: attribution}
  multi: null,            // the whole /api/analyze/multi response
  activeSample: 0,
  gridView: "heatmap",
  linkedOnly: true,       // grids summarise the evidence, not the whole profile
  refFile: null,           // parsed reference records
  result: null,
  selectedPrediction: null,
  selectedNode: null,
  table: "taxa",
};

/* ------------------------------------------------------------------ utils */

function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}
function setStatus(el, msg, kind = "") {
  if (typeof el === "string") el = $(el);
  el.textContent = msg || "";
  el.className = "status" + (kind ? " " + kind : "");
}
async function api(path, opts = {}) {
  const res = await fetch(path, opts);
  let body = null;
  try { body = await res.json(); } catch { /* non-JSON */ }
  if (!res.ok) throw new Error((body && body.error) || `${res.status} ${res.statusText}`);
  return body;
}
function fmtNum(v) {
  if (v === null || v === undefined || v === "") return "";
  const n = Number(v);
  return Number.isFinite(n) ? n.toLocaleString() : esc(v);
}
function truncate(s, n) { s = String(s ?? ""); return s.length > n ? s.slice(0, n - 1) + "…" : s; }
function download(filename, text, type) {
  const blob = new Blob([text], { type });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

/* External links for identifiers, so users can check evidence at the source. */
function idLink(id) {
  if (!id) return "";
  const s = String(id).trim();
  const u = s.toUpperCase();
  let url = null;
  if (/^IPR\d+$/.test(u)) url = `https://www.ebi.ac.uk/interpro/entry/InterPro/${u}/`;
  else if (/^PF\d+$/.test(u)) url = `https://www.ebi.ac.uk/interpro/entry/pfam/${u}/`;
  else if (/^K\d{5}$/.test(u)) url = `https://www.kegg.jp/entry/${u}`;
  else if (/^GO:\d+$/.test(u)) url = `https://www.ebi.ac.uk/QuickGO/term/${u}`;
  else if (/^RHEA:\d+$/.test(u)) url = `https://www.rhea-db.org/rhea/${u.slice(5)}`;
  else if (/^CHEBI:\d+$/.test(u)) url = `https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${u}`;
  else if (/^(EC:)?\d+(\.(\d+|-)){1,3}$/.test(u)) url = `https://enzyme.expasy.org/EC/${u.replace(/^EC:/, "")}`;
  else if (/^MGYA\d+$/.test(u)) url = `https://www.ebi.ac.uk/metagenomics/analyses/${u}`;
  else if (/^MGYS\d+$/.test(u)) url = `https://www.ebi.ac.uk/metagenomics/studies/${u}`;
  return url ? `<a href="${url}" target="_blank" rel="noopener">${esc(s)}</a>` : esc(s);
}
/** A share as a percentage that stays informative when it is very small:
 *  a rare taxon is "<0.01%", not a misleading "0.00%". */
function sharePct(share) {
  const pct = share * 100;
  if (pct >= 10) return `${pct.toFixed(0)}%`;
  if (pct >= 1) return `${pct.toFixed(1)}%`;
  if (pct >= 0.01) return `${pct.toFixed(2)}%`;
  return "<0.01%";
}

function taxidLink(taxId) {
  const t = String(taxId ?? "").trim();
  if (!t) return "";
  return `<a href="https://www.ebi.ac.uk/ena/browser/view/Taxon:${encodeURIComponent(t)}"
    target="_blank" rel="noopener">taxon ${esc(t)}</a>`;
}
function taxonLink(name) {
  return `<a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term=${encodeURIComponent(name)}" target="_blank" rel="noopener">${esc(name)}</a>`;
}
function matchBadge(type) {
  const weak = type === "name" || type === "within-genus";
  const label = { id: "ID match", name: "name match", exact: "exact taxon", "within-genus": "within genus" }[type] || type;
  return `<span class="match${weak ? " weak" : ""}" title="How this detection was linked to the reference">${esc(label)}</span>`;
}

/* --------------------------------------------------------------- jobs */

/** Run a background job to completion.
 *  Slow, outage-prone work (ChEMBL, bulk MGnify downloads) returns 202 + a job
 *  id instead of blocking; poll it, report progress, and offer a cancel. */
async function runJob(request, { statusSel, label = "Working" } = {}) {
  const started = await request;
  if (!started || !started.id) return started;          // server answered directly
  const el = statusSel ? $(statusSel) : null;
  const render = job => {
    if (!el) return;
    const pct = job.progress >= 0 ? Math.round(job.progress * 100) : null;
    el.innerHTML = `<span class="job">
      <span>${esc(job.message || label)}</span>
      <span class="track${pct === null ? " indeterminate" : ""}"><span style="width:${pct ?? 100}%"></span></span>
      <button type="button" data-cancel="${esc(job.id)}">cancel</button></span>`;
    el.className = "status";
  };
  render(started);
  for (let waited = 0; ; waited += 400) {
    await new Promise(r => setTimeout(r, 400));
    let job;
    try {
      job = await api(`/api/jobs/${encodeURIComponent(started.id)}`);
    } catch (err) {
      throw new Error(`lost track of the job: ${err.message}`);
    }
    render(job);
    if (job.state === "done") return job.result;
    if (job.state === "failed") throw new Error(job.error || "job failed");
    if (job.state === "cancelled") throw new Error("cancelled");
  }
}

document.addEventListener("click", e => {
  const btn = e.target.closest("[data-cancel]");
  if (btn) api(`/api/jobs/${encodeURIComponent(btn.dataset.cancel)}/cancel`, { method: "POST" }).catch(() => {});
});

/* ---------------------------------------------------------- segmented UI */

function wireSegmented(group, onChange) {
  const root = $(`[data-group="${group}"]`);
  root.addEventListener("click", e => {
    const b = e.target.closest("button[data-value]");
    if (!b) return;
    $$("button", root).forEach(x => x.classList.toggle("active", x === b));
    onChange(b.dataset.value);
  });
}

wireSegmented("source", v => {
  state.source = v;
  $$("[data-pane]").forEach(p => (p.hidden = p.dataset.pane !== v));
  if (v === "mgnify") { loadDefaultStudy(); refreshMgnifyLibrary(); }
});
wireSegmented("ref", v => {
  state.refSource = v;
  $$("[data-pane-ref]").forEach(p => (p.hidden = p.dataset.paneRef !== v));
  if (v === "chembl") refreshLibrary({ load: !state.chembl.built });
  updateRefStatus();
});
wireSegmented("tables", v => { state.table = v; renderTable(); });

/* ------------------------------------------------------------- examples */

async function loadExample() {
  if (state.example) return state.example;
  state.example = await api("/api/examples");
  return state.example;
}

function describeRefs(records) {
  const g = records.filter(r => r.entity_type === "gene").length;
  const m = records.length - g;
  const rc = new Set(records.map(r => r.reaction_class)).size;
  return `${records.length} statements (${g} gene, ${m} microbe) across ${rc} reaction classes`;
}

async function updateRefStatus() {
  const el = $("#ref-status");
  if (state.refSource === "example") {
    try { setStatus(el, describeRefs((await loadExample()).references), "ok"); }
    catch (e) { setStatus(el, "Could not load example: " + e.message, "err"); }
  } else if (state.refSource === "chembl") {
    const built = state.chembl.built;
    built ? setStatus(el, describeRefs(built.records.filter(r => r.observed)), "ok")
          : setStatus(el, "No ChEMBL assays built yet.");
  } else {
    state.refFile ? setStatus(el, describeRefs(state.refFile), "ok") : setStatus(el, "No file loaded yet.");
  }
}

/* ------------------------------------------------------------ file inputs */

async function uploadParse(endpoint, file, extra = {}) {
  const fd = new FormData();
  fd.append("file", file);
  Object.entries(extra).forEach(([k, v]) => fd.append(k, v));
  return api(endpoint, { method: "POST", body: fd });
}

for (const kind of ["taxa", "functions"]) {
  $(`#file-${kind}`).addEventListener("change", async e => {
    const file = e.target.files[0];
    const status = $(`#file-${kind}-status`);
    state.files[kind] = null;
    if (!file) return setStatus(status, "");
    setStatus(status, "Reading…");
    try {
      const res = await uploadParse("/api/parse/detections", file, { kind });
      state.files[kind] = res.records;
      setStatus(status, `✓ ${res.records.length.toLocaleString()} records from ${file.name}`, "ok");
      if (!$("#file-sample-id").value) $("#file-sample-id").value = file.name.replace(/\.[^.]+$/, "");
    } catch (err) {
      setStatus(status, err.message, "err");
    }
  });
}

$("#file-ref").addEventListener("change", async e => {
  const file = e.target.files[0];
  state.refFile = null;
  if (!file) return updateRefStatus();
  setStatus("#ref-status", "Reading…");
  try {
    state.refFile = (await uploadParse("/api/parse/reference", file)).records;
    updateRefStatus();
  } catch (err) {
    setStatus("#ref-status", err.message, "err");
  }
});

/* ---------------------------------------------------------------- MGnify */

// A pipeline-v6 study known to work end to end; loaded when the MGnify tab is first opened.
const DEFAULT_STUDY = "MGYS00010462";
const MGNIFY_TAXONOMY = "contigs";
$("#mg-keyword").value = DEFAULT_STUDY;

async function loadDefaultStudy() {
  if ($("#mg-study").options.length) return;
  await searchStudies(DEFAULT_STUDY);
  const sel = $("#mg-study");
  if ([...sel.options].some(o => o.value === DEFAULT_STUDY)) {
    sel.value = DEFAULT_STUDY;
    sel.dispatchEvent(new Event("change"));
  }
}

$("#mg-search-form").addEventListener("submit", e => {
  e.preventDefault();
  const q = $("#mg-keyword").value.trim();
  if (q) searchStudies(q);
});

async function searchStudies(q) {
  const sel = $("#mg-study");
  sel.innerHTML = ""; $("#mg-analysis").innerHTML = "";
  state.mgnify = { study: null, analyses: [], data: {} };
  $("#mg-load").disabled = true;
  setStatus("#mg-status", "Searching MGnify…");
  try {
    const { studies } = await api(`/api/mgnify/studies?search=${encodeURIComponent(q)}`);
    studies.forEach(s => {
      const o = new Option(`${s.accession} · ${truncate(s.name || "", 60)}`, s.accession);
      o.title = [s.name, s.biome, (s.ena_accessions || []).join(", ")].filter(Boolean).join("\n");
      sel.add(o);
    });
    setStatus("#mg-status", studies.length ? `${studies.length} pipeline-v6 studies found.` : "No pipeline-v6 studies found.",
      studies.length ? "ok" : "");
  } catch (err) {
    setStatus("#mg-status", err.message, "err");
  }
}

$("#mg-study").addEventListener("change", async e => {
  const acc = e.target.value;
  state.mgnify = { study: acc, analyses: [], data: {} };
  const sel = $("#mg-analysis");
  sel.innerHTML = "";
  $("#mg-load").disabled = true;
  setStatus("#mg-status", `Listing analyses for ${acc}…`);
  try {
    const { analyses, downloaded } = await api(`/api/mgnify/studies/${encodeURIComponent(acc)}/analyses`);
    sel.innerHTML = analyses.map(a => {
      const title = [a.sample_title, a.run && `run ${a.run}`, a.assembly && `assembly ${a.assembly}`,
                     a.downloaded && "already downloaded"].filter(Boolean).join("\n");
      return `<label class="pick" title="${esc(title)}">
        <input type="checkbox" value="${esc(a.accession)}">
        <span><span class="nm">${esc(a.accession)}</span>${
          a.downloaded ? '<span class="flag" title="already downloaded — loads instantly">cached</span>' : ""}${
          a.enzyme_taxonomy ? '<span class="flag" title="enzyme-to-organism links already built">linked</span>' : ""}
          <span class="sub">${esc([a.experiment_type, a.pipeline_version, a.sample].filter(Boolean).join(" · "))}</span>
        </span></label>`;
    }).join("");
    sel.hidden = !analyses.length;
    $("#mg-analysis-tools").hidden = !analyses.length;
    updateAnalysisSelection();
    setStatus("#mg-status", analyses.length
      ? `${analyses.length} analyses in ${acc}${downloaded ? ` · ${downloaded} already downloaded` : ""}.`
      : `No analyses in ${acc} yet.`, analyses.length ? "ok" : "");
  } catch (err) {
    setStatus("#mg-status", err.message, "err");
  }
});

function updateAnalysisSelection() {
  state.mgnify.analyses = [...$$("#mg-analysis input:checked")].map(b => b.value);
  const n = state.mgnify.analyses.length;
  $("#mg-load").disabled = !n;
  $("#mg-analysis-count").textContent = n ? `${n} selected` : "none selected";
  $("#mg-load").textContent = n > 1 ? `Load ${n} analyses` : "Load analysis";
}

$("#mg-analysis").addEventListener("change", e => {
  if (e.target.type === "checkbox") { updateAnalysisSelection(); setStatus("#mg-load-status", ""); }
});
$("#mg-analysis-tools").addEventListener("click", e => {
  const how = e.target.closest("[data-pick]")?.dataset.pick;
  if (!how) return;
  $$("#mg-analysis input[type=checkbox]").forEach(b => (b.checked = how === "all"));
  updateAnalysisSelection();
  setStatus("#mg-load-status", "");
});

/** Download every selected analysis as one background job. Each analysis is
 *  cached server-side, so re-running a comparison is instant and a failure
 *  part-way through keeps whatever already arrived. */
async function loadMgnifyAnalyses() {
  const accs = state.mgnify.analyses;
  if (!accs.length) throw new Error("Tick at least one analysis first.");
  // contig taxonomy only: it is the assignment that ties each enzyme to the
  // organism carrying it, and v6 assemblies have no rRNA profile anyway
  const tax = MGNIFY_TAXONOMY, fun = $("#mg-func-type").value;
  const need = accs.filter(a => !state.mgnify.data[`${a}|${tax}|${fun}`]);

  if (need.length) {
    const out = await runJob(
      api("/api/mgnify/analyses", {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ accessions: need, taxonomy: tax, functional: fun }),
      }),
      { statusSel: "#mg-load-status", label: "Downloading from MGnify" });
    (out.analyses || []).forEach(d => { state.mgnify.data[`${d.accession}|${tax}|${fun}`] = d; });
    state.mgnify.failed = out.failed || [];
  }

  const loaded = accs.map(a => state.mgnify.data[`${a}|${tax}|${fun}`]).filter(Boolean);
  if (!loaded.length) throw new Error((state.mgnify.failed || []).join(" ") || "Nothing loaded.");
  const taxa = loaded.reduce((n, d) => n + d.taxa.length, 0);
  const funcs = loaded.reduce((n, d) => n + d.functions.length, 0);
  const warn = [...(state.mgnify.failed || []), ...loaded.flatMap(d => d.warnings || [])].join(" ");
  setStatus("#mg-load-status",
    `✓ ${loaded.length} analys${loaded.length === 1 ? "is" : "es"}: ` +
    `${taxa.toLocaleString()} taxa, ${funcs.toLocaleString()} functions.${warn ? " " + warn : ""}`,
    warn ? "" : "ok");
  refreshMgnifyLibrary();
  return loaded;
}
$("#mg-load").addEventListener("click", () => loadMgnifyAnalyses().catch(err => setStatus("#mg-load-status", err.message, "err")));

/* The downloaded-analysis library: what is already on disk, so a comparison
   can be re-run with MGnify unreachable. */

async function refreshMgnifyLibrary() {
  let lib;
  try {
    lib = await api("/api/mgnify/library");
  } catch { return; }
  const box = $("#mg-library"), list = $("#mg-library-list");
  box.hidden = !lib.entries.length;
  if (!lib.entries.length) return;
  $("#mg-library-count").textContent =
    `${lib.entries.length} cached${lib.persistent ? "" : " · not saved (no cache)"}`;
  list.innerHTML = lib.entries.map(e => `
    <label class="lib-row" title="${esc((e.warnings || []).join("\n"))}">
      <span class="dot-spacer"></span>
      <span>
        <span class="nm">${esc(e.accession)}</span>
        <span class="sub">${esc(e.taxonomy)} · ${esc(e.functional)}</span>
        <span class="sub">${fmtNum(e.taxa)} taxa · ${fmtNum(e.functions)} functions</span>
      </span>
      <button type="button" class="drop" data-drop="${esc(e.accession)}" title="Remove from the cache">×</button>
    </label>`).join("");
}

$("#mg-library-list").addEventListener("click", async e => {
  const accession = e.target.closest("[data-drop]")?.dataset.drop;
  if (!accession) return;
  e.preventDefault();
  try {
    await api(`/api/mgnify/library/${encodeURIComponent(accession)}`, { method: "DELETE" });
    Object.keys(state.mgnify.data)
      .filter(k => k.split("|")[0] === accession)
      .forEach(k => delete state.mgnify.data[k]);
    await refreshMgnifyLibrary();
  } catch (err) {
    setStatus("#mg-load-status", err.message, "err");
  }
});

/* ---------------------------------------------------------------- ChEMBL */

$("#cb-search-form").addEventListener("submit", async e => {
  e.preventDefault();
  const q = $("#cb-keyword").value.trim();
  if (!q) return;
  const box = $("#cb-results");
  box.hidden = true; box.innerHTML = "";
  state.chembl.selected.clear();
  $("#cb-build").disabled = true;
  setStatus("#cb-status", "Searching ChEMBL…");
  try {
    const { assays, in_library } = await runJob(api(`/api/chembl/assays?search=${encodeURIComponent(q)}`),
      { statusSel: "#cb-status", label: "Searching ChEMBL" });
    state.chembl.hits = assays;
    if (!assays.length) { setStatus("#cb-status", "No biotransformation assays found."); return; }
    box.innerHTML = assays.map(a => `
      <label class="pick">
        <input type="checkbox" value="${esc(a.assay_chembl_id)}">
        <span>
          <span class="nm">${esc(a.assay_chembl_id)}</span>
          <span class="tierbadge ${esc(a.tier)}">${esc(a.tier)}</span>${
          a.in_library ? '<span class="flag" title="already built — no fetching needed">cached</span>' : ""}
          <span class="sub">${esc(a.assay_organism || "")}${a.activities ? ` · ${a.activities} activities` : ""}</span>
          <span class="sub">${esc(truncate(a.description || "", 110))}</span>
        </span>
      </label>`).join("");
    box.hidden = false;
    setStatus("#cb-status",
      `${assays.length} assay(s) found${in_library ? ` · ${in_library} already built` : ""}.`, "ok");
  } catch (err) {
    setStatus("#cb-status", err.message, "err");
  }
});

$("#cb-results").addEventListener("change", e => {
  if (e.target.type !== "checkbox") return;
  const sel = state.chembl.selected;
  e.target.checked ? sel.add(e.target.value) : sel.delete(e.target.value);
  $("#cb-build").disabled = !sel.size;
  // say up front how much of the selection needs fetching: requests are paced
  // at one a second, so a large new set takes a while
  const cached = state.chembl.hits.filter(h => h.in_library && sel.has(h.assay_chembl_id)).length;
  const toFetch = sel.size - cached;
  $("#cb-build").textContent = !sel.size ? "Build references"
    : toFetch ? `Build ${sel.size} (fetch ${toFetch}, ~${Math.max(1, Math.round(toFetch * 6 / 60))} min)`
              : `Build ${sel.size} (all cached)`;
});

/* The library: assays already built, kept on disk between runs so a set of
   references can be toggled on and off without rebuilding or re-fetching. */

async function refreshLibrary({ load = false } = {}) {
  let lib;
  try {
    lib = await api("/api/chembl/library");
  } catch { return; }
  state.chembl.library = lib.entries;
  const box = $("#cb-library"), list = $("#cb-library-list");
  box.hidden = !lib.entries.length;
  if (!lib.entries.length) { state.chembl.built = null; updateRefStatus(); return; }

  const on = lib.entries.filter(e => e.selected).length;
  $("#cb-library-count").textContent =
    `${on} of ${lib.entries.length} on${lib.persistent ? "" : " · not saved (no cache)"}`;
  list.innerHTML = lib.entries.map(e => `
    <label class="lib-row" title="${esc((e.warnings || []).join("\n"))}">
      <input type="checkbox" data-assay="${esc(e.assay_chembl_id)}"${e.selected ? " checked" : ""}>
      <span>
        <span class="nm">${esc(e.assay_chembl_id)}</span>
        <span class="tierbadge ${esc(e.tier || "")}">${esc(e.tier || "")}</span>
        <span class="sub">${esc(e.organism || "")}</span>
        <span class="sub">${fmtNum(e.positive)} positive · ${fmtNum(e.negative)} negative ·
          ${fmtNum(e.reaction_classes)} class(es)${e.communities ? ` · ${fmtNum(e.communities)} community` : ""}</span>
      </span>
      <button type="button" class="drop" data-drop="${esc(e.assay_chembl_id)}" title="Remove from the library">×</button>
    </label>`).join("");
  if (load) await loadSelectedLibrary();
}

async function loadSelectedLibrary() {
  try {
    const refs = await api("/api/chembl/library/references");
    state.chembl.built = refs.assay_ids.length ? refs : null;
  } catch {
    state.chembl.built = null;
  }
  renderBuildSummary();
  updateRefStatus();
}

function renderBuildSummary() {
  const out = state.chembl.built;
  const el = $("#cb-summary");
  if (!out) { el.hidden = true; return; }
  const s = out.summary;
  el.hidden = false;
  el.innerHTML = `
    <div><b>${fmtNum(s.positive)}</b> positive association(s) across <b>${fmtNum(s.reaction_classes)}</b> reaction class(es)
      — ${fmtNum(s.genes)} gene, ${fmtNum(s.microbes)} microbe records; <b>${fmtNum(s.negative)}</b> negative kept but not predicted from.</div>
    ${s.communities ? `<div><b>${fmtNum(s.communities)}</b> measured community assay(s), ${fmtNum(s.community_observations)} observation(s).</div>` : ""}
    ${(out.warnings || []).length ? `<div class="muted">${out.warnings.slice(0, 4).map(esc).join("<br>")}</div>` : ""}`;
}

$("#cb-library-list").addEventListener("change", async e => {
  if (e.target.type !== "checkbox") return;
  const ids = [...$$("#cb-library-list input:checked")].map(b => b.dataset.assay);
  try {
    await api("/api/chembl/library/selection", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ assay_ids: ids }),
    });
    // update in place rather than re-rendering: rebuilding the list would
    // detach the checkbox that was just clicked and lose the scroll position
    state.chembl.library.forEach(entry => (entry.selected = ids.includes(entry.assay_chembl_id)));
    $("#cb-library-count").textContent = `${ids.length} of ${state.chembl.library.length} on`;
    await loadSelectedLibrary();
  } catch (err) {
    setStatus("#cb-build-status", err.message, "err");
  }
});

$("#cb-library-list").addEventListener("click", async e => {
  const id = e.target.closest("[data-drop]")?.dataset.drop;
  if (!id) return;
  e.preventDefault();
  try {
    await api(`/api/chembl/library/${encodeURIComponent(id)}`, { method: "DELETE" });
    await refreshLibrary({ load: true });
  } catch (err) {
    setStatus("#cb-build-status", err.message, "err");
  }
});

/** Cross-match an uploaded id list against the library, then fetch only what
 *  is missing -- the point being to not re-download what is already held. */
$("#cb-idfile").addEventListener("change", async e => {
  const file = e.target.files?.[0];
  if (!file) return;
  setStatus("#cb-idfile-status", `Reading ${file.name}…`);
  try {
    const fd = new FormData();
    fd.append("file", file);
    const match = await api("/api/chembl/library/match", { method: "POST", body: fd });
    const { supplied, present, missing } = match;

    if (!missing.length) {
      await api("/api/chembl/library/selection", {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ assay_ids: supplied }),
      });
      await refreshLibrary({ load: true });
      setStatus("#cb-idfile-status",
        `✓ all ${supplied.length} already in the library; selected.`, "ok");
      return;
    }

    setStatus("#cb-idfile-status",
      `${supplied.length} id(s): ${present.length} cached, fetching ${missing.length}…`);
    await runJob(api("/api/chembl/references", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ assay_ids: missing }),
    }), { statusSel: "#cb-idfile-status", label: `Fetching ${missing.length} new assay(s)` });

    await api("/api/chembl/library/selection", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ assay_ids: supplied }),
    });
    await refreshLibrary({ load: true });
    setStatus("#cb-idfile-status",
      `✓ ${supplied.length} selected (${present.length} were cached, ${missing.length} fetched).`, "ok");
  } catch (err) {
    setStatus("#cb-idfile-status", err.message, "err");
  } finally {
    e.target.value = "";
  }
});

$("#cb-build").addEventListener("click", async () => {
  const ids = [...state.chembl.selected];
  if (!ids.length) return;
  $("#cb-build").disabled = true;
  setStatus("#cb-build-status", `Building references from ${ids.length} assay(s)…`);
  try {
    const out = await runJob(api("/api/chembl/references", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ assay_ids: ids }),
    }), { statusSel: "#cb-build-status", label: "Building references" });
    state.chembl.built = out;
    renderBuildSummary();
    await refreshLibrary({ load: true });
    setStatus("#cb-build-status", `✓ ${out.summary.positive} association(s) ready.`, "ok");
  } catch (err) {
    setStatus("#cb-build-status", err.message, "err");
  } finally {
    $("#cb-build").disabled = !state.chembl.selected.size;
  }
});

/* ------------------------------------------------------------- bioSIFTR */

$("#bs-scan").addEventListener("click", async () => {
  const path = $("#bs-path").value.trim();
  if (!path) { setStatus("#bs-status", "Give the path to a bioSIFTR output directory.", "err"); return; }
  const box = $("#bs-samples");
  box.hidden = true; box.innerHTML = "";
  $("#bs-tools").hidden = true; $("#bs-taxonomy").hidden = true;
  state.biosiftr = { samples: {}, selected: new Set(), summary: null };
  setStatus("#bs-status", "Reading…");
  try {
    const out = await runJob(api("/api/biosiftr/scan", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ path, mapper: $("#bs-mapper").value,
                             resolve_ncbi: $("#bs-resolve").value === "1" }),
    }), { statusSel: "#bs-status", label: "Reading bioSIFTR output" });

    state.biosiftr.samples = out.samples;
    state.biosiftr.summary = out.summary;
    box.innerHTML = out.summary.samples.map(s => `
      <label class="pick">
        <input type="checkbox" value="${esc(s.sample_id)}" checked>
        <span><span class="nm">${esc(s.sample_id)}</span>
          <span class="sub">${fmtNum(s.taxa)} taxa · ${fmtNum(s.functions)} functions</span>
        </span></label>`).join("");
    box.hidden = false; $("#bs-tools").hidden = false;
    updateBiosiftrSelection();

    const tax = out.summary.taxonomy;
    if (tax) {
      const ranks = Object.entries(tax.by_rank || {}).map(([r, n]) => `${n} at ${r}`).join(", ");
      $("#bs-taxonomy").hidden = false;
      $("#bs-taxonomy").innerHTML = `
        <div><b>${fmtNum(tax.resolved)}</b> taxa resolved to NCBI${ranks ? ` (${esc(ranks)})` : ""},
          <b>${fmtNum(tax.unresolved)}</b> unresolved.</div>
        ${tax.unresolved_examples?.length ? `<div class="muted">GTDB-only, e.g.
          ${tax.unresolved_examples.slice(0, 4).map(e => esc(String(e).split(";").pop())).join(", ")}</div>` : ""}`;
    }
    const warn = (out.summary.warnings || []).join(" ");
    setStatus("#bs-status", `✓ ${out.summary.samples.length} sample(s) found.${warn ? " " + warn : ""}`,
      warn ? "" : "ok");
  } catch (err) {
    setStatus("#bs-status", err.message, "err");
  }
});

function updateBiosiftrSelection() {
  state.biosiftr.selected = new Set([...$$("#bs-samples input:checked")].map(b => b.value));
  $("#bs-count").textContent = `${state.biosiftr.selected.size} selected`;
}
$("#bs-samples").addEventListener("change", e => {
  if (e.target.type === "checkbox") updateBiosiftrSelection();
});
$("#bs-tools").addEventListener("click", e => {
  const how = e.target.closest("[data-pick]")?.dataset.pick;
  if (!how) return;
  $$("#bs-samples input[type=checkbox]").forEach(b => (b.checked = how === "all"));
  updateBiosiftrSelection();
});

/* ------------------------------------------------------------------- run */

async function gatherInputs() {
  let samples = [];
  if (state.source === "example") {
    const ex = await loadExample();
    samples = [{ sample_id: ex.sample_id, taxa: ex.taxa, functions: ex.functions,
                 sample_meta: { origin: "bundled example" } }];
  } else if (state.source === "files") {
    const taxa = state.files.taxa || [], functions = state.files.functions || [];
    if (!taxa.length && !functions.length) throw new Error("Load a taxonomy and/or functional annotation file first.");
    samples = [{ sample_id: $("#file-sample-id").value.trim() || "uploaded-sample",
                 taxa, functions, sample_meta: { origin: "uploaded files" } }];
  } else if (state.source === "biosiftr") {
    const picked = [...state.biosiftr.selected];
    if (!picked.length) throw new Error("Scan a bioSIFTR directory and tick at least one sample.");
    samples = picked.map(id => ({
      sample_id: id,
      taxa: state.biosiftr.samples[id].taxa,
      functions: state.biosiftr.samples[id].functions,
      sample_meta: { origin: "bioSIFTR", outdir: state.biosiftr.summary?.outdir,
                     mapper: state.biosiftr.summary?.mapper },
    }));
  } else {
    const loaded = await loadMgnifyAnalyses();
    samples = loaded.map(d => ({
      sample_id: d.accession, taxa: d.taxa, functions: d.functions,
      sample_meta: { origin: "MGnify", study: state.mgnify.study, analysis: d.accession },
    }));
  }

  let references, observations = [];
  if (state.refSource === "example") {
    references = (await loadExample()).references;
  } else if (state.refSource === "chembl") {
    if (!state.chembl.built) throw new Error("Search ChEMBL and build references first.");
    // only positives are predicted from; negatives are kept in the build output
    references = state.chembl.built.records.filter(r => r.observed).map(({ observed, ...rest }) => rest);
    observations = state.chembl.built.observations || [];
    if (!references.length) throw new Error("Those assays produced no positive associations to predict from.");
  } else {
    if (!state.refFile) throw new Error("Load a reference file first.");
    references = state.refFile;
  }
  let enzyme_taxonomy = null;
  if ($("#opt-enzyme-taxonomy").checked) {
    const mgnifySamples = samples.filter(s => (s.sample_meta || {}).origin === "MGnify");
    if (!mgnifySamples.length) {
      setStatus("#enzyme-tax-status",
        "Enzyme attribution needs MGnify assembly analyses; skipped for this source.", "err");
    } else {
      const need = mgnifySamples.map(s => s.sample_id).filter(id => !state.enzymeTaxonomy[id]);
      if (need.length) {
        const out = await runJob(api("/api/mgnify/enzyme-taxonomy", {
          method: "POST", headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ accessions: need }),
        }), { statusSel: "#enzyme-tax-status", label: "Linking enzymes to organisms" });
        Object.entries(out.analyses || {}).forEach(([id, v]) => {
          state.enzymeTaxonomy[id] = v.attribution;
        });
        const failed = (out.failed || []).join(" ");
        setStatus("#enzyme-tax-status",
          `✓ linked ${Object.keys(out.analyses || {}).length} analysis(es).${failed ? " " + failed : ""}`,
          failed ? "" : "ok");
      }
      enzyme_taxonomy = {};
      mgnifySamples.forEach(s => {
        if (state.enzymeTaxonomy[s.sample_id]) enzyme_taxonomy[s.sample_id] = state.enzymeTaxonomy[s.sample_id];
      });
    }
  }

  return { samples, references, observations, enzyme_taxonomy,
           linked_only: state.linkedOnly,
           require_enzyme_attribution: $("#opt-require-attribution").checked,
           allow_name_matching: $("#opt-name-matching").checked,
           taxon_rank_floor: $("#opt-taxon-rank").value };
}

async function run() {
  const btn = $("#run");
  btn.disabled = true;
  setStatus("#run-status", "Running…");
  try {
    const body = await gatherInputs();
    const multi = await api("/api/analyze/multi", {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body),
    });
    state.multi = multi;
    setActiveSample(0, { render: false });
    renderResults();
    const n = multi.samples.length;
    setStatus("#run-status",
      n > 1 ? `Done: ${n} analyses, ${multi.summary.reaction_classes} reaction class(es) in total.`
            : `Done: ${multi.samples[0].predictions.length} predicted reaction class(es).`, "ok");
  } catch (err) {
    setStatus("#run-status", err.message, "err");
  } finally {
    btn.disabled = false;
  }
}
$("#run").addEventListener("click", run);

/* ------------------------------------------------------------- rendering */

function setActiveSample(i, { render = true } = {}) {
  const samples = state.multi?.samples || [];
  state.activeSample = Math.max(0, Math.min(i, samples.length - 1));
  state.result = samples[state.activeSample] || null;
  state.selectedPrediction = state.result?.predictions[0]?.reaction_class ?? null;
  state.selectedNode = null;
  if (render) renderResults();
}

$("#sample-select").addEventListener("change", e => setActiveSample(Number(e.target.value)));
wireSegmented("gridview", v => { state.gridView = v; renderGrid(); });

$("#opt-enzyme-taxonomy").addEventListener("change", e => {
  $("#opt-require-wrap").hidden = !e.target.checked;
  if (!e.target.checked) $("#opt-require-attribution").checked = false;
});

/* -- the comparison grid ------------------------------------------------- */

/** Five-step sequential ramp; `null` score (observed only, nothing predicted)
 *  gets the empty cell. Buckets, not a continuous scale, so the legend can
 *  name each step and every step's label ink clears 4.5:1. */
function scoreBucket(score) {
  if (score == null) return 0;
  return Math.min(5, Math.floor(score * 5) + 1);
}

function renderGrid() {
  const grid = state.multi?.grid;
  const multiple = (state.multi?.samples.length || 0) > 1;
  $("#grid-card").hidden = !multiple;
  $("#sample-switch").hidden = !multiple;
  if (!multiple || !grid) return;

  $("#sample-select").innerHTML = state.multi.samples
    .map((r, i) => `<option value="${i}"${i === state.activeSample ? " selected" : ""}>${esc(r.sample_id)}</option>`).join("");

  const taxa = state.multi.taxon_grid;
  const views = {
    heatmap: { title: "Predicted biotransformations across analyses",
               html: () => heatmapHTML(grid), legend: () => legendHTML(grid) },
    taxa: { title: "Microbes across analyses",
            html: () => taxonHeatmapHTML(taxa), legend: () => taxonLegendHTML(taxa) },
    enzymes: { title: "Enzymes across analyses",
               html: () => enzymeHeatmapHTML(state.multi.enzyme_grid),
               legend: () => enzymeLegendHTML(state.multi.enzyme_grid) },
    multiples: { title: "Predicted biotransformations across analyses",
                 html: () => multiplesHTML(), legend: () => "" },
  };
  const view = views[state.gridView] || views.heatmap;
  $("#grid-title").textContent = view.title;
  $("#grid-wrap").innerHTML = view.html();
  $("#grid-legend").innerHTML = view.legend();
}

/** Microbes x analyses.
 *  Abundance is not comparable between samples as reported (fractions from
 *  bioSIFTR, read counts from MGnify, different depths), so cells are coloured
 *  by each taxon's share of its own sample, scaled against the grid's largest
 *  share -- otherwise every cell lands in the lightest bucket, since a single
 *  organism rarely exceeds a few per cent of a community. */
function taxonHeatmapHTML(grid) {
  if (!grid || !grid.rows.length) {
    return `<p class="muted">${grid && grid.linked_only && grid.total_taxa
      ? `None of the ${fmtNum(grid.total_taxa)} detected taxa matched a reference.`
      : "No taxa were detected."}</p>`;
  }
  const scale = grid.max_share || 1;
  const head = grid.sample_ids.map(s => `<th class="col">${idLink(s)}</th>`).join("");

  const rows = grid.rows.map(row => {
    const cells = row.cells.map(c => {
      if (!c) return `<td><div class="cell blank">·</div></td>`;
      const share = c.share ?? 0;
      const b = share > 0 ? Math.max(1, Math.min(5, Math.ceil((share / scale) * 5))) : 0;
      const pct = share > 0 ? sharePct(share) : "–";
      const title = [`${row.taxon} — ${c.sample_id}`,
                     `abundance ${fmtNum(c.abundance)}`,
                     share ? `${sharePct(share)} of this sample` : "",
                     c.linked ? "linked to a reference association" : ""].filter(Boolean).join("\n");
      return `<td><div class="cell s${b}" title="${esc(title)}">${pct}${
        c.linked ? '<span class="dot" title="linked to a reference">●</span>' : ""}</div></td>`;
    }).join("");
    const label = row.gtdb_organism && row.gtdb_organism !== row.taxon
      ? `${row.taxon}` : row.taxon;
    const title = [row.lineage, row.gtdb_organism && `profiler reported: ${row.gtdb_organism}`,
                   row.rank && `matched at ${row.rank}`].filter(Boolean).join("\n");
    const matched = (row.linked_references || [])
      .map(r => String(r).split(":").slice(-1)[0]).filter(Boolean);
    const matchedHTML = matched.length
      ? matched.slice(0, 2).map(m => `<span class="carrier">${esc(m)}</span>`).join(" ") +
        (matched.length > 2 ? ` <span class="muted">+${matched.length - 2}</span>` : "")
      : `<span class="muted">—</span>`;
    return `<tr><th class="row" title="${esc(title)}">${esc(label)}
      ${row.tax_id ? `<span class="taxid">${taxidLink(row.tax_id)}${
        row.rank && row.rank !== "species" ? ` · ${esc(row.rank)}` : ""}</span>` : ""}</th>
      <td class="organism">${matchedHTML}</td>${cells}</tr>`;
  }).join("");

  return `<table class="grid"><thead><tr><th class="row">Microbe</th>
            <th class="organism">Predicts</th>${head}</tr></thead>
          <tbody>${rows}</tbody></table>`;
}

/** Enzymes x analyses. Same share-based colouring as the microbe grid; the
 *  extra signal is *how* each enzyme matched a reference, because a name-only
 *  hit is much weaker than an id hit and should be visible as such. */
function enzymeHeatmapHTML(grid) {
  if (!grid || !grid.rows.length) {
    return `<p class="muted">${grid && grid.linked_only && grid.total_enzymes
      ? `None of the ${fmtNum(grid.total_enzymes)} detected annotations matched a reference.`
      : "No functional annotations were detected."}</p>`;
  }
  const scale = grid.max_share || 1;
  const head = grid.sample_ids.map(s => `<th class="col">${idLink(s)}</th>`).join("");

  const rows = grid.rows.map(row => {
    const cells = row.cells.map(c => {
      if (!c) return `<td><div class="cell blank">·</div></td>`;
      const share = c.share ?? 0;
      const b = share > 0 ? Math.max(1, Math.min(5, Math.ceil((share / scale) * 5))) : 0;
      // a sample has tens of thousands of annotations, so one entry's share
      // rounds to 0.00% and says nothing; the count is the readable number
      const value = c.abundance != null ? fmtNum(c.abundance) : "–";
      const title = [`${row.annotation_id || row.description} — ${c.sample_id}`,
                     `abundance ${fmtNum(c.abundance)}`,
                     share ? `${sharePct(share)} of this sample's annotations` : "",
                     c.linked ? "linked to a gene reference" : ""].filter(Boolean).join("\n");
      return `<td><div class="cell s${b}" title="${esc(title)}">${value}${
        c.linked ? '<span class="dot" title="linked to a reference">●</span>' : ""}</div></td>`;
    }).join("");
    const flag = row.name_only
      ? ` <span class="badge plain" title="matched only by enzyme name, not by id">name only</span>` : "";
    return `<tr><th class="row" title="${esc(row.description || "")}">
      ${row.annotation_id ? idLink(row.annotation_id) : esc(row.description || "")}${flag}
      ${row.description && row.annotation_id
        ? `<span class="taxid">${esc(truncate(row.description, 44))}</span>` : ""}</th>
      <td class="organism">${carriersHTML(row)}</td>${cells}</tr>`;
  }).join("");

  return `<table class="grid"><thead><tr><th class="row">Enzyme / family</th>
            <th class="organism">Found in</th>${head}</tr></thead>
          <tbody>${rows}</tbody></table>`;
}

/** Which organisms carry this enzyme. With contig taxonomy loaded this is the
 *  whole point of the row: a linked enzyme with no organism behind it is much
 *  weaker evidence than one placed on a named species. */
function carriersHTML(row) {
  const taxa = row.taxa || [];
  if (taxa.length) {
    const shown = taxa.slice(0, 3).map(t =>
      `<span class="carrier" title="${esc(`${t.count} occurrence(s) on ${t.rank || "?"}-level contigs`)}">${
        esc(t.scientific_name)}${t.count > 1 ? ` <em>×${t.count}</em>` : ""}</span>`).join(" ");
    const more = taxa.length > 3 ? ` <span class="muted">+${taxa.length - 3}</span>` : "";
    const placed = row.attributed + row.unattributed;
    return `<span title="${esc(`${row.attributed} of ${placed} occurrences placed at genus/species`)}">${shown}${more}</span>`;
  }
  if (row.unattributable) {
    return `<span class="unattributed" title="every occurrence sat on a contig classified no further than domain">no organism below genus</span>`;
  }
  return `<span class="muted" title="enable 'Attribute enzymes to organisms' under Options">not linked to contigs</span>`;
}

function enzymeLegendHTML(grid) {
  const ramp = [1, 2, 3, 4, 5].map(i => `<i style="background:var(--h${i})"></i>`).join("");
  const max = grid.max_share ? `${(grid.max_share * 100).toFixed(1)}%` : "";
  return `
    <span class="key">counts, shaded by share of annotations
      <span class="ramp"><i style="background:var(--h0)"></i>${ramp}</span> 0 → ${esc(max)}</span>
    <span class="key">● linked to a reference (${fmtNum(grid.linked_enzymes)} of ${fmtNum(grid.total_enzymes)})</span>
    ${grid.name_only ? `<span class="key">${fmtNum(grid.name_only)} matched by name only</span>` : ""}
    ${grid.with_taxonomy ? `<span class="key">${fmtNum(grid.with_taxonomy)} placed on a named organism</span>` : ""}
    ${grid.unattributable ? `<span class="key">${fmtNum(grid.unattributable)} on unclassified contigs</span>` : ""}
    <span class="key">· not detected</span>
    ${grid.truncated ? `<span class="key">${fmtNum(grid.truncated)} rarer entries not shown</span>` : ""}
    ${grid.hidden_unlinked ? `<label class="key check small"><input type="checkbox" data-showall="1">
      show ${fmtNum(grid.hidden_unlinked)} unlinked</label>` : ""}`;
}

function taxonLegendHTML(grid) {
  const ramp = [1, 2, 3, 4, 5].map(i => `<i style="background:var(--h${i})"></i>`).join("");
  const max = grid.max_share ? `${(grid.max_share * 100).toFixed(1)}%` : "";
  return `
    <span class="key">share of sample <span class="ramp"><i style="background:var(--h0)"></i>${ramp}</span>
      0 → ${esc(max)}</span>
    <span class="key">● linked to a reference (${fmtNum(grid.linked_taxa)} of ${fmtNum(grid.total_taxa)} detected)</span>
    <span class="key">· not detected</span>
    ${grid.truncated ? `<span class="key">${fmtNum(grid.truncated)} rarer taxa not shown</span>` : ""}
    ${grid.hidden_unlinked ? `<label class="key check small"><input type="checkbox" data-showall="1">
      show ${fmtNum(grid.hidden_unlinked)} unlinked</label>` : ""}`;
}

/** Reaction classes built from one source share a long prefix
 *  ("Biotransformation of ..."), which wraps every row label onto three lines
 *  and strands the cells. Lift the shared prefix into the column header. */
function commonPrefix(labels) {
  if (labels.length < 2) return "";
  const words = labels.map(l => l.split(" "));
  let n = 0;
  while (words[0][n] !== undefined && words.every(w => w[n] === words[0][n]) &&
         words.some(w => w.length > n + 1)) n++;
  const prefix = words[0].slice(0, n).join(" ");
  return prefix.length >= 8 ? prefix : "";
}

function heatmapHTML(grid) {
  const prefix = commonPrefix(grid.rows.map(r => r.reaction_class));
  const head = grid.samples.map(s =>
    `<th class="col" title="${esc(s.sample_id)}${s.has_observations ? " · has measured outcomes" : ""}">
       ${idLink(s.sample_id)}${s.has_observations ? ' <span class="obs" title="measured outcomes available">◆</span>' : ""}
     </th>`).join("");

  const rows = grid.rows.map(row => {
    const cells = row.cells.map(c => {
      if (!c) return `<td><div class="cell blank">·</div></td>`;
      const b = scoreBucket(c.score);
      const tier = c.tier === 1 ? "▲" : c.tier === 2 ? "○" : "";
      const obs = c.observed === true ? '<span class="obs" title="observed in vitro">✓</span>'
                : c.observed === false ? '<span class="obs" title="not observed in vitro">✗</span>' : "";
      const title = [
        `${row.reaction_class} — ${c.sample_id}`,
        c.score != null ? `score ${c.score.toFixed(2)} (${c.confidence_label || ""})` : "not predicted",
        c.tier1_hits ? `${c.tier1_hits} enzyme reference(s)` : "",
        c.tier2_hits ? `${c.tier2_hits} microbe reference(s)` : "",
        c.observed === true ? "observed in vitro" : c.observed === false ? "not observed in vitro" : "",
      ].filter(Boolean).join("\n");
      return `<td><div class="cell s${b}" tabindex="0" role="button"
        data-rc="${esc(row.reaction_class)}" data-sample="${esc(c.sample_id)}" title="${esc(title)}">
        ${c.score != null ? c.score.toFixed(2) : "–"}<span class="tier">${tier}</span>${obs}</div></td>`;
    }).join("");
    const label = prefix ? row.reaction_class.slice(prefix.length).trim() : row.reaction_class;
    return `<tr><th class="row" title="${esc(row.reaction_class)}">${esc(label)}
      ${row.substrate_chebi ? `<span class="sub">${idLink(row.substrate_chebi)}</span>` : ""}</th>${cells}</tr>`;
  }).join("");

  return `<table class="grid"><thead><tr>
            <th class="row">${prefix ? `${esc(prefix)}…` : "Reaction class"}</th>${head}</tr></thead>
          <tbody>${rows}</tbody></table>`;
}

function legendHTML(grid) {
  const ramp = [1, 2, 3, 4, 5].map(i =>
    `<i class="s${i}" style="background:var(--h${i})"></i>`).join("");
  return `
    <span class="key">confidence <span class="ramp"><i style="background:var(--h0)"></i>${ramp}</span> 0 → 1</span>
    <span class="key">▲ enzyme detected (Tier 1)</span>
    <span class="key">○ microbe only (Tier 2)</span>
    <span class="key">· not predicted</span>
    ${grid.observed_samples.length ? `<span class="key">✓ / ✗ measured in vitro</span>` : ""}`;
}

function multiplesHTML() {
  return `<div class="multiples">${state.multi.samples.map(r => {
    const top = r.predictions.slice(0, 8);
    const meta = r.sample_meta || {};
    return `<div class="multiple">
      <h4>${idLink(r.sample_id)}</h4>
      <div class="sub">${r.predictions.length} class(es) · ${r.summary.tier1_predictions} with enzyme evidence
        ${meta.study ? ` · ${esc(meta.study)}` : ""}</div>
      <ol>${top.length ? top.map(p => `
        <li data-rc="${esc(p.reaction_class)}" data-sample="${esc(r.sample_id)}">
          <span class="nm" title="${esc(p.reaction_class)}">${esc(p.reaction_class)}</span>
          <span>${p.score.toFixed(2)}</span>
          <span class="bar"><span class="${p.tier1_hits.length ? "t1" : "t2"}" style="width:${Math.round(p.score * 100)}%"></span></span>
        </li>`).join("") : `<li class="muted">No predictions.</li>`}</ol>
      ${r.predictions.length > top.length ? `<div class="sub">+${r.predictions.length - top.length} more</div>` : ""}
    </div>`;
  }).join("")}</div>`;
}

$("#grid-legend").addEventListener("change", async e => {
  if (!e.target.matches("[data-showall]")) return;
  state.linkedOnly = !e.target.checked;
  await run();            // the filter is applied server-side when building
});

// clicking a cell (or a small-multiple row) jumps to that sample's evidence
$("#grid-wrap").addEventListener("click", e => {
  const el = e.target.closest("[data-sample][data-rc]");
  if (!el) return;
  const i = state.multi.samples.findIndex(r => r.sample_id === el.dataset.sample);
  if (i < 0) return;
  setActiveSample(i, { render: false });
  state.selectedPrediction = el.dataset.rc;
  renderResults();
  $("#detail").scrollIntoView({ behavior: "smooth", block: "nearest" });
});
$("#grid-wrap").addEventListener("keydown", e => {
  if (e.key === "Enter" && e.target.classList.contains("cell")) e.target.click();
});

function renderResults() {
  const r = state.result;
  $("#empty").hidden = true;
  $("#results").hidden = false;
  renderGrid();

  const samples = state.multi?.samples || [];
  const meta = r.sample_meta || {};
  $("#res-title").innerHTML = samples.length > 1
    ? `${samples.length} analyses compared`
    : `Sample ${idLink(r.sample_id)}`;
  $("#res-sub").innerHTML = [
    samples.length > 1 ? `showing evidence for ${idLink(r.sample_id)}` : (meta.origin && esc(meta.origin)),
    meta.study && `study ${idLink(meta.study)}`,
    r.settings.allow_name_matching ? "enzyme name matching on" : "ID matching only",
  ].filter(Boolean).join(" · ");

  const s = r.summary;
  $("#tiles").innerHTML = [
    [s.predictions, "predicted reaction classes", `${s.tier1_predictions} with enzyme evidence`],
    [s.detections_linked, "detections linked to references",
     `of ${(s.taxa_detected + s.functions_detected).toLocaleString()} (${s.taxa_detected.toLocaleString()} taxa, ${s.functions_detected.toLocaleString()} functions)`],
    [s.references_matched, "reference statements supported", `of ${s.references_total}`],
  ].map(([v, l, sub]) => `<div class="tile"><div class="v">${fmtNum(v)}</div><div class="l">${esc(l)}</div><div class="l">${esc(sub)}</div></div>`).join("");

  renderPredictionList();
  renderDetail();
  renderGraph();
  renderTable();
}

function renderPredictionList() {
  const list = $("#pred-list");
  const preds = state.result.predictions;
  if (!preds.length) {
    list.innerHTML = `<li class="muted">No reference statement is supported by this sample's detections.</li>`;
    return;
  }
  list.innerHTML = preds.map(p => {
    const t = p.tier1_hits.length ? "t1" : "t2";
    return `<li class="pred${p.reaction_class === state.selectedPrediction ? " selected" : ""}" data-rc="${esc(p.reaction_class)}" tabindex="0">
      <span class="name">${esc(p.reaction_class)}</span>
      <span class="score">${p.score.toFixed(2)}</span>
      <span class="bar"><span class="${t}" style="width:${Math.round(p.score * 100)}%"></span></span>
      <span class="meta">
        ${p.tier1_hits.length ? `<span class="badge t1">${p.tier1_hits.length} enzyme</span>` : ""}
        ${p.tier2_hits.length ? `<span class="badge t2">${p.tier2_hits.length} microbe</span>` : ""}
        ${p.substrate_chebi ? `<span class="badge plain">${esc(p.substrate_chebi)}</span>` : ""}
      </span>
    </li>`;
  }).join("");
}
$("#pred-list").addEventListener("click", e => {
  const li = e.target.closest(".pred");
  if (!li) return;
  selectPrediction(li.dataset.rc);
});
$("#pred-list").addEventListener("keydown", e => {
  if (e.key === "Enter" && e.target.classList.contains("pred")) selectPrediction(e.target.dataset.rc);
});

function selectPrediction(rc) {
  state.selectedPrediction = rc;
  state.selectedNode = null;
  $$(".pred").forEach(li => li.classList.toggle("selected", li.dataset.rc === rc));
  renderDetail();
  applyHighlight();
  $("#node-detail").hidden = true;
}

function renderDetail() {
  const el = $("#detail");
  const p = state.result.predictions.find(x => x.reaction_class === state.selectedPrediction);
  if (!p) { el.innerHTML = `<div class="muted">Select a prediction to see its evidence.</div>`; return; }

  const hitBlock = (h, tier) => {
    const dets = h.detections.map(d => `<li>${tier === 1 ? idLink(d.annotation_id) + " " : ""}${esc(tier === 1 ? d.label : "")}${tier === 2 ? taxonLink(d.label) : ""}
        ${matchBadge(d.match_type)}${d.abundance != null ? ` <span class="muted">· abundance ${fmtNum(d.abundance)}</span>` : ""}</li>`).join("");
    const idList = [h.entity_id, ...(h.alt_ids || [])].filter(x => x && x !== h.entity_name);
    const ids = idList.map(idLink).join(", ");
    return `<div class="ev ${tier === 1 ? "t1" : "t2"}">
      <div class="ev-title"><span class="badge ${tier === 1 ? "t1" : "t2"}">Tier ${tier}</span> ${esc(h.entity_name)}
        <span class="muted">${ids}</span></div>
      <div class="ev-body">${esc(h.evidence || "")}${h.source ? ` <em>(${esc(h.source)})</em>` : ""}</div>
      <ul>${dets}</ul>
    </div>`;
  };

  el.innerHTML = `
    <div class="detail-head">
      <div>
        <h3>${esc(p.reaction_class)}</h3>
        <div class="detail-sub">${esc(p.confidence_label)}${p.substrate_chebi ? ` · substrate ${idLink(p.substrate_chebi)}` : ""}</div>
      </div>
      <div class="detail-score"><div class="v">${p.score.toFixed(2)}</div><div class="muted">confidence</div></div>
    </div>
    ${p.tier1_hits.length ? `<div class="ev-group"><h4>Enzyme evidence (detected in functional annotation)</h4>${p.tier1_hits.map(h => hitBlock(h, 1)).join("")}</div>` : ""}
    ${p.tier2_hits.length ? `<div class="ev-group"><h4>Microbe evidence (present in taxonomic profile)</h4>${p.tier2_hits.map(h => hitBlock(h, 2)).join("")}</div>` : ""}
    ${!p.tier1_hits.length ? `<p class="muted">No enzyme for this reaction was detected directly; the prediction rests on microbe presence only.</p>` : ""}
  `;
}

/* ---------------------------------------------------------------- graph */

const NODE_W = 200, NODE_H = 34, COL_GAP = 50, ROW_GAP = 12, TOP = 36, LEFT = 16;
const COLUMNS = [
  { key: "detections", title: "Detected in sample" },
  { key: "references", title: "Reference statements" },
  { key: "reactions", title: "Reaction classes" },
  { key: "compounds", title: "Substrates" },
];
const TYPE_COLOR = {
  taxon: "--c-taxon", functional_annotation: "--c-func", reaction_class: "--c-reaction", compound: "--c-compound",
};

function nodeColor(n) {
  if (n.type === "group") return n.member_type === "taxon" ? "--c-taxon" : "--c-func";
  if (n.type === "reference") return n.entity_type === "gene" ? "--c-ref-gene" : "--c-ref-microbe";
  return TYPE_COLOR[n.type] || "--muted";
}
function nodeLabel(n) {
  switch (n.type) {
    case "taxon": return [n.organism, n.abundance != null ? `abundance ${fmtNum(n.abundance)}` : "taxon"];
    case "functional_annotation": return [n.description || n.annotation_id, n.annotation_id || "function"];
    case "reference": return [n.entity_name, `${n.entity_type === "gene" ? "Tier 1 · gene" : "Tier 2 · microbe"}${n.entity_id !== n.entity_name ? " · " + n.entity_id : ""}`];
    case "reaction_class": {
      const p = state.result.predictions.find(x => x.reaction_class === n.name);
      return [n.name, p ? `score ${p.score.toFixed(2)}` : "not supported"];
    }
    case "group": return [n.member_type === "taxon" ? `${n.members.length} ${n.ref_name} taxa` : `${n.members.length} detected functions`,
      `total abundance ${fmtNum(n.abundance)} · click to list`];
    case "compound": return [n.name || n.chebi_id, n.name && n.name !== n.chebi_id ? n.chebi_id : "substrate (ChEBI)"];
    default: return [n.id, n.type];
  }
}

const GROUP_THRESHOLD = 3;  // more detections than this on one reference -> collapse into a group node

/* Banded layout: each reaction class gets a horizontal band holding its
   reference statements, the detections supporting them and its substrates,
   so evidence reads left-to-right on the same rows. */
function buildLayout() {
  const r = state.result;
  const showUnmatched = $("#opt-unmatched").checked;
  const nodes = new Map(r.graph.nodes.map(n => [n.id, n]));
  const allEdges = r.graph.edges.filter(e => nodes.has(e.source_node) && nodes.has(e.target_node));
  const matchedRefs = new Set(r.references.filter(x => x.matched).map(x => x.node));
  const isDet = n => n && (n.type === "taxon" || n.type === "functional_annotation");

  const refToRc = new Map(), rcToComps = new Map(), refToDets = new Map(), detToRefs = new Map();
  allEdges.forEach(e => {
    const s = nodes.get(e.source_node), t = nodes.get(e.target_node);
    if (s.type === "reference" && t.type === "reaction_class") refToRc.set(s.id, t.id);
    else if (s.type === "reaction_class" && t.type === "compound") (rcToComps.get(s.id) || rcToComps.set(s.id, []).get(s.id)).push(t.id);
    else if (isDet(s) && t.type === "reference") {
      (refToDets.get(t.id) || refToDets.set(t.id, []).get(t.id)).push(e);
      (detToRefs.get(s.id) || detToRefs.set(s.id, []).get(s.id)).push(t.id);
    }
  });

  const refs = [...nodes.values()].filter(n => n.type === "reference" && (showUnmatched || matchedRefs.has(n.id)));
  const predRank = new Map(r.predictions.map((p, i) => [`reaction_class:${p.reaction_class}`, i]));
  const rcIds = [...new Set(refs.map(n => refToRc.get(n.id)))]
    .sort((a, b) => (predRank.get(a) ?? 1e6) - (predRank.get(b) ?? 1e6) || a.localeCompare(b));

  const layoutNodes = new Map(nodes);   // real nodes + synthetic group nodes
  const edges = [];
  const pos = new Map();
  const cols = [[], [], [], []];
  const place = (id, col, row) => { pos.set(id, { col, row }); cols[col].push(id); };
  const placedDet = new Set(), placedComp = new Set();
  let row = 0;

  for (const rcId of rcIds) {
    const bandStart = row;
    const bandRefs = refs.filter(n => refToRc.get(n.id) === rcId)
      .sort((a, b) => a.tier - b.tier || String(a.entity_name).localeCompare(b.entity_name));
    for (const ref of bandRefs) {
      const refRow = row;
      place(ref.id, 1, refRow);
      edges.push(...allEdges.filter(e => e.source_node === ref.id && e.target_node === rcId));
      const detEdges = refToDets.get(ref.id) || [];
      // detections that only support this reference can be grouped when there are many
      const exclusive = detEdges.filter(e => (detToRefs.get(e.source_node) || []).length === 1);
      const shared = detEdges.filter(e => !exclusive.includes(e));
      let items = [];
      if (exclusive.length > GROUP_THRESHOLD) {
        const members = exclusive.map(e => nodes.get(e.source_node));
        const gid = `group:${ref.id}`;
        const types = new Set(exclusive.map(e => e.match_type));
        const kind = members[0].type;
        layoutNodes.set(gid, {
          id: gid, type: "group", member_type: kind, members, ref_name: ref.entity_name,
          abundance: members.reduce((a, m) => a + (Number(m.abundance) || 0), 0),
        });
        items.push({ id: gid, edge: { source_node: gid, target_node: ref.id, source: "integration",
          match_type: types.size === 1 ? [...types][0] : "mixed",
          evidence: `${members.length} detections linked to ${ref.entity_name}` } });
      } else {
        items.push(...exclusive.map(e => ({ id: e.source_node, edge: e })));
      }
      shared.forEach(e => {
        if (placedDet.has(e.source_node)) edges.push(e);      // already drawn in another band
        else items.push({ id: e.source_node, edge: e });
      });
      items.forEach((it, i) => {
        if (!placedDet.has(it.id)) { place(it.id, 0, refRow + i); placedDet.add(it.id); }
        edges.push(it.edge);
      });
      row = refRow + Math.max(1, items.length);
    }
    const bandRows = Math.max(1, row - bandStart);
    const comps = (rcToComps.get(rcId) || []);
    place(rcId, 2, bandStart + (bandRows - 1) / 2);
    comps.forEach((cid, i) => {
      edges.push(...allEdges.filter(e => e.source_node === rcId && e.target_node === cid));
      if (!placedComp.has(cid)) { place(cid, 3, bandStart + (bandRows - comps.length) / 2 + i); placedComp.add(cid); }
    });
    row = bandStart + Math.max(bandRows, comps.length) + 0.5;   // half-row gap between bands
  }

  const rowH = NODE_H + ROW_GAP;
  pos.forEach(p => { p.x = LEFT + p.col * (NODE_W + COL_GAP); p.y = TOP + p.row * rowH; });
  const height = TOP + Math.max(1, row) * rowH;
  const width = LEFT * 2 + COLUMNS.length * NODE_W + (COLUMNS.length - 1) * COL_GAP;
  return { nodes: layoutNodes, cols, pos, edges, width, height };
}

function renderGraph() {
  const wrap = $("#graph-wrap");
  const L = buildLayout();
  state.layout = L;

  $("#legend").innerHTML = [
    ["--c-taxon", "Detected taxon"], ["--c-func", "Detected function"],
    ["--c-ref-gene", "Gene reference (Tier 1)"], ["--c-ref-microbe", "Microbe reference (Tier 2)", true],
    ["--c-reaction", "Reaction class"], ["--c-compound", "Substrate"],
  ].map(([c, l, dash]) => `<span><span class="sw${dash ? " dash" : ""}" style="border-color:var(${c});background:color-mix(in srgb, var(${c}) 15%, transparent)"></span>${l}</span>`).join("")
    + `<span><span class="ln"></span>ID / exact match</span><span><span class="ln dot"></span>name / genus match (weaker)</span>`;

  if (!L.cols[1].length) {
    wrap.innerHTML = `<div class="graph-empty">Nothing to draw: no reference statement is supported by this sample. Tick “Show unmatched references” to see the reference set.</div>`;
    return;
  }

  const edgeSvg = L.edges.map((e, i) => {
    const a = L.pos.get(e.source_node), b = L.pos.get(e.target_node);
    const x1 = a.x + NODE_W, y1 = a.y + NODE_H / 2, x2 = b.x, y2 = b.y + NODE_H / 2, mx = (x1 + x2) / 2;
    return `<path class="gedge ${esc(e.match_type || "")}" data-i="${i}" data-s="${esc(e.source_node)}" data-t="${esc(e.target_node)}"
      d="M${x1},${y1} C${mx},${y1} ${mx},${y2} ${x2},${y2}"><title>${esc(e.evidence || "")}</title></path>`;
  }).join("");

  const nodeSvg = [...L.pos.entries()].map(([id, p]) => {
    const n = L.nodes.get(id);
    const [label, sub] = nodeLabel(n);
    const c = nodeColor(n);
    const dash = n.type === "reference" && n.entity_type === "microbe" ? ` stroke-dasharray="5 3"` : "";
    const sw = n.type === "reference" && n.entity_type === "gene" ? ` stroke-width="2.5"` : "";
    const unmatched = n.type === "reference" && !state.result.references.find(x => x.node === id)?.matched;
    return `<g class="gnode${id === state.selectedNode ? " selected" : ""}" data-id="${esc(id)}" transform="translate(${p.x},${p.y})"${unmatched ? ` opacity=".45"` : ""}>
      <title>${esc(label)}${sub ? "\n" + esc(sub) : ""}</title>
      <rect width="${NODE_W}" height="${NODE_H}" rx="6" fill="color-mix(in srgb, var(${c}) 14%, var(--surface))" stroke="var(${c})"${dash}${sw}></rect>
      <text x="9" y="14">${esc(truncate(label, 30))}</text>
      <text class="sub" x="9" y="27">${esc(truncate(sub, 36))}</text>
    </g>`;
  }).join("");

  const count = i => i === 0
    ? L.cols[0].reduce((acc, id) => acc + (L.nodes.get(id).members?.length || 1), 0)
    : L.cols[i].length;
  const titles = COLUMNS.map((c, i) => `<text class="gcol-title" x="${LEFT + i * (NODE_W + COL_GAP)}" y="20">${c.title} (${count(i)})</text>`).join("");

  wrap.innerHTML = `<svg width="${L.width}" height="${L.height}" viewBox="0 0 ${L.width} ${L.height}" role="img" aria-label="Evidence graph">
    ${titles}<g class="edges">${edgeSvg}</g><g class="nodes">${nodeSvg}</g></svg>`;
  applyHighlight();
}

/* The evidence behind the selected prediction (or node): walk edges both ways. */
function connectedSet(startIds) {
  const L = state.layout;
  const out = new Set(startIds);
  const fwd = new Map(), back = new Map();
  L.edges.forEach(e => {
    (fwd.get(e.source_node) || fwd.set(e.source_node, []).get(e.source_node)).push(e.target_node);
    (back.get(e.target_node) || back.set(e.target_node, []).get(e.target_node)).push(e.source_node);
  });
  const walk = (map, ids) => {
    const stack = [...ids];
    while (stack.length) {
      const id = stack.pop();
      (map.get(id) || []).forEach(n => { if (!out.has(n)) { out.add(n); stack.push(n); } });
    }
  };
  // a reaction class: follow back to references and detections, forward to substrates
  walk(back, startIds); walk(fwd, startIds);
  return out;
}

function applyHighlight() {
  const svg = $("#graph-wrap svg");
  if (!svg || !state.layout) return;
  let ids = null;
  if (state.selectedNode) ids = connectedSet([state.selectedNode]);
  else if (state.selectedPrediction) ids = connectedSet([`reaction_class:${state.selectedPrediction}`]);
  svg.classList.toggle("dimmed", !!ids);
  $$(".gnode", svg).forEach(g => {
    g.classList.toggle("hl", !!ids && ids.has(g.dataset.id));
    g.classList.toggle("selected", g.dataset.id === state.selectedNode);
  });
  $$(".gedge", svg).forEach(p => p.classList.toggle("hl", !!ids && ids.has(p.dataset.s) && ids.has(p.dataset.t)));
}

$("#graph-wrap").addEventListener("click", e => {
  const g = e.target.closest(".gnode");
  if (!g) return;
  state.selectedNode = g.dataset.id;
  const n = state.layout.nodes.get(g.dataset.id);
  if (n.type === "reaction_class" && state.result.predictions.some(p => p.reaction_class === n.name)) {
    selectPrediction(n.name);
    state.selectedNode = g.dataset.id;
  }
  applyHighlight();
  renderNodeDetail(n);
});
$("#graph-clear").addEventListener("click", () => {
  state.selectedNode = null; state.selectedPrediction = null;
  $$(".pred").forEach(li => li.classList.remove("selected"));
  renderDetail(); applyHighlight(); $("#node-detail").hidden = true;
});
$("#opt-unmatched").addEventListener("change", () => state.result && renderGraph());

function renderNodeDetail(n) {
  const el = $("#node-detail");
  const rows = [];
  const add = (k, v) => { if (v !== undefined && v !== null && v !== "" && !(Array.isArray(v) && !v.length)) rows.push([k, v]); };
  add("Type", esc(n.type.replace("_", " ")));
  if (n.type === "group") {
    add("Members", n.members.map(m => m.type === "taxon"
      ? `${taxonLink(m.organism)}${m.abundance != null ? ` <span class="muted">(${fmtNum(m.abundance)})</span>` : ""}`
      : `${idLink(m.annotation_id)} ${esc(m.description || "")}`).join("<br>"));
  }
  if (n.type === "taxon") { add("Organism", taxonLink(n.organism)); add("Lineage", esc(n.lineage)); add("Rank", esc(n.rank)); add("Abundance", fmtNum(n.abundance)); }
  if (n.type === "functional_annotation") { add("Annotation", idLink(n.annotation_id)); add("Description", esc(n.description)); add("Abundance", fmtNum(n.abundance)); }
  if (n.type === "reference") {
    add("Kind", n.entity_type === "gene" ? "Gene / enzyme (Tier 1 when detected)" : "Microbe (Tier 2)");
    add("Entity", n.entity_type === "microbe" ? taxonLink(n.entity_id) : idLink(n.entity_id));
    add("Other IDs", (n.alt_ids || []).map(idLink).join(", "));
    add("Name", esc(n.entity_name)); add("Reaction class", esc(n.reaction_class));
    add("Substrate", idLink(n.chebi_substrate)); add("Source", esc(n.source)); add("Evidence", esc(n.evidence));
  }
  if (n.type === "reaction_class") add("Reaction class", esc(n.name));
  if (n.type === "compound") add("ChEBI", idLink(n.chebi_id));
  add("Evidence source", n.source && n.type !== "reference" ? esc(n.source) : null);
  add("Detection evidence", n.type !== "reference" ? esc(n.evidence) : null);
  el.innerHTML = `<strong>${esc(nodeLabel(n)[0])}</strong><dl>${rows.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v}</dd>`).join("")}</dl>`;
  el.hidden = false;
}

/* ---------------------------------------------------------------- tables */

const TABLE_LIMIT = 500;

function renderTable() {
  const r = state.result;
  if (!r) return;
  const q = $("#table-filter").value.trim().toLowerCase();
  const linkedOnly = $("#opt-linked-only").checked;
  let rows, head, rowHtml;

  if (state.table === "taxa") {
    rows = r.taxa;
    head = ["Organism", "Rank", "Abundance", "Linked to"];
    rowHtml = t => `<tr class="${t.linked_references.length ? "linked" : ""}"><td>${taxonLink(t.organism)}</td><td>${esc(t.rank || "")}</td>
      <td class="num">${fmtNum(t.abundance)}</td><td>${linkedList(t.linked_references)}</td></tr>`;
  } else if (state.table === "functions") {
    rows = r.functions;
    head = ["ID", "Description", "Abundance", "Linked to"];
    rowHtml = f => `<tr class="${f.linked_references.length ? "linked" : ""}"><td><span class="id">${idLink(f.annotation_id)}</span></td><td>${esc(f.description || "")}</td>
      <td class="num">${fmtNum(f.abundance)}</td><td>${linkedList(f.linked_references)}</td></tr>`;
  } else {
    rows = r.references;
    head = ["Tier", "Entity", "Name", "Reaction class", "Supported"];
    rowHtml = x => `<tr class="${x.matched ? "linked" : ""}"><td><span class="badge ${x.entity_type === "gene" ? "t1" : "t2"}">${x.entity_type === "gene" ? "1 · gene" : "2 · microbe"}</span></td>
      <td>${x.entity_type === "microbe" ? taxonLink(x.entity_id) : `<span class="id">${idLink(x.entity_id)}</span>`}</td>
      <td>${esc(x.entity_name)}</td><td>${esc(x.reaction_class)}</td><td>${x.matched ? "yes" : `<span class="muted">no</span>`}</td></tr>`;
  }

  const isLinked = x => (x.linked_references ? x.linked_references.length > 0 : x.matched);
  const filtered = rows.filter(x => (!linkedOnly || isLinked(x)) &&
    (!q || JSON.stringify(x).toLowerCase().includes(q)));
  const shown = filtered.slice(0, TABLE_LIMIT);
  $("#table-wrap").innerHTML = shown.length
    ? `<table><thead><tr>${head.map(h => `<th>${h}</th>`).join("")}</tr></thead><tbody>${shown.map(rowHtml).join("")}</tbody></table>`
    : `<div class="graph-empty">No rows match.</div>`;
  $("#table-note").textContent = filtered.length > TABLE_LIMIT
    ? `Showing ${TABLE_LIMIT} of ${filtered.length.toLocaleString()} rows; use the filter to narrow down.`
    : `${filtered.length.toLocaleString()} of ${rows.length.toLocaleString()} rows.`;
}
function linkedList(refNodes) {
  if (!refNodes.length) return `<span class="muted">—</span>`;
  const refs = new Map(state.result.references.map(x => [x.node, x]));
  return refNodes.map(n => { const x = refs.get(n); return x ? `${esc(x.reaction_class)} <span class="muted">(${esc(x.entity_name)})</span>` : esc(n); }).join("<br>");
}
$("#table-filter").addEventListener("input", renderTable);
$("#opt-linked-only").addEventListener("change", renderTable);

/* ---------------------------------------------------------------- export */

$("#export-json").addEventListener("click", () => {
  if (!state.multi) return;
  const many = state.multi.samples.length > 1;
  const name = many ? `expose_${state.multi.samples.length}_analyses` : `expose_${state.result.sample_id}`;
  download(`${name}.json`, JSON.stringify(many ? state.multi : state.result, null, 2), "application/json");
});
$("#export-csv").addEventListener("click", () => {
  if (!state.multi) return;
  const samples = state.multi.samples;
  const q = v => `"${String(v ?? "").replace(/"/g, '""')}"`;
  const lines = [["sample_id", "reaction_class", "score", "confidence", "substrate_chebi", "tier", "reference_entity", "reference_name", "detected_as", "match_type", "abundance", "reference_evidence"].join(",")];
  samples.forEach(r => r.predictions.forEach(p => {
    [[1, p.tier1_hits], [2, p.tier2_hits]].forEach(([tier, hits]) => hits.forEach(h => h.detections.forEach(d => {
      lines.push([r.sample_id, p.reaction_class, p.score, p.confidence_label, p.substrate_chebi, tier, h.entity_id, h.entity_name,
        d.annotation_id ? `${d.annotation_id} ${d.label}` : d.label, d.match_type, d.abundance, h.evidence].map(q).join(","));
    })));
  }));
  const name = samples.length > 1 ? `expose_${samples.length}_analyses` : `expose_${samples[0].sample_id}`;
  download(`${name}_predictions.csv`, lines.join("\n"), "text/csv");
});

/** Export whichever grid is on screen, as a matrix with analyses across.
 *  Small multiples show the reaction data, so they export that. */
$("#export-grid").addEventListener("click", () => {
  if (!state.multi) return;
  const q = v => `"${String(v ?? "").replace(/"/g, '""')}"`;
  const view = state.gridView === "multiples" ? "heatmap" : state.gridView;

  const spec = {
    heatmap: () => {
      const grid = state.multi.grid;
      return {
        grid, name: "reactions",
        header: ["reaction_class", "substrate_chebi",
                 ...grid.sample_ids.flatMap(s => [`${s} score`, `${s} tier`, `${s} observed`])],
        row: r => [r.reaction_class, r.substrate_chebi,
                   ...r.cells.flatMap(c => c ? [c.score ?? "", c.tier ?? "", c.observed ?? ""] : ["", "", ""])],
      };
    },
    taxa: () => {
      const grid = state.multi.taxon_grid;
      return {
        grid, name: "microbes",
        header: ["taxon", "tax_id", "rank", "linked_references",
                 ...grid.sample_ids.flatMap(s => [`${s} abundance`, `${s} share`])],
        row: r => [r.taxon, r.tax_id, r.rank, (r.linked_references || []).length,
                   ...r.cells.flatMap(c => c ? [c.abundance ?? "", c.share ?? ""] : ["", ""])],
      };
    },
    enzymes: () => {
      const grid = state.multi.enzyme_grid;
      return {
        grid, name: "enzymes",
        header: ["annotation_id", "description", "match_types", "linked_references",
                 ...grid.sample_ids.flatMap(s => [`${s} abundance`, `${s} share`])],
        row: r => [r.annotation_id, r.description, (r.match_types || []).join(";"),
                   (r.linked_references || []).length,
                   ...r.cells.flatMap(c => c ? [c.abundance ?? "", c.share ?? ""] : ["", ""])],
      };
    },
  };

  const { grid, name, header, row } = (spec[view] || spec.heatmap)();
  if (!grid || !grid.rows.length) { setStatus("#run-status", "Nothing to export in this view.", "err"); return; }
  const lines = [header.map(q).join(","), ...grid.rows.map(r => row(r).map(q).join(","))];
  download(`expose_${name}_${grid.sample_ids.length}_analyses.csv`, lines.join("\n"), "text/csv");
});

/* ------------------------------------------------------------ the cache */

/** Clearing the cache throws away hours of downloads and, while ChEMBL is
 *  down, references that cannot be rebuilt at all -- so it takes a second,
 *  deliberate confirmation that states what is about to go. */
function showCacheSize(stats) {
  const mb = (stats.file_bytes || 0) / 1e6;
  $("#cache-clear").textContent = stats.entries
    ? `Cache ${mb >= 1 ? `${mb.toFixed(0)} MB` : "empty"}` : "Cache empty";
  $("#cache-clear").disabled = !stats.entries;
  return stats;
}

async function refreshCacheBox() {
  try {
    return showCacheSize(await api("/api/cache"));
  } catch {
    $("#cachebox").hidden = true;
    return null;
  }
}

$("#cache-clear").addEventListener("click", async () => {
  const stats = await refreshCacheBox();
  if (!stats || !stats.entries) return;
  const mb = (stats.file_bytes || 0) / 1e6;
  // group by source (http.json + http.text are both "http"), largest first
  const bySource = {};
  Object.entries(stats.namespaces || {}).forEach(([name, info]) => {
    const source = name.split(".")[0];
    bySource[source] = (bySource[source] || 0) + info.entries;
  });
  const biggest = Object.entries(bySource).sort((a, b) => b[1] - a[1]).slice(0, 3)
    .map(([source, n]) => `${fmtNum(n)} ${source}`).join(", ");
  $("#cache-confirm-text").textContent =
    `Delete ${fmtNum(stats.entries)} cached entries (${mb.toFixed(0)} MB${biggest ? `: ${biggest}` : ""})? ` +
    `Everything must be downloaded again.`;
  $("#cache-confirm").hidden = false;
  $("#cache-clear").hidden = true;
});

$("#cache-clear-no").addEventListener("click", () => {
  $("#cache-confirm").hidden = true;
  $("#cache-clear").hidden = false;
});

$("#cache-clear-yes").addEventListener("click", async () => {
  $("#cache-clear-yes").disabled = true;
  try {
    const { removed } = await api("/api/cache", { method: "DELETE" });
    // whatever was read from the cache is now stale
    state.chembl.built = null;
    state.chembl.library = [];
    state.mgnify.data = {};
    state.enzymeTaxonomy = {};
    $("#cb-library").hidden = true;
    $("#mg-library").hidden = true;
    await refreshCacheBox();
    updateRefStatus();
    setStatus("#run-status", `Cache cleared: ${fmtNum(removed)} entries removed.`, "ok");
  } catch (err) {
    setStatus("#run-status", err.message, "err");
  } finally {
    $("#cache-clear-yes").disabled = false;
    $("#cache-confirm").hidden = true;
    $("#cache-clear").hidden = false;
  }
});

/* ------------------------------------------------------------------ init */

(async function init() {
  try {
    const h = await api("/api/health");
    $("#version").textContent = `v${h.version}`;
  } catch { /* server not reachable: errors surface on use */ }
  updateRefStatus();
  refreshCacheBox();
})();
