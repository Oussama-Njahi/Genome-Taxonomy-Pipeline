// app.js - Frontend logic: upload, poll job status, render results as plates.
// No build step, no framework: plain fetch() calls against the FastAPI backend.
import {
  METRICS, appendName, clusterOrder, contestedPairs, drawLegend, drawMatrix, drawTree, formatValue,
  layoutTree, parseAAI, parseANI, parseMatrix, parseNewick, verdict,
} from "./plates.js";

const $ = (id) => document.getElementById(id);
const dropzone = $("dropzone");
const fileInput = $("fileInput");
const fileListEl = $("fileList");
const fileCountEl = $("fileCount");
const analyzeBtn = $("analyzeBtn");
const messageBox = $("messageBox");

const progressCard = $("progressCard");
const progressLabel = $("progressLabel");
const progressPct = $("progressPct");
const progressFill = $("progressFill");
const pressSteps = $("pressSteps");
const logBox = $("logBox");
const cancelBtn = $("cancelBtn");

const resultsSection = $("resultsSection");
const newAnalysisBtn = $("newAnalysisBtn");
const stepCheckboxes = document.querySelectorAll(".step-checkbox");

const ALLOWED_EXT = [".fas", ".fasta", ".fna"];
const MIN_GENOMES = 5;

// Same order as STEPS_KEYS in scripts/pipeline_web.py
const STEPS = [
  { key: "qc", num: "I", title: "Quality control" },
  { key: "checkm2", num: "II", title: "Completeness & contamination" },
  { key: "proteines", num: "III", title: "Protein prediction" },
  { key: "ani", num: "IV", title: "Average nucleotide identity" },
  { key: "dddh", num: "V", title: "Digital DNA–DNA hybridization" },
  { key: "aai", num: "VI", title: "Average amino-acid identity" },
  { key: "pocp", num: "VII", title: "Conserved proteins" },
  { key: "arbre", num: "VIII", title: "Phylogenomic tree" },
  { key: "figures", num: "IX", title: "Printed figures" },
];

let selectedFiles = [];
let pollTimer = null;
let currentJobId = null;

// ---------- motion: ink drawing the tree plate (the one signature animation) ----------

const reduceMotion = matchMedia("(prefers-reduced-motion: reduce)").matches;
const gsap = window.gsap;
const canInk = Boolean(gsap && window.DrawSVGPlugin) && !reduceMotion;
if (canInk) gsap.registerPlugin(window.DrawSVGPlugin);

function inkTree(paths, labels) {
  if (!canInk) return;
  gsap.set(paths, { drawSVG: "0% 0%" });
  gsap.set(labels, { opacity: 0 });
  const tl = gsap.timeline({ defaults: { ease: "none" } });
  paths.forEach((p) => tl.to(p, { drawSVG: "0% 100%", duration: 0.35 }, Number(p.dataset.depth) * 1.1));
  tl.to(labels, { opacity: 1, duration: 0.5, ease: "power1.out" }, ">-0.2");
}

// ---------- upload ----------

function getSelectedSteps() {
  return Array.from(stepCheckboxes).filter((cb) => cb.checked).map((cb) => cb.value);
}

function updateAnalyzeButtonState() {
  const enoughFiles = selectedFiles.length >= MIN_GENOMES;
  const steps = getSelectedSteps().length;
  analyzeBtn.disabled = !(enoughFiles && steps > 0);
  analyzeBtn.textContent = enoughFiles
    ? `Analyze ${selectedFiles.length} genomes`
    : "Analyze";
}

stepCheckboxes.forEach((cb) => cb.addEventListener("change", updateAnalyzeButtonState));

function hasAllowedExtension(name) {
  const lower = name.toLowerCase();
  return ALLOWED_EXT.some((ext) => lower.endsWith(ext));
}

function formatSize(bytes) {
  if (bytes < 1024) return bytes + " B";
  if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + " KB";
  return (bytes / (1024 * 1024)).toFixed(1) + " MB";
}

function showMessage(text, kind) {
  messageBox.textContent = text;
  messageBox.className = "message-box " + kind;
}

function clearMessage() {
  messageBox.className = "message-box hidden";
}

function renderFileList() {
  fileListEl.textContent = "";
  selectedFiles.forEach((file, index) => {
    const li = document.createElement("li");
    const name = document.createElement("span");
    name.className = "acc-name";
    appendName(name, file.name.replace(/\.(fas|fasta|fna)$/i, ""));
    const size = document.createElement("span");
    size.className = "acc-size";
    size.textContent = formatSize(file.size);
    const removeBtn = document.createElement("button");
    removeBtn.className = "acc-remove";
    removeBtn.type = "button";
    removeBtn.textContent = "remove";
    removeBtn.setAttribute("aria-label", "Remove " + file.name);
    removeBtn.addEventListener("click", () => {
      selectedFiles.splice(index, 1);
      renderFileList();
    });
    li.append(name, size, removeBtn);
    fileListEl.appendChild(li);
  });

  // empty numbered slots up to the minimum, so the five-genome rule is visible
  for (let k = selectedFiles.length; k < MIN_GENOMES; k++) {
    const li = document.createElement("li");
    li.className = "awaiting";
    li.setAttribute("aria-hidden", "true");
    li.textContent = "awaiting a genome";
    fileListEl.appendChild(li);
  }

  const count = selectedFiles.length;
  fileCountEl.classList.toggle("ready", count >= MIN_GENOMES);
  if (count === 0) fileCountEl.textContent = "No files selected";
  else if (count < MIN_GENOMES) fileCountEl.textContent =
    `${count} of the ${MIN_GENOMES} genomes needed at least — add ${MIN_GENOMES - count} more.`;
  else fileCountEl.textContent = `${count} genomes, ready.`;
  updateAnalyzeButtonState();
}

function addFiles(fileArray) {
  clearMessage();
  const rejected = [];
  fileArray.forEach((file) => {
    if (!hasAllowedExtension(file.name)) {
      rejected.push(file.name);
      return;
    }
    if (!selectedFiles.some((f) => f.name === file.name && f.size === file.size)) {
      selectedFiles.push(file);
    }
  });
  if (rejected.length) {
    showMessage("Ignored unsupported file(s): " + rejected.join(", ") + ". Allowed: .fas, .fasta, .fna", "error");
  }
  renderFileList();
}

dropzone.addEventListener("click", () => fileInput.click());
dropzone.addEventListener("keydown", (e) => {
  if (e.key === "Enter" || e.key === " ") { e.preventDefault(); fileInput.click(); }
});
fileInput.addEventListener("change", (e) => { addFiles(Array.from(e.target.files)); fileInput.value = ""; });
["dragenter", "dragover"].forEach((evt) =>
  dropzone.addEventListener(evt, (e) => { e.preventDefault(); dropzone.classList.add("dragover"); }));
["dragleave", "drop"].forEach((evt) =>
  dropzone.addEventListener(evt, (e) => { e.preventDefault(); dropzone.classList.remove("dragover"); }));
dropzone.addEventListener("drop", (e) => addFiles(Array.from(e.dataTransfer.files)));

analyzeBtn.addEventListener("click", startAnalysis);
newAnalysisBtn.addEventListener("click", resetToUpload);
cancelBtn.addEventListener("click", cancelAnalysis);

// ---------- job lifecycle ----------

function setJobInUrl(jobId) {
  history.replaceState(null, "", jobId ? "#job=" + jobId : location.pathname + location.search);
}

async function startAnalysis() {
  clearMessage();
  analyzeBtn.disabled = true;

  const chosenSteps = getSelectedSteps();
  const formData = new FormData();
  selectedFiles.forEach((file) => formData.append("files", file));
  formData.append("steps", chosenSteps.length === stepCheckboxes.length ? "all" : chosenSteps.join(","));

  let response;
  try {
    response = await fetch("/api/jobs", { method: "POST", body: formData });
  } catch (err) {
    showMessage("Could not reach the server. Is the backend running?", "error");
    updateAnalyzeButtonState();
    return;
  }
  if (!response.ok) {
    showMessage(await safeDetail(response), "error");
    updateAnalyzeButtonState();
    return;
  }
  const data = await response.json();
  beginWatching(data.job_id, chosenSteps.length === stepCheckboxes.length ? "all" : chosenSteps.join(","));
}

function beginWatching(jobId, steps) {
  currentJobId = jobId;
  setJobInUrl(jobId);
  resultsSection.classList.add("hidden");
  $("navPlates").classList.add("hidden");
  progressCard.classList.remove("hidden");
  cancelBtn.disabled = false;
  logBox.textContent = "";
  renderProgress({ status: "queued", step: 0, total: 0, steps, log: [] });
  pollStatus(jobId);
}

async function cancelAnalysis() {
  if (!currentJobId) return;
  cancelBtn.disabled = true;
  try {
    await fetch(`/api/jobs/${currentJobId}/cancel`, { method: "POST" });
  } catch (err) {
    // network hiccup: the job keeps running, the next status poll will reflect it
  }
  backToUploadState("Analysis cancelled. Adjust your selection and click Analyze again.", "info");
}

function backToUploadState(message, kind) {
  if (pollTimer) clearInterval(pollTimer);
  setJobInUrl(null);
  progressCard.classList.add("hidden");
  if (message) showMessage(message, kind); else clearMessage();
  updateAnalyzeButtonState();
}

async function safeDetail(response) {
  try {
    const data = await response.json();
    return data.detail || "Something went wrong.";
  } catch (e) {
    return "Something went wrong (HTTP " + response.status + ").";
  }
}

function handleState(jobId, state) {
  renderProgress(state);
  if (state.status === "done") {
    clearInterval(pollTimer);
    loadResults(jobId);
  } else if (state.status === "error") {
    backToUploadState(state.error || "The pipeline failed. See the run log for details.", "error");
    progressCard.classList.remove("hidden");   // keep the log readable after a failure
    cancelBtn.disabled = true;
  } else if (state.status === "cancelled") {
    backToUploadState("Analysis cancelled. Adjust your selection and click Analyze again.", "info");
  }
}

function pollStatus(jobId) {
  if (pollTimer) clearInterval(pollTimer);
  const tick = async () => {
    let response;
    try {
      response = await fetch(`/api/jobs/${jobId}/status`);
    } catch (err) {
      return; // transient network hiccup, try again next tick
    }
    if (response.status === 404) {
      backToUploadState("This analysis is no longer known to the server: job status is kept in memory, so a server restart forgets it.", "warn");
      return;
    }
    if (!response.ok) return;
    handleState(jobId, await response.json());
  };
  tick();
  pollTimer = setInterval(tick, 1500);
}

function stepsOf(state) {
  const raw = state.steps || "all";
  const chosen = raw === "all" ? STEPS.map((s) => s.key) : raw.split(",").map((s) => s.trim());
  return STEPS.filter((s) => chosen.includes(s.key));
}

function renderProgress(state) {
  const steps = stepsOf(state);
  const total = steps.length;
  const done = state.status === "done";
  const running = state.status === "running";
  // PROGRESS n/total: chapter n is in the press, so n-1 are finished
  const finished = done ? total : Math.max(0, (state.step || 0) - 1);
  const pct = total ? Math.round((finished / total) * 100) : 0;
  progressFill.style.width = pct + "%";
  progressPct.textContent = pct + "%";
  progressLabel.textContent =
    state.status === "queued" ? "Waiting for the press…"
      : done ? "Finished"
      : state.step ? `Chapter ${state.step} of ${total} — ${state.label}` : "Starting…";

  pressSteps.textContent = "";
  steps.forEach((s, k) => {
    const li = document.createElement("li");
    const isDone = done || k < finished;
    const isCurrent = running && k === finished && state.step > 0;
    li.className = isDone ? "done" : isCurrent ? "current" : "";
    const num = document.createElement("span"); num.className = "ps-num"; num.textContent = s.num + ".";
    const title = document.createElement("span"); title.className = "ps-title"; title.textContent = s.title;
    const st = document.createElement("span"); st.className = "ps-state";
    st.textContent = isDone ? "set" : isCurrent ? "in the press" : "waiting";
    li.append(num, title, st);
    pressSteps.appendChild(li);
  });

  const log = state.log || [];
  const atBottom = logBox.scrollHeight - logBox.scrollTop - logBox.clientHeight < 24;
  logBox.textContent = log.join("\n");
  if (atBottom) logBox.scrollTop = logBox.scrollHeight;
}

// ---------- results ----------

// Every step is optional: a plate is shown only if its step produced its
// data file, and the printed PNG only if the Figures step drew it.
const RESULT_TABS = [
  { tab: "tab-qc", data: "qc" },
  { tab: "tab-ani", data: "ani" },
  { tab: "tab-dddh", data: "dddh" },
  { tab: "tab-aai", data: "aai" },
  { tab: "tab-pocp", data: "pocp" },
  { tab: "tab-tree", data: "tree_newick" },
];

const DOWNLOAD_LINKS = {
  "dl-qc": "qc",
  "dl-ani": "ani", "dl-ani-png": "heatmap_ani",
  "dl-dddh": "dddh", "dl-dddh-png": "heatmap_dddh",
  "dl-aai": "aai", "dl-aai-png": "heatmap_aai",
  "dl-pocp": "pocp", "dl-pocp-png": "heatmap_pocp",
  "dl-tree-newick": "tree_newick", "dl-tree-png": "tree_png",
};

const PARSERS = { ani: parseANI, dddh: parseMatrix, aai: parseAAI, pocp: parseMatrix };

// data loaded for the current job, and the pair being compared across plates
let results = null;
let selectedPair = null;
const drawn = new Map();       // tab id -> plate handle, drawn lazily on first view

async function fetchText(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return r.text();
}

async function loadResults(jobId) {
  const response = await fetch(`/api/jobs/${jobId}/results`);
  if (!response.ok) {
    showMessage("Analysis finished but results could not be loaded.", "error");
    return;
  }
  const data = await response.json();
  const files = data.files;

  results = { metrics: {}, layout: null, order: [], files };
  selectedPair = null;
  drawn.clear();
  const problems = [];
  await Promise.all([
    ...Object.keys(PARSERS).filter((k) => files[k]).map(async (k) => {
      try { results.metrics[k] = PARSERS[k](await fetchText(files[k])); }
      catch (e) { problems.push(METRICS[k].label); }
    }),
    (async () => {
      if (!files.tree_newick) return;
      try { results.layout = layoutTree(parseNewick(await fetchText(files.tree_newick))); }
      catch (e) { problems.push("tree"); }
    })(),
  ]);
  await document.fonts.ready;

  results.order = plateOrder();
  renderQcTable(data.qc);

  Object.entries(DOWNLOAD_LINKS).forEach(([id, key]) => {
    const link = $(id);
    link.parentElement.classList.toggle("hidden", !files[key]);
    if (files[key]) link.href = files[key];
    else link.removeAttribute("href");
  });

  let firstTab = null;
  RESULT_TABS.forEach(({ tab, data: key }) => {
    const ran = Boolean(files[key]);
    document.querySelector(`.tab-btn[data-tab="${tab}"]`).classList.toggle("hidden", !ran);
    if (ran && !firstTab && tab !== "tab-qc") firstTab = tab;
  });
  if (!firstTab && files.qc) firstTab = "tab-qc";

  progressCard.classList.add("hidden");
  if (!firstTab) {
    showMessage("Analysis finished. The selected steps produce no table or figure to display here.", "info");
    return;
  }
  if (problems.length) showMessage("Some results could not be read: " + problems.join(", ") + ". They can still be downloaded.", "warn");

  const shown = RESULT_TABS.filter(({ data: key }) => files[key]).length;
  $("platesLede").textContent =
    `${data.n_genomes || results.order.length} genomes · ${shown} plate${shown > 1 ? "s" : ""}. ` +
    "Point at a cell to compare that pair across every measure; arrow keys work too.";
  resultsSection.classList.remove("hidden");
  $("navPlates").classList.remove("hidden");
  activateTab(firstTab);
  resultsSection.scrollIntoView({ behavior: reduceMotion ? "auto" : "smooth", block: "start" });
}

// Rows follow the tree when there is one, so clades read as blocks;
// otherwise genomes are grouped by genus.
function plateOrder() {
  const names = new Set();
  Object.values(results.metrics).forEach((p) => p.names.forEach((n) => names.add(n)));
  if (results.layout) {
    const tips = results.layout.leaves.map((l) => l.name);
    const inTree = tips.filter((t) => names.has(t) || !names.size);
    const rest = [...names].filter((n) => !tips.includes(n)).sort();
    return inTree.concat(rest);
  }
  const sorted = [...names].sort();
  const key = ["aai", "pocp", "dddh", "ani"].find((k) => results.metrics[k]);
  return key ? clusterOrder(sorted, METRICS[key], results.metrics[key]) : sorted;
}

function treeMatchesOrder(order) {
  const layout = results.layout;
  if (!layout || layout.leaves.length !== order.length) return null;
  return layout.leaves.every((l, i) => l.name === order[i]) ? layout : null;
}

function drawPlate(tab) {
  if (drawn.has(tab)) return drawn.get(tab);
  const panel = $(tab);
  let handle = null;
  const key = panel.dataset.metric;
  if (key && results.metrics[key]) {
    const metric = METRICS[key];
    const pairs = results.metrics[key];
    const order = results.order.filter((n) => pairs.names.has(n));
    const plate = drawMatrix($("plot-" + key), metric, pairs, order, treeMatchesOrder(order),
      (a, b) => selectPair(a, b, key));
    const margin = panel.querySelector(".margin");
    margin.textContent = "";
    const legend = section(margin, "Scale");
    legend.classList.add("legend");
    drawLegend(legend, metric, plate.colour);
    const insp = section(margin, "The pair");
    insp.classList.add("inspector");
    renderInspector(insp, null, null, key);
    const cont = section(margin, `Near the ${metric.rank} line`);
    renderContested(cont, metric, pairs, order);
    handle = { plate, inspector: insp, key };
  } else if (tab === "tab-tree" && results.layout) {
    const tree = drawTree($("plot-tree"), results.layout);
    inkTree(tree.paths, tree.labels);
    handle = { tree };
  }
  drawn.set(tab, handle);
  return handle;
}

function section(parent, title) {
  const s = document.createElement("section");
  const h = document.createElement("h4");
  h.textContent = title;
  s.appendChild(h);
  parent.appendChild(s);
  return s;
}

function selectPair(a, b, sourceKey) {
  selectedPair = [a, b];
  const handle = drawn.get("tab-" + sourceKey);
  if (handle) renderInspector(handle.inspector, a, b, sourceKey);
}

function renderInspector(box, a, b, currentKey) {
  box.querySelectorAll(":scope > :not(h4)").forEach((n) => n.remove());
  if (!a) {
    const p = document.createElement("p");
    p.className = "inspector-empty";
    p.textContent = "Point at a cell to read one pair across every measure.";
    box.appendChild(p);
    return;
  }
  const pair = document.createElement("p");
  pair.className = "inspector-pair";
  appendName(pair, a);
  const vs = document.createElement("span");
  vs.className = "vs";
  vs.textContent = a === b ? "against itself" : "compared with";
  pair.appendChild(vs);
  if (a !== b) appendName(pair, b);
  box.appendChild(pair);
  if (a === b) return;

  const table = document.createElement("table");
  table.className = "inspector-table";
  Object.values(METRICS).forEach((m) => {
    const pairs = results.metrics[m.key];
    if (!pairs) return;
    const tr = document.createElement("tr");
    if (m.key === currentKey) tr.className = "current";
    const th = document.createElement("th");
    if (m.key === "dddh") {
      const d = document.createElement("span");
      d.className = "lc";
      d.textContent = "d";
      th.append(d, "DDH");
    } else th.textContent = m.label;
    const td = document.createElement("td");
    const v = pairs.get(a, b);
    const w = pairs.get(b, a);
    // ANI is asymmetric: show both directions when they differ
    td.textContent = m.key === "ani" && w !== undefined && w !== null && v !== null && v !== undefined &&
      Math.abs(v - w) >= 0.05
      ? `${formatValue(m, v)} / ${formatValue(m, w)}`
      : formatValue(m, v === undefined ? null : v);
    const vd = document.createElement("td");
    vd.className = "verdict";
    vd.textContent = verdict(m, v === undefined ? null : v);
    tr.append(th, td, vd);
    table.appendChild(tr);
  });
  box.appendChild(table);
}

function renderContested(box, metric, pairs, order) {
  const list = contestedPairs(metric, pairs, order);
  if (!list.length) {
    const p = document.createElement("p");
    p.className = "contested-none";
    p.textContent = "No pair sits close to the threshold.";
    box.appendChild(p);
    return;
  }
  const ul = document.createElement("ul");
  ul.className = "contested";
  list.forEach(({ a, b, v }) => {
    const li = document.createElement("li");
    const btn = document.createElement("button");
    btn.type = "button";
    const names = document.createElement("span");
    names.textContent = `${order.indexOf(a) + 1} × ${order.indexOf(b) + 1}`;
    names.title = `${a} × ${b}`;
    const val = document.createElement("span");
    val.className = "val";
    val.textContent = formatValue(metric, v) + (metric.unit ? " " + metric.unit : "");
    btn.append(names, val);
    btn.addEventListener("click", () => {
      const handle = drawn.get("tab-" + metric.key);
      handle.plate.select(a, b);
      selectPair(a, b, metric.key);
    });
    li.appendChild(btn);
    ul.appendChild(li);
  });
  box.appendChild(ul);
}

function renderQcTable(rows) {
  const table = $("qcTable");
  table.textContent = "";
  if (!rows || !rows.length) return;

  const COLS = [
    { key: "Genome", label: "Genome" },
    { key: "Taille(pb)", label: "Size (bp)", fmt: (v) => Number(v).toLocaleString("en-US") },
    { key: "GC(%)", label: "GC %" },
    { key: "Contigs", label: "Contigs" },
    { key: "Completeness(%)", label: "Complete %" },
    { key: "Contamination(%)", label: "Contam. %" },
  ];
  const known = new Set(COLS.map((c) => c.key));
  Object.keys(rows[0]).filter((k) => !known.has(k)).forEach((k) => COLS.push({ key: k, label: k }));

  const thead = document.createElement("thead");
  const headRow = document.createElement("tr");
  ["", ...COLS.map((c) => c.label), "Remarks"].forEach((label) => {
    const th = document.createElement("th");
    th.textContent = label;
    if (label === "Genome") th.className = "name";
    headRow.appendChild(th);
  });
  thead.appendChild(headRow);
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  rows.forEach((row, i) => {
    const tr = document.createElement("tr");
    const completeness = parseFloat(row["Completeness(%)"]);
    const contamination = parseFloat(row["Contamination(%)"]);
    const lowCompleteness = !isNaN(completeness) && completeness < 90;
    const highContamination = !isNaN(contamination) && contamination > 5;

    const no = document.createElement("td");
    no.className = "no";
    no.textContent = i + 1;
    tr.appendChild(no);
    COLS.forEach((c) => {
      const td = document.createElement("td");
      const value = row[c.key];
      if (c.key === "Genome") { appendName(td, value); td.className = "name"; }
      else if (value === "NA") { td.textContent = "n.a."; td.className = "na"; td.title = "Run CheckM2 to fill this in"; }
      else td.textContent = c.fmt ? c.fmt(value) : value;
      if ((c.key === "Completeness(%)" && lowCompleteness) || (c.key === "Contamination(%)" && highContamination)) {
        td.classList.add("bad");
      }
      tr.appendChild(td);
    });
    const remark = document.createElement("td");
    remark.className = "remark";
    remark.textContent = [lowCompleteness && "low completeness", highContamination && "contaminated"]
      .filter(Boolean).join("; ");
    tr.appendChild(remark);
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
}

// Tabs
function activateTab(tabId) {
  document.querySelectorAll(".tab-btn").forEach((b) => {
    const on = b.dataset.tab === tabId;
    b.classList.toggle("active", on);
    b.setAttribute("aria-selected", on);
  });
  document.querySelectorAll(".tab-panel").forEach((p) => p.classList.toggle("active", p.id === tabId));
  const handle = drawPlate(tabId);
  // carry the pair being compared over to the new plate
  if (handle && handle.plate && selectedPair) {
    handle.plate.select(selectedPair[0], selectedPair[1]);
    renderInspector(handle.inspector, selectedPair[0], selectedPair[1], handle.key);
  }
}

document.querySelectorAll(".tab-btn").forEach((btn) => {
  btn.addEventListener("click", () => activateTab(btn.dataset.tab));
});

function resetToUpload() {
  selectedFiles = [];
  renderFileList();
  clearMessage();
  setJobInUrl(null);
  progressCard.classList.add("hidden");
  resultsSection.classList.add("hidden");
  $("navPlates").classList.add("hidden");
  analyzeBtn.disabled = true;
  $("analyze").scrollIntoView({ behavior: reduceMotion ? "auto" : "smooth", block: "start" });
}

// A reload keeps the job: its id lives in the URL (#job=<id>)
function resumeFromUrl() {
  const m = /job=([A-Za-z0-9_-]+)/.exec(location.hash);
  if (!m) return;
  currentJobId = m[1];
  progressCard.classList.remove("hidden");
  pollStatus(m[1]);
}

renderFileList();
resumeFromUrl();
