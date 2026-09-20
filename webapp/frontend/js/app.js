// app.js - Frontend logic: upload, poll job status, render results.
// No build step, no framework: plain fetch() calls against the FastAPI backend.

const dropzone = document.getElementById("dropzone");
const fileInput = document.getElementById("fileInput");
const fileListEl = document.getElementById("fileList");
const fileCountEl = document.getElementById("fileCount");
const analyzeBtn = document.getElementById("analyzeBtn");
const messageBox = document.getElementById("messageBox");

const progressCard = document.getElementById("progressCard");
const progressLabel = document.getElementById("progressLabel");
const progressPct = document.getElementById("progressPct");
const progressFill = document.getElementById("progressFill");
const logBox = document.getElementById("logBox");
const cancelBtn = document.getElementById("cancelBtn");

const resultsSection = document.getElementById("resultsSection");
const newAnalysisBtn = document.getElementById("newAnalysisBtn");

const stepCheckboxes = document.querySelectorAll(".step-checkbox");

const ALLOWED_EXT = [".fas", ".fasta", ".fna"];
let selectedFiles = [];
let pollTimer = null;
let currentJobId = null;

function getSelectedSteps() {
  return Array.from(stepCheckboxes)
    .filter((cb) => cb.checked)
    .map((cb) => cb.value);
}

function updateAnalyzeButtonState() {
  const enoughFiles = selectedFiles.length >= 5;
  const atLeastOneStep = getSelectedSteps().length > 0;
  analyzeBtn.disabled = !(enoughFiles && atLeastOneStep);
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
  fileListEl.innerHTML = "";
  selectedFiles.forEach((file, index) => {
    const li = document.createElement("li");
    const name = document.createElement("span");
    name.textContent = file.name;
    const meta = document.createElement("span");
    meta.className = "file-size";
    meta.textContent = formatSize(file.size);
    const removeBtn = document.createElement("button");
    removeBtn.className = "file-remove";
    removeBtn.innerHTML = "&times;";
    removeBtn.title = "Remove";
    removeBtn.addEventListener("click", () => {
      selectedFiles.splice(index, 1);
      renderFileList();
    });
    const right = document.createElement("span");
    right.style.display = "flex";
    right.style.alignItems = "center";
    right.style.gap = "10px";
    right.appendChild(meta);
    right.appendChild(removeBtn);
    li.appendChild(name);
    li.appendChild(right);
    fileListEl.appendChild(li);
  });

  const count = selectedFiles.length;
  fileCountEl.textContent =
    count === 0 ? "No files selected" : count + " file" + (count > 1 ? "s" : "") + " selected (minimum 5 required)";
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
    showMessage(
      "Ignored unsupported file(s): " + rejected.join(", ") + ". Allowed: .fas, .fasta, .fna",
      "error"
    );
  }
  renderFileList();
}

dropzone.addEventListener("click", () => fileInput.click());
fileInput.addEventListener("change", (e) => addFiles(Array.from(e.target.files)));

["dragenter", "dragover"].forEach((evt) =>
  dropzone.addEventListener(evt, (e) => {
    e.preventDefault();
    dropzone.classList.add("dragover");
  })
);
["dragleave", "drop"].forEach((evt) =>
  dropzone.addEventListener(evt, (e) => {
    e.preventDefault();
    dropzone.classList.remove("dragover");
  })
);
dropzone.addEventListener("drop", (e) => {
  addFiles(Array.from(e.dataTransfer.files));
});

analyzeBtn.addEventListener("click", startAnalysis);
newAnalysisBtn.addEventListener("click", resetToUpload);
cancelBtn.addEventListener("click", cancelAnalysis);

async function startAnalysis() {
  clearMessage();
  analyzeBtn.disabled = true;

  const chosenSteps = getSelectedSteps();
  const allSteps = stepCheckboxes.length;

  const formData = new FormData();
  selectedFiles.forEach((file) => formData.append("files", file));
  formData.append("steps", chosenSteps.length === allSteps ? "all" : chosenSteps.join(","));

  let response;
  try {
    response = await fetch("/api/jobs", { method: "POST", body: formData });
  } catch (err) {
    showMessage("Could not reach the server. Is the backend running?", "error");
    analyzeBtn.disabled = false;
    return;
  }

  if (!response.ok) {
    const detail = await safeDetail(response);
    showMessage(detail, "error");
    analyzeBtn.disabled = false;
    return;
  }

  const data = await response.json();
  currentJobId = data.job_id;
  progressCard.classList.remove("hidden");
  cancelBtn.disabled = false;
  logBox.textContent = "";
  pollStatus(data.job_id);
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

function pollStatus(jobId) {
  if (pollTimer) clearInterval(pollTimer);
  pollTimer = setInterval(async () => {
    let response;
    try {
      response = await fetch(`/api/jobs/${jobId}/status`);
    } catch (err) {
      return; // transient network hiccup, try again next tick
    }
    if (!response.ok) return;
    const state = await response.json();
    renderProgress(state);

    if (state.status === "done") {
      clearInterval(pollTimer);
      loadResults(jobId);
    } else if (state.status === "error") {
      backToUploadState(state.error || "The pipeline failed. See the run log for details.", "error");
    } else if (state.status === "cancelled") {
      backToUploadState("Analysis cancelled. Adjust your selection and click Analyze again.", "info");
    }
  }, 1500);
}

function renderProgress(state) {
  const pct = Math.round((state.step / state.total) * 100);
  progressFill.style.width = pct + "%";
  progressPct.textContent = pct + "%";
  progressLabel.textContent =
    state.status === "queued" ? "Waiting in queue..." : `Step ${state.step}/${state.total} - ${state.label}`;
  logBox.textContent = (state.log || []).join("\n");
  logBox.scrollTop = logBox.scrollHeight;
}

async function loadResults(jobId) {
  const response = await fetch(`/api/jobs/${jobId}/results`);
  if (!response.ok) {
    showMessage("Analysis finished but results could not be loaded.", "error");
    return;
  }
  const data = await response.json();

  renderQcTable(data.qc);

  document.getElementById("img-ani").src = data.files.heatmap_ani;
  document.getElementById("img-dddh").src = data.files.heatmap_dddh;
  document.getElementById("img-aai").src = data.files.heatmap_aai;
  document.getElementById("img-pocp").src = data.files.heatmap_pocp;
  document.getElementById("img-tree").src = data.files.tree_png;

  document.getElementById("dl-qc").href = data.files.qc;
  document.getElementById("dl-ani").href = data.files.ani;
  document.getElementById("dl-dddh").href = data.files.dddh;
  document.getElementById("dl-aai").href = data.files.aai;
  document.getElementById("dl-pocp").href = data.files.pocp;
  document.getElementById("dl-tree-png").href = data.files.tree_png;
  document.getElementById("dl-tree-newick").href = data.files.tree_newick;

  progressCard.classList.add("hidden");
  resultsSection.classList.remove("hidden");
  resultsSection.scrollIntoView({ behavior: "smooth", block: "start" });
}

function renderQcTable(rows) {
  const table = document.getElementById("qcTable");
  table.innerHTML = "";
  if (!rows || !rows.length) return;

  const thead = document.createElement("thead");
  const headRow = document.createElement("tr");
  Object.keys(rows[0]).forEach((col) => {
    const th = document.createElement("th");
    th.textContent = col;
    headRow.appendChild(th);
  });
  thead.appendChild(headRow);
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  rows.forEach((row) => {
    const tr = document.createElement("tr");

    const completeness = parseFloat(row["Completeness(%)"]);
    const contamination = parseFloat(row["Contamination(%)"]);
    const lowCompleteness = !isNaN(completeness) && completeness < 90;
    const highContamination = !isNaN(contamination) && contamination > 5;

    Object.entries(row).forEach(([col, value]) => {
      const td = document.createElement("td");
      td.textContent = value;
      if ((col === "Completeness(%)" && lowCompleteness) ||
          (col === "Contamination(%)" && highContamination)) {
        td.classList.add("value-bad");
      }
      tr.appendChild(td);
    });

    if (lowCompleteness || highContamination) {
      tr.classList.add("row-low-quality");
    }
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
}

// Tabs
document.querySelectorAll(".tab-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tab-btn").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".tab-panel").forEach((p) => p.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById(btn.dataset.tab).classList.add("active");
  });
});

function resetToUpload() {
  selectedFiles = [];
  renderFileList();
  clearMessage();
  progressCard.classList.add("hidden");
  resultsSection.classList.add("hidden");
  analyzeBtn.disabled = true;
  document.getElementById("analyze").scrollIntoView({ behavior: "smooth", block: "start" });
}
if (location.search.includes("test=qc")) {
  document.getElementById("resultsSection").classList.remove("hidden");
  renderQcTable([
    {"Genome": "Good-Genome-1", "Taille(pb)": "3000000", "GC(%)": "45.0", "Contigs": "1", "Completeness(%)": "98.5", "Contamination(%)": "1.2"},
    {"Genome": "Bad-Completeness", "Taille(pb)": "2500000", "GC(%)": "42.0", "Contigs": "5", "Completeness(%)": "85.0", "Contamination(%)": "2.0"},
    {"Genome": "Bad-Contamination", "Taille(pb)": "3200000", "GC(%)": "50.0", "Contigs": "2", "Completeness(%)": "97.0", "Contamination(%)": "8.5"}
  ]);
}