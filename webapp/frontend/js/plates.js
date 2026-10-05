// plates.js - read the pipeline's raw result files and draw them as plates:
// pairwise matrices (ANI, dDDH, AAI, POCP) beside the tree, and the tree itself.
// Everything is plain SVG built from the same files the user can download.

const SVGNS = "http://www.w3.org/2000/svg";

// ---------- names ----------

export function cleanName(path) {
  return path.trim().replace(/^.*\//, "").replace(/\.(fas|fasta|fna|faa)$/i, "");
}

// "Bacillus-subtilis-168" -> genus "Bacillus", species "subtilis", strain "168".
// Names that do not look like a binomial are shown as they are.
export function splitName(name) {
  const m = /^([A-Z][a-z]+)[-_ ]([a-z][a-z.]*)(?:[-_ ](.*))?$/.exec(name);
  if (!m) return { raw: name };
  return { genus: m[1], species: m[2], strain: m[3] || "" };
}

export function genusOf(name) {
  const parts = splitName(name);
  return parts.genus || name.split(/[-_ ]/)[0];
}

// Append a taxonomic name to an HTML or SVG element: genus and species in
// italic, strain in roman, optionally with the genus abbreviated.
export function appendName(el, name, { abbreviate = false, svg = false } = {}) {
  const parts = splitName(name);
  const make = (text, italic) => {
    const node = svg ? document.createElementNS(SVGNS, "tspan") : document.createElement(italic ? "i" : "span");
    if (svg && italic) node.setAttribute("class", "sp");
    if (!svg && italic) node.className = "sp";
    node.textContent = text;
    return node;
  };
  if (parts.raw !== undefined) {
    el.appendChild(make(parts.raw, false));
    return el;
  }
  const genus = abbreviate ? parts.genus[0] + "." : parts.genus;
  el.appendChild(make(genus + " " + parts.species, true));
  if (parts.strain) el.appendChild(make(" " + parts.strain, false));
  return el;
}

function nameText(name, abbreviate) {
  const parts = splitName(name);
  if (parts.raw !== undefined) return parts.raw;
  return (abbreviate ? parts.genus[0] + "." : parts.genus) + " " + parts.species + (parts.strain ? " " + parts.strain : "");
}

// ---------- metrics ----------

export const METRICS = {
  ani:  { key: "ani",  label: "ANI",  unit: "%", threshold: 95,     kind: "sim",  decimals: 1, rank: "species", band: 1 },
  dddh: { key: "dddh", label: "dDDH", unit: "",  threshold: 0.0412, kind: "dist", decimals: 3, rank: "species", band: 0.008 },
  aai:  { key: "aai",  label: "AAI",  unit: "%", threshold: 65,     kind: "sim",  decimals: 1, rank: "genus",   band: 3 },
  pocp: { key: "pocp", label: "POCP", unit: "%", threshold: 50,     kind: "sim",  decimals: 1, rank: "genus",   band: 3 },
};

// true when the value puts the pair on the same side of the threshold as
// identical genomes (same species / same genus)
export function sameSide(metric, v) {
  return metric.kind === "sim" ? v >= metric.threshold : v <= metric.threshold;
}

export function verdict(metric, v) {
  if (v === null || v === undefined) return "no value";
  const rank = metric.rank;
  if (sameSide(metric, v)) return "same " + rank;
  return rank === "species" ? "distinct species" : "other genus";
}

export function formatValue(metric, v, compact = false) {
  if (v === null || v === undefined) return "n.a.";
  const s = v.toFixed(metric.decimals);
  return compact && metric.kind === "dist" ? s.replace(/^0\./, ".") : s;
}

// A pairwise table: names + get(a, b) -> number | null (no value) | undefined (absent)
class Pairs {
  constructor() { this.values = new Map(); this.names = new Set(); }
  set(a, b, v) { this.values.set(a + "\u0000" + b, v); this.names.add(a); this.names.add(b); }
  get(a, b) { return this.values.get(a + "\u0000" + b); }
}

// fastANI: query, reference, ANI, mapped fragments, total fragments (no header)
export function parseANI(text) {
  const p = new Pairs();
  for (const line of text.split("\n")) {
    const c = line.split("\t");
    if (c.length < 3) continue;
    const v = parseFloat(c[2]);
    if (!isNaN(v)) p.set(cleanName(c[0]), cleanName(c[1]), v);
  }
  // fastANI omits pairs it cannot align (below ~80 %): record them as "no value"
  for (const a of p.names) for (const b of p.names) if (p.get(a, b) === undefined) p.set(a, b, null);
  return p;
}

// EzAAI: header row, labels in "Label 1"/"Label 2", value in "AAI"
export function parseAAI(text) {
  const p = new Pairs();
  const lines = text.split("\n").filter((l) => l.trim());
  const head = lines.shift().split("\t").map((h) => h.trim().toLowerCase());
  const ia = head.indexOf("label 1") >= 0 ? head.indexOf("label 1") : 2;
  const ib = head.indexOf("label 2") >= 0 ? head.indexOf("label 2") : 3;
  const iv = head.indexOf("aai") >= 0 ? head.indexOf("aai") : 4;
  for (const line of lines) {
    const c = line.split("\t");
    const v = parseFloat(c[iv]);
    p.set(cleanName(c[ia]), cleanName(c[ib]), isNaN(v) ? null : v);
  }
  return p;
}

// square matrix with a header row and row names (dDDH distance, POCP)
export function parseMatrix(text) {
  const p = new Pairs();
  const lines = text.split("\n").filter((l) => l.trim());
  const cols = lines.shift().split("\t").slice(1).map(cleanName);
  for (const line of lines) {
    const c = line.split("\t");
    const a = cleanName(c[0]);
    cols.forEach((b, j) => {
      const raw = (c[j + 1] || "").trim();
      const v = parseFloat(raw);
      p.set(a, b, raw === "" || raw.toUpperCase() === "NA" || isNaN(v) ? null : v);
    });
  }
  return p;
}

// ---------- Newick ----------

export function parseNewick(text) {
  const s = text.trim().replace(/;\s*$/, "");
  let i = 0;
  const label = () => {
    if (s[i] === "'") {
      let out = "";
      i++;
      while (i < s.length) {
        if (s[i] === "'" && s[i + 1] === "'") { out += "'"; i += 2; continue; }
        if (s[i] === "'") { i++; break; }
        out += s[i++];
      }
      return out;
    }
    const m = /^[^:,();[]*/.exec(s.slice(i));
    i += m[0].length;
    return m[0].trim();
  };
  const node = () => {
    const n = { name: "", length: 0, children: [] };
    if (s[i] === "(") {
      i++;
      for (;;) {
        n.children.push(node());
        if (s[i] === ",") { i++; continue; }
        if (s[i] === ")") { i++; break; }
        throw new Error("Malformed Newick tree near character " + i);
      }
    }
    const lab = label();
    if (s[i] === "[") i = s.indexOf("]", i) + 1;          // skip [comments]
    if (s[i] === ":") {
      i++;
      const m = /^[-+0-9.eE]+/.exec(s.slice(i)) || [""];
      n.length = parseFloat(m[0]) || 0;
      i += m[0].length;
    }
    if (n.children.length) {
      const sup = parseFloat(lab);
      if (lab !== "" && !isNaN(sup)) n.support = sup;
      else n.name = lab;
    } else {
      n.name = cleanName(lab);
    }
    return n;
  };
  return node();
}

// Rectangular phylogram layout: x = distance from the root, y = leaf rank.
export function layoutTree(root) {
  const leaves = [];
  let maxX = 0;
  let maxDepth = 0;
  (function walk(n, x, depth, parent) {
    n.x = x; n.depth = depth; n.parent = parent;
    maxX = Math.max(maxX, x); maxDepth = Math.max(maxDepth, depth);
    if (!n.children.length) { n.leafIndex = leaves.length; leaves.push(n); }
    n.children.forEach((c) => walk(c, x + (c.length || 0), depth + 1, n));
  })(root, 0, 0, null);
  // a tree without branch lengths is drawn as a cladogram
  const cladogram = maxX === 0;
  if (cladogram) {
    (function walk(n) { n.x = n.depth; n.children.forEach(walk); })(root);
    leaves.forEach((l) => { l.x = maxDepth; });
    maxX = maxDepth || 1;
  }
  (function setY(n) {
    if (!n.children.length) { n.y = n.leafIndex; return; }
    n.children.forEach(setY);
    n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
  })(root);
  const nodes = [];
  (function all(n) { nodes.push(n); n.children.forEach(all); })(root);
  return { root, leaves, nodes, maxX, cladogram };
}

// ---------- colour ----------

// Two inks diverging exactly at the threshold: Prussian blue on the side of
// "same taxon", umber on the other; lighter near the line, darker away from it.
const SAME_RAMP = ["#DCE2E4", "#93A6B6", "#4F6A82", "#2F4A63"];
const DIFF_RAMP = ["#EFE1CA", "#D2B48D", "#A47E57", "#7A5636"];

function hex(c) { return [1, 3, 5].map((k) => parseInt(c.slice(k, k + 2), 16)); }
function ramp(stops, t) {
  t = Math.max(0, Math.min(1, t));
  const x = t * (stops.length - 1);
  const k = Math.min(stops.length - 2, Math.floor(x));
  const a = hex(stops[k]); const b = hex(stops[k + 1]); const f = x - k;
  return "#" + a.map((v, i) => Math.round(v + (b[i] - v) * f).toString(16).padStart(2, "0")).join("");
}
export function luminance(c) {
  const [r, g, b] = hex(c).map((v) => { v /= 255; return v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4; });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

export function colourScale(metric, values) {
  const thr = metric.threshold;
  const finite = values.filter((v) => v !== null && v !== undefined);
  let lo = Math.min(...finite, thr);
  let hi = Math.max(...finite, thr);
  // similarities always run up to identity, so plates compare across runs
  if (metric.kind === "sim") { hi = 100; lo = Math.min(lo, thr - 1); }
  else { lo = 0; hi = Math.max(hi, thr * 1.5); }
  const fn = (v) => {
    if (metric.kind === "sim") {
      return v >= thr ? ramp(SAME_RAMP, ((v - thr) / (hi - thr)) ** 0.8) : ramp(DIFF_RAMP, ((thr - v) / (thr - lo)) ** 0.8);
    }
    return v <= thr ? ramp(SAME_RAMP, ((thr - v) / (thr - lo)) ** 0.8) : ramp(DIFF_RAMP, ((v - thr) / (hi - thr)) ** 0.8);
  };
  fn.lo = lo; fn.hi = hi;
  return fn;
}

// ---------- SVG helpers ----------

function el(name, attrs = {}, parent) {
  const node = document.createElementNS(SVGNS, name);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (parent) parent.appendChild(node);
  return node;
}

let measureCtx = null;
function textWidth(text, font) {
  measureCtx = measureCtx || document.createElement("canvas").getContext("2d");
  measureCtx.font = font;
  return measureCtx.measureText(text).width;
}

function hatchDefs(svg, id) {
  const defs = el("defs", {}, svg);
  const pat = el("pattern", { id, width: 6, height: 6, patternUnits: "userSpaceOnUse", patternTransform: "rotate(45)" }, defs);
  el("rect", { width: 6, height: 6, fill: "#F6EEE3" }, pat);
  el("line", { x1: 0, y1: 0, x2: 0, y2: 6, stroke: "#C2B29B", "stroke-width": 0.8 }, pat);
  return defs;
}

let uid = 0;

// Branch paths of a tree, root to tips, each an elbow from its parent.
function drawBranches(g, layout, sx, sy, stub = 8) {
  const paths = [];
  el("path", { class: "tr-branch", d: `M${sx(0) - stub} ${sy(layout.root.y)}H${sx(0)}` }, g);
  for (const n of layout.nodes) {
    if (!n.parent) continue;
    const p = el("path", {
      class: "tr-branch",
      d: `M${sx(n.parent.x)} ${sy(n.parent.y)}V${sy(n.y)}H${sx(n.x)}`,
    }, g);
    p.dataset.depth = n.parent.x / layout.maxX;
    paths.push(p);
  }
  return paths;
}

// ---------- matrix plate ----------

// order: genome names, top to bottom. layout: tree layout or null.
// Returns { svg, select(i, j) }; onSelect(a, b) fires on hover, click or keys.
export function drawMatrix(container, metric, pairs, order, layout, onSelect) {
  container.textContent = "";
  const n = order.length;
  const labelFont = "13px Spectral";
  const nameW = Math.max(...order.map((nm) => textWidth(nameText(nm, true), "italic " + labelFont)));
  const numW = 26;
  const labelW = numW + 8 + nameW + 12;
  const treeW = layout ? Math.min(170, Math.max(90, n * 7)) : 0;
  // cells fill the plate's width; past ~45 genomes the plate scrolls instead
  const avail = (container.clientWidth || 700) - treeW - labelW - 4;
  // never below a readable size for small sets: on a phone the plate scrolls sideways instead
  const cell = Math.max(n > 40 ? 15 : 24, Math.min(42, Math.floor(avail / n)));
  const showValues = cell >= 22;
  const compact = cell < 34;
  const headH = 26;
  const x0 = treeW + labelW;
  const y0 = headH;
  const W = Math.ceil(x0 + n * cell + 2);
  const H = Math.ceil(y0 + n * cell + 2);

  const svg = el("svg", {
    viewBox: `0 0 ${W} ${H}`, width: W, height: H, tabindex: 0, role: "grid",
    "aria-label": `${metric.label} for ${n} genomes. Use the arrow keys to move between pairs.`,
  });
  const hatchId = "hatch" + (++uid);
  hatchDefs(svg, hatchId);

  const values = [];
  order.forEach((a) => order.forEach((b) => { if (a !== b) values.push(pairs.get(a, b)); }));
  const colour = colourScale(metric, values);

  // tree, rows aligned with the matrix
  if (layout) {
    const g = el("g", {}, svg);
    const sx = (x) => 10 + (x / layout.maxX) * (treeW - 20);
    const rowOf = new Map(order.map((nm, i) => [nm, i]));
    const sy = (y) => y0 + y * cell + cell / 2;
    // leaves follow `order`, which is the tree's own leaf order
    layout.leaves.forEach((leaf) => { if (rowOf.has(leaf.name)) leaf.y = rowOf.get(leaf.name); });
    drawBranches(g, layout, sx, sy);
    layout.leaves.forEach((leaf) => {
      el("line", { class: "tr-leader", x1: sx(leaf.x) + 2, x2: treeW - 2, y1: sy(leaf.y), y2: sy(leaf.y) }, g);
    });
  }

  // row labels and column numerals
  const rowLabels = []; const colLabels = [];
  order.forEach((nm, i) => {
    const y = y0 + i * cell + cell / 2 + 4.5;
    const num = el("text", { class: "hm-num", x: treeW + numW, y, "text-anchor": "end" }, svg);
    num.textContent = i + 1;
    const t = el("text", { class: "hm-label", x: treeW + numW + 8, y }, svg);
    appendName(t, nm, { abbreviate: true, svg: true });
    el("title", {}, t).textContent = nameText(nm, false);
    rowLabels.push([num, t]);
    const c = el("text", { class: "hm-num", x: x0 + i * cell + cell / 2, y: headH - 9, "text-anchor": "middle" }, svg);
    c.textContent = i + 1;
    colLabels.push(c);
  });

  // cells
  const cells = el("g", {}, svg);
  order.forEach((a, i) => order.forEach((b, j) => {
    const x = x0 + j * cell; const y = y0 + i * cell;
    const v = pairs.get(a, b);
    let fill;
    if (a === b) fill = "#E6DCCB";
    else if (v === null || v === undefined) fill = `url(#${hatchId})`;
    else fill = colour(v);
    const r = el("rect", { class: "hm-cell", x: x + 0.5, y: y + 0.5, width: cell - 1, height: cell - 1, fill }, cells);
    r.dataset.i = i; r.dataset.j = j;
    if (a === b) {
      el("line", { x1: x + 4, y1: y + cell - 4, x2: x + cell - 4, y2: y + 4, stroke: "#A8957D", "stroke-width": 1 }, cells);
    } else if (showValues && v !== null && v !== undefined) {
      const dark = luminance(fill) < 0.3;
      const t = el("text", {
        class: "hm-val" + (sameSide(metric, v) ? " strong" : ""), x: x + cell / 2, y: y + cell / 2 + 3.6,
        "text-anchor": "middle", fill: dark ? "#F6EEE3" : "#2A211B",
      }, cells);
      t.textContent = formatValue(metric, v, compact);
    }
  }));

  const focus = el("rect", { class: "hm-focus", width: cell + 1, height: cell + 1, visibility: "hidden" }, svg);
  let cur = null;

  function select(i, j, fire = true) {
    if (cur) { rowLabels[cur[0]].forEach((t) => t.style.fill = ""); colLabels[cur[1]].style.fill = ""; }
    cur = [i, j];
    focus.setAttribute("x", x0 + j * cell - 0.5);
    focus.setAttribute("y", y0 + i * cell - 0.5);
    focus.setAttribute("visibility", "visible");
    rowLabels[i].forEach((t) => t.style.fill = "#9B2F3A");
    colLabels[j].style.fill = "#9B2F3A";
    if (fire) onSelect(order[i], order[j]);
  }

  cells.addEventListener("click", (e) => {
    const r = e.target.closest("rect");
    if (r) select(+r.dataset.i, +r.dataset.j);
  });
  cells.addEventListener("mouseover", (e) => {
    const r = e.target.closest("rect");
    if (r) select(+r.dataset.i, +r.dataset.j);
  });
  svg.addEventListener("keydown", (e) => {
    const moves = { ArrowUp: [-1, 0], ArrowDown: [1, 0], ArrowLeft: [0, -1], ArrowRight: [0, 1] };
    if (!moves[e.key]) return;
    e.preventDefault();
    const [i, j] = cur || [0, 1];
    const [di, dj] = cur ? moves[e.key] : [0, 0];
    select(Math.max(0, Math.min(n - 1, i + di)), Math.max(0, Math.min(n - 1, j + dj)));
  });
  svg.addEventListener("focus", () => { if (!cur && n > 1) select(0, 1); });

  container.appendChild(svg);
  return { svg, colour, select: (a, b) => {
    const i = order.indexOf(a); const j = order.indexOf(b);
    if (i >= 0 && j >= 0) select(i, j, false);
  } };
}

// ---------- legend ----------

export function drawLegend(container, metric, colour) {
  container.textContent = "";
  const W = 260; const H = 54; const barY = 8; const barH = 12;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
    "aria-label": `Colour scale for ${metric.label}, threshold ${metric.threshold}` });
  const id = "lg" + (++uid);
  const grad = el("linearGradient", { id }, el("defs", {}, svg));
  const lo = colour.lo; const hi = colour.hi; const thr = metric.threshold;
  // diverging scale: the threshold sits mid-bar and each side spans its own range
  const valueAt = (f) => (f <= 0.5 ? lo + (thr - lo) * (f / 0.5) : thr + (hi - thr) * ((f - 0.5) / 0.5));
  for (let k = 0; k <= 24; k++) {
    el("stop", { offset: (k / 24 * 100).toFixed(1) + "%", "stop-color": colour(valueAt(k / 24)) }, grad);
  }
  el("rect", { x: 0, y: barY, width: W, height: barH, fill: `url(#${id})`, stroke: "#5C4A3D", "stroke-width": 0.8 }, svg);
  const tx = W / 2;
  el("line", { x1: tx, x2: tx, y1: barY - 5, y2: barY + barH + 5, stroke: "#9B2F3A", "stroke-width": 1.6 }, svg);
  const label = (x, text, anchor, fill = "#5C4A3D") => {
    const t = el("text", { x, y: barY + barH + 20, "text-anchor": anchor, fill,
      "font-family": "Spectral", "font-size": 12.5, style: "font-variant-numeric: lining-nums" }, svg);
    t.textContent = text;
  };
  label(0, formatValue(metric, lo), "start");
  label(W, formatValue(metric, hi), "end");
  label(tx, formatValue(metric, thr), "middle", "#9B2F3A");
  container.appendChild(svg);

  const left = metric.kind === "sim" ? (metric.rank === "species" ? "distinct species" : "other genus") : "same " + metric.rank;
  const right = metric.kind === "sim" ? "same " + metric.rank : "distinct species";
  const sides = document.createElement("p");
  sides.className = "legend-sides";
  sides.style.cssText = "display:flex;justify-content:space-between;font-style:italic;color:#5C4A3D;font-size:14px;margin-top:2px";
  const l = document.createElement("span"); l.textContent = "← " + left;
  const r = document.createElement("span"); r.textContent = right + " →";
  sides.append(l, r);
  container.appendChild(sides);

  const key = document.createElement("ul");
  key.className = "legend-key";
  const swatch = (draw, text) => {
    const li = document.createElement("li");
    const s = el("svg", { viewBox: "0 0 16 16", "aria-hidden": "true" });
    draw(s);
    const span = document.createElement("span");
    span.textContent = text;
    li.append(s, span);
    key.appendChild(li);
  };
  swatch((s) => {
    const hid = "lk" + (++uid);
    hatchDefs(s, hid);
    el("rect", { x: 0.5, y: 0.5, width: 15, height: 15, fill: `url(#${hid})`, stroke: "#5C4A3D", "stroke-width": 0.6 }, s);
  }, metric.key === "ani" ? "no value: too distant for FastANI" : "no value: no significant alignment");
  swatch((s) => {
    el("rect", { x: 0.5, y: 0.5, width: 15, height: 15, fill: "#E6DCCB", stroke: "#5C4A3D", "stroke-width": 0.6 }, s);
    el("line", { x1: 3, y1: 13, x2: 13, y2: 3, stroke: "#A8957D" }, s);
  }, "a genome against itself");
  const li = document.createElement("li");
  li.innerHTML = '<span style="width:16px;text-align:center;font-weight:600;font-size:12px">9</span>';
  const span = document.createElement("span");
  span.textContent = "bold figures: same " + metric.rank;
  li.appendChild(span);
  key.appendChild(li);
  container.appendChild(key);
}

// ---------- pairs near the threshold ----------

export function contestedPairs(metric, pairs, order, limit = 6) {
  const out = [];
  order.forEach((a, i) => order.forEach((b, j) => {
    if (j <= i) return;
    const vs = [pairs.get(a, b), pairs.get(b, a)].filter((v) => v !== null && v !== undefined);
    if (!vs.length) return;
    const v = vs.reduce((s, x) => s + x, 0) / vs.length;
    if (Math.abs(v - metric.threshold) <= metric.band) out.push({ a, b, v });
  }));
  return out.sort((p, q) => Math.abs(p.v - metric.threshold) - Math.abs(q.v - metric.threshold)).slice(0, limit);
}

// ---------- tree plate ----------

function niceStep(x) {
  const p = 10 ** Math.floor(Math.log10(x));
  const m = x / p;
  return (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p;
}

// Returns the branch paths (for the ink animation) and the labels group.
export function drawTree(container, layout) {
  container.textContent = "";
  const n = layout.leaves.length;
  const rowH = n <= 30 ? 30 : 20;
  const tipFont = n <= 30 ? 15 : 12.5;
  const labelW = Math.max(...layout.leaves.map((l) => textWidth(nameText(l.name, false), `italic ${tipFont}px Spectral`)));
  const avail = Math.max(container.clientWidth || 900, 640);
  const genusW = 150;
  const treeW = Math.max(260, Math.min(620, avail - labelW - genusW - 60));
  const top = 16; const bottom = 52;
  const W = Math.ceil(16 + treeW + 18 + labelW + 24 + genusW);
  const H = top + n * rowH + bottom;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
    "aria-label": `Phylogenomic tree of ${n} genomes` });
  const sx = (x) => 16 + (x / layout.maxX) * treeW;
  const sy = (y) => top + y * rowH + rowH / 2;
  const labelX = 16 + treeW + 18;

  const branches = el("g", {}, svg);
  const paths = drawBranches(branches, layout, sx, sy);

  const labels = el("g", {}, svg);
  layout.leaves.forEach((leaf) => {
    el("line", { class: "tr-leader", x1: sx(leaf.x) + 3, x2: labelX - 6, y1: sy(leaf.y), y2: sy(leaf.y) }, labels);
    const t = el("text", { class: "tr-tip", x: labelX, y: sy(leaf.y) + tipFont * 0.33, "font-size": tipFont }, labels);
    appendName(t, leaf.name, { svg: true });
  });

  // support values, when the tree has them (FastTree: 0-1)
  layout.nodes.forEach((nd) => {
    if (!nd.children.length || !nd.parent || nd.support === undefined) return;
    const val = nd.support <= 1 ? Math.round(nd.support * 100) : Math.round(nd.support);
    const t = el("text", { class: "tr-support", x: sx(nd.x) - 4, y: sy(nd.y) - 5, "text-anchor": "end" }, labels);
    t.textContent = val;
  });

  // genus brackets on the right margin
  const bx = labelX + labelW + 18;
  let start = 0;
  const leaves = layout.leaves;
  for (let k = 1; k <= leaves.length; k++) {
    if (k < leaves.length && genusOf(leaves[k].name) === genusOf(leaves[start].name)) continue;
    const y1 = sy(leaves[start].y) - rowH * 0.32; const y2 = sy(leaves[k - 1].y) + rowH * 0.32;
    el("path", { class: "tr-bracket", d: `M${bx - 5} ${y1}H${bx}V${y2}H${bx - 5}` }, labels);
    const t = el("text", { class: "tr-genus", x: bx + 10, y: (y1 + y2) / 2 + 5 }, labels);
    t.textContent = genusOf(leaves[start].name);
    start = k;
  }

  // scale bar
  if (!layout.cladogram) {
    const step = niceStep(layout.maxX / 4);
    const y = top + n * rowH + 26;
    el("path", { class: "tr-scale", d: `M${sx(0)} ${y - 4}V${y}H${sx(step)}V${y - 4}` }, labels);
    const t = el("text", { class: "tr-scale-label", x: sx(step) + 8, y: y + 4 }, labels);
    t.textContent = `${step} substitutions per site`;
  }

  container.appendChild(svg);
  return { svg, paths, labels };
}

// ---------- order without a tree ----------

// Average-linkage clustering on one pairwise measure, so that related
// genomes still sit next to each other as blocks when no tree was built.
export function clusterOrder(names, metric, pairs) {
  const far = metric.kind === "sim" ? 100 : 1;
  const d = (a, b) => {
    const vs = [pairs.get(a, b), pairs.get(b, a)].filter((v) => v !== null && v !== undefined);
    if (!vs.length) return far;
    const v = vs.reduce((s, x) => s + x, 0) / vs.length;
    return metric.kind === "sim" ? 100 - v : v;
  };
  let clusters = names.map((n) => ({ members: [n] }));
  const dist = (p, q) => {
    let s = 0;
    p.members.forEach((a) => q.members.forEach((b) => { s += d(a, b); }));
    return s / (p.members.length * q.members.length);
  };
  while (clusters.length > 1) {
    let best = [0, 1, Infinity];
    for (let i = 0; i < clusters.length; i++) {
      for (let j = i + 1; j < clusters.length; j++) {
        const x = dist(clusters[i], clusters[j]);
        if (x < best[2]) best = [i, j, x];
      }
    }
    const [i, j] = best;
    const merged = { members: clusters[i].members.concat(clusters[j].members) };
    clusters = clusters.filter((_, k) => k !== i && k !== j).concat([merged]);
  }
  return clusters.length ? clusters[0].members : [];
}
