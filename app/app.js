const state = {
  tool: "atom",
  weightMode: "weighted",
  latticeShape: "TLSG",
  columns: 9,
  rows: 7,
  nodes: new Map([
    ["3:3", { id: "3:3", q: 3, r: 3, weight: 1 }],
    ["4:3", { id: "4:3", q: 4, r: 3, weight: 1 }],
    ["5:3", { id: "5:3", q: 5, r: 3, weight: 1 }],
  ]),
  pins: ["3:3", "5:3"],
  selected: null,
  result: null,
};

const svg = document.querySelector("#lattice");
const resultSection = document.querySelector("#result-section");
const weightInput = document.querySelector("#weight-input");
const selectionEmpty = document.querySelector("#selection-empty");
const selectionFields = document.querySelector("#selection-fields");
const weightField = document.querySelector("#weight-field");
const fixedWeightNote = document.querySelector("#fixed-weight-note");
const fileInput = document.querySelector("#file-input");
const columnsInput = document.querySelector("#lattice-columns");
const rowsInput = document.querySelector("#lattice-rows");
const resizeLatticeButton = document.querySelector("#resize-lattice-button");
const computeButton = document.querySelector("#compute-button");
const computeLabel = document.querySelector("#compute-label");
const SVG_NS = "http://www.w3.org/2000/svg";

function latticeLayout() {
  const width = svg.clientWidth || 680;
  const height = svg.clientHeight || 500;
  const horizontalPadding = 44;
  const topPadding = 36;
  const bottomPadding = 72;
  const horizontalUnits = state.latticeShape === "KSG"
    ? state.columns - 1
    : state.columns - 0.5;
  const verticalUnits = state.latticeShape === "KSG"
    ? state.rows - 1
    : (state.rows - 1) * Math.sqrt(3) / 2;
  const availableHeight = height - topPadding - bottomPadding;
  const step = Math.min(
    (width - horizontalPadding * 2) / horizontalUnits,
    availableHeight / verticalUnits,
  );
  const latticeWidth = horizontalUnits * step;
  const latticeHeight = verticalUnits * step;
  return {
    step,
    originX: (width - latticeWidth) / 2,
    originY: topPadding + (availableHeight - latticeHeight) / 2,
  };
}

function point(q, r, layout) {
  const physical = physicalPoint(q, r);
  return {
    x: layout.originX + physical.x * layout.step,
    y: layout.originY + physical.y * layout.step,
  };
}

function physicalPoint(q, r) {
  if (state.latticeShape === "KSG") return { x: q, y: r };
  return { x: q + (r % 2 ? 0.5 : 0), y: r * Math.sqrt(3) / 2 };
}

function element(name, attributes = {}) {
  const item = document.createElementNS(SVG_NS, name);
  Object.entries(attributes).forEach(([key, value]) => item.setAttribute(key, value));
  return item;
}

function neighbors(a, b) {
  const pa = physicalPoint(a.q, a.r);
  const pb = physicalPoint(b.q, b.r);
  const radiusSquared = state.latticeShape === "KSG" ? 1.5 ** 2 : 1.1 ** 2;
  return (pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 < radiusSquared;
}

function edges() {
  const nodes = [...state.nodes.values()];
  const result = [];
  for (let sourceIndex = 0; sourceIndex < nodes.length; sourceIndex += 1) {
    for (let targetIndex = sourceIndex + 1; targetIndex < nodes.length; targetIndex += 1) {
      const source = nodes[sourceIndex];
      const target = nodes[targetIndex];
      if (neighbors(source, target)) result.push({ source, target });
    }
  }
  return result;
}

function drawLine(group, source, target, className, layout) {
  const a = point(source.q, source.r, layout);
  const b = point(target.q, target.r, layout);
  group.append(element("line", { x1: a.x, y1: a.y, x2: b.x, y2: b.y, class: className }));
}

function renderCanvas() {
  svg.replaceChildren();
  const layout = latticeLayout();
  const graphEdges = edges();
  const gridGroup = element("g");
  const edgeGroup = element("g");
  const pointGroup = element("g");
  const nodeGroup = element("g");
  const nodeRadius = Math.min(11, layout.step * 0.34);
  const haloRadius = Math.min(18, layout.step * 0.48);
  const labelOffset = nodeRadius + Math.min(7, layout.step * 0.12);
  const nodeStroke = Math.min(3, nodeRadius * 0.28);
  const selectedStroke = Math.min(4, nodeRadius * 0.36);
  const nodeLabelSize = Math.min(10, layout.step * 0.34);
  const pinLabelSize = Math.min(9, layout.step * 0.32);

  for (let r = 0; r < state.rows; r += 1) {
    for (let q = 0; q < state.columns; q += 1) {
      const current = { q, r };
      const forwardNeighbors = state.latticeShape === "KSG"
        ? [[q + 1, r], [q - 1, r + 1], [q, r + 1], [q + 1, r + 1]]
        : [[q + 1, r], [q, r + 1], [q + (r % 2 ? 1 : -1), r + 1]];
      forwardNeighbors.forEach(([nq, nr]) => {
        if (nq >= 0 && nq < state.columns && nr < state.rows) {
          drawLine(gridGroup, current, { q: nq, r: nr }, "grid-line", layout);
        }
      });
    }
  }

  graphEdges.forEach((edge) => drawLine(edgeGroup, edge.source, edge.target, "edge", layout));

  for (let r = 0; r < state.rows; r += 1) {
    for (let q = 0; q < state.columns; q += 1) {
      const id = `${q}:${r}`;
      const location = point(q, r, layout);
      if (!state.nodes.has(id)) {
        const dot = element("circle", { cx: location.x, cy: location.y, r: 3, class: "grid-point", tabindex: "0" });
        dot.addEventListener("click", () => interact(id, q, r));
        dot.addEventListener("keydown", (event) => {
          if (event.key === "Enter" || event.key === " ") interact(id, q, r);
        });
        pointGroup.append(dot);
      }
    }
  }

  state.nodes.forEach((node) => {
    const location = point(node.q, node.r, layout);
    const group = element("g");
    const pinIndex = state.pins.indexOf(node.id);
    group.style.setProperty("--node-stroke", `${nodeStroke}px`);
    group.style.setProperty("--selected-stroke", `${selectedStroke}px`);
    group.style.setProperty("--node-label-size", `${nodeLabelSize}px`);
    group.style.setProperty("--pin-label-size", `${pinLabelSize}px`);
    group.append(element("circle", { cx: location.x, cy: location.y, r: haloRadius, class: "node-halo" }));
    const circle = element("circle", {
      cx: location.x,
      cy: location.y,
      r: nodeRadius,
      class: `node${pinIndex >= 0 ? " pin" : ""}${state.selected === node.id ? " selected" : ""}`,
      tabindex: "0",
    });
    circle.addEventListener("click", () => interact(node.id, node.q, node.r));
    circle.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") interact(node.id, node.q, node.r);
    });
    group.append(circle);

    if (pinIndex >= 0) {
      const label = element("text", { x: location.x, y: location.y + 0.5, class: "pin-label" });
      label.textContent = `P${pinIndex + 1}`;
      group.append(label);
    }
    if (state.weightMode === "weighted") {
      const label = element("text", { x: location.x, y: location.y - labelOffset, class: "node-label" });
      label.textContent = node.weight;
      group.append(label);
    }
    nodeGroup.append(group);
  });

  svg.append(gridGroup, edgeGroup, pointGroup, nodeGroup);
  document.querySelector("#graph-summary").textContent =
    `${state.nodes.size} vertices · ${graphEdges.length} edges · ${state.pins.length} pins · ${state.latticeShape} · ${state.weightMode === "weighted" ? "Weighted" : "Unweighted"}`;
  document.querySelector("#lattice-dimensions").textContent =
    `${state.latticeShape} / ${String(state.columns).padStart(2, "0")} × ${String(state.rows).padStart(2, "0")}`;
  renderPins();
  renderInspector();
}

function interact(id, q, r) {
  if (state.tool === "erase") {
    state.nodes.delete(id);
    state.pins = state.pins.filter((pin) => pin !== id);
    if (state.selected === id) state.selected = null;
  } else if (state.tool === "pin") {
    if (!state.nodes.has(id)) state.nodes.set(id, { id, q, r, weight: 1 });
    togglePin(id);
    state.selected = id;
  } else if (!state.nodes.has(id)) {
    state.nodes.set(id, { id, q, r, weight: 1 });
    state.selected = id;
  } else {
    state.selected = id;
  }
  invalidateResult();
  renderCanvas();
}

function togglePin(id) {
  const index = state.pins.indexOf(id);
  if (index >= 0) state.pins.splice(index, 1);
  else state.pins.push(id);
}

function renderPins() {
  const container = document.querySelector("#pin-order");
  container.replaceChildren();
  state.pins.forEach((id, index) => {
    const chip = document.createElement("span");
    chip.className = "pin-chip";
    chip.textContent = `P${index + 1} · (${id.replace(":", ", ")})`;
    container.append(chip);
  });
}

function renderInspector() {
  const node = state.nodes.get(state.selected);
  selectionEmpty.hidden = Boolean(node);
  selectionFields.hidden = !node;
  if (!node) return;
  document.querySelector("#selected-coordinate").textContent = `q ${node.q}  ·  r ${node.r}`;
  weightInput.value = node.weight;
  weightField.hidden = state.weightMode === "unweighted";
  fixedWeightNote.hidden = state.weightMode === "weighted";
  const pinIndex = state.pins.indexOf(node.id);
  document.querySelector("#toggle-pin-button").textContent =
    pinIndex >= 0 ? `Remove P${pinIndex + 1}` : "Set as next pin";
}

function setTool(tool) {
  state.tool = tool;
  document.querySelectorAll(".tool").forEach((button) => {
    const active = button.dataset.tool === tool;
    button.classList.toggle("active", active);
    button.setAttribute("aria-pressed", active);
  });
}

function setWeightMode(weightMode) {
  state.weightMode = weightMode;
  document.querySelectorAll("[data-weight-mode]").forEach((button) => {
    const active = button.dataset.weightMode === weightMode;
    button.classList.toggle("active", active);
    button.setAttribute("aria-pressed", active);
  });
  computeLabel.textContent =
    weightMode === "weighted" ? "Compute Ground States" : "Compute Reduced Alpha Tensor";
  document.querySelector("#computation-label").textContent =
    weightMode === "weighted" ? "Ground-state computation" : "Alpha-tensor computation";
  document.querySelector("#computation-title").textContent =
    weightMode === "weighted" ? "Direct MIS solver" : "Reduced alpha tensor";
  const computationHelper = document.querySelector("#computation-helper");
  if (weightMode === "weighted") {
    computationHelper.textContent = "Computes weighted ground states. Pin projections are reported in P1 → Pn order.";
  } else {
    computationHelper.innerHTML = `Treats pins as open vertices and computes
      <math class="inline-math" aria-label="alpha tilde of R">
        <mover accent="true"><mi>α</mi><mo>~</mo></mover><mo>(</mo><mi>R</mi><mo>)</mo>
      </math>. Vertex weights are ignored.`;
  }
  state.result = null;
  renderIdle();
  renderCanvas();
}

function setLatticeShape(latticeShape) {
  state.latticeShape = latticeShape;
  document.querySelectorAll("[data-lattice-shape]").forEach((button) => {
    const active = button.dataset.latticeShape === latticeShape;
    button.classList.toggle("active", active);
    button.setAttribute("aria-pressed", active);
  });
  invalidateResult();
  renderCanvas();
}

function setLatticeSize() {
  const columns = Number(columnsInput.value);
  const rows = Number(rowsInput.value);
  if (columns < 2 || columns > 20 || rows < 2 || rows > 20) {
    columnsInput.value = state.columns;
    rowsInput.value = state.rows;
    showToast("Lattice dimensions must be between 2 and 20");
    return;
  }
  const outside = [...state.nodes.values()].find((node) => node.q >= columns || node.r >= rows);
  if (outside) {
    columnsInput.value = state.columns;
    rowsInput.value = state.rows;
    showToast(`Remove the vertex at (${outside.q}, ${outside.r}) before shrinking the lattice`);
    return;
  }
  state.columns = columns;
  state.rows = rows;
  state.result = null;
  renderIdle();
  renderCanvas();
}

function renderIdle() {
  const description = state.weightMode === "weighted"
    ? "The solver enumerates maximal independent sets and returns the maximum-energy states."
    : `The solver contracts the independent-set tensor network and compactifies
      <math class="inline-math" aria-label="alpha tilde of R">
        <mover accent="true"><mi>α</mi><mo>~</mo></mover><mo>(</mo><mi>R</mi><mo>)</mo>
      </math>.`;
  resultSection.innerHTML = `
    <div class="result-idle">
      <svg class="solver-mark" viewBox="0 0 96 64" aria-hidden="true">
        <path d="M18 32 33 14h30l15 18-15 18H33Z M18 32h60 M33 14l30 36 M63 14 33 50"/>
        <circle cx="18" cy="32" r="4"/>
        <circle cx="33" cy="14" r="4"/>
        <circle cx="63" cy="14" r="4"/>
        <circle cx="78" cy="32" r="4"/>
        <circle cx="63" cy="50" r="4"/>
        <circle cx="33" cy="50" r="4"/>
        <circle class="solver-mark-focus" cx="48" cy="32" r="8"/>
        <circle class="solver-mark-core" cx="48" cy="32" r="3"/>
      </svg>
      <h3>Ready to compute</h3>
      <p>${description}</p>
    </div>
  `;
}

function invalidateResult() {
  state.result = null;
  renderIdle();
}

function payload() {
  return {
    model: "rydberg",
    weight_mode: state.weightMode,
    lattice: {
      shape: state.latticeShape,
      columns: state.columns,
      rows: state.rows,
    },
    nodes: [...state.nodes.values()],
    pins: state.pins,
  };
}

async function compute() {
  if (computeButton.disabled) return;
  computeButton.disabled = true;
  const progress = state.weightMode === "weighted"
    ? "Enumerating maximal independent sets…"
    : "Contracting the alpha tensor…";
  resultSection.innerHTML = `<div class="result-idle"><h3>Computing…</h3><p>${progress}</p></div>`;
  try {
    const response = await fetch("/api/compute", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload()),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.error);
    state.result = data;
    renderResult(data);
  } catch (error) {
    resultSection.innerHTML = `<div class="error-card"><strong>Could not compute</strong><p>${escapeHtml(error.message)}</p></div>`;
  } finally {
    computeButton.disabled = false;
  }
}

function renderResult(data) {
  if (data.operation === "reduced_alpha_tensor") {
    resultSection.innerHTML = `
      <div class="verdict valid">
        <strong>Reduced alpha tensor computed</strong>
        <span>${data.boundary_count} open vertices · ${data.tensor.length} boundary configurations</span>
      </div>
      <div class="metric-grid">
        <div class="metric"><span>Vertices</span><strong>${data.vertex_count}</strong></div>
        <div class="metric"><span>Edges</span><strong>${data.edge_count}</strong></div>
      </div>
      <p class="state-heading">
        <math class="inline-math" aria-label="alpha tilde of R">
          <mover accent="true"><mi>α</mi><mo>~</mo></mover><mo>(</mo><mi>R</mi><mo>)</mo>
        </math>
      </p>
      <table class="tensor-table">
        <thead><tr><th>Boundary state</th><th>Value</th></tr></thead>
        <tbody>${data.tensor.map((entry) => `
          <tr><td>${entry.configuration || "∅"}</td><td>${entry.value}</td></tr>
        `).join("")}</tbody>
      </table>
    `;
    return;
  }

  const pinProjection = state.pins.length > 0
    ? `Pin projections: ${data.observed.join(", ")}`
    : "No pins marked";
  resultSection.innerHTML = `
    <div class="verdict valid">
      <strong>Ground states computed</strong><span>${pinProjection}</span>
    </div>
    <div class="metric-grid">
      <div class="metric"><span>Maximum energy</span><strong>${formatNumber(data.max_energy)}</strong></div>
      <div class="metric"><span>State space</span><strong>${data.state_count}</strong></div>
      <div class="metric"><span>Vertices</span><strong>${data.vertex_count}</strong></div>
      <div class="metric"><span>Degeneracy</span><strong>${data.ground_states.length}</strong></div>
    </div>
    <p class="state-heading">Ground states</p>
    ${data.ground_states.map((item) => `
      <div class="state-card">
        <strong>pins · ${item.pins || "—"}</strong>
        <small>σ = ${item.configuration}</small>
      </div>
    `).join("")}
  `;
}

function formatNumber(value) {
  return Number.isInteger(value) ? String(value) : Number(value).toFixed(4).replace(/0+$/, "").replace(/\.$/, "");
}

function escapeHtml(value) {
  const item = document.createElement("div");
  item.textContent = value;
  return item.innerHTML;
}

function showToast(message) {
  const toast = document.querySelector("#toast");
  toast.textContent = message;
  toast.classList.add("show");
  window.setTimeout(() => toast.classList.remove("show"), 1600);
}

function exportJson() {
  const blob = new Blob([JSON.stringify(payload(), null, 2)], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `${state.latticeShape.toLowerCase()}-gadget.json`;
  link.click();
  URL.revokeObjectURL(link.href);
  showToast("Gadget exported");
}

function importJson(file) {
  const reader = new FileReader();
  reader.addEventListener("load", () => {
    const data = JSON.parse(reader.result);
    state.columns = data.lattice.columns;
    state.rows = data.lattice.rows;
    columnsInput.value = state.columns;
    rowsInput.value = state.rows;
    state.nodes = new Map(data.nodes.map((node) => [node.id, node]));
    state.pins = data.pins;
    state.selected = null;
    setLatticeShape(data.lattice.shape);
    setWeightMode(data.weight_mode);
    showToast("Gadget imported");
  });
  reader.readAsText(file);
}

document.querySelectorAll(".tool").forEach((button) => {
  button.addEventListener("click", () => setTool(button.dataset.tool));
});
document.querySelectorAll("[data-weight-mode]").forEach((button) => {
  button.addEventListener("click", () => setWeightMode(button.dataset.weightMode));
});
document.querySelectorAll("[data-lattice-shape]").forEach((button) => {
  button.addEventListener("click", () => setLatticeShape(button.dataset.latticeShape));
});
resizeLatticeButton.addEventListener("click", setLatticeSize);
weightInput.addEventListener("change", () => {
  state.nodes.get(state.selected).weight = Number(weightInput.value);
  invalidateResult();
  renderCanvas();
});
document.querySelector("#toggle-pin-button").addEventListener("click", () => {
  togglePin(state.selected);
  invalidateResult();
  renderCanvas();
});
computeButton.addEventListener("click", compute);
document.querySelector("#clear-button").addEventListener("click", () => {
  state.nodes.clear();
  state.pins = [];
  state.selected = null;
  invalidateResult();
  renderCanvas();
});
document.querySelector("#export-button").addEventListener("click", exportJson);
document.querySelector("#import-button").addEventListener("click", () => fileInput.click());
fileInput.addEventListener("change", () => {
  if (fileInput.files[0]) importJson(fileInput.files[0]);
  fileInput.value = "";
});
window.addEventListener("resize", renderCanvas);
window.addEventListener("keydown", (event) => {
  if (event.target.matches("input, textarea")) return;
  if (event.key.toLowerCase() === "a") setTool("atom");
  if (event.key.toLowerCase() === "p") setTool("pin");
  if (event.key.toLowerCase() === "e") setTool("erase");
  if (event.key === "Enter" && (event.metaKey || event.ctrlKey)) compute();
});

renderCanvas();
