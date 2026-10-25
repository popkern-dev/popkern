// The Simulate tab (PK-948): hold faces and load faces picked in the
// viewer, a material, a force (or, for heat, a held temperature and a
// heat input), Solve. The panel writes an analysis setup file
// (docs/architecture/analysis-integration.md) that names each picked
// face by its stable name, so the same setup solves again after a
// parameter edit, and the worker solves it on the live model. The
// result is a summary card and a field overlay on the part.

import { FIELDS, FieldOverlay, rampGradient } from "./simulate-overlay.js";

const $ = (id) => document.getElementById(id);

/** A TOML string literal. */
function quoted(text) {
  return `"${String(text).replace(/\\/g, "\\\\").replace(/"/g, '\\"')}"`;
}

/** The setup file for the panel's state: `{analysis, material, hold,
 *  load, force, direction, temperature, power}`; `hold` and `load` are
 *  stable face names. */
export function setupText(s) {
  const lines = [];
  if (s.analysis === "thermal") lines.push('analysis = "thermal"');
  lines.push('units = "mm"', `material = ${quoted(s.material)}`, "", "[faces]");
  s.hold.forEach((name, i) => lines.push(`hold_${i + 1} = { name = ${quoted(name)}, expect = { at_least = 1 } }`));
  s.load.forEach((name, i) => lines.push(`load_${i + 1} = { name = ${quoted(name)}, expect = { at_least = 1 } }`));
  const list = (prefix, names) => `[${names.map((_, i) => quoted(`${prefix}_${i + 1}`)).join(", ")}]`;
  if (s.analysis === "thermal") {
    lines.push("", "[[fixed_temperature]]", `faces = ${list("hold", s.hold)}`, `temperature_c = ${Number(s.temperature)}`);
    lines.push("", "[[heat_flux]]", `faces = ${list("load", s.load)}`, `power_w = ${Number(s.power)}`);
  } else {
    lines.push("", "[[support]]", `faces = ${list("hold", s.hold)}`, 'kind = "fixed"');
    const direction = s.direction === "face_normal" ? '"face_normal"' : quoted(s.direction);
    lines.push("", "[[load]]", `faces = ${list("load", s.load)}`, `force_n = ${Number(s.force)}`, `direction = ${direction}`);
  }
  return `${lines.join("\n")}\n`;
}

/** The project file a panel setup is saved as. */
export function setupPath(analysis) {
  return `simulate/${analysis === "thermal" ? "thermal" : "structural"}.setup.toml`;
}

// The saved file's second line: the panel's own state, so reopening the
// project restores the picks and inputs (the TOML below it is what
// solves).
const PANEL_PREFIX = "# panel: ";

/** The saved setup file: a header, the panel state, the setup. */
export function savedText(state) {
  const panel = { analysis: state.analysis, material: state.material, hold: state.hold, load: state.load, force: state.force, direction: state.direction, temperature: state.temperature, power: state.power };
  return `# popkern live Simulate setup (PK-948): faces by their stable popkern names, so it solves this project's part after an edit.\n${PANEL_PREFIX}${JSON.stringify(panel)}\n${setupText(state)}`;
}

/** The panel state saved in `text`, or null. */
export function panelState(text) {
  const line = String(text ?? "").split("\n").find((l) => l.startsWith(PANEL_PREFIX));
  if (!line) return null;
  try {
    const state = JSON.parse(line.slice(PANEL_PREFIX.length));
    return Array.isArray(state.hold) && Array.isArray(state.load) ? state : null;
  } catch {
    return null;
  }
}

/** Number text for the summary card. */
function n(x, digits = 3) {
  if (x == null || !Number.isFinite(x)) return "–";
  return Math.abs(x) >= 1e-3 || x === 0 ? x.toPrecision(digits) : x.toExponential(2);
}

/** The summary card rows of a solved `simulation`. */
export function summaryRows(simulation) {
  const r = simulation.result ?? {};
  const rows = [];
  if (simulation.analysis === "thermal") {
    rows.push(["hottest point", `${n(r.max_temperature_c)} °C`]);
    rows.push(["heat in", `${n(r.flux_in_w)} W`]);
  } else {
    rows.push(["max deflection", `${n(r.max_displacement_mm)} mm`]);
    // The worst peak over the solids, with its face's stable name.
    const peaks = (r.assessment?.peaks ?? []).map((p) => p.von_mises).filter(Boolean);
    const peak = peaks.reduce((a, b) => (b.value > (a?.value ?? -Infinity) ? b : a), null);
    if (peak) rows.push(["peak stress", `${n(peak.value)} ± ${n(peak.value * peak.error, 2)} MPa${peak.face_label ? ` on ${peak.face_label}` : ""}`]);
    // The governing check sets the factor; its band covers the mesh
    // error estimate and the material's scatter.
    const safety = r.assessment?.safety;
    const check = safety?.checks?.[safety.governing];
    if (check) {
      const band = check.band ? ` (${n(check.band[0], 2)}–${n(check.band[1], 2)})` : "";
      rows.push(["safety factor", `${n(check.factor, 2)}${band}, ${check.check.replace(/([a-z])([A-Z])/g, "$1 $2").toLowerCase()}`]);
    }
    if (r.global_error != null) rows.push(["estimated error", `${n(100 * r.global_error, 2)}%`]);
  }
  rows.push(["material", r.material ?? "–"]);
  for (const w of simulation.solve?.warnings ?? []) {
    const kind = typeof w === "string" ? w : Object.keys(w)[0];
    rows.push(["note", kind.replace(/([a-z])([A-Z])/g, "$1 $2").toLowerCase()]);
  }
  return rows;
}

export class SimulatePanel {
  constructor(page) {
    this.page = page;
    this.overlay = new FieldOverlay(page.viewer);
    this.hold = [];
    this.load = [];
    this.result = null;
    this.picking = null;
    $("sim-pick-hold").addEventListener("click", () => this.pick("hold"));
    $("sim-pick-load").addEventListener("click", () => this.pick("load"));
    $("sim-analysis").addEventListener("change", () => this.renderMode());
    $("sim-solve").addEventListener("click", () => this.solve());
    $("sim-clear").addEventListener("click", () => this.clearResult());
    $("sim-field").addEventListener("change", () => this.draw());
    $("sim-deform").addEventListener("input", () => this.draw());
    $("sim-legend-bar").style.background = rampGradient();
    this.renderMode();
  }

  /** Called when a new model is on screen: a field solved for the old
   *  one must not be drawn on it (the picked names are kept, so the
   *  same setup can be solved again). */
  modelChanged() {
    this.clearResult();
  }

  async loadMaterials() {
    if (this.materials) return;
    const reply = await this.page.host.materials();
    if (!reply?.ok) return;
    this.materials = reply.materials;
    const select = $("sim-material");
    select.innerHTML = "";
    for (const m of reply.materials) {
      const option = document.createElement("option");
      option.value = m.name;
      option.textContent = m.name;
      select.append(option);
    }
    if (this.pendingMaterial) select.value = this.pendingMaterial;
  }

  renderMode() {
    const thermal = $("sim-analysis").value === "thermal";
    $("sim-hold-title").textContent = thermal ? "Held temperature on" : "Held faces";
    $("sim-load-title").textContent = thermal ? "Heat in through" : "Loaded faces";
    $("sim-structural").hidden = thermal;
    $("sim-thermal").hidden = !thermal;
    this.renderFaces();
  }

  renderFaces() {
    for (const [role, names] of [["hold", this.hold], ["load", this.load]]) {
      const list = $(`sim-${role}`);
      list.innerHTML = "";
      for (const name of names) {
        const row = document.createElement("li");
        row.className = "tree-row";
        const label = document.createElement("span");
        label.className = "tree-name";
        label.textContent = name;
        const remove = document.createElement("button");
        remove.type = "button";
        remove.className = "sim-remove";
        remove.textContent = "×";
        remove.title = `Remove ${name}`;
        remove.setAttribute("aria-label", `Remove ${name}`);
        remove.addEventListener("click", () => {
          names.splice(names.indexOf(name), 1);
          this.renderFaces();
        });
        row.append(label, remove);
        list.append(row);
      }
      $(`sim-pick-${role}`).classList.toggle("armed", this.picking === role);
    }
  }

  /** Arm the viewer: each face clicked joins `role` until the button is
   *  pressed again. A face without a stable name is refused (a raw id
   *  would point at another face after an edit). */
  async pick(role) {
    if (this.picking === role) {
      this.picking = null;
      this.page.viewer.pickHook = null;
      this.renderFaces();
      return;
    }
    if (!(await this.page.ensureLive())) return;
    this.clearResult();
    this.picking = role;
    const arm = () => {
      this.page.viewer.pickHook = (hit) => {
        const label = hit.mesh.userData.body.face_labels?.[hit.face];
        const names = role === "hold" ? this.hold : this.load;
        if (!label || label === "unnamed") {
          this.page.status("simulate: that face has no stable name, so a setup cannot keep it; pick another", true);
        } else if (!names.includes(label)) {
          names.push(label);
        }
        this.renderFaces();
        if (this.picking === role) arm();
      };
    };
    arm();
    this.renderFaces();
    this.page.status(`simulate: click the faces to ${role === "hold" ? "hold" : "load"}; press the button again when done`);
  }

  state() {
    return {
      analysis: $("sim-analysis").value,
      material: $("sim-material").value,
      hold: this.hold,
      load: this.load,
      force: $("sim-force").value,
      direction: $("sim-direction").value,
      temperature: $("sim-temperature").value,
      power: $("sim-power").value,
    };
  }

  /** Restore the panel from a project's saved setup (the structural
   *  one first), or reset it for a project without one. */
  restore(files) {
    const saved = panelState(files?.[setupPath("structural")]) ?? panelState(files?.[setupPath("thermal")]);
    this.clearResult();
    this.hold = saved?.hold?.slice() ?? [];
    this.load = saved?.load?.slice() ?? [];
    if (saved) {
      $("sim-analysis").value = saved.analysis === "thermal" ? "thermal" : "structural";
      for (const [id, key] of [["sim-force", "force"], ["sim-direction", "direction"], ["sim-temperature", "temperature"], ["sim-power", "power"]]) {
        if (saved[key] != null) $(id).value = saved[key];
      }
      this.pendingMaterial = saved.material;
      if (this.materials && saved.material) $("sim-material").value = saved.material;
    }
    this.renderMode();
  }

  /** Save the panel's setup with the open project (not a read-only
   *  example). */
  persist(state) {
    const projects = this.page.projects;
    const p = projects?.current;
    if (!p || p.readOnly) return;
    const path = setupPath(state.analysis);
    const text = savedText(state);
    if (p.files[path] === text) return;
    p.files = { ...p.files, [path]: text };
    p.unsaved = true;
    projects.save();
  }

  async solve() {
    if (!this.hold.length || !this.load.length) {
      $("sim-summary").textContent = "Pick at least one held face and one loaded face.";
      return null;
    }
    if (!(await this.page.ensureLive())) return null;
    this.picking = null;
    this.page.viewer.pickHook = null;
    this.renderFaces();
    const good = this.page.lastGood;
    const state = this.state();
    const setup = setupText(state);
    this.persist(state);
    $("sim-summary").textContent = "solving…";
    this.page.status("simulate: solving…");
    const started = performance.now();
    const reply = await this.page.host.simulate(setup);
    // A newer model on screen drops this reply.
    if (this.page.lastGood !== good) return null;
    this.render(reply, performance.now() - started);
    return reply;
  }

  render(reply, ms) {
    const card = $("sim-card");
    card.innerHTML = "";
    if (!reply?.ok) {
      $("sim-summary").textContent = reply?.outcome?.detail ?? "the solve failed";
      this.page.status(`simulate: ${$("sim-summary").textContent}`, true);
      return;
    }
    const s = reply.simulation;
    if (s.status !== "solved") {
      $("sim-summary").textContent = `Declined (${s.error.kind}): ${s.error.message}`;
      this.page.status(`simulate: declined, ${s.error.kind}`, true);
      return;
    }
    this.result = { meshes: reply.meshes, mmPerUnit: s.mm_per_unit ?? 1, analysis: s.analysis };
    $("sim-summary").textContent = `Solved in ${(ms / 1000).toFixed(1)} s on ${s.solve?.cells_active ?? "?"} cells.`;
    for (const [what, value] of summaryRows(s)) {
      const row = document.createElement("li");
      const label = document.createElement("span");
      label.className = "sim-what";
      label.textContent = what;
      const text = document.createElement("span");
      text.className = "sim-value";
      text.textContent = value;
      row.append(label, text);
      card.append(row);
    }
    const select = $("sim-field");
    select.innerHTML = "";
    for (const name of reply.meshes[0]?.field_names ?? []) {
      const option = document.createElement("option");
      option.value = name;
      option.textContent = FIELDS[name]?.label ?? name;
      select.append(option);
    }
    // Stress first for strength (where it breaks), else the one field.
    const names = reply.meshes[0]?.field_names ?? [];
    select.value = names.includes("von_mises_mpa") ? "von_mises_mpa" : names[0] ?? "";
    $("sim-deform-row").hidden = s.analysis === "thermal";
    $("sim-view").hidden = false;
    this.draw();
    this.page.status(`simulate: solved in ${(ms / 1000).toFixed(1)} s`);
  }

  /** Draw the selected field (and deformation) of the last result. */
  draw() {
    if (!this.result) return;
    const field = $("sim-field").value;
    const deform = this.result.analysis === "thermal" ? 0 : Number($("sim-deform").value);
    const [lo, hi] = this.overlay.show(this.result.meshes, field, deform, this.result.mmPerUnit);
    const unit = FIELDS[field]?.unit ?? "";
    $("sim-legend-title").textContent = FIELDS[field]?.label ?? field;
    $("sim-legend-min").textContent = `${n(lo)} ${unit}`;
    $("sim-legend-max").textContent = `${n(hi)} ${unit}`;
    $("sim-legend").hidden = false;
    document.body.classList.add("sim-shown");
  }

  clearResult() {
    this.result = null;
    this.overlay.clear();
    $("sim-legend").hidden = true;
    document.body.classList.remove("sim-shown");
    $("sim-view").hidden = true;
    $("sim-card").innerHTML = "";
    $("sim-summary").textContent = "";
  }
}
