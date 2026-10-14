// The inspector's model summary: on the Parts tab, a part document's
// bodies (an assembly's tree is assembly-panel.js) with the mesh
// closure chip, and the clearance check for assemblies; on the Checks
// tab, one row per check of the model on screen — the features served,
// mesh closure, blend continuity against what was asked, deferred seam
// certification, report warnings — each ok, warn or bad. The first
// model of a kind opens the tab that suits it (Parts for an assembly,
// Features for a part).
import { showInspectorTab } from "./layout.js";
import { deferredSeams } from "./report-table.js";
import { paletteCss } from "./shading.js";

const $ = (id) => document.getElementById(id);

// The kind of the model last shown ("assembly" | "part").
let shownKind = null;

/** The check rows for `envelope`: `{ level, title, detail }` with
 *  level "ok" | "warn" | "bad". */
export function modelChecks(envelope) {
  const rows = [];
  const report = envelope.report;
  const features = report?.features ?? [];
  if (report) {
    const failed = report.failures ?? 0;
    rows.push(failed
      ? { level: "bad", title: `${failed} feature(s) failed`, detail: "the Features tab names them" }
      : { level: "ok", title: `${features.length} feature(s) served`, detail: envelope.outcome === "served" ? "the model replays cleanly" : `outcome: ${String(envelope.outcome).replace("_", " ")}` });
  }
  const bodies = envelope.bodies ?? [];
  if (envelope.outcome === "served" && bodies.length) {
    const open = bodies.reduce((n, b) => n + (b.open_edges ?? 0), 0);
    rows.push(open
      ? { level: "warn", title: `${open} open mesh edge(s)`, detail: "the display mesh is not closed" }
      : { level: "ok", title: "Watertight", detail: `${bodies.length} body(ies), closed display meshes` });
  }
  const blends = features.filter((f) => f.blend?.continuity?.requested);
  if (blends.length) {
    const downgraded = report.downgraded_blends ?? 0;
    rows.push(downgraded
      ? { level: "warn", title: `${downgraded} blend(s) downgraded`, detail: "achieved less continuity than requested" }
      : { level: "ok", title: "Blend continuity as requested", detail: `${blends.length} blend(s)` });
    const pending = deferredSeams(report);
    rows.push(pending
      ? { level: "warn", title: "Certifying continuity…", detail: `${pending} seam(s) measured from samples so far` }
      : { level: "ok", title: "Continuity certified", detail: "every blend seam carries a certificate" });
  }
  if (report && !report.is_clean && !(report.failures ?? 0)) {
    const warned = features.filter((f) => Array.isArray(f.warnings) && f.warnings.length).length;
    rows.push({ level: "warn", title: "Built with warnings", detail: warned ? `${warned} feature(s) warn: the Features tab lists them` : "see the Features tab" });
  }
  return rows;
}

function renderChecks(envelope) {
  const list = $("checks-list");
  list.replaceChildren();
  for (const row of modelChecks(envelope)) {
    const li = document.createElement("li");
    li.className = `check-row check-${row.level}`;
    const mark = document.createElement("span");
    mark.className = "check-mark";
    mark.setAttribute("aria-label", row.level);
    const title = document.createElement("span");
    title.className = "check-title";
    title.textContent = row.title;
    const detail = document.createElement("span");
    detail.className = "check-detail";
    detail.textContent = row.detail;
    li.append(mark, title, detail);
    list.append(li);
  }
}

/** The closure chips' text and class for a served envelope. */
function closure(envelope) {
  if (envelope.outcome !== "served" || !envelope.bodies?.length) return null;
  const open = envelope.bodies.reduce((n, b) => n + (b.open_edges ?? 0), 0);
  return open ? { text: `${open} open edges`, css: "chip chip-warn closure-chip" } : { text: "watertight", css: "chip chip-ok closure-chip" };
}

function renderBodies(envelope) {
  const section = $("bodies-section");
  const assembly = envelope.assembly != null;
  section.hidden = assembly;
  if (assembly) return;
  const bodies = envelope.bodies ?? [];
  $("bodies-summary").textContent = bodies.length ? `${bodies.length} body(ies) · ${envelope.features ?? 0} feature(s)` : "no bodies";
  const list = $("bodies-list");
  list.replaceChildren();
  bodies.forEach((body, i) => {
    const row = document.createElement("li");
    row.className = "tree-row";
    const dot = document.createElement("span");
    dot.className = "tree-dot";
    dot.style.background = paletteCss(i);
    const name = document.createElement("span");
    name.className = "tree-name";
    name.textContent = body.label ?? `body ${i + 1}`;
    const kind = document.createElement("span");
    kind.className = "tree-kind";
    kind.textContent = `${body.face_count} faces`;
    row.append(dot, name, kind);
    list.append(row);
  });
}

/** Render the model summary for `envelope` (the model on screen). */
export function renderModelPanel(envelope) {
  renderBodies(envelope);
  const chip = closure(envelope);
  for (const el of document.querySelectorAll(".closure-chip")) {
    el.hidden = !chip;
    if (chip) {
      el.textContent = chip.text;
      el.className = chip.css;
    }
  }
  $("clearance-row").hidden = envelope.assembly == null;
  renderChecks(envelope);
  const kind = envelope.assembly ? "assembly" : "part";
  if (envelope.outcome === "served" && kind !== shownKind) {
    shownKind = kind;
    showInspectorTab(kind === "assembly" ? "parts" : "features");
  }
}
