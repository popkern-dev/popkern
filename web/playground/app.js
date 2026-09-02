// The playground page: editor → worker → viewer, with the report,
// parameters, downloads and show-as-code panels. Everything the user
// types stays in this tab (design §7: client-side only; the worker is
// terminated, never trusted, when a script defeats its budgets).
import { Editor } from "./editor.js";
import { Viewer } from "./viewer.js";
import { WorkerHost } from "./worker-host.js";

const BUDGET = { max_operations: 20_000_000, timeout_ms: 20_000 };
const IDLE_RUN_MS = 900;
const CANCEL_GRACE_MS = 4_000;
const DRAFT_KEY = "popkern-playground-draft";

const DEFAULT_SCRIPT = `// A boss plate: plate, boss, union, rounds.
let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch")
             .extrude(blind(0.4));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
             .extrude(blind(0.8)));
plate = plate.union(boss);
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{ continuity: "g1", label: "rounds" });
m.assert_watertight();
print(m.feature("rounds").describe());
m
`;

const $ = (id) => document.getElementById(id);

function storage(action) {
  try {
    return action(window.localStorage);
  } catch {
    return null;
  }
}

class Page {
  constructor() {
    this.host = new WorkerHost(() => new Worker("./worker.js", { type: "module" }), BUDGET);
    this.host.onCrash = (detail) => this.status(`worker error: ${detail}`, true);
    this.viewer = new Viewer($("viewer"), (pick) => this.showPick(pick));
    this.editor = new Editor($("editor"), {
      onChange: (text) => this.changed(text),
      onRun: () => this.run(),
    });
    this.idleTimer = null;
    this.cancelTimer = null;
    this.lastGood = null;
    this.running = false;
    this.wire();
    this.editor.value = storage((s) => s.getItem(DRAFT_KEY)) || DEFAULT_SCRIPT;
    this.loadExamples();
    this.status("loading the kernel…");
    this.host
      .start()
      .then((info) => {
        this.kernelReadyMs = performance.now();
        this.status(`kernel ready in ${this.kernelReadyMs.toFixed(0)} ms (protocol ${info.protocol_version}${info.clock_available ? "" : ", no clock: wall budget off"})`);
        this.run();
      })
      .catch((error) => this.status(`kernel failed to load: ${error.message}`, true));
  }

  wire() {
    $("run").addEventListener("click", () => this.run());
    $("cancel").addEventListener("click", () => this.cancel("cancelled by the user"));
    $("fit").addEventListener("click", () => this.viewer.fitAgain());
    $("edges").addEventListener("change", (e) => {
      this.viewer.toggleEdges(e.target.checked);
      if (this.lastGood) this.viewer.show(this.lastGood.envelope.bodies, this.lastGood.buffers);
    });
    $("download-step").addEventListener("click", () => this.download("step"));
    $("download-pkn").addEventListener("click", () => this.download("pkn"));
    $("show-code").addEventListener("click", () => this.showCode());
    $("code-close").addEventListener("click", () => ($("code-dialog").hidden = true));
    $("code-use").addEventListener("click", () => {
      this.editor.value = $("code-text").textContent;
      $("code-dialog").hidden = true;
      this.run();
    });
    $("examples").addEventListener("change", (e) => this.loadExample(e.target.value));
    $("reset").addEventListener("click", () => {
      this.editor.value = DEFAULT_SCRIPT;
      this.run();
    });
  }

  status(text, failed = false) {
    const el = $("status");
    el.textContent = text;
    el.classList.toggle("status-failed", failed);
  }

  changed(text) {
    storage((s) => s.setItem(DRAFT_KEY, text));
    clearTimeout(this.idleTimer);
    this.idleTimer = setTimeout(() => this.run(), IDLE_RUN_MS);
  }

  async run() {
    clearTimeout(this.idleTimer);
    if (!this.host.info) return;
    if (this.running) {
      this.rerun = true;
      return;
    }
    this.running = true;
    this.rerun = false;
    $("run").disabled = true;
    this.status("building…");
    this.cancelTimer = setTimeout(() => this.cancel(`no answer within ${(BUDGET.timeout_ms + CANCEL_GRACE_MS) / 1000}s: the worker was stopped`), BUDGET.timeout_ms + CANCEL_GRACE_MS);
    const reply = await this.host.run(this.editor.value, "playground.rhai");
    clearTimeout(this.cancelTimer);
    this.running = false;
    $("run").disabled = false;
    this.receive(reply);
    if (this.rerun) this.run();
    return reply;
  }

  async cancel(reason) {
    clearTimeout(this.cancelTimer);
    this.status("stopping the worker…", true);
    await this.host.cancel(reason);
    this.running = false;
    $("run").disabled = false;
    this.status(`${reason} — kernel restarted`, true);
  }

  receive(reply) {
    if (!reply.ok) {
      this.showOutcome(reply.outcome);
      return;
    }
    const { envelope, buffers } = reply;
    // A run that raised carries no report: the panels keep describing
    // the last good model, like the viewer does.
    if (envelope.report) {
      this.renderReport(envelope);
      this.renderParams(envelope.params);
    }
    if (envelope.outcome === "served") {
      this.lastGood = { envelope, buffers };
      this.viewer.show(envelope.bodies, buffers);
      this.editor.mark(null);
      $("overlay").hidden = true;
      const ms = envelope.elapsed_ms == null ? "" : ` in ${envelope.elapsed_ms.toFixed(0)} ms`;
      const open = envelope.bodies.reduce((n, b) => n + b.open_edges, 0);
      const first = this.firstRenderMs == null ? ` — first render ${(this.firstRenderMs = performance.now()).toFixed(0)} ms after page start` : "";
      this.status(`served: ${envelope.features} feature(s), ${envelope.bodies.length} body(ies)${ms}${open ? `, ${open} open mesh edges` : ", watertight"}${first}`);
    } else {
      this.showOutcome(envelope);
    }
  }

  showOutcome(outcome) {
    const overlay = $("overlay");
    overlay.hidden = false;
    const lines = [];
    if (outcome.outcome === "gate_failure") {
      lines.push(`${outcome.step} failed: ${outcome.detail}`);
    } else {
      lines.push(outcome.rendered ?? `${outcome.kind}`);
    }
    $("overlay-text").textContent = lines.join("\n");
    $("overlay-class").textContent = outcome.outcome.replace("_", " ");
    const position = /:(\d+):(\d+):/.exec(outcome.rendered ?? "");
    if (position) {
      this.editor.mark(Number(position[1]));
      $("overlay-goto").hidden = false;
      $("overlay-goto").onclick = () => this.editor.goTo(Number(position[1]), Number(position[2]));
    } else {
      this.editor.mark(null);
      $("overlay-goto").hidden = true;
    }
    this.status(`${outcome.outcome.replace("_", " ")}: ${outcome.kind ?? outcome.step}${this.lastGood ? " (showing the last good model)" : ""}`, true);
  }

  renderReport(envelope) {
    const table = $("report");
    table.innerHTML = "";
    const features = envelope.report?.features ?? [];
    for (const f of features) {
      const row = document.createElement("tr");
      const label = f.feature?.label ?? f.label ?? "?";
      const operation = f.feature?.operation ?? f.operation ?? "?";
      const outcome = typeof f.outcome === "string" ? f.outcome : Object.keys(f.outcome ?? {})[0] ?? "?";
      const elapsed = typeof f.elapsed_ms === "number" ? `${f.elapsed_ms.toFixed(1)} ms` : "";
      const downgraded = f.blend?.downgraded === true;
      const continuity = f.blend?.continuity ?? null;
      let verdict = outcome;
      if (continuity && continuity.requested) {
        // A blend row always states what was measured against what
        // was asked — a downgrade is never a bare success mark.
        verdict = downgraded
          ? `downgraded: ${continuity.achieved} measured, ${continuity.requested} requested`
          : `${outcome} · ${continuity.achieved} measured`;
      }
      row.innerHTML = `<td>${label}</td><td>${operation}</td><td class="${downgraded ? "cell-warn" : outcome === "Ok" || outcome === "ok" ? "cell-ok" : "cell-bad"}">${verdict}</td><td>${elapsed}</td>`;
      if (Array.isArray(f.warnings) && f.warnings.length) {
        const warn = document.createElement("tr");
        warn.className = "row-warn";
        warn.innerHTML = `<td colspan="4">${f.warnings.map((w) => (typeof w === "string" ? w : JSON.stringify(w))).join("; ")}</td>`;
        table.append(row, warn);
      } else {
        table.append(row);
      }
    }
    const summary = envelope.report
      ? `${envelope.report.failures} failure(s), ${envelope.report.downgraded_blends} downgraded blend(s)${envelope.report.is_clean ? ", clean" : ", with warnings"}`
      : "no report";
    $("report-summary").textContent = summary;
  }

  renderParams(params) {
    const list = $("params");
    list.innerHTML = "";
    for (const p of params) {
      if (p.read_only || p.expression) continue;
      const row = document.createElement("label");
      row.className = "param";
      const bounds = p.bounds ?? [p.value > 0 ? p.value / 4 : p.value * 2 - 1, p.value > 0 ? p.value * 3 : p.value / 2 + 1];
      const lo = Math.min(bounds[0], p.value);
      const hi = Math.max(bounds[1], p.value);
      row.innerHTML = `<span class="param-name">${p.name}</span>
        <input type="range" min="${lo}" max="${hi}" step="${(hi - lo) / 200}" value="${p.value}" aria-label="${p.name} slider">
        <input type="number" value="${p.value}" step="any" aria-label="${p.name} value">`;
      const [slider, number] = row.querySelectorAll("input");
      const apply = (value) => {
        if (!Number.isFinite(value)) return;
        number.value = value;
        slider.value = value;
        this.setParam(p.name, value);
      };
      slider.addEventListener("change", () => apply(Number(slider.value)));
      number.addEventListener("change", () => apply(Number(number.value)));
      list.append(row);
    }
    $("params-empty").hidden = list.children.length > 0;
  }

  async setParam(name, value) {
    if (this.running) return { ok: false, outcome: { outcome: "gate_failure", step: "busy", detail: "a run is in progress" } };
    this.running = true;
    this.status(`replaying with ${name} = ${value}…`);
    const reply = await this.host.setParam(name, value);
    this.running = false;
    this.receive(reply);
    return reply;
  }

  showPick(pick) {
    const el = $("pick");
    if (!pick) {
      el.textContent = "click a face to see its stable name";
      return;
    }
    el.textContent = `${pick.body.label ? `${pick.body.label} · ` : ""}face ${pick.face}: ${pick.label}`;
  }

  async download(kind) {
    const reply = kind === "step" ? await this.host.exportStep("model") : await this.host.exportPkn();
    if (!reply.ok) {
      this.showOutcome(reply.outcome);
      return;
    }
    const blob = new Blob([reply.bytes], { type: "application/octet-stream" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `model.${kind}`;
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2_000);
  }

  async showCode() {
    const reply = await this.host.toRhai();
    if (!reply.ok) {
      this.showOutcome(reply.outcome);
      return;
    }
    $("code-text").textContent = reply.code;
    $("code-dialog").hidden = false;
  }

  async loadExamples() {
    try {
      const response = await fetch("./examples.json");
      if (!response.ok) return;
      const manifest = await response.json();
      const select = $("examples");
      for (const entry of manifest.examples ?? []) {
        const option = document.createElement("option");
        option.value = entry.file;
        option.textContent = entry.title;
        select.append(option);
      }
      $("examples-row").hidden = select.options.length <= 1;
    } catch {
      // No gallery beside the page: the example picker stays hidden.
    }
  }

  async loadExample(file) {
    if (!file) return;
    const response = await fetch(`./examples/${file}`);
    if (!response.ok) return;
    this.editor.value = await response.text();
    this.viewer.fitted = false;
    this.run();
  }
}

window.addEventListener("error", (event) => {
  const status = $("status");
  if (status) status.textContent = `page error: ${event.message}`;
});
window.addEventListener("unhandledrejection", (event) => {
  const status = $("status");
  if (status) status.textContent = `page error: ${event.reason?.message ?? event.reason}`;
});
window.addEventListener("DOMContentLoaded", () => {
  window.playground = new Page();
});
