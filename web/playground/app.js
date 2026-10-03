// popkern live: editor → worker → viewer, with the feature tree,
// parameters, exports and show-as-code. Everything the user
// types stays in this tab (design §7: client-side only; the worker is
// terminated, never trusted, when a script defeats its budgets).
import { AssemblyPanel } from "./assembly-panel.js";
import { Bakes } from "./bakes.js";
import { Editor } from "./editor.js";
import { Viewer } from "./viewer.js";
import { WorkerHost } from "./worker-host.js";

// Fillets serve on sampled continuity evidence; the page certifies
// their seams afterwards in short slices (PROTOCOL.md, "Deferred
// continuity certification"), so a G2 part is on screen long before
// its continuity is proven.
const CERTIFY_SLICE_MS = 250;
const IDLE_RUN_MS = 900;
const DRAFT_KEY = "popkern-playground-draft";
// An edited script from an earlier visit, held aside while the page
// opens on the pre-built default (restore puts it back).
const SAVED_DRAFT_KEY = "popkern-playground-saved-draft";
// The document's title names its tab and its downloads (see fileStem):
// an example's gallery title, "Boss plate" for the default script, and
// "Document" for a script that is neither.
const DEFAULT_TITLE = "Boss plate";
const UNTITLED = "Document";

/** A title as a file stem: lowercase, spaces as underscores, and no
 *  characters a file name cannot carry ("Planetary gearbox" ->
 *  planetary_gearbox). */
function fileStem(title) {
  const stem = title
    .trim()
    .toLowerCase()
    .replace(/\s+/g, "_")
    .replace(/[\/\\:*?"<>|\x00-\x1f]/g, "");
  return stem || "document";
}
const TIME_LIMIT_KEY = "popkern-playground-time-limit";

// The build time limit the user picks, in seconds (null = none). The
// kernel's wall budget stops the script at the limit less a grace, and
// the page terminates the worker at the limit itself (the backstop for
// a long kernel call the script cannot interrupt). The operation
// budget, a guard against runaway script loops, scales with the time.
const TIME_LIMITS = [10, 24, 60, 120, 300, 900, null];
const DEFAULT_TIME_LIMIT = 24;
const OPERATIONS_PER_SECOND = 1_000_000;

function budgetFor(limit) {
  if (limit == null) return { max_operations: null, timeout_ms: null, continuity_certification: "Deferred" };
  const kernel = limit - Math.min(4, limit / 6);
  return { max_operations: Math.round(kernel * OPERATIONS_PER_SECOND), timeout_ms: kernel * 1000, continuity_certification: "Deferred" };
}

function describeLimit(limit) {
  if (limit == null) return "no limit";
  return limit < 60 ? `${limit} s` : `${limit / 60} min`;
}

function storedLimit() {
  const stored = storage((s) => s.getItem(TIME_LIMIT_KEY));
  if (stored === "none") return null;
  const value = Number(stored);
  return TIME_LIMITS.includes(value) ? value : DEFAULT_TIME_LIMIT;
}

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
    this.timeLimit = storedLimit();
    this.host = new WorkerHost(() => new Worker("./worker.js", { type: "module" }), budgetFor(this.timeLimit));
    this.bakes = new Bakes();
    // Whether the worker's model is the one on screen. A pre-built
    // example is shown without the kernel; anything that needs the
    // live model (a slider, an export, show-as-code, zebra) builds it
    // first.
    this.live = false;
    // Bumped by every run request: a reply or bake lookup from an
    // older request never replaces a newer view.
    this.viewToken = 0;
    this.runDone = Promise.resolve();
    this.host.onCrash = (detail) => this.status(`worker error: ${detail}`, true);
    this.viewer = new Viewer($("viewer"), (pick) => this.showPick(pick));
    // The assembly tree: selection and visibility are view state here
    // (a rebuild keeps them); "open" scopes the inspector to one part.
    this.assembly = new AssemblyPanel({
      select: (path) => this.viewer.select(path),
      visible: (path, visible) => this.viewer.setVisible(path, visible),
      open: () => this.renderInspector(),
    });
    this.editor = new Editor($("editor"), {
      onChange: (text) => this.changed(text),
      onRun: () => this.run(),
    });
    this.idleTimer = null;
    this.cancelTimer = null;
    this.lastGood = null;
    this.running = false;
    // Bumped by anything that replaces the report (run, parameter
    // edit, stop): a certify loop from an older generation stops.
    this.certifyGeneration = 0;
    this.wire();
    this.renderTimeLimit();
    this.examplesReady = this.loadExamples();
    this.status("loading the kernel…");
    this.openFirstScript();
    this.host
      .start()
      .then((info) => {
        this.kernelReadyMs = performance.now();
        if (!this.lastGood) this.status(`kernel ready in ${this.kernelReadyMs.toFixed(0)} ms (protocol ${info.protocol_version}${info.clock_available ? "" : ", no clock: wall budget off"})`);
        if (this.pendingLive) this.run({ live: true });
      })
      .catch((error) => this.status(`kernel failed to load: ${error.message}`, true));
  }

  wire() {
    $("run").addEventListener("click", () => this.run({ live: true }));
    $("time-limit").addEventListener("change", (e) => this.setTimeLimit(e.target.value === "none" ? null : Number(e.target.value)));
    $("cancel").addEventListener("click", () => this.cancel("cancelled by the user"));
    $("fit").addEventListener("click", () => this.viewer.fitAgain());
    $("view-mode").addEventListener("change", async (e) => {
      const mode = e.target.value;
      if (mode === "zebra") await this.ensureAnalysisMesh();
      this.viewer.setMode(mode);
    });
    $("edges").addEventListener("change", (e) => this.viewer.toggleEdges(e.target.checked));
    $("explode").addEventListener("input", (e) => this.viewer.setExplode(Number(e.target.value)));
    $("clearance-check").addEventListener("click", () => this.checkClearance());
    $("open-step").addEventListener("click", () => $("open-step-file").click());
    $("open-step-file").addEventListener("change", async (e) => {
      const file = e.target.files?.[0];
      if (file) await this.openStep(file.name, await file.arrayBuffer());
      e.target.value = "";
    });
    $("download-step").addEventListener("click", () => this.download("step"));
    $("download-pkn").addEventListener("click", () => this.download("pkn"));
    $("show-code").addEventListener("click", () => this.showCode());
    $("code-close").addEventListener("click", () => ($("code-dialog").hidden = true));
    $("code-use").addEventListener("click", () => {
      this.editor.value = $("code-text").textContent;
      this.showDocument();
      $("code-dialog").hidden = true;
      this.run();
    });
    $("examples").addEventListener("change", (e) => this.loadExample(e.target.value));
    $("draft-restore").addEventListener("click", () => this.restoreDraft());
    $("about").addEventListener("click", () => ($("about-dialog").hidden = false));
    $("about-close").addEventListener("click", () => ($("about-dialog").hidden = true));
    for (const dialog of ["about-dialog", "code-dialog"]) {
      // A click on the backdrop (not the box) or Escape closes it.
      $(dialog).addEventListener("click", (e) => {
        if (e.target === e.currentTarget) e.currentTarget.hidden = true;
      });
    }
    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") for (const dialog of ["about-dialog", "code-dialog"]) $(dialog).hidden = true;
    });
    $("draft-dismiss").addEventListener("click", () => {
      storage((s) => s.removeItem(SAVED_DRAFT_KEY));
      this.showDraftNotice();
    });
    $("reset").addEventListener("click", async () => {
      this.openDocument(DEFAULT_TITLE, await this.defaultScript());
      storage((s) => s.setItem(DRAFT_KEY, this.editor.value));
      this.run();
    });
  }

  /** The page opens on something pre-built, so a visitor never waits
   *  on the kernel: the last script when it has a bake (the default or
   *  an example), else the default — and an edited script from last
   *  time is held aside with a notice offering to restore it. */
  async openFirstScript() {
    const draft = storage((s) => s.getItem(DRAFT_KEY));
    const fallback = await this.defaultScript();
    const prebuilt = draft ? await this.bakes.find(draft) : null;
    if (draft && !prebuilt && draft !== fallback) {
      storage((s) => s.setItem(SAVED_DRAFT_KEY, draft));
    }
    const name = prebuilt ? await this.bakes.nameOf(draft) : null;
    const titles = await this.examplesReady;
    this.openDocument(titles?.get(name) ?? DEFAULT_TITLE, prebuilt ? draft : fallback);
    this.showDraftNotice();
    this.run();
  }

  /** Put `text` in the editor as the document titled `title`. */
  openDocument(title, text) {
    this.docTitle = title;
    this.docSource = text;
    this.editor.value = text;
    this.showDocument();
  }

  /** The tab: the document's name, marked once it is edited. */
  showDocument() {
    const edited = this.docSource != null && this.editor.value !== this.docSource;
    $("doc-name").textContent = `${fileStem(this.docTitle ?? UNTITLED)}.rhai${edited ? " ●" : ""}`;
    $("doc-name").title = edited ? "edited: builds live" : "";
  }

  showDraftNotice() {
    $("draft-notice").hidden = !storage((s) => s.getItem(SAVED_DRAFT_KEY));
  }

  restoreDraft() {
    const saved = storage((s) => s.getItem(SAVED_DRAFT_KEY));
    storage((s) => s.removeItem(SAVED_DRAFT_KEY));
    this.showDraftNotice();
    if (!saved) return;
    this.openDocument(UNTITLED, saved);
    storage((s) => s.setItem(DRAFT_KEY, saved));
    this.viewer.fitted = false;
    this.run();
  }

  /** The script a fresh page (and reset) shows, fetched once. It
   *  lives in default.rhai so tools/bake-playground.sh pre-builds it
   *  with the gallery: a refresh that has not edited it never waits on
   *  the kernel. */
  defaultScript() {
    this.defaultText ??= fetch("./default.rhai")
      .then((r) => (r.ok ? r.text() : Promise.reject(new Error(`HTTP ${r.status}`))))
      .catch((error) => `// The default script could not be loaded (${error.message}).\n// Pick an example above, or write a script here.\nlet m = model();\nm\n`);
    return this.defaultText;
  }

  status(text, failed = false) {
    const el = $("status");
    el.textContent = text;
    el.classList.toggle("status-failed", failed);
  }

  changed(text) {
    storage((s) => s.setItem(DRAFT_KEY, text));
    this.showDocument();
    clearTimeout(this.idleTimer);
    this.idleTimer = setTimeout(() => this.run(), IDLE_RUN_MS);
  }

  /** Show the editor's script: its pre-built bake when one exists
   *  (unless `live`), else a kernel run. */
  async run({ live = false } = {}) {
    clearTimeout(this.idleTimer);
    const source = this.editor.value;
    // An imported STEP document has no script to build: keep it.
    if (this.importedSource != null && source === this.importedSource) {
      this.status("an imported STEP assembly is view only: pick an example or write a script to build");
      return null;
    }
    const token = ++this.viewToken;
    if (!live) {
      const bake = await this.bakes.find(source);
      if (token !== this.viewToken) return null;
      if (bake) {
        if (this.lastGood?.bakedSource !== source || this.live) this.showBake(bake, source);
        return { ok: true, baked: true };
      }
    }
    if (!this.host.info) {
      this.pendingLive = true;
      return null;
    }
    this.pendingLive = false;
    if (this.running) {
      this.rerun = true;
      return null;
    }
    this.running = true;
    this.rerun = false;
    this.certifyGeneration += 1;
    $("run").disabled = true;
    $("time-limit").disabled = true;
    document.body.classList.add("busy");
    let finished;
    this.runDone = new Promise((resolve) => (finished = resolve));
    const started = performance.now();
    const tick = () => this.status(`building… ${this.clock(started)}`);
    tick();
    this.clockTimer = setInterval(tick, 1000);
    if (this.timeLimit != null) {
      this.cancelTimer = setTimeout(() => this.cancel(`no answer within the ${describeLimit(this.timeLimit)} time limit: the worker was stopped`), this.timeLimit * 1000);
    }
    const reply = await this.host.run(source, "playground.rhai");
    clearTimeout(this.cancelTimer);
    clearInterval(this.clockTimer);
    this.running = false;
    $("run").disabled = false;
    $("time-limit").disabled = false;
    document.body.classList.remove("busy");
    if (token === this.viewToken) {
      // A described model (a report) is now the worker's; it is the one
      // on screen only when it served. A run that raised kept the
      // previous model, so `live` stands.
      if (reply.ok && reply.envelope?.report) this.live = reply.envelope.outcome === "served";
      await this.receive(reply);
    } else {
      // Superseded (by an edit queued as a rerun, or a pre-built
      // example now on screen): its model is in the worker but not on
      // screen, and its reply is dropped.
      this.live = false;
    }
    finished();
    if (this.rerun) this.run();
    return reply;
  }

  /** Elapsed build time against the limit, for the status line. */
  clock(started) {
    const seconds = Math.floor((performance.now() - started) / 1000);
    return this.timeLimit == null ? `${seconds} s (no time limit; stop cancels)` : `${seconds} s of ${describeLimit(this.timeLimit)}`;
  }

  /** Make the worker's model the one on screen (a pre-built example
   *  is rebuilt live). Resolves to whether it now is. */
  async ensureLive() {
    await this.runDone;
    if (this.live) return true;
    if (!this.host.info) {
      this.status("the kernel is still loading", true);
      return false;
    }
    await this.run({ live: true });
    return this.live;
  }

  /** Show a pre-built example: the certified report, parameters and
   *  display meshes exactly as the native build left them. */
  showBake(bake, source) {
    this.live = false;
    this.certifyGeneration += 1;
    const { envelope, buffers } = bake;
    this.lastGood = { envelope, buffers, bakedSource: source };
    this.assembly.render(envelope.assembly ?? null);
    this.renderInspector();
    this.viewer.setAnalysis(null);
    this.showModel(envelope, buffers);
    this.editor.mark(null);
    $("overlay").hidden = true;
    const open = envelope.bodies.reduce((n, b) => n + b.open_edges, 0);
    const native = ((bake.native?.run_ms ?? 0) + (bake.native?.certify_ms ?? 0)) / 1000;
    const first = this.firstRenderMs == null ? ` — first render ${(this.firstRenderMs = performance.now()).toFixed(0)} ms after page start` : "";
    this.servedText = `served (pre-built): ${envelope.features} feature(s), ${contents(envelope)}${open ? `, ${open} open mesh edges` : ", watertight"}${first}`;
    this.status(`${this.servedText} — continuity certified; builds natively in ${native.toFixed(1)} s, and edits, sliders and exports rebuild it live`);
    if (this.viewer.mode === "zebra") this.ensureAnalysisMesh();
  }

  renderTimeLimit() {
    const select = $("time-limit");
    select.innerHTML = "";
    for (const limit of TIME_LIMITS) {
      const option = document.createElement("option");
      option.value = limit == null ? "none" : String(limit);
      option.textContent = describeLimit(limit) + (limit === DEFAULT_TIME_LIMIT ? " (default)" : "");
      option.selected = limit === this.timeLimit;
      select.append(option);
    }
  }

  /** Change the build time limit for later runs (remembered here). */
  async setTimeLimit(limit) {
    this.timeLimit = limit;
    storage((s) => s.setItem(TIME_LIMIT_KEY, limit == null ? "none" : String(limit)));
    const budget = budgetFor(limit);
    // A restarted worker (stop, crash) starts under the new limit too.
    this.host.options = budget;
    if (this.host.info) await this.host.setBudget(budget);
  }

  async cancel(reason) {
    clearTimeout(this.cancelTimer);
    clearInterval(this.clockTimer);
    this.live = false;
    this.certifyGeneration += 1;
    this.status("stopping the worker…", true);
    await this.host.cancel(reason);
    this.running = false;
    $("run").disabled = false;
    $("time-limit").disabled = false;
    document.body.classList.remove("busy");
    this.status(`${reason} — kernel restarted`, true);
  }

  async receive(reply) {
    if (!reply.ok) {
      this.showOutcome(reply.outcome);
      return;
    }
    const { envelope, buffers } = reply;
    // A run that raised carries no report: the panels keep describing
    // the last good model, like the viewer does.
    if (envelope.report) {
      this.reportEnvelope = envelope;
      if (envelope.outcome === "served") this.assembly.render(envelope.assembly ?? null);
      this.renderInspector(envelope);
    }
    if (envelope.outcome === "served") {
      this.lastGood = { envelope, buffers };
      this.viewer.setAnalysis(null);
      if (this.viewer.mode === "zebra") await this.ensureAnalysisMesh();
      this.showModel(envelope, buffers);
      this.editor.mark(null);
      $("overlay").hidden = true;
      const ms = envelope.elapsed_ms == null ? "" : ` in ${envelope.elapsed_ms.toFixed(0)} ms`;
      const open = envelope.bodies.reduce((n, b) => n + b.open_edges, 0);
      const first = this.firstRenderMs == null ? ` — first render ${(this.firstRenderMs = performance.now()).toFixed(0)} ms after page start` : "";
      this.servedText = `served: ${envelope.features} feature(s), ${contents(envelope)}${ms}${open ? `, ${open} open mesh edges` : ", watertight"}${first}`;
      this.status(this.servedText);
      this.certify(envelope);
    } else {
      this.showOutcome(envelope);
    }
  }

  showOutcome(outcome) {
    const overlay = $("overlay");
    overlay.hidden = false;
    const lines = [];
    if (outcome.outcome === "gate_failure") {
      lines.push(outcome.step === "cancelled" ? outcome.detail : `${outcome.step} failed: ${outcome.detail}`);
    } else {
      lines.push(outcome.rendered ?? `${outcome.kind}`);
    }
    if (outcome.kind === "BudgetExceeded" || outcome.step === "cancelled") {
      lines.push("", this.timeLimit == null ? "Press run to try again." : `The time limit is ${describeLimit(this.timeLimit)}: raise it in the top bar and run again.`);
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

  /** The zebra view renders from a finer ANALYSIS mesh (6°, PK-577):
   *  fetched once per served model, kept until the next run. On a
   *  failure the zebra view falls back to the display mesh. */
  async ensureAnalysisMesh() {
    if (!this.lastGood || this.lastGood.analysis) return;
    if (!(await this.ensureLive())) return;
    const good = this.lastGood;
    if (good.analysis) return;
    const reply = await this.host.analysisMesh();
    if (this.lastGood !== good) return;
    if (!reply?.ok) {
      const why = reply?.outcome?.detail ?? "unavailable";
      this.status(`zebra: ${why}`, reply?.outcome?.step !== "analysis_mesh");
      return;
    }
    good.analysis = reply.buffers;
    this.viewer.setAnalysis(reply.buffers);
    const { triangles, mesh_ms, facet_deg } = reply.stats;
    const ms = mesh_ms == null ? "" : ` in ${mesh_ms.toFixed(0)} ms`;
    this.status(`zebra: ${triangles} triangles at ${facet_deg}°${ms}`);
  }

  /** Certify the served model's deferred seams slice by slice, then
   *  swap in the upgraded report. Abandoned when a newer run,
   *  parameter edit or stop bumps the generation. */
  async certify(envelope) {
    const total = deferredSeams(envelope.report);
    if (!total) return;
    const generation = ++this.certifyGeneration;
    const started = performance.now();
    let remaining = total;
    while (remaining > 0) {
      this.status(`${this.servedText} — certifying continuity: ${total - remaining}/${total} seams`);
      const reply = await this.host.certify(CERTIFY_SLICE_MS);
      if (generation !== this.certifyGeneration) return;
      const step = reply.ok ? reply.step : { failure: reply.outcome };
      if (step.failure) {
        this.status(`${this.servedText} — continuity certification stopped: ${step.failure.detail ?? step.failure.kind ?? "failed"}`, true);
        return;
      }
      remaining = step.remaining;
      if (step.report) {
        envelope.report = step.report;
        this.renderInspector(envelope);
      }
    }
    this.status(`${this.servedText} — continuity certified (${total} seams, ${((performance.now() - started) / 1000).toFixed(1)} s)`);
  }

  /** Draw a served model, re-applying the tree's visibility and
   *  selection (a rebuild keeps what the user hid). */
  showModel(envelope, buffers) {
    // The explode slider belongs to assemblies; a new document starts
    // assembled.
    const exploding = envelope.assembly != null;
    if (!exploding || $("explode-row").hidden) {
      $("explode").value = "0";
      this.viewer.explode = 0;
    }
    $("explode-row").hidden = !exploding;
    this.viewer.show(envelope.bodies, buffers, envelope.assembly ?? null);
    for (const path of this.assembly.hidden) this.viewer.setVisible(path, false);
    if (this.assembly.selected != null) this.viewer.select(this.assembly.selected);
  }

  /** Open a STEP file as a view-only assembly (PK-631): the product
   *  tree, names and placements, with no script behind it. The file is
   *  read here and handed to the worker; it never leaves the tab. */
  async openStep(fileName, bytes) {
    if (!this.host.info) {
      this.status("the kernel is still loading", true);
      return null;
    }
    await this.runDone;
    const stem = fileName.replace(/\.(step|stp)$/i, "");
    this.status(`opening ${fileName}…`);
    const token = ++this.viewToken;
    this.certifyGeneration += 1;
    const reply = await this.host.openStep(bytes, stem);
    if (token !== this.viewToken) return reply;
    if (reply.ok && reply.envelope?.outcome === "served") {
      this.importedSource = `// ${fileName}: an imported STEP assembly — view only (no script).\n// Pick an example or write a script to build.\n`;
      this.openDocument(stem, this.importedSource);
      this.viewer.fitted = false;
      this.live = true;
    }
    await this.receive(reply);
    if (reply.ok && reply.envelope?.printed?.length) this.status(`${reply.envelope.printed[0]} — ${this.servedText ?? ""}`);
    return reply;
  }

  /** Interference and clearance between the assembly's parts (PK-630):
   *  certified clearances, contact and overlap volumes, worst first. A
   *  row selects its two occurrences. */
  async checkClearance() {
    if (!(await this.ensureLive())) return null;
    const clearance = Math.max(0, Number($("clearance").value) || 0);
    this.status(`checking clearance within ${clearance}…`);
    const reply = await this.host.interference(clearance);
    const list = $("clearance-list");
    list.innerHTML = "";
    list.hidden = false;
    $("clearance-summary").hidden = false;
    if (!reply.ok) {
      $("clearance-summary").textContent = reply.outcome?.detail ?? reply.outcome?.rendered ?? "the check failed";
      return reply;
    }
    const { pairs, occurrences, candidate_pairs } = reply.report;
    const clashes = pairs.filter((p) => p.status === "interfering").length;
    $("clearance-summary").textContent = `${occurrences} occurrence(s), ${candidate_pairs} pair(s) near: ${clashes} interfering, ${pairs.length - clashes} within ${clearance}`;
    for (const pair of pairs) {
      const row = document.createElement("li");
      row.className = `tree-row ${pair.status === "interfering" ? "clash" : pair.status === "touching" ? "contact" : ""}`;
      const name = document.createElement("span");
      name.className = "tree-name";
      name.textContent = `${pair.a} ↔ ${pair.b}`;
      const what = document.createElement("span");
      what.className = "tree-kind";
      what.textContent = clearanceText(pair);
      row.title = (pair.faces ?? []).join("  ↔  ") || what.textContent;
      row.append(name, what);
      row.addEventListener("click", () => this.viewer.select([pair.a, pair.b]));
      list.append(row);
    }
    this.status(`clearance checked: ${$("clearance-summary").textContent}`);
    return reply;
  }

  /** The feature tree and parameters for the model on screen, scoped
   *  to the opened part of an assembly. */
  renderInspector(envelope = this.lastGood?.envelope) {
    if (!envelope) return;
    this.renderReport(envelope);
    this.renderParams(envelope.params);
  }

  renderReport(envelope) {
    const table = $("report");
    table.innerHTML = "";
    const part = this.assembly?.openPart ?? null;
    const scoped = envelope.assembly != null;
    const features = (envelope.report?.features ?? []).filter((f) => part == null || f.part === part);
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
        // was asked — a downgrade is never a bare success mark — and
        // how: certified over every seam, or from samples (pending
        // certification, or where a certificate could not decide).
        const how = continuity.deferred_seams > 0 ? "sampled, certifying…" : continuity.certified_achieved ? "certified" : "measured";
        verdict = downgraded
          ? `downgraded: ${continuity.achieved} ${how}, ${continuity.requested} requested`
          : `${outcome} · ${continuity.achieved} ${how}`;
      }
      const shown = scoped && part == null && f.part ? `${f.part}: ${label}` : label;
      row.innerHTML = `<td>${shown}</td><td>${operation}</td><td class="${downgraded ? "cell-warn" : outcome === "Served" ? "cell-ok" : "cell-bad"}">${verdict}</td><td>${elapsed}</td>`;
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
    const part = this.assembly?.openPart ?? null;
    for (const p of params) {
      if (p.read_only || p.expression) continue;
      // An assembly's parameters are qualified "part:name"; an opened
      // part shows only its own.
      if (part != null && !p.name.startsWith(`${part}:`)) continue;
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
    if (!this.live && !(await this.ensureLive())) return { ok: false, outcome: { outcome: "gate_failure", step: "live", detail: "the live model could not be built" } };
    if (this.running) return { ok: false, outcome: { outcome: "gate_failure", step: "busy", detail: "a run is in progress" } };
    this.running = true;
    this.certifyGeneration += 1;
    this.status(`replaying with ${name} = ${value}…`);
    document.body.classList.add("busy");
    const reply = await this.host.setParam(name, value);
    document.body.classList.remove("busy");
    this.running = false;
    await this.receive(reply);
    return reply;
  }

  showPick(pick) {
    const el = $("pick");
    el.textContent = pickText(pick);
    // One line (the full name is the tooltip), so a long pick never
    // reflows the viewer under the pointer.
    el.title = el.textContent;
  }

  async download(kind) {
    if (!(await this.ensureLive())) return;
    // Named after the document: "Planetary gearbox" downloads as
    // planetary_gearbox.step (the STEP product carries the name too).
    const name = fileStem(this.docTitle ?? UNTITLED);
    // Visibility set in the tree is document data: the saved .pkn
    // carries it (the live model was built with the script's own).
    if (kind === "pkn" && this.lastGood?.envelope.assembly) {
      for (const instance of this.lastGood.envelope.assembly.instances) {
        const hidden = this.assembly.hidden.has(instance.path);
        const scripted = this.lastGood.envelope.assembly.instances.find((i) => i.path === instance.path);
        if (hidden || !scripted.visible) await this.host.setVisible(instance.path, !hidden);
      }
    }
    const reply = kind === "step" ? await this.host.exportStep(name) : await this.host.exportPkn();
    if (!reply.ok) {
      this.showOutcome(reply.outcome);
      return;
    }
    const blob = new Blob([reply.bytes], { type: "application/octet-stream" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${name}.${kind}`;
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2_000);
  }

  async showCode() {
    if (!(await this.ensureLive())) return;
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
      const titles = new Map();
      for (const entry of manifest.examples ?? []) {
        titles.set(entry.file, entry.title);
        const option = document.createElement("option");
        option.value = entry.file;
        option.textContent = entry.title;
        select.append(option);
      }
      $("examples-row").hidden = select.options.length <= 1;
      return titles;
    } catch {
      // No gallery beside the page: the example picker stays hidden.
      return null;
    }
  }

  async loadExample(file) {
    if (!file) return;
    const response = await fetch(`./examples/${file}`);
    if (!response.ok) return;
    const titles = await this.examplesReady;
    this.openDocument(titles?.get(file) ?? file.replace(/\.rhai$/, ""), await response.text());
    // A refresh reopens the chosen example (pre-built, so instant).
    storage((s) => s.setItem(DRAFT_KEY, this.editor.value));
    this.viewer.fitted = false;
    this.run();
  }
}

/** Deferred seams across a report's blend rows. */
/** The viewer bar's text for a face or edge pick (null: the hint). */
function pickText(pick) {
  if (!pick) return "click a face or edge to see its stable name";
  if (pick.instance) {
    // An assembly pick is an entity path: the occurrence, then the
    // part's own stable name (design §4): "planet_1:planet/…".
    const at = (label) => `${pick.instance}:${label}`;
    if (pick.edge) {
      const [la, lb] = pick.edge.labels;
      return pick.edge.faces[1] === null ? `open edge of ${at(la)}` : `edge of ${at(la)} | ${at(lb)}`;
    }
    return `face ${at(pick.label)}`;
  }
  const body = pick.body.label ? `${pick.body.label} · ` : "";
  if (pick.edge) {
    const [a, b] = pick.edge.faces;
    const [la, lb] = pick.edge.labels;
    return b === null ? `${body}open edge of face ${a}: ${la}` : `${body}edge of faces ${a} | ${b}: ${la} | ${lb}`;
  }
  return `${body}face ${pick.face}: ${pick.label}`;
}

/** A clearance row's verdict. */
function clearanceText(pair) {
  const n = (x) => (Math.abs(x) >= 1e-3 || x === 0 ? x.toPrecision(4) : x.toExponential(3));
  if (pair.status === "interfering") return pair.volume != null ? `interfering: ${n(pair.volume)} overlap` : `interfering (volume unknown: ${pair.decline ?? "declined"})`;
  if (pair.status === "touching") return "touching";
  return `clear: ${n(pair.distance)}${pair.refined ? "" : " (mesh)"}, ≥ ${n(pair.lower_bound)}`;
}

/** "3 body(ies)" for a part document, "4 part(s), 9 occurrence(s)" for an assembly. */
function contents(envelope) {
  const view = envelope.assembly;
  if (!view) return `${envelope.bodies.length} body(ies)`;
  const parts = view.definitions.filter((d) => d.kind === "part").length;
  return `${parts} part(s), ${view.instances.filter((i) => i.is_part).length} occurrence(s)`;
}

function deferredSeams(report) {
  return (report?.features ?? []).reduce((n, f) => n + (f.blend?.continuity?.deferred_seams ?? 0), 0);
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
