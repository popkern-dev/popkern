// popkern live: editor → worker → viewer, with the feature tree,
// parameters, exports and show-as-code. Everything the user
// types stays in this tab (design §7: client-side only; the worker is
// terminated, never trusted, when a script defeats its budgets).
//
// The page controller: it owns the build lifecycle (a run, its bake
// shortcut, cancel, deferred certification, parameter replays) and the
// state the panels share, and delegates rendering, persistence and the
// document flows to the modules imported below. `window.popkernLive`
// is this Page.
import { ensureAnalysisMesh } from "./analysis-mesh.js";
import { AssemblyPanel } from "./assembly-panel.js";
import { Bakes } from "./bakes.js";
import { checkClearance } from "./clearance-panel.js";
import { certifyContinuity } from "./continuity-certification.js";
import { wireControls } from "./controls.js";
import { DocumentTab } from "./document-tab.js";
import * as drafts from "./drafts.js";
import { Editor } from "./editor.js";
import { downloadModel, showModelCode } from "./exports.js";
import { fetchDefaultScript } from "./gallery.js";
import { adoptPlaygroundKeys } from "./local-store.js";
import { drawModel, showBake, showPlaced, showServed } from "./model-display.js";
import { showOutcome } from "./outcome-overlay.js";
import { renderParams, wireParamFilter } from "./params-panel.js";
import { hostMode } from "./host-mode.js";
import { Projects } from "./projects.js";
import { SocketWorker } from "./socket-worker.js";
import { Scrubber } from "./scrubber.js";
import { MigrationsRail } from "./migrations-rail.js";
import { wireProjects } from "./projects-dialog.js";
import { renderReport } from "./report-table.js";
import { reportPageErrors, setStatus, showPick } from "./status-line.js";
import { openStepAssembly } from "./step-import.js";
import { applyTimeLimit, budgetFor, describeLimit, renderTimeLimitOptions, storedLimit } from "./time-limit.js";
import { Viewer } from "./viewer.js";
import { WorkerHost } from "./worker-host.js";

// An edit rebuilds once typing pauses this long.
const IDLE_RUN_MS = 900;

const $ = (id) => document.getElementById(id);

class Page {
  constructor() {
    adoptPlaygroundKeys();
    this.timeLimit = storedLimit();
    // Where the kernel runs and projects live (host-mode.js): the
    // public page has no host.json; `popkern serve` names server
    // storage and (by default) native compute.
    this.modeReady = hostMode();
    this.host = new WorkerHost(
      () => (this.mode?.compute === "native" ? new SocketWorker() : new Worker("./worker.js", { type: "module" })),
      budgetFor(this.timeLimit),
    );
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
      open: (part) => {
        this.renderInspector();
        this.scrubber.scopePart(part);
      },
    });
    this.editor = new Editor($("editor"), {
      onChange: (text) => this.changed(text),
      onRun: () => this.run(),
    });
    this.tab = new DocumentTab(this.editor);
    // The open project (projects.js): autosaved in this browser.
    this.projects = new Projects(this);
    this.projectsReady = this.projects.init();
    this.idleTimer = null;
    this.cancelTimer = null;
    this.lastGood = null;
    this.running = false;
    // Bumped by anything that replaces the report (run, parameter
    // edit, stop): a certify loop from an older generation stops.
    this.certifyGeneration = 0;
    // The history scrubber under the viewport (live models only).
    this.scrubber = new Scrubber(this);
    // A project's migration files beside the editor.
    this.rail = new MigrationsRail(this);
    wireControls(this);
    wireProjects(this);
    wireParamFilter();
    renderTimeLimitOptions($("time-limit"), this.timeLimit);
    this.status("loading the kernel…");
    drafts.openFirstScript(this);
    this.modeReady
      .then((mode) => {
        this.mode = mode;
        document.body.dataset.host = `${mode.storage}-${mode.compute}`;
        if (mode.storage === "server") {
          // Self-hosted: say where things actually go.
          const note = $("privacy-note");
          note.textContent = mode.compute === "native" ? "your server" : "projects on your server";
          note.title = mode.compute === "native"
            ? "popkern serve: projects are folders on the machine serving this page, and the kernel runs there too. Nothing goes anywhere else."
            : "popkern serve: projects are folders on the machine serving this page; the kernel runs in this tab. Nothing goes anywhere else.";
        }
        return this.host.start();
      })
      .then((info) => {
        this.kernelReadyMs = performance.now();
        if (!this.lastGood) this.status(`kernel ready in ${this.kernelReadyMs.toFixed(0)} ms (protocol ${info.protocol_version}${info.clock_available ? "" : ", no clock: wall budget off"})`);
        if (this.pendingLive) this.run({ live: true });
      })
      .catch((error) => this.status(`kernel failed to load: ${error.message}`, true));
  }

  /** The script a fresh page (and reset) shows, fetched once. */
  defaultScript() {
    this.defaultText ??= fetchDefaultScript();
    return this.defaultText;
  }

  status(text, failed = false) {
    setStatus(text, failed);
  }

  changed(text) {
    // An edit always builds at the end of the history.
    this.scrubber.toEnd();
    this.tab.show();
    this.projects.edited(text);
    clearTimeout(this.idleTimer);
    this.idleTimer = setTimeout(() => this.run(), IDLE_RUN_MS);
  }

  /** Show the editor's script: its pre-built bake when one exists
   *  (unless `live`), else a kernel run. */
  async run({ live = false } = {}) {
    clearTimeout(this.idleTimer);
    const source = this.editor.value;
    // An imported STEP document has no script to build: keep it.
    if (this.tab.isImported(source)) {
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
    this.setBusy(true);
    let finished;
    this.runDone = new Promise((resolve) => (finished = resolve));
    const started = performance.now();
    const tick = () => this.status(`building… ${this.clock(started)}`);
    tick();
    this.clockTimer = setInterval(tick, 1000);
    if (this.timeLimit != null) {
      this.cancelTimer = setTimeout(() => this.cancel(`no answer within the ${describeLimit(this.timeLimit)} time limit: the worker was stopped`), this.timeLimit * 1000);
    }
    const request = this.projects.buildRequest(source);
    const reply = request?.kind === "project"
      ? await this.host.request({ type: "run_project", files: request.files })
      : await this.host.run(source, request?.name ?? "script.rhai");
    clearTimeout(this.cancelTimer);
    clearInterval(this.clockTimer);
    this.running = false;
    this.setBusy(false);
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

  /** A build is (not) in progress: the build controls lock meanwhile. */
  setBusy(busy) {
    $("run").disabled = busy;
    $("time-limit").disabled = busy;
    document.body.classList.toggle("busy", busy);
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

  async cancel(reason) {
    clearTimeout(this.cancelTimer);
    clearInterval(this.clockTimer);
    this.certifyGeneration += 1;
    this.status("stopping the build…", true);
    const result = await this.host.cancel(reason);
    this.running = false;
    this.setBusy(false);
    if (result?.in_place) {
      // A native session stopped the script and kept its model.
      this.status(`${reason} — the session kept its model`, true);
      return;
    }
    this.live = false;
    this.status(`${reason} — kernel restarted`, true);
  }

  async receive(reply) {
    if (!reply.ok) {
      this.showOutcome(reply.outcome);
      return;
    }
    const { envelope, buffers } = reply;
    // A joint moved placements only: the bodies and report on screen
    // still hold.
    if (envelope.bodies_unchanged && this.lastGood != null) {
      showPlaced(this, envelope);
      return;
    }
    if (envelope.project) this.rail.setVersions(envelope);
    // A run that raised carries no report: the panels keep describing
    // the last good model, like the viewer does.
    if (envelope.report) {
      this.reportEnvelope = envelope;
      if (envelope.outcome === "served") this.assembly.render(envelope.assembly ?? null);
      this.renderInspector(envelope);
    }
    if (envelope.outcome !== "served") {
      this.showOutcome(envelope);
      return;
    }
    await showServed(this, envelope, buffers);
  }

  /** The overlay for a build that did not serve (outcome-overlay.js). */
  showOutcome(outcome) {
    showOutcome(outcome, {
      editor: this.editor,
      timeLimit: this.timeLimit,
      keptModel: this.lastGood != null,
      // A project failure names its migration: show that file first.
      openScript: (file) => {
        const path = `migrations/${file}`;
        if (this.projects.current?.files?.[path] != null && this.projects.path !== path) this.projects.select(path);
      },
    });
  }

  /** The feature tree and parameters for the model on screen, scoped
   *  to the opened part of an assembly. */
  renderInspector(envelope = this.lastGood?.envelope) {
    if (!envelope) return;
    this.renderReport(envelope);
    this.renderParams(envelope.params);
  }

  renderReport(envelope) {
    renderReport(envelope, this.assembly?.openPart ?? null);
  }

  renderParams(params) {
    renderParams(params, this.assembly?.openPart ?? null, (name, value) => this.setParam(name, value));
  }

  async setParam(name, value) {
    if (!this.live && !(await this.ensureLive())) return { ok: false, outcome: { outcome: "gate_failure", step: "live", detail: "the live model could not be built" } };
    if (this.running) return { ok: false, outcome: { outcome: "gate_failure", step: "busy", detail: "a run is in progress" } };
    await this.scrubber.toEnd();
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

  // The page's other actions stay on the Page (the harness and the
  // console reach them here); each lives with the module that owns it.
  setTimeLimit(limit) { return applyTimeLimit(this, limit); }
  showBake(bake, source) { showBake(this, bake, source); }
  showModel(envelope, buffers) { drawModel(this, envelope, buffers); }
  ensureAnalysisMesh() { return ensureAnalysisMesh(this); }
  certify(envelope) { return certifyContinuity(this, envelope); }
  openStep(fileName, bytes) { return openStepAssembly(this, fileName, bytes); }
  checkClearance() { return checkClearance(this); }
  download(kind) { return downloadModel(this, kind); }
  showCode() { return showModelCode(this); }
  showPick(pick) { showPick(pick); }
  openDocument(title, text) { this.tab.open(title, text); }
  showDocument() { this.tab.show(); }
  openFirstScript() { return drafts.openFirstScript(this); }
  showDraftNotice() { drafts.showDraftNotice(); }
  restoreDraft() { drafts.restoreDraft(this); }
  loadExample(file) { return drafts.openExample(this, file); }
  copyToEdit() { return this.projects.copyToEdit(); }
}

reportPageErrors();
window.addEventListener("DOMContentLoaded", () => {
  window.popkernLive = new Page();
});
