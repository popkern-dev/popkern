// The page's side of PROTOCOL.md: owns the worker, matches replies to
// requests, and terminates + restarts the worker when a run must be
// abandoned (the backstop for a script that defeats the engine's
// budgets). Transport-agnostic: the page passes a `spawn` that returns
// a browser `Worker`; the node harness passes worker_threads.

// How long an in-place cancel may take before the session is replaced.
const CANCEL_GRACE_MS = 5000;

export class WorkerHost {
  /**
   * @param {() => {postMessage, terminate, onmessage|on}} spawn creates a worker
   * @param {{max_operations: number|null, timeout_ms: number|null}} options budgets
   */
  constructor(spawn, options) {
    this.spawn = spawn;
    this.options = options;
    this.worker = null;
    this.pending = new Map();
    this.next = 1;
    this.ready = null;
    this.generation = 0;
    this.restarts = 0;
    this.info = null;
    // Solved generative designs from bakes (PK-956): sent to every
    // worker this host starts, so a replaced worker keeps them.
    this.seeds = [];
    // Design histories from bakes (PK-959), likewise.
    this.historySeeds = [];
  }

  /** Start (or restart) the worker and wait for `init`. */
  start() {
    this.generation += 1;
    const worker = this.spawn();
    this.worker = worker;
    const onMessage = (reply) => this.#deliver(reply);
    const onError = (error) => this.#crashed(error);
    if (typeof worker.addEventListener === "function") {
      worker.addEventListener("message", (event) => onMessage(event.data));
      worker.addEventListener("error", (event) => onError(event.message ?? String(event)));
      worker.addEventListener("messageerror", () => onError("a reply could not be deserialised"));
    } else {
      worker.on("message", onMessage);
      worker.on("error", (error) => onError(String(error?.message ?? error)));
    }
    this.ready = this.request({ type: "init", options: this.options }).then((reply) => {
      if (!reply.ok) throw new Error(reply.outcome?.detail ?? "init failed");
      this.info = reply;
      for (const designs of this.seeds) this.request({ type: "seed_generative", designs });
      for (const bytes of this.historySeeds) this.request({ type: "seed_design_history", bytes });
      return reply;
    });
    return this.ready;
  }

  /** The worker itself failed (a load error, an uncaught exception):
   *  every pending request answers with one gate failure. */
  #crashed(detail) {
    const pending = [...this.pending.values()];
    this.pending.clear();
    for (const entry of pending) {
      entry.resolve({ ok: false, outcome: { outcome: "gate_failure", step: "worker", detail: `worker error: ${detail}` } });
    }
    if (this.onCrash) this.onCrash(detail);
  }

  #deliver(reply) {
    // A build's generative progress (PK-956) is not a reply.
    if (reply?.type === "progress") {
      if (this.onProgress) this.onProgress(reply.progress, reply.id);
      return;
    }
    const entry = this.pending.get(reply.id);
    if (!entry) return;
    this.pending.delete(reply.id);
    entry.resolve(reply);
  }

  /** Send one request; resolves with its reply. */
  request(message) {
    const id = this.next++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ ...message, id });
    });
  }

  /** Whether a request is in flight. */
  get busy() {
    return this.pending.size > 0;
  }

  /**
   * Abandon whatever the worker is doing: terminate it, fail every
   * pending request with a `gate_failure`, and start a fresh worker.
   * Resolves when the replacement has answered `init`.
   */
  async cancel(reason = "cancelled by the page") {
    const worker = this.worker;
    // A native session (socket-worker.js) stops the running script and
    // keeps its document; only a session that does not answer within
    // the grace is replaced.
    if (worker?.cancelRun && this.pending.size > 0) {
      worker.cancelRun();
      for (let waited = 0; waited < CANCEL_GRACE_MS && this.pending.size > 0; waited += 50) {
        await new Promise((resolve) => setTimeout(resolve, 50));
      }
      if (this.pending.size === 0) return { ok: true, in_place: true };
    }
    const pending = [...this.pending.values()];
    this.pending.clear();
    if (worker) {
      await Promise.resolve(worker.terminate());
    }
    for (const entry of pending) {
      entry.resolve({ ok: false, outcome: { outcome: "gate_failure", step: "cancelled", detail: reason } });
    }
    this.restarts += 1;
    return this.start();
  }

  /** Hand the worker a bake's solved generative designs (`designs`:
   *  the bake's `generative`, as JSON text) so builds serve them
   *  without solving; kept for any worker started later. */
  seedGenerative(designs) {
    if (this.seeds.includes(designs)) return Promise.resolve({ ok: true, added: 0 });
    this.seeds.push(designs);
    return this.info ? this.request({ type: "seed_generative", designs }) : Promise.resolve({ ok: true, added: 0 });
  }

  /** Hand the worker a bake's design histories (`bytes`: an
   *  ArrayBuffer of `<stem>.designs.bin`, PK-959); kept for any worker
   *  started later. */
  seedDesignHistory(bytes) {
    if (!this.historySeeds.includes(bytes)) this.historySeeds.push(bytes);
    return this.request({ type: "seed_design_history", bytes });
  }

  /** Design `design`'s recorded designs: `{ ok, history: { design,
   *  updates, cell, frames } }`, or `ok: false` when the worker holds
   *  no history for it. */
  designHistory(design) {
    return this.request({ type: "design_history", design });
  }

  /** Design `design` at position `at` among its records (or its design
   *  space with `space`): `{ ok, mesh: { positions, normals, indices } }`. */
  designMesh(design, at, space = false) {
    return this.request({ type: "design_mesh", design, at, space });
  }

  /** Replace the operation and wall-clock budgets for later runs. */
  setBudget(options) {
    return this.request({ type: "set_budget", options: { max_operations: options.max_operations, timeout_ms: options.timeout_ms } });
  }

  /** Build `source`; `options.generative === "solve"` lets this one
   *  build solve generative designs (the page's re-solve, PK-956). */
  run(source, name, options = {}) {
    return this.request({ type: "run", source, name, ...(options.generative ? { generative: options.generative } : {}) });
  }

  /** Open a STEP file's bytes (an ArrayBuffer) as a view-only assembly. */
  openStep(bytes, name) {
    return this.request({ type: "open_step", bytes, name });
  }

  setParam(name, value) {
    return this.request({ type: "set_param", name, value });
  }

  /** Set several joints at once (`{ j1: 30, … }`): one solve, one
   *  placement envelope (PK-853). */
  setJoints(values) {
    return this.request({ type: "set_joints", values });
  }

  /** Certify deferred blend seams for about `sliceMs` (see PROTOCOL.md). */
  certify(sliceMs) {
    return this.request({ type: "certify", slice_ms: sliceMs });
  }

  toRhai() {
    return this.request({ type: "to_rhai" });
  }

  exportStep(name) {
    return this.request({ type: "export_step", name });
  }

  exportPkn() {
    return this.request({ type: "export_pkn" });
  }

  /** The exact section of the built document by a clipping plane
   *  `{ normal, offset }` (world SI): `{ section: { svg, dxf, bodies,
   *  frame } }` (PK-856). */
  section(plane) {
    return this.request({ type: "section", normal: plane.normal, offset: plane.offset });
  }

  /** Interference and clearance between assembly occurrences. */
  interference(clearance) {
    return this.request({ type: "interference", clearance });
  }

  /** Show or hide an assembly occurrence in the worker's document (so
   *  a saved .pkn carries it); no geometry changes. */
  setVisible(path, visible) {
    return this.request({ type: "set_visible", path, visible });
  }

  /** Solve an analysis setup (setup file text) on the current part
   *  (PK-948): `{simulation, meshes}`. */
  simulate(setup) {
    return this.request({ type: "simulate", setup });
  }

  /** The material presets a simulation setup may name. */
  materials() {
    return this.request({ type: "materials" });
  }

  /** The zebra analysis mesh of the current model (PK-577). */
  analysisMesh() {
    return this.request({ type: "analysis_mesh" });
  }

  stats() {
    return this.request({ type: "stats" });
  }
}
