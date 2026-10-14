// The history scrubber (live-projects-and-history.md §6.3, §12): the
// history dock with one card per timeline feature (sketches and
// references minor, user-labelled features as named marks; drawn by
// history-strip.js), playback and step controls, and, for a project,
// one chapter per migration. The script's gutter marks each step's
// line (a click scrubs there). Moving the playhead shows
// the document after that step — bodies fetched once per key and
// cached, untouched bodies never re-sent — with the faces the step
// created glowing and fading, the step's script line marked in the
// editor and its report row highlighted. The camera is never refitted
// while scrubbing. A scrubbed-back view is read-only: any edit, slider
// or build returns the playhead to the end. `#step=N` (and `#v=K`, the
// K-th migration's last step) deep-links a step. An assembly scrubs its
// construction steps instead (PK-666): each place, move and mate, the
// occurrences placed by then at their placements then, the step's
// occurrences glowing — placements only, the definition meshes are the
// run's. A part opened in the assembly panel scrubs its own features
// (a step before the end shows the part's bodies alone). On a part, "insert a feature here" adds a
// migration that inserts one feature after the step on screen (inside
// `a.edit_part` for an assembly's part).

import { HistoryStrip } from "./history-strip.js";
import { showSheet } from "./layout.js";

const $ = (id) => document.getElementById(id);
// Play advances one step this often.
const PLAY_MS = 650;
// What the dock says while there is no timeline to step through.
const NO_HISTORY = "The model's history appears here after a build.";
const PREBUILT = "A pre-built model: build it live to step through its history.";
const ONE_STEP = "One step: nothing to step through.";

export class Scrubber {
  constructor(page) {
    this.page = page;
    this.steps = [];
    this.at = 0;
    this.cache = new Map();
    // Rows of steps already shown, by index (redraws need no worker).
    this.visited = new Map();
    this.playing = null;
    this.drawToken = 0;
    // The assembly timeline (`{ steps, recorded }`) when scrubbing one.
    this.assembly = null;
    // The assembly part whose features are scrubbed (opened in the
    // assembly panel), or null.
    this.part = null;
    this.strip = new HistoryStrip((index) => this.go(index));
    $("scrubber-play").addEventListener("click", () => (this.playing ? this.stop() : this.play()));
    $("scrubber-prev").addEventListener("click", () => this.go(this.at - 1));
    $("scrubber-next").addEventListener("click", () => this.go(this.at + 1));
    // The phone's step card jumps to the step's script line.
    $("step-card-line").addEventListener("click", (e) => {
      const line = Number(e.currentTarget.dataset.line);
      if (!line) return;
      showSheet("code");
      this.page.editor.goTo(line, 1);
    });
    // A feature row in the inspector scrubs to its step.
    $("report").addEventListener("click", (e) => {
      const row = e.target.closest("tr[data-feature]");
      const step = row && this.steps.find((x) => x.label === row.dataset.feature);
      if (step) this.go(step.index);
    });
    $("scrubber-range").addEventListener("input", (e) => this.go(Number(e.target.value)));
    $("scrubber").addEventListener("keydown", (e) => {
      if (e.key === "ArrowLeft") this.go(this.at - 1);
      else if (e.key === "ArrowRight") this.go(this.at + 1);
      else if (e.key === "Home") this.go(1);
      else if (e.key === "End") this.go(this.steps.length);
      else return;
      e.preventDefault();
    });
    $("scrubber-insert").addEventListener("click", () => {
      const step = this.steps[this.at - 1];
      // Kept so a caller (the harness) can await the new migration.
      if (step && !this.assembly) this.inserting = this.page.rail.insertAt(step, this.steps, this.part);
    });
  }

  /** Whether a step before the end is on screen. */
  get scrubbed() {
    return this.steps.length > 0 && this.at < this.steps.length;
  }

  /** Forget the timeline (a new model is coming). */
  reset() {
    this.stop();
    this.steps = [];
    this.at = 0;
    this.assembly = null;
    this.cache.clear();
    this.visited.clear();
    this.page.editor.setStepMarks?.(new Map(), null);
    this.strip.empty(this.page.lastGood && !this.page.live ? PREBUILT : NO_HISTORY, this.page.lastGood != null && !this.page.live);
    document.body.classList.remove("scrubbed");
    $("scrubber-insert").hidden = true;
  }

  /** A joint moved occurrences (a placement-only change): the steps
   *  and every part's meshes stand, but each visited assembly step's
   *  placements were solved at the old value. Nothing is requested. */
  placementsMoved() {
    if (this.assembly) this.visited.clear();
  }

  /** Load the live model's timeline — a part's features, or an
   *  assembly's construction steps — and put the playhead at the end,
   *  or at the deep-linked step. */
  load(envelope) {
    // Kept so a caller (the harness) can await the timeline.
    this.loading = this.loadTimeline(envelope);
    return this.loading;
  }

  async loadTimeline(envelope) {
    this.reset();
    if (!this.page.live) return;
    const reply = await this.page.host.request({ type: "timeline" });
    if (!reply.ok) return this.strip.empty(`No history: ${reply.outcome?.detail ?? "the timeline is unavailable"}`);
    const scoped = envelope.assembly && this.part ? reply.parts?.find((p) => p.part === this.part) : null;
    if (scoped) {
      this.steps = scoped.steps;
    } else if (envelope.assembly) {
      this.part = null;
      if (!reply.assembly) return;
      this.assembly = reply.assembly;
      this.steps = reply.assembly.steps.map((s) => ({ ...s, operation: s.kind }));
    } else {
      this.part = null;
      this.steps = reply.parts?.[0]?.steps ?? [];
    }
    if (this.steps.length < 2) {
      this.steps = [];
      return this.strip.empty(ONE_STEP);
    }
    this.render();
    const target = this.deepLink();
    this.at = this.steps.length;
    this.showLabel();
    if (target != null && target !== this.steps.length) await this.go(target);
  }

  /** Scrub assembly part `name`'s own features (null: back to the
   *  assembly's construction steps). The playhead starts at the end,
   *  which is the live assembly: opening a part never changes the
   *  view until a step is chosen. */
  async scopePart(name) {
    const live = this.page.lastGood;
    if (!live?.envelope.assembly || name === this.part) return;
    const wasScrubbed = this.scrubbed;
    this.part = name;
    await this.load(live.envelope);
    if (wasScrubbed) this.page.showModel(live.envelope, live.buffers);
  }

  /** The step a `#step=N` / `#v=K` link asks for. */
  deepLink() {
    const params = new URLSearchParams(location.hash.slice(1));
    if (params.has("step")) return Math.max(1, Math.min(this.steps.length, Number(params.get("step")) || this.steps.length));
    if (params.has("v")) {
      const chapters = this.chapters();
      const k = Number(params.get("v"));
      if (k >= 1 && k <= chapters.length) return chapters[k - 1].last;
    }
    return null;
  }

  /** Migration chapters: `{ name, first, last }` in order. */
  chapters() {
    const out = [];
    for (const s of this.steps) {
      const name = s.migration ?? null;
      if (!name) continue;
      const last = out.at(-1);
      if (last?.name === name) last.last = s.index;
      else out.push({ name, first: s.index, last: s.index });
    }
    return out;
  }

  render() {
    this.strip.render(this.steps, {
      flags: this.flags(),
      assembly: this.assembly != null,
      chapters: this.chapters(),
      goChapter: (c, k) => this.go(c.last, { link: `v=${k + 1}` }),
    });
    this.markLines();
  }

  /** Feature labels the build report flags: "failed" or "warn"
   *  (warnings, a downgraded blend). */
  flags() {
    const out = new Map();
    for (const f of this.page.lastGood?.envelope?.report?.features ?? []) {
      const label = f.feature?.label ?? f.label;
      if (!label) continue;
      const outcome = typeof f.outcome === "string" ? f.outcome : Object.keys(f.outcome ?? {})[0];
      if (outcome && outcome !== "Served") out.set(label, "failed");
      else if (f.blend?.downgraded || (Array.isArray(f.warnings) && f.warnings.length)) out.set(label, "warn");
    }
    return out;
  }

  /** The open script's line of `step`, or null when it is elsewhere. */
  lineOf(step) {
    const open = this.page.projects?.path?.split("/").pop() ?? null;
    return step?.source && (!open || step.source.script === open) ? step.source.line : null;
  }

  /** Gutter marks: each line of the open script that records a step
   *  (its last step, when a loop records several). */
  markLines() {
    const marks = new Map();
    for (const step of this.steps) {
      const line = this.lineOf(step);
      if (line != null) marks.set(line, { index: step.index, operation: step.operation, label: step.label });
    }
    this.page.editor.setStepMarks?.(marks, (index) => this.go(index));
  }

  showLabel() {
    const s = this.steps[this.at - 1];
    if (!s) return;
    const opened = this.assembly && !this.assembly.recorded ? " · opened file: placements only" : "";
    this.strip.paint(this.steps, this.at, { line: this.lineOf(s), note: opened });
    document.body.classList.toggle("scrubbed", this.scrubbed);
    const editable = this.page.projects?.current && !this.page.projects.current.readOnly;
    $("scrubber-insert").hidden = !(this.scrubbed && !this.assembly && editable);
  }

  /** Show the document after step `index`. Resolves to the redraw time
   *  in ms (measured from request to drawn). */
  async go(index, { link = null } = {}) {
    if (!this.steps.length) return null;
    const n = this.steps.length;
    index = Math.max(1, Math.min(n, index));
    const token = ++this.drawToken;
    const started = performance.now();
    this.at = index;
    this.showLabel();
    this.updateLink(link ?? (index === n ? null : `step=${index}`));
    if (this.assembly) return this.goAssembly(index, token, started);
    const rows = await this.stepRows(index);
    if (token !== this.drawToken || !rows) return null;
    const live = this.page.lastGood;
    if (index === n && live) {
      // The end is the live model itself (its analysis mesh, picks).
      this.page.showModel(live.envelope, live.buffers);
      if (this.page.viewer.mode === "zebra") this.page.ensureAnalysisMesh();
    } else {
      // An intermediate step has no analysis mesh: zebra reads the
      // display normals until the end is back on screen.
      this.page.viewer.setAnalysis(null);
      const bodies = rows.map((r) => this.cache.get(r.key).mesh);
      const buffers = rows.map((r) => this.cache.get(r.key).buffers);
      this.page.viewer.show(bodies, buffers, null);
    }
    const created = new Map();
    rows.forEach((r, i) => {
      if (r.created_faces.length) created.set(i, new Set(r.created_faces));
    });
    this.page.viewer.glowFaces(created);
    this.markSource(this.steps[index - 1]);
    if (index < n) this.warmNext(index + 1);
    return performance.now() - started;
  }

  /** An assembly after construction step `index`: the run's meshes at
   *  the step's occurrences, the step's own occurrences glowing. */
  async goAssembly(index, token, started) {
    const live = this.page.lastGood;
    let instances = this.visited.get(index);
    if (!instances) {
      const reply = await this.page.host.request({ type: "assembly_step", index });
      if (!reply.ok) {
        this.page.status(`step ${index}: ${reply.outcome?.detail ?? "unavailable"}`, true);
        return null;
      }
      instances = reply.instances;
      this.visited.set(index, instances);
    }
    if (token !== this.drawToken || !live?.envelope.assembly) return null;
    if (index === this.steps.length) this.page.showModel(live.envelope, live.buffers);
    else this.page.viewer.show(live.envelope.bodies, live.buffers, { ...live.envelope.assembly, instances });
    const step = this.steps[index - 1];
    this.page.viewer.glowFaces(new Map(), step.instances);
    this.markSource(step);
    return performance.now() - started;
  }

  /** The rows of step `index`, buffers cached by key; a visited step
   *  answers from the page without asking the worker. */
  async stepRows(index) {
    if (this.visited.has(index)) return this.visited.get(index);
    const known = [...this.cache.keys()];
    const reply = await this.page.host.request({ type: "step_meshes", part: this.part ?? "", index, known });
    if (!reply.ok) {
      this.page.status(`step ${index}: ${reply.outcome?.detail ?? "unavailable"}`, true);
      return null;
    }
    for (const r of reply.rows) {
      if (r.buffers) {
        this.cache.set(r.key, {
          mesh: r.mesh,
          buffers: { positions: r.positions, normals: r.normals, indices: r.indices, face_ids: r.face_ids },
        });
      }
    }
    const rows = reply.rows.filter((r) => this.cache.has(r.key));
    this.visited.set(index, rows);
    return rows;
  }

  /** Load a neighbour's new bodies while idle. */
  warmNext(index) {
    clearTimeout(this.warmTimer);
    this.warmTimer = setTimeout(() => {
      if (!this.page.running) this.stepRows(index);
    }, 120);
  }

  /** The step's line in the editor and row in the feature tree. */
  markSource(step) {
    this.page.editor.markStep?.(this.lineOf(step));
    for (const row of document.querySelectorAll("#report tr")) {
      row.classList.toggle("step-current", row.dataset.feature === step?.label);
    }
  }

  updateLink(fragment) {
    const url = new URL(location.href);
    url.hash = fragment ?? "";
    history.replaceState(null, "", url);
  }

  play() {
    if (!this.steps.length) return;
    if (this.at >= this.steps.length) this.go(1);
    const button = $("scrubber-play");
    button.classList.add("playing");
    button.setAttribute("aria-label", "Pause");
    this.playing = setInterval(async () => {
      if (this.at >= this.steps.length) return this.stop();
      await this.go(this.at + 1);
    }, PLAY_MS);
  }

  stop() {
    clearInterval(this.playing);
    this.playing = null;
    const button = $("scrubber-play");
    if (!button) return;
    button.classList.remove("playing");
    button.setAttribute("aria-label", "Play the history");
  }

  /** Back to the end: the live model (edits always build there). */
  async toEnd() {
    this.stop();
    if (!this.scrubbed) return;
    await this.go(this.steps.length);
  }
}
