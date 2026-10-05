// The history scrubber (live-projects-and-history.md §6.3, §12): a
// strip under the viewport with one tick per timeline feature (sketches
// and references minor, user-labelled features as named marks) and,
// for a project, one chapter per migration. Moving the playhead shows
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

const $ = (id) => document.getElementById(id);
// Play advances one step this often.
const PLAY_MS = 650;
// Operations drawn as minor ticks.
const MINOR = new Set(["sketch", "reference"]);
// Auto labels (`extrude-1`, `boolean-3`) are not chapter marks.
const AUTO_LABEL = /^[a-z_]+-\d+$/;

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
    this.strip = $("scrubber-ticks");
    $("scrubber-play").addEventListener("click", () => (this.playing ? this.stop() : this.play()));
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
    $("scrubber").hidden = true;
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
    if (!reply.ok) return;
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
    if (this.steps.length < 2) return;
    this.render();
    $("scrubber").hidden = false;
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
    const n = this.steps.length;
    const range = $("scrubber-range");
    range.min = "1";
    range.max = String(n);
    range.value = String(n);
    this.strip.replaceChildren();
    for (const s of this.steps) {
      const tick = document.createElement("span");
      tick.className = "tick";
      if (MINOR.has(s.operation)) tick.classList.add("minor");
      if (s.failed) tick.classList.add("failed");
      if (s.suppressed) tick.classList.add("suppressed");
      tick.style.left = `${((s.index - 1) / Math.max(1, n - 1)) * 100}%`;
      tick.title = `${s.index}. ${s.label} (${s.operation})`;
      if (this.assembly ? s.kind === "mate" : !AUTO_LABEL.test(s.label) && !MINOR.has(s.operation)) {
        tick.classList.add("named");
        tick.dataset.label = this.assembly ? "mate" : s.label;
      }
      this.strip.append(tick);
    }
    const chapters = this.chapters();
    const bar = $("scrubber-chapters");
    bar.replaceChildren();
    bar.hidden = chapters.length < 2;
    for (const [k, c] of chapters.entries()) {
      const chip = document.createElement("button");
      chip.type = "button";
      chip.className = "chapter";
      chip.textContent = c.name;
      chip.title = `${c.name}: steps ${c.first}–${c.last}`;
      chip.style.left = `${((c.first - 1) / Math.max(1, n - 1)) * 100}%`;
      chip.addEventListener("click", () => this.go(c.last, { link: `v=${k + 1}` }));
      bar.append(chip);
    }
  }

  showLabel() {
    const s = this.steps[this.at - 1];
    if (!s) return;
    const where = s.migration ? ` · ${s.migration}` : "";
    const what = this.assembly ? s.label : `${s.label} (${s.operation})`;
    const opened = this.assembly && !this.assembly.recorded ? " · opened file: placements only" : "";
    $("scrubber-label").textContent = `${this.at}/${this.steps.length} · ${what}${where}${opened}`;
    $("scrubber-range").value = String(this.at);
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
    const editor = this.page.editor;
    const open = this.page.projects?.path?.split("/").pop() ?? null;
    const line = step?.source && (!open || step.source.script === open) ? step.source.line : null;
    editor.markStep?.(line);
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
    $("scrubber-play").textContent = "pause";
    this.playing = setInterval(async () => {
      if (this.at >= this.steps.length) return this.stop();
      await this.go(this.at + 1);
    }, PLAY_MS);
  }

  stop() {
    clearInterval(this.playing);
    this.playing = null;
    const button = $("scrubber-play");
    if (button) button.textContent = "play";
  }

  /** Back to the end: the live model (edits always build there). */
  async toEnd() {
    this.stop();
    if (!this.scrubbed) return;
    await this.go(this.steps.length);
  }
}
