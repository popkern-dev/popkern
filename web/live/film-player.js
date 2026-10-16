// The film player (docs/architecture/presentation-films.md §6): runs a
// film script in the worker against the document on screen, then
// applies the film's state at a time to the viewer — occurrence
// matrices from their placements, explode, fades, visibility,
// highlights, the section, joints (placement-only, via the worker) and
// the camera — and plays it in real time. The export (film-export.js)
// drives the same `apply`.
import * as THREE from "./vendor/three.module.min.js";
import { cameraAt, sectionAt, stateAt } from "./film-state.js";
import { FilmLook } from "./film-look.js";

const HIGHLIGHT = new THREE.Color(0x3b82f6);

/** The film selector of a drawn mesh: its occurrence path, or
 *  `body N` in a part document. */
function selectorOf(mesh) {
  return mesh.userData.instance ?? `body ${mesh.userData.index}`;
}

export class FilmPlayer {
  constructor(page) {
    this.page = page;
    this.film = null;
    this.look = new FilmLook(page.viewer);
    this.time = 0;
    this.playing = false;
    this.frame = 0;
    this.jointsApplied = {};
    this.jointQueue = Promise.resolve();
    this.listeners = new Set();
    // Opening a film (prepare): what it waits on, for the stage to show
    // (`{ stage, started, quiet }`), and the request it belongs to.
    this.preparing = null;
    this.ticket = 0;
  }

  get active() {
    return this.film != null;
  }

  onChange(listener) {
    this.listeners.add(listener);
  }

  notify() {
    for (const listener of this.listeners) listener();
  }

  /** Open the project's film `path` from wherever the page is: the
   *  stage shows at once with what it waits on — the kernel loading,
   *  the model building (live when the shown model cannot carry a film),
   *  the film script running — and each wait is retried until the film
   *  opens, the build fails (its outcome shows like any build's), or
   *  another request (an edit, another film, closing) supersedes it.
   *  A film already open (an edit to its script) keeps playing on
   *  screen meanwhile. Resolves to the film reply, or `{ ok: false }`. */
  async prepare(path, { live = false } = {}) {
    const page = this.page;
    const ticket = ++this.ticket;
    const current = () => ticket === this.ticket && page.projects.filmOpen && page.projects.path === path;
    const superseded = { ok: false, superseded: true };
    const started = performance.now();
    const stage = (text) => {
      if (!current()) return;
      this.preparing = { stage: text, started, quiet: this.active };
      this.notify();
    };
    const pause = () => new Promise((r) => setTimeout(r, 120));
    try {
      for (let attempt = 0; attempt < 6; attempt++) {
        if (!page.host.info) {
          stage("Loading the kernel");
          while (!page.host.info) {
            await pause();
            if (!current()) return superseded;
          }
        }
        if (page.running) {
          stage("Building the model");
          await page.runDone;
          if (!current()) return superseded;
          continue;
        }
        const geometry = page.projects.geometrySource();
        const context = page.lastGood?.envelope?.film_context;
        const ready = !live && page.lastGood && page.geometryKey === geometry && (page.live || context);
        if (!ready) {
          stage("Building the model");
          // A pre-built model from before films carries no film context:
          // build it live.
          const built = await page.run({ live: live || (page.geometryKey === geometry && !context), forFilm: true });
          if (!current()) return superseded;
          // Not run (the kernel was loading, or another build superseded
          // this one): go round.
          if (built == null) continue;
          live = false;
          if (!built.ok || (!built.baked && built.envelope?.outcome !== "served")) return { ok: false, outcome: built.outcome ?? null };
          continue;
        }
        stage("Running the film script");
        const reply = await this.run(page.editor.value, path);
        return current() ? reply : superseded;
      }
      page.status("the film could not start: the model kept changing while it was prepared", true);
      return { ok: false };
    } finally {
      if (ticket === this.ticket) {
        this.preparing = null;
        this.notify();
      }
    }
  }

  /** Run the film script `source` (project file `path`) and open it.
   *  Resolves to the worker's reply (`ok` false: the outcome was shown
   *  on the editor like a build failure). */
  async run(source, path) {
    const page = this.page;
    const context = page.lastGood?.envelope?.film_context ?? null;
    if (!page.live && !context && !this.ensuring) {
      // A pre-built model from before films: build it live (that run
      // plays the film when it serves).
      this.ensuring = true;
      try {
        if (!(await page.ensureLive())) return { ok: false, outcome: { outcome: "gate_failure", step: "film", detail: "the model did not build" } };
      } finally {
        this.ensuring = false;
      }
      return { ok: this.active, film: this.film, context: this.context };
    }
    const reply = await page.host.request({ type: "film", source, name: path, context });
    if (!reply.ok) {
      page.showOutcome(reply.outcome);
      return reply;
    }
    if (reply.film.tracks.some((t) => t.channel.kind === "joint") && !page.live) {
      // Joints move occurrences in the worker's document: build it live.
      page.status("building the assembly live for the film's joints…");
      if (!(await page.ensureLive())) return { ok: false };
    }
    this.open(reply.film, reply.context);
    page.status(`film ${reply.film.name}: ${reply.film.steps.length} step(s), ${reply.film.duration.toFixed(1)} s`);
    return reply;
  }

  /** Show `film`: capture the placements the film animates from, put
   *  the look on, and show its first frame. */
  open(film, context) {
    const viewer = this.page.viewer;
    const reopening = this.active;
    this.pause();
    this.film = film;
    this.context = context;
    if (!reopening) {
      this.saved = {
        explode: viewer.explode ?? 0,
        section: viewer.section.state(),
        handles: viewer.section.gizmo.enabled,
        camera: { target: viewer.orbit.target.clone(), spherical: viewer.orbit.spherical.clone(), fov: viewer.camera.fov },
      };
    }
    viewer.section.gizmo.setEnabled(false);
    viewer.section.linkFrozen = true;
    viewer.setExplode(0);
    this.look.apply(film, viewer.visibleBounds());
    this.capture();
    this.jointsApplied = {};
    this.time = Math.min(this.time, film.duration);
    document.body.classList.add("film-open");
    this.apply(this.time);
    this.notify();
  }

  /** The animated meshes and their placements before any film offset,
   *  and the bounds the camera and section fit (static: the camera does
   *  not chase moving parts). */
  capture() {
    const viewer = this.page.viewer;
    this.look.castShadows();
    this.objects = [...viewer.bodies, ...viewer.edgeLines].map((object) => ({
      object,
      // The film owns the matrix while it plays: an auto-updating one
      // (a pre-built part's meshes) is rebuilt from position/rotation
      // at every render, which silently dropped every move and turn.
      autoUpdate: (object.userData.filmAutoUpdate ??= object.matrixAutoUpdate),
      mesh: object.isLineSegments ? viewer.bodies.find((m) => m.userData.lines === object) : object,
      base: (object.userData.base ?? object.matrix).clone(),
      // The viewer's own visibility, recorded once: a frame the film
      // already applied hides parts it has not shown yet.
      visible: (object.userData.filmVisible ??= object.visible),
    }));
    for (const { object } of this.objects) object.matrixAutoUpdate = false;
    this.selectors = [...new Set(viewer.bodies.map(selectorOf))];
    this.bounds = viewer.visibleBounds();
    this.targetBounds = new Map();
    for (const mesh of viewer.bodies) {
      const [lo, hi] = mesh.userData.body.bounds;
      const box = new THREE.Box3(new THREE.Vector3(...lo), new THREE.Vector3(...hi)).applyMatrix4(mesh.matrixWorld);
      const path = selectorOf(mesh);
      // Every prefix of the path (a sub-assembly target) gets the box.
      const names = path.startsWith("body ") ? [path] : path.split("/").map((_, i, a) => a.slice(0, i + 1).join("/"));
      for (const name of [...names, "all"]) {
        const b = this.targetBounds.get(name) ?? new THREE.Box3();
        this.targetBounds.set(name, b.union(box));
      }
    }
  }

  /** Close the film and put the viewer back as it was. */
  close() {
    // A film still being prepared stops too.
    this.ticket += 1;
    if (this.preparing) {
      this.preparing = null;
      this.notify();
    }
    if (!this.active) return;
    this.pause();
    const viewer = this.page.viewer;
    for (const { object, base, visible, autoUpdate } of this.objects) {
      object.matrix.copy(base);
      object.matrixAutoUpdate = autoUpdate;
      delete object.userData.filmAutoUpdate;
      object.matrixWorldNeedsUpdate = true;
      object.visible = visible;
      delete object.userData.filmVisible;
      this.setOpacity(object, 1);
      if (object.material.emissive) object.material.emissive.setHex(0);
    }
    this.film = null;
    document.body.classList.remove("film-open");
    this.look.remove();
    viewer.section.linkFrozen = false;
    viewer.section.setPlanes(this.saved.section.planes, this.saved.section);
    viewer.section.gizmo.setEnabled(this.saved.handles);
    viewer.setExplode(this.saved.explode);
    viewer.camera.fov = this.saved.camera.fov;
    viewer.camera.updateProjectionMatrix();
    viewer.orbit.target.copy(this.saved.camera.target);
    viewer.orbit.spherical.copy(this.saved.camera.spherical);
    viewer.orbit.update();
    this.notify();
  }

  /** After the viewer rebuilt its meshes (a mode change or a rebuild),
   *  re-capture and re-apply. */
  rebound() {
    if (!this.active) return;
    this.capture();
    this.apply(this.time);
  }

  setOpacity(object, opacity) {
    const m = object.material;
    const see = opacity < 0.999;
    if (m.transparent !== see || (see && m.opacity !== opacity)) {
      m.transparent = see;
      m.opacity = opacity;
      m.depthWrite = !see;
      m.needsUpdate = true;
    }
  }

  /** Apply the film at time `t` to the scene, `camera` (default the
   *  viewer's) framed for `aspect` (default the film's). Joints are
   *  awaited. */
  async apply(t, { camera = this.page.viewer.camera, aspect = null, render = true } = {}) {
    if (!this.active) return;
    const film = this.film;
    const viewer = this.page.viewer;
    this.time = Math.max(0, Math.min(film.duration, t));
    const state = stateAt(film, this.time, this.selectors);
    await this.applyJoints(state.joints, render);
    const explodeDirections = viewer.explodePlan?.directions;
    const m = new THREE.Matrix4();
    const turn = new THREE.Matrix4();
    for (const entry of this.objects) {
      const { object, mesh, base } = entry;
      const item = state.items.get(selectorOf(mesh));
      if (!item) continue;
      m.identity();
      for (const r of item.turns) {
        const o = new THREE.Vector3(...r.origin);
        turn.makeTranslation(-o.x, -o.y, -o.z);
        const rot = new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(...r.direction), THREE.MathUtils.degToRad(r.degrees));
        const back = new THREE.Matrix4().makeTranslation(o.x, o.y, o.z);
        m.premultiply(turn).premultiply(rot).premultiply(back);
      }
      const d = explodeDirections?.get(mesh.userData.instance?.split("/")[0]);
      const ex = d ? [d.x * state.explode, d.y * state.explode, d.z * state.explode] : [0, 0, 0];
      m.premultiply(new THREE.Matrix4().makeTranslation(item.offset[0] + ex[0], item.offset[1] + ex[1], item.offset[2] + ex[2]));
      object.matrix.multiplyMatrices(m, base);
      object.matrixWorldNeedsUpdate = true;
      object.visible = entry.visible && item.visible && item.opacity > 0.001;
      this.setOpacity(object, item.opacity);
      if (object.material.emissive) object.material.emissive.copy(HIGHLIGHT).multiplyScalar(0.45 * item.highlight);
    }
    viewer.group.updateMatrixWorld(true);
    this.applySection(state.section);
    this.caption = state.caption;
    this.frameCamera(camera, aspect);
    if (render) {
      viewer.render();
      this.notify();
    }
  }

  /** Joints (placement-only assembly parameters) in the worker; each
   *  value is asked for once, and the occurrences' placements re-read. */
  async applyJoints(joints, render = true) {
    const names = Object.keys(joints).filter((n) => Math.abs((this.jointsApplied[n] ?? NaN) - joints[n]) > 1e-9 || !(n in this.jointsApplied));
    if (!names.length) return;
    for (const name of names) {
      const reply = await this.page.host.setParam(name, joints[name]);
      this.jointsApplied[name] = joints[name];
      const assembly = reply.ok ? reply.envelope?.assembly : null;
      if (assembly) {
        this.page.viewer.place(assembly, { render });
        // The placements moved under the film: re-read its bases.
        for (const entry of this.objects) entry.base.copy(entry.object.userData.base ?? entry.object.matrix);
      }
    }
  }

  applySection(section) {
    const s = this.page.viewer.section;
    const plane = section && this.bounds ? sectionAt(section, this.bounds) : null;
    if (!plane) {
      if (s.active) s.setPlanes([]);
      return;
    }
    if (s.planes.length === 1 && s.combine === section.combine) {
      s.planes[0] = { normal: plane.normal, offset: plane.offset };
      s.update();
    } else {
      s.setPlanes([plane], { combine: section.combine });
    }
  }

  /** Frame `camera` on the film's camera pose: the target's bounding
   *  sphere fits inside the narrower field of view of `aspect`. */
  frameCamera(camera, aspect = null) {
    const film = this.film;
    const viewer = this.page.viewer;
    const pose = cameraAt(film, this.time);
    const centre = new THREE.Vector3();
    let radius = 0;
    let weight = 0;
    for (const [name, w] of Object.entries(pose.targets)) {
      const box = this.targetBounds.get(name) ?? this.targetBounds.get("all");
      if (!box || box.isEmpty()) continue;
      centre.addScaledVector(box.getCenter(new THREE.Vector3()), w);
      radius += (w * box.getSize(new THREE.Vector3()).length()) / 2;
      weight += w;
    }
    if (weight <= 0) return;
    centre.multiplyScalar(1 / weight);
    radius = Math.max(1e-9, radius / weight);
    const filmAspect = film.size[0] / film.size[1];
    const frameAspect = aspect ?? filmAspect;
    // The film's vertical field of view is 30 degrees; a pane wider or
    // narrower than the film shows the film's frame inside it.
    const base = THREE.MathUtils.degToRad(30);
    let vfov = base;
    if (aspect == null && camera === viewer.camera) {
      const paneAspect = camera.aspect;
      if (paneAspect < filmAspect) vfov = 2 * Math.atan(Math.tan(base / 2) * (filmAspect / paneAspect));
    }
    const hfov = 2 * Math.atan(Math.tan(base / 2) * frameAspect);
    const fit = Math.min(base, hfov);
    const distance = radius / Math.sin(fit / 2) / pose.zoom;
    const az = THREE.MathUtils.degToRad(pose.azimuth);
    const el = THREE.MathUtils.degToRad(pose.elevation);
    camera.fov = THREE.MathUtils.radToDeg(vfov);
    if (aspect != null) camera.aspect = aspect;
    camera.position.set(
      centre.x + distance * Math.cos(el) * Math.cos(az),
      centre.y + distance * Math.cos(el) * Math.sin(az),
      centre.z + distance * Math.sin(el),
    );
    camera.up.set(0, 0, 1);
    camera.lookAt(centre);
    camera.near = Math.max(distance * 1e-3, distance - 4 * radius);
    camera.far = distance + 6 * radius;
    camera.updateProjectionMatrix();
    if (camera === viewer.camera) {
      // Keep the orbit control on the same pose, so a drag after
      // pausing starts from the film's view.
      viewer.orbit.target.copy(centre);
      viewer.orbit.spherical.set(distance, Math.PI / 2 - el, az);
    }
  }

  play() {
    if (!this.active || this.playing) return;
    if (this.time >= this.film.duration - 1e-6) this.time = 0;
    this.playing = true;
    let last = performance.now();
    let busy = false;
    const tick = async (now) => {
      if (!this.playing) return;
      const dt = (now - last) / 1000;
      last = now;
      if (!busy) {
        busy = true;
        await this.apply(this.time + dt);
        busy = false;
        if (this.time >= this.film.duration) {
          this.playing = false;
          this.notify();
          return;
        }
      }
      this.frame = requestAnimationFrame(tick);
    };
    this.frame = requestAnimationFrame(tick);
    this.notify();
  }

  pause() {
    this.playing = false;
    cancelAnimationFrame(this.frame);
    this.notify();
  }

  seek(t) {
    this.pause();
    return this.apply(t);
  }
}
