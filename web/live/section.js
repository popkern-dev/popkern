// The clipping view (docs/architecture/clipping-view.md; six-axis-arm
// §4.6): up to three display-only planes through the model, each solid's
// cut capped so a section reads as solid material, not a hollow shell,
// and outlined. Nothing is cut in the kernel: the planes clip every
// body's material, and per (body, plane) a stencil pass (back faces
// count up, front faces count down — the parity of a closed mesh) marks
// where that plane lies inside that body, where its cap — the body's
// colour, hatched — is drawn, clipped by the other planes.
import * as THREE from "./vendor/three.module.min.js";
import { fatLines } from "./edges.js";
import { SectionGizmo } from "./section-gizmo.js";
import {
  MAX_PLANES,
  axisPlane,
  cutOutline,
  flipped,
  fromLink,
  keeps,
  offsetAt,
  planeThrough,
  sliderPosition,
  straddles,
  toLink,
  unit,
  viewPlane,
} from "./section-geometry.js";

const OUTLINE = 0x111318;
const OUTLINE_PX = 2;

/** A diagonal hatch, multiplied over the cap colour. */
function hatchTexture() {
  const size = 64;
  const canvas = document.createElement("canvas");
  canvas.width = canvas.height = size;
  const ctx = canvas.getContext("2d");
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, size, size);
  ctx.strokeStyle = "#8a8f99";
  ctx.lineWidth = 6;
  for (let k = -size; k < 2 * size; k += size / 2) {
    ctx.beginPath();
    ctx.moveTo(k, 0);
    ctx.lineTo(k + size, size);
    ctx.stroke();
  }
  const texture = new THREE.CanvasTexture(canvas);
  texture.wrapS = texture.wrapT = THREE.RepeatWrapping;
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

/** A mesh's world box from its body's local bounds. */
function worldBox(mesh) {
  const [lo, hi] = mesh.userData.body.bounds;
  mesh.updateMatrixWorld();
  const box = new THREE.Box3(new THREE.Vector3(...lo), new THREE.Vector3(...hi)).applyMatrix4(mesh.matrixWorld);
  return [box.min.toArray(), box.max.toArray()];
}

/** three.js clips points with n·p + c < 0; a section plane removes
 *  normal·p > offset, so its three.js plane is (−normal, offset). */
function threePlane(plane, target = new THREE.Plane()) {
  const [x, y, z] = plane.normal;
  return target.set(new THREE.Vector3(-x, -y, -z), plane.offset);
}

export class Section {
  constructor(viewer) {
    this.viewer = viewer;
    this.planes = [];
    this.combine = "any";
    this.exclude = new Set();
    // The plane the gizmo and the slider act on.
    this.current = 0;
    // Axis-preset state kept for the old `set(axis, position)` API.
    this.axis = null;
    this.position = 0.5;
    this.extent = null;
    this.three = [];
    this.caps = new THREE.Group();
    this.caps.name = "section-caps";
    this.outlines = new THREE.Group();
    this.outlines.name = "section-outlines";
    this.hatch = null;
    this.capGeometry = null;
    this.showOutline = true;
    this.lastOutlineMs = 0;
    this.outlineFrame = 0;
    this.listeners = new Set();
    viewer.renderer.localClippingEnabled = true;
    this.gizmo = new SectionGizmo(this);
    this.pendingLink = new URLSearchParams(location.hash.slice(1)).get("section");
  }

  /** Whether a section is shown. */
  get active() {
    return this.planes.length > 0;
  }

  /** Notified (with no arguments) after every change of the cut. */
  onChange(listener) {
    this.listeners.add(listener);
    return () => this.listeners.delete(listener);
  }

  /** The cut as plain data: `{ planes, combine, exclude }`. */
  state() {
    return {
      planes: this.planes.map((p) => ({ normal: [...p.normal], offset: p.offset })),
      combine: this.combine,
      exclude: [...this.exclude],
    };
  }

  /** Replace the whole cut. `planes`: `[{ normal, offset }]` (world
   *  space; each removes the side its normal points to), at most three;
   *  `combine`: "any" (remove beyond any plane) or "all" (only where
   *  every plane removes: the corner cut-away); `exclude`: occurrence
   *  paths drawn whole. */
  setPlanes(planes, { combine = this.combine, exclude = [...this.exclude], current = this.current } = {}) {
    this.planes = planes
      .map((p) => ({ normal: unit(p.normal), offset: Number(p.offset) }))
      .filter((p) => p.normal && Number.isFinite(p.offset))
      .slice(0, MAX_PLANES);
    this.combine = combine === "all" ? "all" : "any";
    this.exclude = new Set(exclude);
    this.current = Math.max(0, Math.min(this.planes.length - 1, current));
    if (!this.active) this.axis = null;
    this.extent ??= this.viewer.visibleBounds();
    if (!this.active) this.extent = null;
    this.apply();
  }

  clear() {
    this.setPlanes([]);
  }

  /** Cut along `axis` (x | y | z, or null to clear) at `position`
   *  (0..1 across the model's bounds when the cut was first made),
   *  removing the axis's high side: the single-plane axis preset. */
  set(axis, position) {
    const changedAxis = axis !== this.axis;
    this.position = Math.min(1, Math.max(0, position));
    if (!["x", "y", "z"].includes(axis)) {
      this.axis = null;
      this.extent = null;
      return this.setPlanes([]);
    }
    if (changedAxis || !this.extent) this.extent = this.viewer.visibleBounds();
    this.axis = axis;
    const plane = axisPlane(axis, this.extent, this.position);
    if (!plane) return this.setPlanes([]);
    this.setPlanes([plane], { current: 0 });
  }

  /** Slide the current plane to `position` (0..1 across the bounds
   *  along its normal) without rebuilding caps. */
  move(position) {
    if (!this.active) return;
    this.position = Math.min(1, Math.max(0, position));
    const bounds = this.extent ?? this.viewer.visibleBounds();
    if (!bounds) return;
    const plane = this.planes[this.current];
    this.movePlane(this.current, offsetAt(plane.normal, bounds, this.position));
  }

  /** The current plane's 0..1 slider position. */
  sliderValue() {
    const bounds = this.extent ?? this.viewer.visibleBounds();
    const plane = this.planes[this.current];
    return plane && bounds ? sliderPosition(plane, bounds) : 0.5;
  }

  /** Add a plane from a preset: "x" | "y" | "z" | "view", or
   *  `{ point, normal }`. Replaces the current plane when `replace`. */
  addPlane(preset, { replace = false } = {}) {
    const bounds = this.extent ?? this.viewer.visibleBounds();
    if (!bounds) return null;
    this.extent = bounds;
    let plane = null;
    if (typeof preset === "string" && "xyz".includes(preset) && preset.length === 1) {
      plane = axisPlane(preset, bounds, 0.5);
    } else if (preset === "view") {
      const forward = new THREE.Vector3();
      this.viewer.camera.getWorldDirection(forward);
      plane = viewPlane(forward.toArray(), this.viewer.orbit.target.toArray());
    } else if (preset && preset.point && preset.normal) {
      plane = planeThrough(preset.point, preset.normal);
    }
    if (!plane) return null;
    const planes = this.planes.map((p) => ({ ...p }));
    if (replace && planes.length) planes[this.current] = plane;
    else if (planes.length < MAX_PLANES) planes.push(plane);
    else return null;
    const current = replace && this.planes.length ? this.current : planes.length - 1;
    this.axis = null;
    this.setPlanes(planes, { current });
    return plane;
  }

  removePlane(index = this.current) {
    const planes = this.planes.filter((_, i) => i !== index);
    this.setPlanes(planes, { current: Math.min(this.current, planes.length - 1) });
  }

  /** Remove the other side of plane `index`. */
  flip(index = this.current) {
    if (!this.planes[index]) return;
    this.planes[index] = flipped(this.planes[index]);
    this.apply();
  }

  /** Move plane `index` along its normal (offset only: caps and
   *  clipping follow by reference; the outline recomputes once per
   *  frame). */
  movePlane(index, offset) {
    const plane = this.planes[index];
    if (!plane || !Number.isFinite(offset)) return;
    plane.offset = offset;
    this.update();
  }

  /** Turn plane `index` to `normal`, keeping it through `pivot` (a
   *  world point; default the plane's anchor). */
  orientPlane(index, normal, pivot = null) {
    const plane = this.planes[index];
    const n = unit(normal);
    if (!plane || !n) return;
    const through = pivot ?? this.anchor(index);
    plane.normal = n;
    plane.offset = n[0] * through[0] + n[1] * through[1] + n[2] * through[2];
    this.update();
  }

  setCombine(combine) {
    this.combine = combine === "all" ? "all" : "any";
    this.apply();
  }

  /** Draw the occurrence at `path` (and what lies under it) whole. */
  setExcluded(path, excluded) {
    if (excluded) this.exclude.add(path);
    else this.exclude.delete(path);
    this.apply();
  }

  excluded(mesh) {
    const path = mesh.userData.instance;
    if (path == null) return false;
    for (const p of this.exclude) if (path === p || path.startsWith(`${p}/`)) return true;
    return false;
  }

  /** The point on plane `index` nearest the visible bounds' centre: the
   *  gizmo's anchor and the pivot for tilting. */
  anchor(index = this.current) {
    const plane = this.planes[index];
    const bounds = this.viewer.visibleBounds() ?? this.extent;
    if (!plane || !bounds) return [0, 0, 0];
    const [lo, hi] = bounds;
    const c = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
    const n = plane.normal;
    const d = n[0] * c[0] + n[1] * c[1] + n[2] * c[2] - plane.offset;
    return [c[0] - d * n[0], c[1] - d * n[1], c[2] - d * n[2]];
  }

  /** Whether the cut keeps world point `p` (a ray hit there is drawn). */
  keepsPoint(p) {
    return keeps(this.planes, this.combine, p, -1, 1e-9 * (this.capSpan ?? 1));
  }

  /** Clip `material` like the bodies (hover highlights, glows). */
  clipMaterial(material) {
    material.clippingPlanes = this.active ? this.three : null;
    material.clipIntersection = this.combine === "all";
    material.needsUpdate = true;
  }

  /** Re-apply to the current scene: called after every redraw. */
  apply() {
    if (this.pendingLink != null && this.viewer.bodies.length) {
      const link = fromLink(this.pendingLink);
      this.pendingLink = null;
      if (link.planes.length) {
        this.planes = link.planes;
        this.combine = link.combine;
        this.exclude = new Set(link.exclude);
        this.current = 0;
        this.extent = this.viewer.visibleBounds();
      }
    }
    const viewer = this.viewer;
    this.clearCaps();
    this.three = this.planes.map((p) => threePlane(p));
    const intersect = this.combine === "all";
    for (const object of [...viewer.bodies, ...viewer.edgeLines]) {
      const mesh = object.isLineSegments ? viewer.bodies.find((m) => m.userData.lines === object) : object;
      const cut = this.active && !(mesh && this.excluded(mesh));
      object.material.clippingPlanes = cut ? this.three : null;
      object.material.clipIntersection = intersect;
      object.material.needsUpdate = true;
    }
    if (this.active) this.buildCaps();
    this.gizmo.sync();
    this.writeLink();
    for (const listener of this.listeners) listener();
    viewer.render();
  }

  /** The capped (body, plane) pairs: visible solid bodies the cut
   *  applies to whose world box the plane passes through. A body wholly
   *  on one side is either drawn whole or removed whole: no cap, no
   *  outline, no stencil passes (the cycloid's 307 occurrences). */
  capPairs() {
    const pairs = [];
    for (const mesh of this.viewer.bodies) {
      if (!mesh.visible || mesh.userData.sheet || mesh.userData.ghost || this.excluded(mesh)) continue;
      const box = worldBox(mesh);
      this.planes.forEach((plane, i) => {
        if (straddles(plane, box)) pairs.push({ mesh, plane: i, box });
      });
    }
    return pairs;
  }

  /** A cheap re-apply after a plane moved or turned: the three.js
   *  planes are shared by reference, so only they, the caps' placement,
   *  the gizmo and (once per frame) the outline need updating — unless
   *  the plane now crosses a different set of bodies, which rebuilds. */
  update() {
    if (this.three.length !== this.planes.length) return this.apply();
    const pairs = this.capPairs();
    const key = pairs.map((p) => `${p.mesh.id}:${p.plane}`).join(",");
    if (key !== this.pairKey) return this.apply();
    this.planes.forEach((p, i) => threePlane(p, this.three[i]));
    for (const cap of this.caps.children) {
      // Negated copies (combine "all") follow their source plane.
      for (const [i, plane] of cap.userData.negated ?? []) plane.copy(this.three[i]).negate();
      this.placeCap(cap);
    }
    this.gizmo.sync();
    this.scheduleOutline();
    this.writeLink();
    for (const listener of this.listeners) listener();
    this.viewer.render();
  }

  buildCaps() {
    const viewer = this.viewer;
    const bounds = viewer.visibleBounds() ?? this.extent;
    if (!bounds) return;
    const [lo, hi] = bounds;
    const span = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
    this.hatch ??= hatchTexture();
    this.capGeometry = new THREE.PlaneGeometry(1, 1);
    this.capSpan = span;
    const intersect = this.combine === "all";
    let order = 10;
    const pairs = this.capPairs();
    this.pairKey = pairs.map((p) => `${p.mesh.id}:${p.plane}`).join(",");
    // One hatch density across the model, whatever each cap's size.
    const hatchPerUnit = 24 / (4 * span);
    for (const { mesh, plane: i, box } of pairs) {
      {
        for (const [side, op] of [
          [THREE.BackSide, THREE.IncrementWrapStencilOp],
          [THREE.FrontSide, THREE.DecrementWrapStencilOp],
        ]) {
          const material = new THREE.MeshBasicMaterial({
            side,
            clippingPlanes: [this.three[i]],
            depthWrite: false,
            depthTest: false,
            colorWrite: false,
            stencilWrite: true,
            stencilFunc: THREE.AlwaysStencilFunc,
            stencilFail: op,
            stencilZFail: op,
            stencilZPass: op,
          });
          const stencil = new THREE.Mesh(mesh.geometry, material);
          stencil.renderOrder = order;
          stencil.userData.sectionStencil = true;
          mesh.add(stencil);
        }
        // The cap on plane i, limited by the other planes: "any" keeps
        // what none of them removes; "all" bounds the removed corner,
        // so it lies where every other plane removes (their negations
        // clip the rest).
        const others = [];
        const negated = [];
        this.three.forEach((plane, j) => {
          if (j === i) return;
          if (intersect) {
            const copy = plane.clone().negate();
            negated.push([j, copy]);
            others.push(copy);
          } else {
            others.push(plane);
          }
        });
        const color = mesh.material.color ? mesh.material.color.clone().multiplyScalar(0.85) : new THREE.Color(0x9aa3ad);
        // A cap as large as its body (not the model): each cap pass fills
        // only the body's footprint on screen.
        const size = 1.5 * Math.hypot(box[1][0] - box[0][0], box[1][1] - box[0][1], box[1][2] - box[0][2]) || span;
        const hatch = this.hatch.clone();
        hatch.repeat.set(size * hatchPerUnit, size * hatchPerUnit);
        hatch.needsUpdate = true;
        const cap = new THREE.Mesh(
          this.capGeometry,
          new THREE.MeshStandardMaterial({
            color,
            map: hatch,
            metalness: 0.0,
            roughness: 0.8,
            side: THREE.DoubleSide,
            clippingPlanes: others.length ? others : null,
            clipIntersection: false,
            stencilWrite: true,
            stencilRef: 0,
            stencilFunc: THREE.NotEqualStencilFunc,
            stencilFail: THREE.ReplaceStencilOp,
            stencilZFail: THREE.ReplaceStencilOp,
            stencilZPass: THREE.ReplaceStencilOp,
          }),
        );
        cap.renderOrder = order + 1;
        cap.onAfterRender = (renderer) => renderer.clearStencil();
        cap.userData = { plane: i, negated, mesh, size };
        cap.scale.set(size, size, 1);
        this.placeCap(cap);
        this.caps.add(cap);
        order += 2;
      }
    }
    viewer.scene.add(this.caps);
    this.outlineNow();
    viewer.scene.add(this.outlines);
  }

  /** Put a cap quad on its plane over its body, facing along the
   *  plane's normal. */
  placeCap(cap) {
    const plane = this.planes[cap.userData.plane];
    const box = worldBox(cap.userData.mesh);
    const c = [(box[0][0] + box[1][0]) / 2, (box[0][1] + box[1][1]) / 2, (box[0][2] + box[1][2]) / 2];
    const n = plane.normal;
    const d = n[0] * c[0] + n[1] * c[1] + n[2] * c[2] - plane.offset;
    cap.position.set(c[0] - d * n[0], c[1] - d * n[1], c[2] - d * n[2]);
    cap.lookAt(cap.position.clone().add(new THREE.Vector3(...plane.normal)));
  }

  /** Recompute the cut outline at most once per animation frame. */
  scheduleOutline() {
    if (this.outlineFrame) return;
    this.outlineFrame = requestAnimationFrame(() => {
      this.outlineFrame = 0;
      this.outlineNow();
      this.viewer.render();
    });
  }

  /** The cut outline: every capped body's display mesh intersected with
   *  each plane, the other planes' removed parts dropped. */
  outlineNow() {
    for (const line of this.outlines.children) {
      line.geometry.dispose();
      line.material.dispose();
    }
    this.outlines.clear();
    if (!this.active || !this.showOutline) return;
    const started = performance.now();
    const ratio = this.viewer.renderer.getPixelRatio();
    // Every crossing body's segments in one fat-line draw.
    const parts = [];
    let total = 0;
    for (const { mesh, plane } of this.capPairs()) {
      const position = mesh.geometry.getAttribute("position").array;
      const index = mesh.geometry.index.array;
      mesh.updateMatrixWorld();
      const segments = cutOutline(position, index, mesh.matrixWorld.elements, this.planes, plane, this.combine);
      if (!segments.length) continue;
      parts.push(segments);
      total += segments.length;
    }
    if (total) {
      const all = new Float32Array(total);
      let at = 0;
      for (const part of parts) {
        all.set(part, at);
        at += part.length;
      }
      const line = fatLines(THREE, all, OUTLINE, OUTLINE_PX * ratio);
      line.renderOrder = 1000;
      line.userData.segments = total / 6;
      this.outlines.add(line);
    }
    this.lastOutlineMs = performance.now() - started;
  }

  /** Fat-line overlays whose viewport uniform the viewer refreshes. */
  overlays() {
    return this.outlines.children;
  }

  clearCaps() {
    for (const mesh of this.viewer.bodies) {
      for (const child of [...mesh.children]) {
        if (child.userData.sectionStencil) {
          mesh.remove(child);
          child.material.dispose();
        }
      }
    }
    for (const cap of this.caps.children) {
      cap.material.map?.dispose();
      cap.material.dispose();
    }
    this.capGeometry?.dispose();
    this.capGeometry = null;
    this.caps.clear();
    this.viewer.scene.remove(this.caps);
    for (const line of this.outlines.children) {
      line.geometry.dispose();
      line.material.dispose();
    }
    this.outlines.clear();
    this.viewer.scene.remove(this.outlines);
  }

  /** Mirror the cut into the page link (`#section=…`), keeping the
   *  link's other parameters. */
  writeLink() {
    if (this.pendingLink != null || this.linkFrozen) return;
    const params = new URLSearchParams(location.hash.slice(1));
    const text = toLink(this.state());
    if ((params.get("section") ?? "") === text) return;
    if (text) params.set("section", text);
    else params.delete("section");
    const hash = params.toString();
    history.replaceState(null, "", hash ? `#${hash}` : location.pathname + location.search);
  }

  /** What the tests and the diagnostics read. */
  describe() {
    return {
      ...this.state(),
      current: this.current,
      caps: this.caps.children.length,
      outline_segments: this.outlines.children.reduce((n, l) => n + (l.userData.segments ?? 0), 0),
      outline_ms: this.lastOutlineMs,
      link: toLink(this.state()),
    };
  }
}
