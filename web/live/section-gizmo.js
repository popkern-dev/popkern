// The clipping view's in-viewport handle (docs/architecture/clipping-view.md
// §5): at the current plane's anchor, a translucent square showing the
// plane, an arrow along its normal (drag: slide the plane) and two rings
// about its in-plane axes (drag: tilt it about that axis). Constant size
// on screen. A grab on a handle takes the pointer before the orbit
// control (capture-phase pointerdown), so dragging a handle never orbits.
import * as THREE from "./vendor/three.module.min.js";
import { inPlaneAxes, rotated } from "./section-geometry.js";

const PLANE_COLOR = 0x3b82f6;
const ARROW_COLOR = 0xf59e0b;
const RING_COLORS = [0xef4444, 0x22c55e];
const HOVER_COLOR = 0xffffff;
// Gizmo size as a fraction of the anchor's distance from the camera.
const SCREEN_SIZE = 0.085;

export class SectionGizmo {
  constructor(section) {
    this.section = section;
    this.viewer = section.viewer;
    this.enabled = true;
    this.group = new THREE.Group();
    this.group.name = "section-gizmo";
    this.group.renderOrder = 2000;
    this.dragging = null;
    this.justDragged = false;
    this.build();
    const canvas = this.viewer.canvas;
    canvas.addEventListener("pointerdown", (e) => this.down(e), { capture: true });
    canvas.addEventListener("pointermove", (e) => this.moveTo(e), { capture: true });
    canvas.addEventListener("pointerup", (e) => this.up(e), { capture: true });
    canvas.addEventListener("pointercancel", (e) => this.up(e), { capture: true });
    // A click that ends a handle drag must not also pick a face.
    canvas.addEventListener(
      "click",
      (e) => {
        if (this.justDragged) {
          this.justDragged = false;
          e.stopImmediatePropagation();
        }
      },
      { capture: true },
    );
  }

  build() {
    const overlay = (color, opacity = 1) =>
      new THREE.MeshBasicMaterial({ color, transparent: true, opacity, depthTest: false, depthWrite: false, side: THREE.DoubleSide });
    // The plane: a unit square in the gizmo's local xy (normal +z).
    this.square = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), overlay(PLANE_COLOR, 0.12));
    const border = new THREE.LineLoop(
      new THREE.BufferGeometry().setFromPoints([
        new THREE.Vector3(-1, -1, 0),
        new THREE.Vector3(1, -1, 0),
        new THREE.Vector3(1, 1, 0),
        new THREE.Vector3(-1, 1, 0),
      ]),
      new THREE.LineBasicMaterial({ color: PLANE_COLOR, transparent: true, opacity: 0.8, depthTest: false }),
    );
    // The arrow along +z (the plane's normal: the removed side).
    const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.03, 0.03, 1.1, 12), overlay(ARROW_COLOR));
    shaft.rotation.x = Math.PI / 2;
    shaft.position.z = 0.55;
    const head = new THREE.Mesh(new THREE.ConeGeometry(0.1, 0.3, 16), overlay(ARROW_COLOR));
    head.rotation.x = Math.PI / 2;
    head.position.z = 1.25;
    const arrowHit = new THREE.Mesh(new THREE.CylinderGeometry(0.16, 0.16, 1.6, 8), overlay(0, 0));
    arrowHit.rotation.x = Math.PI / 2;
    arrowHit.position.z = 0.7;
    this.arrow = new THREE.Group();
    this.arrow.add(shaft, head, arrowHit);
    this.arrow.userData = { handle: "slide", parts: [shaft, head], color: ARROW_COLOR };
    arrowHit.userData.handle = this.arrow;
    // Rings about local x (u) and local y (v).
    this.rings = RING_COLORS.map((color, k) => {
      const ring = new THREE.Mesh(new THREE.TorusGeometry(1.05, 0.02, 8, 64), overlay(color, 0.9));
      const hit = new THREE.Mesh(new THREE.TorusGeometry(1.05, 0.09, 6, 48), overlay(0, 0));
      const group = new THREE.Group();
      group.add(ring, hit);
      if (k === 0) group.rotation.y = Math.PI / 2;
      else group.rotation.x = Math.PI / 2;
      group.userData = { handle: k === 0 ? "tilt-u" : "tilt-v", parts: [ring], color };
      hit.userData.handle = group;
      return group;
    });
    // The square itself slides the plane too: the arrow is foreshortened
    // to a dot when the plane faces the viewer.
    this.square.userData = { handle: "slide", parts: [], color: PLANE_COLOR };
    this.square.userData.handle = this.square;
    this.square.userData.slide = true;
    this.handles = [this.arrow, ...this.rings];
    this.hitTargets = [arrowHit, ...this.rings.map((r) => r.children[1]), this.square];
    this.group.add(this.square, border, this.arrow, ...this.rings);
    for (const object of this.group.children) object.renderOrder = 2000;
    this.group.traverse((o) => (o.renderOrder = 2000));
    this.group.visible = false;
  }

  /** Follow the current plane (or hide); called after every change. */
  sync() {
    const section = this.section;
    const plane = section.planes[section.current];
    const show = this.enabled && !!plane && this.viewer.bodies.length > 0;
    if (!show) {
      this.group.visible = false;
      this.viewer.scene.remove(this.group);
      return;
    }
    if (!this.group.parent) this.viewer.scene.add(this.group);
    this.group.visible = true;
    const n = plane.normal;
    const [u, v] = inPlaneAxes(n);
    this.axes = { n, u, v };
    const basis = new THREE.Matrix4().makeBasis(new THREE.Vector3(...u), new THREE.Vector3(...v), new THREE.Vector3(...n));
    this.group.quaternion.setFromRotationMatrix(basis);
    this.group.position.set(...section.anchor());
    this.resize();
  }

  /** Keep a constant on-screen size: called before every render. */
  resize() {
    if (!this.group.visible) return;
    const distance = this.viewer.camera.position.distanceTo(this.group.position);
    this.group.scale.setScalar(Math.max(1e-9, distance * SCREEN_SIZE));
  }

  setEnabled(on) {
    this.enabled = on;
    this.sync();
    this.viewer.render();
  }

  /** The handle under a client point, or null. */
  handleAt(clientX, clientY) {
    if (!this.group.visible) return null;
    const rect = this.viewer.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
    const ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, this.viewer.camera);
    this.group.updateMatrixWorld(true);
    const hits = ray.intersectObjects(this.hitTargets, false);
    // The arrow and the rings win over the square they stand on.
    const hit = hits.find((h) => h.object !== this.square) ?? hits[0];
    return hit ? hit.object.userData.handle : null;
  }

  /** Screen position (CSS pixels, canvas-relative) of a world point. */
  screen(point) {
    const rect = this.viewer.canvas.getBoundingClientRect();
    const p = new THREE.Vector3(...point).project(this.viewer.camera);
    return [((p.x + 1) / 2) * rect.width + rect.left, ((1 - p.y) / 2) * rect.height + rect.top];
  }

  down(e) {
    if (e.button !== 0 || this.dragging) return;
    const handle = this.handleAt(e.clientX, e.clientY);
    if (!handle) return;
    e.stopImmediatePropagation();
    e.preventDefault();
    // Capture keeps the drag when the pointer leaves the canvas; it
    // throws for a pointer the browser no longer tracks (released mid
    // event, or synthetic), which must not cancel the grab.
    try {
      this.viewer.canvas.setPointerCapture(e.pointerId);
    } catch {
      // The drag still follows moves over the canvas.
    }
    const section = this.section;
    const plane = section.planes[section.current];
    const anchor = section.anchor();
    this.dragging = {
      pointer: e.pointerId,
      kind: handle === this.square ? "slide" : handle.userData.handle,
      start: [e.clientX, e.clientY],
      offset: plane.offset,
      normal: [...plane.normal],
      anchor,
      axes: { ...this.axes },
      centre: this.screen(anchor),
      moved: 0,
    };
    this.highlight(handle);
  }

  moveTo(e) {
    const drag = this.dragging;
    if (!drag || e.pointerId !== drag.pointer) {
      if (!drag && e.pointerType === "mouse" && !e.buttons) this.highlight(this.handleAt(e.clientX, e.clientY));
      return;
    }
    e.stopImmediatePropagation();
    drag.moved = Math.max(drag.moved, Math.hypot(e.clientX - drag.start[0], e.clientY - drag.start[1]));
    const section = this.section;
    if (drag.kind === "slide") {
      // The normal's direction on screen and its pixels per world unit.
      const a = this.screen(drag.anchor);
      const tip = drag.anchor.map((c, i) => c + drag.normal[i]);
      const b = this.screen(tip);
      const sx = b[0] - a[0];
      const sy = b[1] - a[1];
      const pixels = Math.hypot(sx, sy);
      // World size of a pixel at the anchor's depth.
      const camera = this.viewer.camera;
      const rect = this.viewer.canvas.getBoundingClientRect();
      const depth = camera.position.distanceTo(new THREE.Vector3(...drag.anchor));
      const worldPerPixel = (2 * depth * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / Math.max(1, rect.height);
      if (pixels * worldPerPixel < 0.3) {
        // The normal faces the viewer (foreshortened under 0.3 of its
        // length): dragging up pushes the plane away from the camera.
        const toCamera = new THREE.Vector3().subVectors(camera.position, new THREE.Vector3(...drag.anchor)).normalize();
        const facing = drag.normal[0] * toCamera.x + drag.normal[1] * toCamera.y + drag.normal[2] * toCamera.z;
        const up = drag.start[1] - e.clientY;
        section.movePlane(section.current, drag.offset - Math.sign(facing || 1) * up * worldPerPixel);
        return;
      }
      const along = ((e.clientX - drag.start[0]) * sx + (e.clientY - drag.start[1]) * sy) / pixels;
      section.movePlane(section.current, drag.offset + along / pixels);
      return;
    }
    // Tilt: the pointer's angle about the anchor on screen, applied about
    // the ring's axis (sign: the axis's sense toward the viewer).
    const angle = (x, y) => Math.atan2(y - drag.centre[1], x - drag.centre[0]);
    let delta = angle(e.clientX, e.clientY) - angle(...drag.start);
    const axis = drag.kind === "tilt-u" ? drag.axes.u : drag.axes.v;
    const forward = new THREE.Vector3();
    this.viewer.camera.getWorldDirection(forward);
    const toward = axis[0] * forward.x + axis[1] * forward.y + axis[2] * forward.z;
    if (toward > 0) delta = -delta;
    const normal = rotated(drag.normal, axis, -delta);
    if (normal) section.orientPlane(section.current, normal, drag.anchor);
  }

  up(e) {
    const drag = this.dragging;
    if (!drag || e.pointerId !== drag.pointer) return;
    e.stopImmediatePropagation();
    this.dragging = null;
    this.justDragged = drag.moved > 2;
    this.highlight(null);
    // A full re-apply settles the caps and the outline at the final pose.
    this.section.apply();
  }

  highlight(handle) {
    if (handle === this.hovered) return;
    for (const h of this.handles) {
      for (const part of h.userData.parts) part.material.color.setHex(h === handle ? HOVER_COLOR : h.userData.color);
    }
    this.hovered = handle;
    this.viewer.canvas.style.cursor = handle ? "grab" : "";
    this.viewer.render();
  }
}
