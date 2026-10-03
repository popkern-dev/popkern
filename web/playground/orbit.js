// A small orbit control (no addon): one pointer orbits (or pans with
// the right button / shift), two pinch-zoom and pan, the wheel zooms,
// arrow keys orbit. Z is up.
import * as THREE from "./vendor/three.module.min.js";

export class Orbit {
  constructor(camera, element, onChange) {
    this.camera = camera;
    this.element = element;
    this.onChange = onChange;
    this.target = new THREE.Vector3();
    this.spherical = new THREE.Spherical(5, Math.PI / 3, Math.PI / 4);
    this.dragging = null;
    // Active pointers by id: one drives orbit (or pan with the right
    // button / shift), two drive pinch-zoom plus two-finger pan.
    this.pointers = new Map();
    this.pinch = null;
    // Path length of the current gesture, so a drag does not also pick.
    this.travel = 0;
    element.addEventListener("pointerdown", (e) => {
      if (this.pointers.size === 0) this.travel = 0;
      this.pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
      element.setPointerCapture(e.pointerId);
      if (this.pointers.size === 2) {
        this.dragging = null;
        this.pinch = this.pinchState();
      } else if (this.pointers.size === 1) {
        this.dragging = { x: e.clientX, y: e.clientY, button: e.button, shift: e.shiftKey };
      }
    });
    element.addEventListener("pointermove", (e) => {
      const p = this.pointers.get(e.pointerId);
      if (!p) return;
      this.travel += Math.hypot(e.clientX - p.x, e.clientY - p.y);
      p.x = e.clientX;
      p.y = e.clientY;
      if (this.pinch) {
        const next = this.pinchState();
        if (!next) return;
        if (next.distance > 0 && this.pinch.distance > 0) this.spherical.radius *= this.pinch.distance / next.distance;
        this.pan(next.x - this.pinch.x, next.y - this.pinch.y);
        this.pinch = next;
        this.update();
        return;
      }
      if (!this.dragging) return;
      const dx = e.clientX - this.dragging.x;
      const dy = e.clientY - this.dragging.y;
      this.dragging.x = e.clientX;
      this.dragging.y = e.clientY;
      if (this.dragging.button === 2 || this.dragging.shift) {
        this.pan(dx, dy);
      } else {
        this.spherical.theta -= dx * 0.008;
        this.spherical.phi = Math.min(Math.PI - 0.05, Math.max(0.05, this.spherical.phi - dy * 0.008));
      }
      this.update();
    });
    const release = (e) => {
      this.pointers.delete(e.pointerId);
      this.pinch = this.pointers.size >= 2 ? this.pinchState() : null;
      // Lifting one finger of a pinch hands orbit to the other from
      // where it now is, so the model does not jump.
      const rest = this.pointers.size === 1 ? [...this.pointers.values()][0] : null;
      this.dragging = rest ? { x: rest.x, y: rest.y, button: 0, shift: false } : null;
    };
    element.addEventListener("pointerup", release);
    element.addEventListener("pointercancel", release);
    element.addEventListener("contextmenu", (e) => e.preventDefault());
    element.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        this.spherical.radius *= Math.exp(e.deltaY * 0.001);
        this.update();
      },
      { passive: false },
    );
    element.addEventListener("keydown", (e) => {
      const step = 0.05;
      if (e.key === "ArrowLeft") this.spherical.theta += step;
      else if (e.key === "ArrowRight") this.spherical.theta -= step;
      else if (e.key === "ArrowUp") this.spherical.phi = Math.max(0.05, this.spherical.phi - step);
      else if (e.key === "ArrowDown") this.spherical.phi = Math.min(Math.PI - 0.05, this.spherical.phi + step);
      else if (e.key === "+" || e.key === "=") this.spherical.radius *= 0.9;
      else if (e.key === "-") this.spherical.radius *= 1.1;
      else return;
      e.preventDefault();
      this.update();
    });
  }

  /** Midpoint and separation of the first two active pointers. */
  pinchState() {
    const [a, b] = [...this.pointers.values()];
    if (!a || !b) return null;
    return { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, distance: Math.hypot(a.x - b.x, a.y - b.y) };
  }

  /** Whether the last gesture moved far enough to count as a drag. */
  wasDrag() {
    return this.travel > 6;
  }

  pan(dx, dy) {
    const scale = this.spherical.radius * 0.0015;
    const right = new THREE.Vector3();
    const up = new THREE.Vector3(0, 0, 1);
    this.camera.getWorldDirection(right);
    right.cross(up).normalize();
    const trueUp = right.clone().cross(new THREE.Vector3().subVectors(this.target, this.camera.position).normalize()).normalize();
    this.target.addScaledVector(right, -dx * scale);
    this.target.addScaledVector(trueUp, dy * scale);
  }

  fit(bounds) {
    const min = new THREE.Vector3(...bounds[0]);
    const max = new THREE.Vector3(...bounds[1]);
    this.target = min.clone().add(max).multiplyScalar(0.5);
    const size = max.clone().sub(min).length() || 1;
    this.spherical.radius = size * 1.4;
    this.update();
  }

  update() {
    // Z is up: spherical about the Z axis.
    const s = this.spherical;
    const x = s.radius * Math.sin(s.phi) * Math.cos(s.theta);
    const y = s.radius * Math.sin(s.phi) * Math.sin(s.theta);
    const z = s.radius * Math.cos(s.phi);
    this.camera.position.set(this.target.x + x, this.target.y + y, this.target.z + z);
    this.camera.up.set(0, 0, 1);
    this.camera.lookAt(this.target);
    this.onChange();
  }
}
