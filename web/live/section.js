// The section view (six-axis-arm §4.6): a display-only clipping plane
// through the model, with each solid's cut capped so a section reads
// as solid material, not a hollow shell. Nothing is cut in the kernel:
// the plane clips every body's material, and a per-body stencil pass
// (back faces count up, front faces count down — the parity of a
// closed mesh) marks where the plane lies inside that body, where its
// cap — the body's colour, hatched — is drawn.
import * as THREE from "./vendor/three.module.min.js";

const AXES = { x: new THREE.Vector3(1, 0, 0), y: new THREE.Vector3(0, 1, 0), z: new THREE.Vector3(0, 0, 1) };

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

export class Section {
  constructor(viewer) {
    this.viewer = viewer;
    this.axis = null;
    this.position = 0.5;
    this.plane = new THREE.Plane(new THREE.Vector3(-1, 0, 0), 0);
    this.caps = new THREE.Group();
    this.caps.name = "section-caps";
    this.hatch = null;
    this.extent = null;
    viewer.renderer.localClippingEnabled = true;
  }

  /** Whether a section is shown. */
  get active() {
    return this.axis != null;
  }

  /** Cut along `axis` (x | y | z, or null to clear) at `position`
   *  (0..1 across the model's bounds when the cut was first made). The
   *  material on the axis's positive side is removed. */
  set(axis, position) {
    const changedAxis = axis !== this.axis;
    this.axis = AXES[axis] ? axis : null;
    this.position = Math.min(1, Math.max(0, position));
    if (!this.active) {
      this.extent = null;
      this.apply();
      return;
    }
    if (changedAxis || !this.extent) this.extent = this.viewer.visibleBounds();
    this.apply();
  }

  /** The plane for the current axis and position (keeps x ≤ x0). */
  updatePlane() {
    const [lo, hi] = this.extent ?? [[0, 0, 0], [0, 0, 0]];
    const i = "xyz".indexOf(this.axis);
    const at = lo[i] + (hi[i] - lo[i]) * this.position;
    const normal = AXES[this.axis].clone().negate();
    this.plane.set(normal, at);
    return at;
  }

  /** Re-apply to the current scene: called after every redraw. */
  apply() {
    const viewer = this.viewer;
    this.clearCaps();
    const planes = this.active ? [this.plane] : null;
    if (this.active) this.updatePlane();
    for (const object of [...viewer.bodies, ...viewer.edgeLines]) {
      object.material.clippingPlanes = planes;
      object.material.needsUpdate = true;
    }
    if (!this.active || !this.extent) return viewer.render();
    const [lo, hi] = this.extent;
    const span = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
    this.hatch ??= hatchTexture();
    this.hatch.repeat.set(24, 24);
    const geometry = new THREE.PlaneGeometry(2 * span, 2 * span);
    viewer.bodies.forEach((mesh, k) => {
      if (!mesh.visible) return;
      const order = 10 + 2 * k;
      for (const [side, op] of [
        [THREE.BackSide, THREE.IncrementWrapStencilOp],
        [THREE.FrontSide, THREE.DecrementWrapStencilOp],
      ]) {
        const material = new THREE.MeshBasicMaterial({
          side,
          clippingPlanes: [this.plane],
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
      const color = mesh.material.color ? mesh.material.color.clone().multiplyScalar(0.85) : new THREE.Color(0x9aa3ad);
      const cap = new THREE.Mesh(
        geometry,
        new THREE.MeshStandardMaterial({
          color,
          map: this.hatch,
          metalness: 0.0,
          roughness: 0.8,
          side: THREE.DoubleSide,
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
      this.placeCap(cap);
      this.caps.add(cap);
    });
    viewer.scene.add(this.caps);
    viewer.render();
  }

  /** Put a cap quad on the plane, facing along its normal. */
  placeCap(cap) {
    const centre = new THREE.Vector3();
    const [lo, hi] = this.extent;
    centre.set((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2);
    this.plane.projectPoint(centre, cap.position);
    cap.lookAt(cap.position.clone().add(this.plane.normal));
  }

  /** Move the cut without rebuilding: the plane is shared by reference. */
  move(position) {
    if (!this.active) return;
    this.position = Math.min(1, Math.max(0, position));
    this.updatePlane();
    for (const cap of this.caps.children) this.placeCap(cap);
    this.viewer.render();
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
    for (const cap of this.caps.children) cap.material.dispose();
    if (this.caps.children[0]) this.caps.children[0].geometry.dispose();
    this.caps.clear();
    this.viewer.scene.remove(this.caps);
  }
}
