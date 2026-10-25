// The simulation overlay (PK-948): a solved field drawn on the analysis
// surface it was sampled on, with a colour ramp and a legend, and the
// part optionally drawn in its deformed shape. The overlay lives beside
// the viewer's body group (a rebuild clears the bodies, never this);
// the page clears it whenever the model on screen changes, so a field
// is never drawn on a part it was not solved for.
import * as THREE from "./vendor/three.module.min.js";

/** The fields a solve returns, with their legend text. */
export const FIELDS = {
  von_mises_mpa: { label: "von Mises stress", unit: "MPa" },
  max_principal_mpa: { label: "max principal stress", unit: "MPa" },
  displacement_mm: { label: "displacement", unit: "mm" },
  zz_error_mpa: { label: "estimated stress error", unit: "MPa" },
  temperature_c: { label: "temperature", unit: "°C" },
};

// Blue -> cyan -> green -> yellow -> red: low values cool, high hot.
const RAMP = [
  [0.19, 0.3, 0.86],
  [0.13, 0.72, 0.9],
  [0.3, 0.78, 0.36],
  [0.98, 0.8, 0.2],
  [0.86, 0.2, 0.17],
];

/** The ramp colour at `t` in [0, 1]. */
export function rampColor(t) {
  const x = Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0)) * (RAMP.length - 1);
  const i = Math.min(RAMP.length - 2, Math.floor(x));
  const f = x - i;
  return RAMP[i].map((c, k) => c + f * (RAMP[i + 1][k] - c));
}

/** The CSS gradient of the ramp, for the legend bar. */
export function rampGradient() {
  const stops = RAMP.map((c, i) => `rgb(${c.map((v) => Math.round(v * 255)).join(",")}) ${(100 * i) / (RAMP.length - 1)}%`);
  return `linear-gradient(to right, ${stops.join(", ")})`;
}

/** The [min, max] of `field` over every mesh (negative ZZ entries mark
 *  "no estimate" and are left out). */
export function fieldRange(meshes, field) {
  let lo = Infinity;
  let hi = -Infinity;
  for (const mesh of meshes) {
    for (const v of mesh.fields[field] ?? []) {
      if (field === "zz_error_mpa" && v < 0) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  }
  return Number.isFinite(lo) ? [lo, hi] : [0, 0];
}

/** Deformation scale at "full": the largest displacement drawn as a
 *  tenth of the field meshes' bounding-box diagonal (model units). */
export function autoDeformScale(meshes, mmPerUnit) {
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  let most = 0;
  for (const mesh of meshes) {
    const p = mesh.positions;
    for (let i = 0; i < p.length; i += 3) {
      for (let k = 0; k < 3; k++) {
        lo[k] = Math.min(lo[k], p[i + k]);
        hi[k] = Math.max(hi[k], p[i + k]);
      }
    }
    const d = mesh.displacement;
    for (let i = 0; i < d.length; i += 3) most = Math.max(most, Math.hypot(d[i], d[i + 1], d[i + 2]) / mmPerUnit);
  }
  const diagonal = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  return most > 0 && Number.isFinite(diagonal) ? (0.1 * diagonal) / most : 0;
}

export class FieldOverlay {
  constructor(viewer) {
    this.viewer = viewer;
    this.group = new THREE.Group();
    viewer.scene.add(this.group);
    this.meshes = null;
  }

  /** Draw `field` of `meshes` (a simulate reply's), displaced by
   *  `deform` (0 = undeformed, 1 = the auto scale). Hides the bodies. */
  show(meshes, field, deform, mmPerUnit) {
    this.disposeMeshes();
    this.meshes = meshes;
    const [lo, hi] = fieldRange(meshes, field);
    const span = hi - lo || 1;
    const scale = deform * autoDeformScale(meshes, mmPerUnit);
    for (const mesh of meshes) {
      const positions = Float32Array.from(mesh.positions);
      if (scale > 0 && mesh.displacement.length === positions.length) {
        for (let i = 0; i < positions.length; i++) positions[i] += (scale * mesh.displacement[i]) / mmPerUnit;
      }
      const values = mesh.fields[field] ?? [];
      const colors = new Float32Array(positions.length);
      for (let v = 0; v < positions.length / 3; v++) {
        const value = values[v];
        const c = field === "zz_error_mpa" && value < 0 ? [0.6, 0.6, 0.6] : rampColor((value - lo) / span);
        colors.set(c, 3 * v);
      }
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
      geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
      geometry.setIndex(new THREE.BufferAttribute(Uint32Array.from(mesh.indices), 1));
      geometry.computeVertexNormals();
      const material = new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide });
      this.group.add(new THREE.Mesh(geometry, material));
    }
    this.setBodiesVisible(false);
    this.viewer.render();
    return [lo, hi];
  }

  /** Remove the overlay and show the bodies again. */
  clear() {
    if (!this.meshes) return;
    this.disposeMeshes();
    this.meshes = null;
    this.setBodiesVisible(true);
    this.viewer.render();
  }

  get active() {
    return this.meshes !== null;
  }

  disposeMeshes() {
    for (const m of this.group.children) {
      m.geometry.dispose();
      m.material.dispose();
    }
    this.group.clear();
  }

  setBodiesVisible(visible) {
    // The viewer's own per-body visibility (assembly occurrences hidden
    // in the tree) is kept: only hide what is shown, and restore that.
    if (!visible) {
      this.hidden = [...this.viewer.bodies, ...this.viewer.edgeLines].filter((m) => m.visible);
      for (const m of this.hidden) m.visible = false;
    } else {
      for (const m of this.hidden ?? []) m.visible = true;
      this.hidden = null;
    }
  }
}
