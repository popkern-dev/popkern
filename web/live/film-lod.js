// Screen-space detail for film export (PK-853). A display mesh is built
// for close inspection: an M4 screw's real thread is 20k triangles, the
// six-axis arm 4.3M, and headless Chrome's software GL spends ~16 s per
// 720p frame on them (vertex work: half the resolution or a flat
// material barely changes it). In a 720p frame most parts cover fewer
// pixels than they have triangles, so each export frame draws every
// occurrence at a vertex-clustered level whose cells are about a pixel
// and a half at its depth, keeping the full mesh wherever clustering would not
// save much. Export only: the viewer keeps full meshes for picking.
import * as THREE from "./vendor/three.module.min.js";

/** Cell size in pixels. Measured on the arm's 720p joints film: 1.5 px
 *  cells draw 0.48M of 4.3M triangles (4.3 s -> 1.4 s a frame) and
 *  change 153 of 921,600 pixels by more than 16 levels (2 by more than
 *  64); 0.75 px cells draw 0.98M and change 103. */
const CELL_PIXELS = 1.5;
/** Coarsest level: the body's bounding diagonal over 2^MIN_LEVEL. */
const MIN_LEVEL = 3;
/** Finest level tried. */
const MAX_LEVEL = 12;
/** A level is used only when it keeps under this share of triangles. */
const WORTH = 0.7;

/** The dominant axis and sign of a normal: clusters never merge across
 *  a crease (a screw head's flat top and its side stay separate). */
function facing(nx, ny, nz) {
  const ax = Math.abs(nx);
  const ay = Math.abs(ny);
  const az = Math.abs(nz);
  if (ax >= ay && ax >= az) return nx >= 0 ? 0 : 1;
  if (ay >= az) return ny >= 0 ? 2 : 3;
  return nz >= 0 ? 4 : 5;
}

/** Vertex-clustered copy of an indexed `geometry` with cells of side
 *  `cell` (its own units): each cluster's vertex is the mean of its
 *  members, and triangles collapsing onto fewer than three clusters go. */
export function clusterGeometry(geometry, cell) {
  const pos = geometry.attributes.position.array;
  const nor = geometry.attributes.normal?.array;
  const col = geometry.attributes.color;
  const colArr = col?.array;
  const colSize = col?.itemSize ?? 0;
  const index = geometry.index.array;
  const n = pos.length / 3;
  const of = new Int32Array(n);
  const keys = new Map();
  const sums = [];
  for (let i = 0; i < n; i++) {
    const x = pos[3 * i];
    const y = pos[3 * i + 1];
    const z = pos[3 * i + 2];
    const f = nor ? facing(nor[3 * i], nor[3 * i + 1], nor[3 * i + 2]) : 0;
    const key = `${Math.floor(x / cell)},${Math.floor(y / cell)},${Math.floor(z / cell)},${f}`;
    let c = keys.get(key);
    if (c === undefined) {
      c = sums.length;
      keys.set(key, c);
      sums.push({ x: 0, y: 0, z: 0, nx: 0, ny: 0, nz: 0, c: colSize ? new Array(colSize).fill(0) : null, k: 0 });
    }
    const s = sums[c];
    s.x += x;
    s.y += y;
    s.z += z;
    if (nor) {
      s.nx += nor[3 * i];
      s.ny += nor[3 * i + 1];
      s.nz += nor[3 * i + 2];
    }
    for (let j = 0; j < colSize; j++) s.c[j] += colArr[colSize * i + j];
    s.k++;
    of[i] = c;
  }
  const outIndex = [];
  for (let t = 0; t < index.length; t += 3) {
    const a = of[index[t]];
    const b = of[index[t + 1]];
    const c = of[index[t + 2]];
    if (a !== b && b !== c && a !== c) outIndex.push(a, b, c);
  }
  const m = sums.length;
  const p = new Float32Array(3 * m);
  const q = new Float32Array(3 * m);
  const k = colSize ? new Float32Array(colSize * m) : null;
  sums.forEach((s, i) => {
    p[3 * i] = s.x / s.k;
    p[3 * i + 1] = s.y / s.k;
    p[3 * i + 2] = s.z / s.k;
    const len = Math.hypot(s.nx, s.ny, s.nz) || 1;
    q[3 * i] = s.nx / len;
    q[3 * i + 1] = s.ny / len;
    q[3 * i + 2] = s.nz / len;
    for (let j = 0; j < colSize; j++) k[colSize * i + j] = s.c[j] / s.k;
  });
  const out = new THREE.BufferGeometry();
  out.setAttribute("position", new THREE.BufferAttribute(p, 3));
  if (nor) out.setAttribute("normal", new THREE.BufferAttribute(q, 3));
  if (k) out.setAttribute("color", new THREE.BufferAttribute(k, colSize));
  out.setIndex(outIndex);
  return out;
}

/** Per-frame level choice for the viewer's bodies during an export. */
export class FilmDetail {
  constructor(viewer) {
    this.viewer = viewer;
    /** geometry uuid -> Map(level -> geometry | null when not worth it) */
    this.levels = new Map();
    this.swapped = new Map();
  }

  level(full, k, diag) {
    let byLevel = this.levels.get(full.uuid);
    if (!byLevel) this.levels.set(full.uuid, (byLevel = new Map()));
    if (!byLevel.has(k)) {
      const lod = clusterGeometry(full, diag / 2 ** k);
      const worth = lod.index.count < WORTH * full.index.count;
      if (!worth) lod.dispose();
      byLevel.set(k, worth ? lod : null);
    }
    return byLevel.get(k);
  }

  /** Give every visible body the level its size on screen needs, for a
   *  frame `height` pixels tall through `camera`. */
  apply(camera, height) {
    // Section caps fill each cut body's closed geometry through the
    // stencil: a clustered level may open it, so a cut frame draws full.
    if (this.viewer.section?.active) return this.restore();
    const centre = new THREE.Vector3();
    const toCamera = new THREE.Vector3();
    const camPos = camera.getWorldPosition(new THREE.Vector3());
    const tanHalf = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    for (const mesh of this.viewer.bodies) {
      const full = this.swapped.get(mesh) ?? mesh.geometry;
      if (!mesh.visible || !full.index) continue;
      if (!full.boundingSphere) full.computeBoundingSphere();
      const radius = full.boundingSphere.radius;
      centre.copy(full.boundingSphere.center).applyMatrix4(mesh.matrixWorld);
      const depth = Math.max(toCamera.subVectors(centre, camPos).length() - radius, 1e-9);
      // World units per pixel at the body's nearest depth.
      const unit = (2 * depth * tanHalf) / height;
      const cell = CELL_PIXELS * unit;
      const diag = 2 * radius;
      let geometry = full;
      for (let k = MIN_LEVEL; k <= MAX_LEVEL; k++) {
        if (diag / 2 ** k > cell) continue;
        geometry = this.level(full, k, diag) ?? full;
        break;
      }
      if (geometry !== mesh.geometry) {
        if (!this.swapped.has(mesh)) this.swapped.set(mesh, full);
        mesh.geometry = geometry;
      }
    }
  }

  /** Put the full meshes back and free the levels. */
  restore() {
    for (const [mesh, full] of this.swapped) mesh.geometry = full;
    this.swapped.clear();
    for (const byLevel of this.levels.values()) for (const g of byLevel.values()) g?.dispose();
    this.levels.clear();
  }
}
