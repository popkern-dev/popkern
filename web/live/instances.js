// Assembly occurrences in the viewer: what to draw per occurrence, how
// to fix an object at a world matrix, and the exploded-view directions.
import * as THREE from "./vendor/three.module.min.js";

/** What to draw: one entry per body for a part document; one per body
 *  of every part occurrence for an assembly, with the occurrence's
 *  world matrix (16 column-major numbers), path and visibility. */
export function drawList(bodies, assembly) {
  if (!assembly) return bodies.map((_, index) => ({ index, matrix: null, instance: null, visible: true, colorKey: index }));
  const draws = [];
  for (const instance of assembly.instances) {
    if (!instance.is_part) continue;
    const definition = assembly.definitions[instance.definition];
    for (const index of definition?.bodies ?? []) {
      draws.push({ index, matrix: instance.world, instance: instance.path, visible: instance.visible, colorKey: instance.definition });
    }
  }
  return draws;
}

/** Fix `object` at the column-major world matrix `m`. */
export function placeAt(object, m) {
  object.matrixAutoUpdate = false;
  if (m.isMatrix4) object.matrix.copy(m);
  else object.matrix.fromArray(m);
  object.matrixWorldNeedsUpdate = true;
  object.updateMatrixWorld(true);
}

/** Explode directions per top-level occurrence: from the assembly's
 *  centre to the occurrence's centre, scaled to the model size; an
 *  occurrence centred on the assembly's axis (a sun, a ring, a
 *  carrier) moves along z instead, ordered by its height, so coaxial
 *  parts separate rather than stay nested. */
export function explodePlan(meshes) {
  const boxes = new Map();
  for (const mesh of meshes) {
    const top = mesh.userData.instance?.split("/")[0];
    if (top == null) continue;
    const [lo, hi] = mesh.userData.body.bounds;
    const box = new THREE.Box3(new THREE.Vector3(...lo), new THREE.Vector3(...hi)).applyMatrix4(mesh.userData.base ?? mesh.matrix);
    boxes.set(top, boxes.has(top) ? boxes.get(top).union(box) : box);
  }
  if (boxes.size < 2) return null;
  const whole = new THREE.Box3();
  for (const box of boxes.values()) whole.union(box);
  const centre = whole.getCenter(new THREE.Vector3());
  const size = whole.getSize(new THREE.Vector3()).length();
  const directions = new Map();
  const axial = [];
  for (const [top, box] of boxes) {
    const c = box.getCenter(new THREE.Vector3());
    const radial = new THREE.Vector3(c.x - centre.x, c.y - centre.y, 0);
    if (radial.length() > 1e-6 * size) {
      directions.set(top, radial.normalize().multiplyScalar(size * 0.6).setZ((c.z - centre.z) * 0.6));
    } else {
      axial.push([top, c.z]);
    }
  }
  axial.sort((a, b) => a[1] - b[1]);
  axial.forEach(([top], k) => {
    const rank = axial.length > 1 ? k / (axial.length - 1) - 0.5 : 0;
    directions.set(top, new THREE.Vector3(0, 0, rank * size * 1.2));
  });
  return { directions, size };
}
