// Screen-space detail for film export (web/live/film-lod.js): the
// clustered level of a dense mesh stays on its surface within a cell,
// keeps its attributes, and does not merge faces across a crease; the
// clustered edge lines stay on their curve, drop collapsed and repeated
// segments, and keep a fine polyline whole.
//   node web/live/test/film-lod.test.mjs  -> prints {"checks": n}
import assert from "node:assert/strict";
import * as THREE from "../vendor/three.module.min.js";
import { clusterGeometry, clusterLines } from "../film-lod.js";

let checks = 0;
const check = (name, fn) => {
  fn();
  checks += 1;
  void name;
};

// A unit sphere, 64 x 64 segments (8k triangles), with vertex colours.
const sphere = new THREE.SphereGeometry(1, 64, 64);
const n = sphere.attributes.position.count;
sphere.setAttribute("color", new THREE.BufferAttribute(new Float32Array(3 * n).fill(0.5), 3));

check("a coarse cell collapses the mesh", () => {
  const lod = clusterGeometry(sphere, 0.25);
  assert.ok(lod.index.count < sphere.index.count / 4, `${lod.index.count} of ${sphere.index.count}`);
  assert.ok(lod.index.count > 0);
});

check("every clustered vertex stays within a cell of the surface", () => {
  const cell = 0.1;
  const lod = clusterGeometry(sphere, cell);
  const p = lod.attributes.position.array;
  for (let i = 0; i < p.length; i += 3) {
    const r = Math.hypot(p[i], p[i + 1], p[i + 2]);
    assert.ok(Math.abs(r - 1) < cell * Math.sqrt(3), `vertex at radius ${r}`);
  }
});

check("normals are unit and colours carried", () => {
  const lod = clusterGeometry(sphere, 0.1);
  const q = lod.attributes.normal.array;
  for (let i = 0; i < q.length; i += 3) assert.ok(Math.abs(Math.hypot(q[i], q[i + 1], q[i + 2]) - 1) < 1e-5);
  const c = lod.attributes.color.array;
  assert.equal(c.length, lod.attributes.position.count * 3);
  assert.ok(c.every((x) => Math.abs(x - 0.5) < 1e-6));
});

check("a fine cell keeps the triangles", () => {
  const lod = clusterGeometry(sphere, 1e-4);
  // Only the poles' fan triangles (already degenerate in chart) go.
  assert.ok(lod.index.count >= 0.95 * sphere.index.count, `${lod.index.count} of ${sphere.index.count}`);
});

check("a crease survives: a cube's faces stay apart", () => {
  const cube = new THREE.BoxGeometry(1, 1, 1, 8, 8, 8);
  const lod = clusterGeometry(cube, 2.0);
  // Cells far larger than the cube: still every cluster belongs to one
  // face (its averaged normal stays axis-aligned), so no face's vertices
  // mix with a neighbour's.
  const q = lod.attributes.normal.array;
  assert.ok(q.length > 0);
  for (let i = 0; i < q.length; i += 3) {
    assert.ok(Math.max(Math.abs(q[i]), Math.abs(q[i + 1]), Math.abs(q[i + 2])) > 1 - 1e-6);
  }
});

// A unit circle as 4096 line segments (an edge sampled for close
// inspection), every segment drawn twice (two faces' outlines).
const ring = [];
for (let i = 0; i < 4096; i++) {
  for (const j of [i, i + 1]) {
    const a = (j / 4096) * 2 * Math.PI;
    ring.push(Math.cos(a), Math.sin(a), 0);
  }
}
const circle = new THREE.BufferGeometry();
circle.setAttribute("position", new THREE.BufferAttribute(new Float32Array([...ring, ...ring]), 3));

check("a coarse cell collapses an edge and drops repeats", () => {
  const lod = clusterLines(circle, 0.1);
  const count = lod.attributes.position.count;
  assert.ok(count > 0 && count % 2 === 0);
  // About 2π/0.1 cells along the curve: far under one 4096-segment copy.
  assert.ok(count / 2 < 200, `${count / 2} segments`);
  const p = lod.attributes.position.array;
  const keys = new Set();
  for (let i = 0; i < p.length; i += 6) {
    const a = `${p[i]},${p[i + 1]},${p[i + 2]}`;
    const b = `${p[i + 3]},${p[i + 4]},${p[i + 5]}`;
    assert.notEqual(a, b, "a collapsed segment survived");
    const key = a < b ? `${a}|${b}` : `${b}|${a}`;
    assert.ok(!keys.has(key), "a repeated segment survived");
    keys.add(key);
  }
});

check("clustered edge ends stay within a cell of the curve", () => {
  const cell = 0.05;
  const p = clusterLines(circle, cell).attributes.position.array;
  for (let i = 0; i < p.length; i += 3) {
    const r = Math.hypot(p[i], p[i + 1], p[i + 2]);
    assert.ok(Math.abs(r - 1) < cell * Math.sqrt(3), `end at radius ${r}`);
  }
});

check("a fine cell keeps the polyline (one copy of it)", () => {
  const lod = clusterLines(circle, 1e-5);
  assert.equal(lod.attributes.position.count / 2, 4096);
});

console.log(JSON.stringify({ checks }));
