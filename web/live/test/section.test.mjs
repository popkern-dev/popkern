// The clipping view's pure geometry (section-geometry.js;
// docs/architecture/clipping-view.md §9), run by
// crates/popkern-live/tests/page_section.rs. Exits non-zero with the
// failing check's message; prints the check count on success.
import assert from "node:assert/strict";
import {
  MAX_PLANES,
  axisPlane,
  cutOutline,
  flipped,
  fromLink,
  keeps,
  keepsCut,
  offsetAt,
  past,
  planeThrough,
  rotated,
  sliderPosition,
  toLink,
  viewPlane,
} from "../section-geometry.js";

let checks = 0;
const check = (name, fn) => {
  try {
    fn();
    checks++;
  } catch (e) {
    console.error(`FAIL ${name}: ${e.message}`);
    process.exit(1);
  }
};
const near = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps, `${a} vs ${b}`);

// A unit cube [0,1]^3 as 12 triangles.
const cube = (() => {
  const v = [];
  for (let k = 0; k < 8; k++) v.push(k & 1, (k >> 1) & 1, (k >> 2) & 1);
  const quads = [
    [0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5],
  ];
  const idx = [];
  for (const [a, b, c, d] of quads) idx.push(a, b, c, a, c, d);
  return { positions: Float32Array.from(v), indices: Uint32Array.from(idx) };
})();
const bounds = [[0, 0, 0], [1, 1, 1]];

/** Segments as [[p, q], ...] from a flat array. */
const segs = (flat) => {
  const out = [];
  for (let i = 0; i < flat.length; i += 6) out.push([[flat[i], flat[i + 1], flat[i + 2]], [flat[i + 3], flat[i + 4], flat[i + 5]]]);
  return out;
};
const length = (s) => s.reduce((n, [p, q]) => n + Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]), 0);

check("axis preset removes the high side at the position", () => {
  const p = axisPlane("x", bounds, 0.25);
  assert.deepEqual(p.normal, [1, 0, 0]);
  near(p.offset, 0.25);
  assert.ok(past(p, [0.9, 0.5, 0.5]) > 0, "high x removed");
  assert.ok(past(p, [0.1, 0.5, 0.5]) < 0, "low x kept");
});

check("slider position and offset invert, 0 removes all and 1 nothing", () => {
  const n = [0.6, 0.8, 0];
  for (const t of [0, 0.3, 1]) near(sliderPosition({ normal: n, offset: offsetAt(n, bounds, t) }, bounds), t);
  const all = { normal: n, offset: offsetAt(n, bounds, 0) };
  const none = { normal: n, offset: offsetAt(n, bounds, 1) };
  for (let k = 0; k < 8; k++) {
    const c = [k & 1, (k >> 1) & 1, (k >> 2) & 1];
    assert.ok(past(all, c) >= -1e-12, "t=0 removes every corner");
    assert.ok(past(none, c) <= 1e-12, "t=1 keeps every corner");
  }
});

check("flip removes the other side", () => {
  const p = planeThrough([0.5, 0.5, 0.5], [0, 0, 1]);
  const q = flipped(p);
  assert.ok(past(p, [0.5, 0.5, 0.9]) > 0 && past(q, [0.5, 0.5, 0.9]) < 0);
});

check("facing the viewer removes the half nearer the camera", () => {
  // Camera on +y looking toward -y at the origin.
  const p = viewPlane([0, -1, 0], [0, 0, 0]);
  assert.ok(past(p, [0, 1, 0]) > 0, "near half removed");
  assert.ok(past(p, [0, -1, 0]) < 0, "far half kept");
});

check("combine: any removes beyond either plane, all only the corner", () => {
  const planes = [axisPlane("x", bounds, 0.5), axisPlane("y", bounds, 0.5)];
  const corner = [0.9, 0.9, 0.5];
  const edge = [0.9, 0.1, 0.5];
  const kept = [0.1, 0.1, 0.5];
  assert.equal(keeps(planes, "any", corner), false);
  assert.equal(keeps(planes, "any", edge), false);
  assert.equal(keeps(planes, "any", kept), true);
  assert.equal(keeps(planes, "all", corner), false);
  assert.equal(keeps(planes, "all", edge), true);
  assert.equal(keeps(planes, "all", kept), true);
});

check("a cube cut at any height outlines exactly its square", () => {
  for (const t of [0.13, 0.5, 0.87]) {
    const planes = [axisPlane("z", bounds, t)];
    const s = segs(cutOutline(cube.positions, cube.indices, null, planes, 0));
    near(length(s), 4, 1e-6);
    for (const [p, q] of s) {
      near(p[2], t, 1e-6);
      near(q[2], t, 1e-6);
    }
  }
});

check("a tilted cut through the cube's centre is a closed loop", () => {
  const planes = [planeThrough([0.5, 0.5, 0.5], [1, 1, 1])];
  const s = segs(cutOutline(cube.positions, cube.indices, null, planes, 0));
  // The regular hexagon of side sqrt(2)/2.
  near(length(s), 6 * Math.SQRT1_2, 1e-6);
  const ends = new Map();
  const key = (p) => p.map((c) => c.toFixed(5)).join(",");
  for (const [p, q] of s) for (const k of [key(p), key(q)]) ends.set(k, (ends.get(k) ?? 0) + 1);
  for (const [k, n] of ends) assert.equal(n % 2, 0, `endpoint ${k} used ${n} times`);
});

check("a cut through mesh vertices and edges is the closed diagonal rectangle", () => {
  // z = x through the edges x=z=0 and x=z=1 (every vertex there lies
  // exactly on the plane).
  const planes = [planeThrough([0, 0, 0], [1, 0, -1])];
  const s = segs(cutOutline(cube.positions, cube.indices, null, planes, 0));
  near(length(s), 2 + 2 * Math.SQRT2, 1e-6);
  const ends = new Map();
  const key = (p) => p.map((c) => c.toFixed(5)).join(",");
  for (const [p, q] of s) for (const k of [key(p), key(q)]) ends.set(k, (ends.get(k) ?? 0) + 1);
  for (const [k, n] of ends) assert.equal(n % 2, 0, `endpoint ${k} used ${n} times`);
});

check("the world matrix moves the outline", () => {
  // Column-major translation by (10, 0, 0).
  const m = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 10, 0, 0, 1];
  const planes = [planeThrough([10.5, 0, 0], [1, 0, 0])];
  const s = segs(cutOutline(cube.positions, cube.indices, m, planes, 0));
  near(length(s), 4, 1e-6);
  for (const [p] of s) near(p[0], 10.5, 1e-6);
});

check("other planes trim the outline: any keeps what none removes, all the corner face", () => {
  const planes = [axisPlane("z", bounds, 0.5), axisPlane("x", bounds, 0.5)];
  const any = segs(cutOutline(cube.positions, cube.indices, null, planes, 0, "any"));
  near(length(any), 2, 1e-6); // x <= 0.5 half of the z square: 0.5 + 1 + 0.5
  const all = segs(cutOutline(cube.positions, cube.indices, null, planes, 0, "all"));
  near(length(all), 2, 1e-6); // the x >= 0.5 half bounds the removed corner
  for (const [p, q] of all) assert.ok((p[0] + q[0]) / 2 >= 0.5 - 1e-9);
  assert.equal(keepsCut(planes, "all", [0.9, 0.5, 0.5], 0), true);
  assert.equal(keepsCut(planes, "any", [0.9, 0.5, 0.5], 0), false);
});

check("a plane missing the mesh outlines nothing", () => {
  assert.equal(cutOutline(cube.positions, cube.indices, null, [axisPlane("z", [[0, 0, 2], [1, 1, 3]], 0.5)], 0).length, 0);
});

check("link round trip keeps planes, combine and exclusions", () => {
  const state = {
    planes: [planeThrough([0.5, 0.25, 0], [0, 0.6, 0.8]), axisPlane("x", bounds, 0.3)],
    combine: "all",
    exclude: ["housing_1/screw 2", "a,b"],
  };
  const back = fromLink(toLink(state));
  assert.equal(back.combine, "all");
  assert.deepEqual(back.exclude, state.exclude);
  assert.equal(back.planes.length, 2);
  back.planes.forEach((p, i) => {
    p.normal.forEach((c, k) => near(c, state.planes[i].normal[k], 1e-6));
    near(p.offset, state.planes[i].offset, 1e-6);
  });
});

check("link: malformed planes dropped, at most three, normals renormalised", () => {
  const text = "0,0,2,4;1,2;0,0,0,1;nan,1,0,0;1,0,0,0;0,1,0,0;0,0,1,0";
  const back = fromLink(text);
  assert.equal(back.planes.length, MAX_PLANES);
  assert.deepEqual(back.planes[0].normal, [0, 0, 1]);
  near(back.planes[0].offset, 2);
  assert.equal(toLink({ planes: [] }), "");
  assert.deepEqual(fromLink("").planes, []);
});

check("rotation about an in-plane axis keeps unit length", () => {
  const n = rotated([0, 0, 1], [1, 0, 0], Math.PI / 2);
  near(Math.hypot(...n), 1);
  near(n[1], -1, 1e-12);
});

console.log(JSON.stringify({ checks }));
