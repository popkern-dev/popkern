// Pure geometry for the clipping view (docs/architecture/clipping-view.md):
// plane presets, the combine rule, the cut outline (triangle × plane
// segments), and the page-link encoding. No three.js, so the node test
// (test/section.test.mjs) exercises it directly.
//
// A plane is `{ normal: [x, y, z], offset }`: the points with
// normal·p = offset, and it REMOVES the side its normal points to
// (normal·p > offset), which is three.js's clipping convention.

/** At most this many planes at once (UI, link, and stencil cost). */
export const MAX_PLANES = 3;

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

/** `v` scaled to unit length, or null for a zero / non-finite vector. */
export function unit(v) {
  const n = Math.hypot(v[0], v[1], v[2]);
  if (!Number.isFinite(n) || n < 1e-12) return null;
  return [v[0] / n, v[1] / n, v[2] / n];
}

/** The plane through `point` with `normal` (removing the side it points to). */
export function planeThrough(point, normal) {
  const n = unit(normal);
  if (!n) return null;
  return { normal: n, offset: dot(n, point) };
}

/** The same cut, removing the other side. */
export function flipped(plane) {
  return { normal: plane.normal.map((c) => -c), offset: -plane.offset };
}

/** Signed distance of `p` past `plane` (> 0: on the removed side). */
export function past(plane, p) {
  return dot(plane.normal, p) - plane.offset;
}

/** Whether the combined planes keep point `p`: with "any" a point is
 *  removed beyond ANY plane; with "all" only where EVERY plane removes
 *  it (the corner cut-away). `skip` leaves one plane out (a cap or an
 *  outline on that plane is judged by the others). `eps` keeps points
 *  on a plane. */
export function keeps(planes, combine, p, skip = -1, eps = 0) {
  let removedAny = false;
  let removedAll = true;
  let counted = 0;
  planes.forEach((plane, i) => {
    if (i === skip) return;
    counted++;
    const removed = past(plane, p) > eps;
    removedAny ||= removed;
    removedAll &&= removed;
  });
  if (counted === 0) return true;
  return combine === "all" ? !removedAll : !removedAny;
}

/** Axis presets, as before: normal +axis at `position` (0..1 across
 *  the bounds), removing the high side and keeping the low. */
export function axisPlane(axis, bounds, position = 0.5) {
  const i = "xyz".indexOf(axis);
  if (i < 0 || !bounds) return null;
  const [lo, hi] = bounds;
  const normal = [0, 0, 0];
  normal[i] = 1;
  return { normal, offset: lo[i] + (hi[i] - lo[i]) * position };
}

/** Facing the viewer: removes the half nearer the camera, through `target`. */
export function viewPlane(cameraForward, target) {
  const n = unit(cameraForward);
  if (!n) return null;
  return planeThrough(target, n.map((c) => -c));
}

/** The 0..1 position of `plane` along its normal across `bounds` — the
 *  slider's value: 0 where the plane removes everything (at the bounds'
 *  lowest corner along the normal), 1 where it removes nothing. */
export function sliderPosition(plane, bounds) {
  const [lo, hi] = extentAlong(plane.normal, bounds);
  if (hi - lo < 1e-12) return 0.5;
  return Math.min(1, Math.max(0, (plane.offset - lo) / (hi - lo)));
}

/** Offset of a plane with `normal` at slider position `t` (inverse of
 *  `sliderPosition`). */
export function offsetAt(normal, bounds, t) {
  const [lo, hi] = extentAlong(normal, bounds);
  return lo + t * (hi - lo);
}

/** Whether `plane` passes through the box `bounds` (strictly: a box
 *  wholly on one side needs no cap or outline). */
export function straddles(plane, bounds) {
  const [lo, hi] = extentAlong(plane.normal, bounds);
  return lo < plane.offset && plane.offset < hi;
}

/** Range of normal·p over the corners of `bounds`. */
export function extentAlong(normal, bounds) {
  const [lo, hi] = bounds;
  let min = Infinity;
  let max = -Infinity;
  for (let k = 0; k < 8; k++) {
    const p = [k & 1 ? hi[0] : lo[0], k & 2 ? hi[1] : lo[1], k & 4 ? hi[2] : lo[2]];
    const d = dot(normal, p);
    min = Math.min(min, d);
    max = Math.max(max, d);
  }
  return [min, max];
}

/** The cut outline of one triangle mesh on `planes[index]`: a flat
 *  `Float32Array` of segment endpoints (x0 y0 z0 x1 y1 z1 …) in WORLD
 *  space. `positions` / `indices` are the mesh buffers, `matrix` a
 *  column-major 4×4 world matrix (or null for identity). Segments the
 *  other planes remove (per `combine`) are dropped, judged at their
 *  midpoints. A triangle with a vertex exactly on the plane contributes
 *  the segment through that vertex once per crossing edge pair, so a
 *  closed mesh's outline is a closed polyline. */
export function cutOutline(positions, indices, matrix, planes, index, combine = "any") {
  const plane = planes[index];
  const n = plane.normal;
  const vertexCount = positions.length / 3;
  const m = matrix;
  // The plane in the mesh's own frame (world = A·local + t, so
  // n·(A l + t) = d is (Aᵀn)·l = d − n·t): distances need no vertex
  // transform, and only crossing points are taken to world space.
  const ln = m ? [n[0] * m[0] + n[1] * m[1] + n[2] * m[2], n[0] * m[4] + n[1] * m[5] + n[2] * m[6], n[0] * m[8] + n[1] * m[9] + n[2] * m[10]] : n;
  const ld = m ? plane.offset - (n[0] * m[12] + n[1] * m[13] + n[2] * m[14]) : plane.offset;
  const dist = new Float64Array(vertexCount);
  for (let v = 0; v < vertexCount; v++) {
    dist[v] = ln[0] * positions[3 * v] + ln[1] * positions[3 * v + 1] + ln[2] * positions[3 * v + 2] - ld;
  }
  const out = [];
  const point = (a, b) => {
    const t = dist[a] / (dist[a] - dist[b]);
    const x = positions[3 * a] + t * (positions[3 * b] - positions[3 * a]);
    const y = positions[3 * a + 1] + t * (positions[3 * b + 1] - positions[3 * a + 1]);
    const z = positions[3 * a + 2] + t * (positions[3 * b + 2] - positions[3 * a + 2]);
    if (!m) return [x, y, z];
    return [m[0] * x + m[4] * y + m[8] * z + m[12], m[1] * x + m[5] * y + m[9] * z + m[13], m[2] * x + m[6] * y + m[10] * z + m[14]];
  };
  const others = planes.length > 1;
  for (let t = 0; t < indices.length; t += 3) {
    const tri = [indices[t], indices[t + 1], indices[t + 2]];
    // Classify with a half-open rule (d ≥ 0 counts as "beyond") so a
    // vertex on the plane belongs to exactly one side: each crossing
    // edge has one endpoint on each side, and a triangle crosses on
    // exactly zero or two edges.
    const hits = [];
    for (let k = 0; k < 3; k++) {
      const a = tri[k];
      const b = tri[(k + 1) % 3];
      if (dist[a] >= 0 !== dist[b] >= 0) hits.push(point(a, b));
    }
    if (hits.length !== 2) continue;
    const [p, q] = hits;
    if (Math.abs(p[0] - q[0]) + Math.abs(p[1] - q[1]) + Math.abs(p[2] - q[2]) < 1e-15) continue;
    if (others) {
      const mid = [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2, (p[2] + q[2]) / 2];
      if (!keepsCut(planes, combine, mid, index)) continue;
    }
    out.push(p[0], p[1], p[2], q[0], q[1], q[2]);
  }
  return Float32Array.from(out);
}

/** Whether a point ON `planes[index]` is on a drawn cap: "any" — no
 *  other plane removes it; "all" — every other plane removes it (cap i
 *  bounds the removed corner). */
export function keepsCut(planes, combine, p, index) {
  for (let i = 0; i < planes.length; i++) {
    if (i === index) continue;
    const removed = past(planes[i], p) > 0;
    if (combine === "all" ? !removed : removed) return false;
  }
  return true;
}

const fmt = (x) => {
  const s = Number(x.toPrecision(6));
  return Object.is(s, -0) ? "0" : String(s);
};

/** `#section=` text for a cut: `nx,ny,nz,d;…` then `~all` when the
 *  combine mode is "all", then `~x:` and the excluded occurrence paths. */
export function toLink({ planes, combine = "any", exclude = [] }) {
  if (!planes.length) return "";
  let text = planes.map((p) => [...p.normal, p.offset].map(fmt).join(",")).join(";");
  if (combine === "all" && planes.length > 1) text += "~all";
  if (exclude.length) text += "~x:" + exclude.map(encodeURIComponent).join(",");
  return text;
}

/** Inverse of `toLink`; malformed planes are dropped, at most
 *  `MAX_PLANES` kept. */
export function fromLink(text) {
  const out = { planes: [], combine: "any", exclude: [] };
  if (!text) return out;
  const [planeText, ...flags] = String(text).split("~");
  for (const part of planeText.split(";")) {
    const values = part.split(",").map(Number);
    if (values.length !== 4 || !values.every(Number.isFinite)) continue;
    const normal = unit(values.slice(0, 3));
    if (!normal) continue;
    // Renormalising rescales the offset with the normal.
    const scale = Math.hypot(values[0], values[1], values[2]);
    out.planes.push({ normal, offset: values[3] / scale });
    if (out.planes.length === MAX_PLANES) break;
  }
  for (const flag of flags) {
    if (flag === "all") out.combine = "all";
    else if (flag.startsWith("x:") && flag.length > 2) out.exclude = flag.slice(2).split(",").map(decodeURIComponent);
  }
  return out;
}

/** Rotate unit vector `v` about unit `axis` by `angle` (Rodrigues). */
export function rotated(v, axis, angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  const k = axis;
  const d = dot(k, v);
  const cross = [k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]];
  return unit([
    v[0] * c + cross[0] * s + k[0] * d * (1 - c),
    v[1] * c + cross[1] * s + k[1] * d * (1 - c),
    v[2] * c + cross[2] * s + k[2] * d * (1 - c),
  ]);
}

/** Two unit vectors spanning the plane with `normal` (a right-handed
 *  frame u, v, normal), stable for the axis presets. */
export function inPlaneAxes(normal) {
  const n = normal;
  const seed = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  const u = unit([seed[1] * n[2] - seed[2] * n[1], seed[2] * n[0] - seed[0] * n[2], seed[0] * n[1] - seed[1] * n[0]]);
  const v = [n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]];
  return [u, v];
}
