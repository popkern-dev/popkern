// Face-boundary edges from a body's display buffers: the mesh edges
// whose two triangles belong to different BRep faces (or that only one
// triangle uses). These are the model's topological edges, tangent
// fillet seams included, where a dihedral-angle crease filter draws
// only the sharp ones. Segments are grouped by the face pair they
// separate, so a hover or pick can name the edge by its two faces.

/** Weld vertices that coincide to within `tolerance` (faces are
 *  tessellated separately, so a shared boundary may repeat vertices)
 *  and answer each vertex's welded id. */
function weldVertices(positions, tolerance) {
  const count = positions.length / 3;
  const ids = new Uint32Array(count);
  const cells = new Map();
  let next = 0;
  for (let v = 0; v < count; v++) {
    const key = `${Math.round(positions[v * 3] / tolerance)},${Math.round(positions[v * 3 + 1] / tolerance)},${Math.round(positions[v * 3 + 2] / tolerance)}`;
    let id = cells.get(key);
    if (id === undefined) {
      id = next++;
      cells.set(key, id);
    }
    ids[v] = id;
  }
  return ids;
}

/** Boundary segments of one body's buffers `b` (positions, indices,
 *  face_ids): `{ positions, segmentEdge, edges }`, where `positions`
 *  holds two points per segment, `segmentEdge[s]` is the index into
 *  `edges` of segment `s`, and each edge is `{ faces: [a, b] }` (b is
 *  `null` on an open boundary). */
export function boundaryEdges(b) {
  const empty = { positions: new Float32Array(0), segmentEdge: new Uint32Array(0), edges: [] };
  if (!b.indices.length) return empty;
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < b.positions.length; i++) {
    lo[i % 3] = Math.min(lo[i % 3], b.positions[i]);
    hi[i % 3] = Math.max(hi[i % 3], b.positions[i]);
  }
  // f32 buffers: a weld well above single-precision noise on the
  // model's scale, well below any feature the kernel would build.
  const diagonal = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
  const weld = weldVertices(b.positions, diagonal * 1e-6);
  let welded = 0;
  for (const id of weld) welded = Math.max(welded, id + 1);
  // Welded edge key -> one of its mesh edges, the first face that uses
  // it, a different face that also does (-1 if none), and its use count.
  const uses = new Map();
  const triangles = b.indices.length / 3;
  for (let t = 0; t < triangles; t++) {
    const face = b.face_ids[t];
    for (let k = 0; k < 3; k++) {
      const p = b.indices[t * 3 + k];
      const q = b.indices[t * 3 + ((k + 1) % 3)];
      const a = weld[p];
      const c = weld[q];
      if (a === c) continue;
      const key = a < c ? a * welded + c : c * welded + a;
      const use = uses.get(key);
      if (!use) uses.set(key, { p, q, face, other: -1, count: 1 });
      else {
        use.count++;
        if (face !== use.face) use.other = face;
      }
    }
  }
  const edgeOf = new Map();
  const edges = [];
  const points = [];
  const segmentEdge = [];
  for (const use of uses.values()) {
    const open = use.count === 1;
    if (!open && use.other < 0) continue;
    const a = Math.min(use.face, use.other < 0 ? use.face : use.other);
    const c = open ? null : Math.max(use.face, use.other);
    const key = `${a}:${c}`;
    let edge = edgeOf.get(key);
    if (edge === undefined) {
      edge = edges.length;
      edgeOf.set(key, edge);
      edges.push({ faces: [a, c] });
    }
    for (const v of [use.p, use.q]) points.push(b.positions[v * 3], b.positions[v * 3 + 1], b.positions[v * 3 + 2]);
    segmentEdge.push(edge);
  }
  return { positions: new Float32Array(points), segmentEdge: Uint32Array.from(segmentEdge), edges };
}

// Fat lines: each segment is a quad widened in screen space, since
// WebGL draws GL lines one pixel wide whatever `linewidth` says. The
// ribbon is pulled slightly towards the camera so the faces it bounds
// do not cut it in half.
const FAT_VERTEX = `
  attribute vec3 other;
  attribute float side;
  uniform vec2 viewport;
  uniform float width;
  void main() {
    vec4 a = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    vec4 b = projectionMatrix * modelViewMatrix * vec4(other, 1.0);
    vec2 d = b.xy / b.w * viewport - a.xy / a.w * viewport;
    vec2 dir = length(d) > 1e-6 ? normalize(d) : vec2(1.0, 0.0);
    vec2 normal = vec2(-dir.y, dir.x);
    a.xy += normal * side * width / viewport * a.w;
    a.z -= 1e-3 * a.w;
    gl_Position = a;
  }
`;
const FAT_FRAGMENT = `
  uniform vec3 color;
  void main() { gl_FragColor = vec4(color, 1.0); }
`;

/** A mesh drawing `segments` (two points per segment) `widthPx`
 *  drawing-buffer pixels wide; set `material.uniforms.viewport` to the
 *  drawing-buffer size before each render. */
export function fatLines(THREE, segments, color, widthPx) {
  const count = segments.length / 6;
  const position = new Float32Array(count * 12);
  const other = new Float32Array(count * 12);
  const side = new Float32Array(count * 4);
  const index = new Uint32Array(count * 6);
  for (let s = 0; s < count; s++) {
    const p = segments.subarray(s * 6, s * 6 + 3);
    const q = segments.subarray(s * 6 + 3, s * 6 + 6);
    // At q the direction reverses, so the side signs swap to keep the
    // quad's corners on the same two sides.
    const corners = [[p, q, 1], [p, q, -1], [q, p, -1], [q, p, 1]];
    corners.forEach(([at, to, sign], c) => {
      position.set(at, (s * 4 + c) * 3);
      other.set(to, (s * 4 + c) * 3);
      side[s * 4 + c] = sign;
    });
    index.set([s * 4, s * 4 + 1, s * 4 + 3, s * 4, s * 4 + 3, s * 4 + 2], s * 6);
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(position, 3));
  geometry.setAttribute("other", new THREE.BufferAttribute(other, 3));
  geometry.setAttribute("side", new THREE.BufferAttribute(side, 1));
  geometry.setIndex(new THREE.BufferAttribute(index, 1));
  const material = new THREE.ShaderMaterial({
    vertexShader: FAT_VERTEX,
    fragmentShader: FAT_FRAGMENT,
    uniforms: {
      viewport: { value: new THREE.Vector2(1, 1) },
      width: { value: widthPx },
      color: { value: new THREE.Color(color) },
    },
    side: THREE.DoubleSide,
  });
  const mesh = new THREE.Mesh(geometry, material);
  // The vertex shader moves the corners off the centre line the
  // bounding sphere is built from.
  mesh.frustumCulled = false;
  return mesh;
}
