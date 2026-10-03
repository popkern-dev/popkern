// The three.js viewer: one mesh per body from the engine's buffers,
// per-face colouring and picking through the face-id attribute, the
// faces' boundary edges, hover pre-highlighting in the CAD view, a
// small orbit control (no addon), fit-to-bounds on the first render.
import * as THREE from "./vendor/three.module.min.js";
import { boundaryEdges, fatLines } from "./edges.js";
import { Orbit } from "./orbit.js";

const PALETTE = [0x6c8ebf, 0x8fb56c, 0xd9a25f, 0xb87ab8, 0x66b3b3, 0xc97a7a, 0x9c9c6c, 0x7a8fc9];

function faceColor(index) {
  return new THREE.Color(PALETTE[index % PALETTE.length]);
}

// The CAD view: uniform grey faces, black edges, and a blue
// pre-highlight on the face or edge under the pointer.
const CAD_FACE = 0xb4b8be;
const CAD_EDGE = 0x000000;
const HOVER = 0x3b82f6;
// Pointer slack for edge hover and pick, and the hovered edge's
// width, in CSS pixels.
const EDGE_PICK_PX = 5;
const HOVER_EDGE_PX = 3;
// Depth offsets so edges (offset 0) draw over a hovered face (1),
// which draws over the face itself (2).
const FACE_OFFSET = 2;
const HOVER_OFFSET = 1;

// Boundary edges per body buffer, computed once per served mesh.
const edgeCache = new WeakMap();
function cachedEdges(b) {
  let edges = edgeCache.get(b);
  if (!edges) {
    edges = boundaryEdges(b);
    edgeCache.set(b, edges);
  }
  return edges;
}

// Reflection-line zebra: per-pixel stripes over the reflected view
// direction — the class-A continuity inspection (a G1 defect kinks
// the stripes across a seam; a G2 defect kinks their derivative).
const ZEBRA_VERTEX = `
  varying vec3 vWorldPos;
  varying vec3 vWorldNormal;
  void main() {
    vec4 world = modelMatrix * vec4(position, 1.0);
    vWorldPos = world.xyz;
    vWorldNormal = normalize(mat3(modelMatrix) * normal);
    gl_Position = projectionMatrix * viewMatrix * world;
  }
`;
const ZEBRA_FRAGMENT = `
  varying vec3 vWorldPos;
  varying vec3 vWorldNormal;
  uniform vec3 stripeAxis;
  uniform float stripeCount;
  void main() {
    vec3 n = normalize(vWorldNormal);
    vec3 view = normalize(vWorldPos - cameraPosition);
    vec3 r = reflect(view, n);
    float phase = (dot(r, stripeAxis) * 0.5 + 0.5) * stripeCount;
    float f = fract(phase);
    float aa = fwidth(phase);
    float band = smoothstep(0.25 - aa, 0.25 + aa, f) - smoothstep(0.75 - aa, 0.75 + aa, f);
    float shade = mix(0.09, 0.93, band);
    gl_FragColor = vec4(vec3(shade), 1.0);
  }
`;

// Diverging blue-white-red map over a symmetric range for curvature.
function divergingColor(value, halfRange) {
  const t = halfRange > 0 ? Math.max(-1, Math.min(1, value / halfRange)) : 0;
  const white = [0.95, 0.95, 0.95];
  const cold = [0.23, 0.30, 0.75];
  const warm = [0.71, 0.02, 0.15];
  const to = t < 0 ? cold : warm;
  const f = Math.abs(t);
  return [white[0] + (to[0] - white[0]) * f, white[1] + (to[1] - white[1]) * f, white[2] + (to[2] - white[2]) * f];
}

// Robust symmetric half-range: 1.25x the median NONZERO magnitude
// across all bodies — flat regions would drag a percentile to zero
// and degenerate spikes would inflate a max-based range; the nonzero
// median saturates the curved features (mirrors popkern-viz).
function robustHalfRange(buffers, key) {
  const mags = [];
  for (const b of buffers) {
    const field = b[key];
    if (!field) continue;
    for (let i = 0; i < field.length; i++) {
      const m = Math.abs(field[i]);
      if (Number.isFinite(m) && m > 1e-12) mags.push(m);
    }
  }
  if (!mags.length) return 0;
  mags.sort((a, b) => a - b);
  return 1.25 * mags[Math.floor(mags.length / 2)];
}

export class Viewer {
  constructor(canvas, onPick) {
    this.canvas = canvas;
    this.onPick = onPick;
    this.mode = "cad";
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x1d2026);
    this.camera = new THREE.PerspectiveCamera(40, 1, 0.01, 1000);
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x334455, 0.9));
    const key = new THREE.DirectionalLight(0xffffff, 0.8);
    key.position.set(3, -4, 6);
    this.scene.add(key);
    this.group = new THREE.Group();
    this.scene.add(this.group);
    this.bodies = [];
    this.edgeLines = [];
    this.highlight = null;
    this.hovered = null;
    this.showEdges = true;
    this.fitted = false;
    this.raycaster = new THREE.Raycaster();
    this.orbit = new Orbit(this.camera, canvas, () => this.render());
    canvas.tabIndex = 0;
    canvas.addEventListener("click", (e) => {
      if (!this.orbit.wasDrag()) this.pick(e);
    });
    // Hover pre-highlight (CAD view, mouse only, not mid-drag), one hit
    // test per frame however fast the pointer moves.
    canvas.addEventListener("pointermove", (e) => {
      if (this.mode !== "cad" || e.pointerType !== "mouse" || e.buttons) return;
      this.pendingHover = { x: e.clientX, y: e.clientY };
      if (this.hoverFrame) return;
      this.hoverFrame = requestAnimationFrame(() => {
        this.hoverFrame = 0;
        const at = this.pendingHover;
        this.setHover(at ? this.hitTest(at.x, at.y) : null);
      });
    });
    canvas.addEventListener("pointerleave", () => {
      this.pendingHover = null;
      this.setHover(null);
    });
    canvas.addEventListener("pointerdown", () => this.setHover(null));
    this.resize();
    // The container, not the window, sets the canvas size: the narrow
    // layout resizes it on rotation and when mobile toolbars collapse.
    if (typeof ResizeObserver === "function") new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
    else window.addEventListener("resize", () => this.resize());
    this.orbit.update();
  }

  resize() {
    const rect = this.canvas.parentElement.getBoundingClientRect();
    const width = Math.max(1, Math.floor(rect.width));
    const height = Math.max(1, Math.floor(rect.height));
    this.renderer.setSize(width, height, false);
    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.render();
  }

  /** Zebra analysis buffers for the shown model (PK-577), or `null`. */
  setAnalysis(buffers) {
    this.analysis = buffers;
  }

  /** Replace the scene with the envelope's bodies and their buffers.
   *  The zebra view draws the analysis mesh when one is loaded. */
  show(bodies, displayBuffers) {
    this.lastShown = { bodies, buffers: displayBuffers };
    const buffers =
      this.mode === "zebra" && this.analysis?.length === bodies.length ? this.analysis : displayBuffers;
    this.group.clear();
    this.bodies = [];
    this.edgeLines = [];
    this.highlight = null;
    this.hovered = null;
    this.zebraMaterials = [];
    let union = null;
    const halfRange = this.mode === "gaussian" || this.mode === "mean" ? robustHalfRange(buffers, this.mode) : 0;
    bodies.forEach((body, i) => {
      const b = buffers[i];
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(b.positions, 3));
      geometry.setAttribute("normal", new THREE.BufferAttribute(b.normals, 3));
      geometry.setIndex(new THREE.BufferAttribute(b.indices, 1));
      geometry.setAttribute("color", this.vertexColors(b, halfRange));
      const mesh = new THREE.Mesh(geometry, this.material());
      mesh.userData = { body, faceIds: b.face_ids, index: i };
      this.group.add(mesh);
      if (this.showEdges) {
        const edges = cachedEdges(b);
        const lineGeometry = new THREE.BufferGeometry();
        lineGeometry.setAttribute("position", new THREE.BufferAttribute(edges.positions, 3));
        const color = this.mode === "cad" ? CAD_EDGE : 0x0f1115;
        const lines = new THREE.LineSegments(lineGeometry, new THREE.LineBasicMaterial({ color }));
        lines.userData = { edges, index: i };
        this.group.add(lines);
        this.edgeLines.push(lines);
      }
      this.bodies.push(mesh);
      const [lo, hi] = body.bounds;
      union = union ? [lo.map((v, k) => Math.min(v, union[0][k])), hi.map((v, k) => Math.max(v, union[1][k]))] : [lo, hi];
    });
    if (union && !this.fitted) {
      this.orbit.fit(union);
      this.fitted = true;
    }
    this.render();
  }

  /** Per-vertex colour attribute for the current mode. */
  vertexColors(b, halfRange) {
    const colors = new Float32Array(b.positions.length);
    if (this.mode === "gaussian" || this.mode === "mean") {
      const field = b[this.mode] ?? new Float32Array(b.positions.length / 3);
      for (let v = 0; v < field.length; v++) {
        const [r, g, bl] = divergingColor(field[v], halfRange);
        colors[v * 3] = r;
        colors[v * 3 + 1] = g;
        colors[v * 3 + 2] = bl;
      }
    } else {
      // Face colour from the first triangle that uses each vertex.
      const seen = new Uint8Array(b.positions.length / 3);
      for (let t = 0; t < b.face_ids.length; t++) {
        const color = faceColor(b.face_ids[t]);
        for (let k = 0; k < 3; k++) {
          const v = b.indices[t * 3 + k];
          if (seen[v]) continue;
          seen[v] = 1;
          colors[v * 3] = color.r;
          colors[v * 3 + 1] = color.g;
          colors[v * 3 + 2] = color.b;
        }
      }
    }
    return new THREE.BufferAttribute(colors, 3);
  }

  material() {
    if (this.mode === "zebra") {
      const material = new THREE.ShaderMaterial({
        vertexShader: ZEBRA_VERTEX,
        fragmentShader: ZEBRA_FRAGMENT,
        uniforms: {
          stripeAxis: { value: new THREE.Vector3(1, 0, 0) },
          stripeCount: { value: 14.0 },
        },
        side: THREE.DoubleSide,
      });
      this.zebraMaterials.push(material);
      return material;
    }
    const offset = { polygonOffset: true, polygonOffsetFactor: FACE_OFFSET, polygonOffsetUnits: FACE_OFFSET };
    if (this.mode === "gaussian" || this.mode === "mean") {
      return new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide, ...offset });
    }
    if (this.mode === "cad") {
      return new THREE.MeshStandardMaterial({ color: CAD_FACE, metalness: 0.05, roughness: 0.55, side: THREE.DoubleSide, ...offset });
    }
    return new THREE.MeshStandardMaterial({ vertexColors: true, metalness: 0.1, roughness: 0.65, side: THREE.DoubleSide, ...offset });
  }

  /** Switch the surface view: cad | shaded | zebra | gaussian | mean. */
  setMode(mode) {
    this.mode = mode;
    if (this.lastShown) this.show(this.lastShown.bodies, this.lastShown.buffers);
  }

  toggleEdges(on) {
    this.showEdges = on;
    if (this.lastShown) this.show(this.lastShown.bodies, this.lastShown.buffers);
  }

  fitAgain() {
    this.fitted = false;
    const bodies = this.bodies.map((m) => m.userData.body);
    if (bodies.length) {
      const lo = bodies[0].bounds[0].slice();
      const hi = bodies[0].bounds[1].slice();
      for (const b of bodies) {
        for (let k = 0; k < 3; k++) {
          lo[k] = Math.min(lo[k], b.bounds[0][k]);
          hi[k] = Math.max(hi[k], b.bounds[1][k]);
        }
      }
      this.orbit.fit([lo, hi]);
      this.fitted = true;
    }
  }

  /** What lies under the client point: `{ kind: "face", mesh, face,
   *  point }`, `{ kind: "edge", mesh, edge, point }` (an edge within
   *  EDGE_PICK_PX of the pointer and not hidden behind a face), or null. */
  hitTest(clientX, clientY) {
    const rect = this.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const faceHit = this.raycaster.intersectObjects(this.bodies, false)[0];
    // World size of the pixel slack at the depth of the hit (or of the
    // orbit target when the ray misses every face).
    const depth = faceHit ? faceHit.distance : this.camera.position.distanceTo(this.orbit.target);
    const perPixel = (2 * depth * Math.tan(THREE.MathUtils.degToRad(this.camera.fov) / 2)) / Math.max(1, rect.height);
    const slack = EDGE_PICK_PX * perPixel;
    this.raycaster.params.Line.threshold = slack;
    const edgeHit = this.edgeLines.length ? this.raycaster.intersectObjects(this.edgeLines, false)[0] : null;
    if (edgeHit && (!faceHit || edgeHit.distance <= faceHit.distance + 2 * slack)) {
      const { edges, index } = edgeHit.object.userData;
      const edge = edges.segmentEdge[Math.floor(edgeHit.index / 2)];
      return { kind: "edge", mesh: this.bodies[index], edge, point: edgeHit.point.toArray() };
    }
    if (!faceHit) return null;
    const face = faceHit.object.userData.faceIds[faceHit.faceIndex];
    return { kind: "face", mesh: faceHit.object, face, point: faceHit.point.toArray() };
  }

  /** Pre-highlight the hovered face or edge (null clears it). */
  setHover(hit) {
    const same =
      hit && this.hovered && hit.kind === this.hovered.kind && hit.mesh === this.hovered.mesh &&
      (hit.kind === "face" ? hit.face === this.hovered.face : hit.edge === this.hovered.edge);
    if (same || (!hit && !this.hovered)) return;
    if (this.highlight) {
      this.group.remove(this.highlight);
      this.highlight.geometry.dispose();
      this.highlight.material.dispose();
      this.highlight = null;
    }
    this.hovered = hit;
    if (hit?.kind === "face") this.highlight = this.faceHighlight(hit.mesh, hit.face);
    else if (hit?.kind === "edge") this.highlight = this.edgeHighlight(hit.mesh.userData.index, hit.edge);
    if (this.highlight) this.group.add(this.highlight);
    this.canvas.style.cursor = hit ? "pointer" : "";
    this.render();
  }

  /** The triangles of `face` on `mesh`, over the mesh's own vertices. */
  faceHighlight(mesh, face) {
    const { faceIds } = mesh.userData;
    const source = mesh.geometry;
    const index = source.index.array;
    const picked = [];
    for (let t = 0; t < faceIds.length; t++) {
      if (faceIds[t] === face) picked.push(index[t * 3], index[t * 3 + 1], index[t * 3 + 2]);
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", source.getAttribute("position"));
    geometry.setAttribute("normal", source.getAttribute("normal"));
    geometry.setIndex(picked);
    const material = new THREE.MeshStandardMaterial({
      color: HOVER,
      metalness: 0.05,
      roughness: 0.55,
      side: THREE.DoubleSide,
      polygonOffset: true,
      polygonOffsetFactor: HOVER_OFFSET,
      polygonOffsetUnits: HOVER_OFFSET,
    });
    return new THREE.Mesh(geometry, material);
  }

  /** The segments of boundary edge `edge` on body `index`. */
  edgeHighlight(index, edge) {
    const lines = this.edgeLines.find((l) => l.userData.index === index);
    if (!lines) return null;
    const { edges } = lines.userData;
    const points = [];
    for (let s = 0; s < edges.segmentEdge.length; s++) {
      if (edges.segmentEdge[s] !== edge) continue;
      for (let k = 0; k < 6; k++) points.push(edges.positions[s * 6 + k]);
    }
    const highlight = fatLines(THREE, new Float32Array(points), HOVER, HOVER_EDGE_PX * this.renderer.getPixelRatio());
    highlight.renderOrder = 1;
    return highlight;
  }

  pick(event) {
    const hit = this.hitTest(event.clientX, event.clientY);
    if (!hit) {
      this.onPick(null);
      return;
    }
    const { body } = hit.mesh.userData;
    const name = (face) => (face === null ? "open boundary" : body.face_labels[face] ?? "unnamed");
    if (hit.kind === "edge") {
      const lines = this.edgeLines.find((l) => l.userData.index === hit.mesh.userData.index);
      const faces = lines.userData.edges.edges[hit.edge].faces;
      this.onPick({ body, edge: { faces, labels: faces.map(name) }, point: hit.point });
      return;
    }
    this.onPick({ body, face: hit.face, label: name(hit.face), point: hit.point });
  }

  render() {
    // Zebra stripes follow the camera's right axis so the family
    // stays perpendicular to the view while orbiting.
    if (this.zebraMaterials?.length) {
      const forward = new THREE.Vector3();
      this.camera.getWorldDirection(forward);
      const right = forward.clone().cross(this.camera.up).normalize();
      for (const m of this.zebraMaterials) m.uniforms.stripeAxis.value.copy(right);
    }
    const viewport = this.highlight?.material.uniforms?.viewport;
    if (viewport) this.renderer.getDrawingBufferSize(viewport.value);
    this.renderer.render(this.scene, this.camera);
  }
}
