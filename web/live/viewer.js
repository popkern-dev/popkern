// The three.js viewer: one mesh per body from the engine's buffers,
// per-face colouring and picking through the face-id attribute, the
// faces' boundary edges, hover pre-highlighting in the CAD view, a
// small orbit control (no addon), fit-to-bounds on the first render.
// An assembly draws every part occurrence as a mesh SHARING its
// definition's geometry under the occurrence's world matrix (design
// docs/architecture/assemblies.md §6.4): GPU memory scales with
// definitions, not instances.
import * as THREE from "./vendor/three.module.min.js";
import { boundaryEdges, enclosesNoVolume, fatLines } from "./edges.js";
import { drawList, explodePlan, placeAt } from "./instances.js";
import { Orbit } from "./orbit.js";
import { Section } from "./section.js";
import { ZEBRA_FRAGMENT, ZEBRA_VERTEX, divergingColor, faceColor, robustHalfRange } from "./shading.js";

// The CAD view: uniform grey faces, black edges, and a blue
// pre-highlight on the face or edge under the pointer.
const CAD_FACE = 0xb4b8be;
const CAD_EDGE = 0x000000;
const HOVER = 0x3b82f6;
// History scrubber: the faces a step created glow, then settle.
const GLOW = 0xf59e0b;
const GLOW_MS = 1200;
// A history step's tool bodies (consumed by a later step) are drawn
// see-through, so a slab or cutter does not hide the part it shapes.
const GHOST_OPACITY = 0.3;
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

export class Viewer {
  constructor(canvas, onPick) {
    this.canvas = canvas;
    this.onPick = onPick;
    this.mode = "cad";
    // The stencil buffer carries the section caps (section.js); three.js
    // allocates none by default since r163.
    // A transparent canvas: the model floats over the page's ground.
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true, stencil: true });
    this.renderer.setClearColor(0x000000, 0);
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.scene = new THREE.Scene();
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
    this.section = new Section(this);
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

  /** Replace the scene with the envelope's bodies and their buffers;
   *  `assembly` (the envelope's assembly view, or null) places shared
   *  definition meshes per occurrence. The zebra view draws the
   *  analysis mesh when one is loaded. */
  /** `ghosts`: indices of bodies drawn see-through (a history step's
   *  tool bodies, which a later step consumes). */
  show(bodies, displayBuffers, assembly = null, ghosts = null) {
    this.lastShown = { bodies, buffers: displayBuffers, assembly, ghosts };
    const buffers =
      this.mode === "zebra" && this.analysis?.length === bodies.length ? this.analysis : displayBuffers;
    this.group.clear();
    this.glows = [];
    this.bodies = [];
    this.edgeLines = [];
    this.highlight = null;
    this.hovered = null;
    this.zebraMaterials = [];
    this.geometries = [];
    const halfRange = this.mode === "gaussian" || this.mode === "mean" ? robustHalfRange(buffers, this.mode) : 0;
    // One geometry (and edge set) per body buffer, shared by every
    // occurrence that draws it.
    const shared = bodies.map((body, i) => {
      const b = buffers[i];
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", new THREE.BufferAttribute(b.positions, 3));
      geometry.setAttribute("normal", new THREE.BufferAttribute(b.normals, 3));
      geometry.setIndex(new THREE.BufferAttribute(b.indices, 1));
      geometry.setAttribute("color", this.vertexColors(b, halfRange));
      this.geometries.push(geometry);
      let lineGeometry = null;
      const edges = this.showEdges ? cachedEdges(b) : null;
      if (edges) {
        lineGeometry = new THREE.BufferGeometry();
        lineGeometry.setAttribute("position", new THREE.BufferAttribute(edges.positions, 3));
      }
      return { geometry, lineGeometry, edges, faceIds: b.face_ids };
    });
    for (const draw of drawList(bodies, assembly)) {
      const { geometry, lineGeometry, edges, faceIds } = shared[draw.index];
      const mesh = new THREE.Mesh(geometry, this.material(draw.colorKey, bodies[draw.index].appearance));
      // A sketch's lamina: each of its coincident faces from its own side.
      const sheet = bodies[draw.index].open_edges === 0 && enclosesNoVolume(buffers[draw.index]);
      if (sheet) mesh.material.side = THREE.FrontSide;
      const ghost = ghosts?.has(draw.index) ?? false;
      if (ghost) Object.assign(mesh.material, { transparent: true, opacity: GHOST_OPACITY, depthWrite: false });
      mesh.userData = { body: bodies[draw.index], faceIds, index: draw.index, instance: draw.instance, sheet };
      if (draw.matrix) placeAt(mesh, draw.matrix);
      mesh.userData.base = mesh.matrix.clone();
      mesh.visible = draw.visible;
      this.group.add(mesh);
      if (lineGeometry) {
        const color = this.mode === "cad" || this.mode === "parts" ? CAD_EDGE : 0x0f1115;
        const lines = new THREE.LineSegments(
          lineGeometry,
          new THREE.LineBasicMaterial(ghost ? { color, transparent: true, opacity: GHOST_OPACITY } : { color }),
        );
        lines.userData = { edges, index: draw.index, instance: draw.instance };
        if (draw.matrix) placeAt(lines, draw.matrix);
        lines.userData.base = lines.matrix.clone();
        lines.visible = draw.visible;
        mesh.userData.lines = lines;
        this.group.add(lines);
        this.edgeLines.push(lines);
      }
      this.bodies.push(mesh);
    }
    this.explodePlan = assembly ? explodePlan(this.bodies) : null;
    if (this.explode && this.explodePlan) this.setExplode(this.explode);
    if (this.section.active) this.section.apply();
    if (!this.fitted) {
      const box = this.visibleBounds();
      if (box) {
        this.orbit.fit(box);
        this.fitted = true;
      }
    }
    this.render();
  }

  /** Re-place the meshes on screen at `assembly`'s occurrence
   *  matrices (a joint moved: the same bodies at new placements).
   *  Nothing rebuilds; hidden occurrences stay hidden. */
  place(assembly) {
    const instances = new Map(assembly.instances.map((i) => [i.path, i]));
    for (const object of [...this.bodies, ...this.edgeLines]) {
      const instance = instances.get(object.userData.instance);
      if (!instance) continue;
      placeAt(object, instance.world);
      object.userData.base = object.matrix.clone();
    }
    if (this.lastShown) this.lastShown = { ...this.lastShown, assembly };
    this.explodePlan = explodePlan(this.bodies);
    if (this.explode && this.explodePlan) this.setExplode(this.explode);
    if (this.section.active) this.section.apply();
    if (this.highlight && this.hovered) this.setHover(null);
    this.render();
  }

  /** World bounds of the drawn (visible) meshes, or null. */
  visibleBounds() {
    const box = new THREE.Box3();
    for (const mesh of this.bodies) {
      if (!mesh.visible) continue;
      const [lo, hi] = mesh.userData.body.bounds;
      const local = new THREE.Box3(new THREE.Vector3(...lo), new THREE.Vector3(...hi));
      box.union(local.applyMatrix4(mesh.matrixWorld));
    }
    return box.isEmpty() ? null : [box.min.toArray(), box.max.toArray()];
  }

  /** Exploded view (PK-630): move every top-level occurrence away from
   *  the assembly centre by `factor` (0 = assembled) times the model
   *  size. Pure transforms on the shared geometry: nothing rebuilds and
   *  nothing re-meshes, so the slider runs at frame rate. */
  setExplode(factor) {
    this.explode = factor;
    if (!this.explodePlan) return;
    const offset = new THREE.Matrix4();
    for (const object of [...this.bodies, ...this.edgeLines]) {
      const base = object.userData.base;
      if (!base) continue;
      const direction = this.explodePlan.directions.get(object.userData.instance?.split("/")[0]);
      if (!direction) continue;
      offset.makeTranslation(direction.x * factor, direction.y * factor, direction.z * factor);
      object.matrix.multiplyMatrices(offset, base);
      object.matrixWorldNeedsUpdate = true;
    }
    this.group.updateMatrixWorld(true);
    if (this.highlight && this.hovered) this.setHover(null);
    this.render();
  }

  /** Show or hide every occurrence at or under `path`. */
  setVisible(path, visible) {
    for (const object of [...this.bodies, ...this.edgeLines]) {
      const instance = object.userData.instance;
      if (instance === path || instance?.startsWith(`${path}/`)) object.visible = visible;
    }
    if (this.hovered && !this.hovered.mesh.visible) this.setHover(null);
    this.render();
  }

  /** Emphasise the occurrences at or under `path` (a path, an array of
   *  paths, or null to clear). */
  select(path) {
    this.selected = path;
    const paths = path == null ? [] : Array.isArray(path) ? path : [path];
    for (const mesh of this.bodies) {
      const instance = mesh.userData.instance;
      const on = paths.some((p) => instance === p || instance?.startsWith(`${p}/`));
      if (mesh.material.emissive) mesh.material.emissive.setHex(on ? 0x1d3d73 : 0x000000);
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

  material(colorKey = 0, appearance = null) {
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
      // A part's own appearance (six-axis-arm §4.5) is its CAD material.
      if (appearance) {
        const [r, g, b] = appearance.color;
        const color = new THREE.Color().setRGB(r, g, b, THREE.SRGBColorSpace);
        return new THREE.MeshStandardMaterial({
          color,
          metalness: appearance.metal,
          roughness: appearance.rough,
          envMap: this.studioEnvironment(),
          envMapIntensity: 0.9,
          side: THREE.DoubleSide,
          ...offset,
        });
      }
      return new THREE.MeshStandardMaterial({ color: CAD_FACE, metalness: 0.05, roughness: 0.55, side: THREE.DoubleSide, ...offset });
    }
    if (this.mode === "parts") {
      // One colour per part definition (per body for a part document).
      return new THREE.MeshStandardMaterial({ color: faceColor(colorKey), metalness: 0.05, roughness: 0.55, side: THREE.DoubleSide, ...offset });
    }
    return new THREE.MeshStandardMaterial({ vertexColors: true, metalness: 0.1, roughness: 0.65, side: THREE.DoubleSide, ...offset });
  }

  /** A soft studio environment (built once): a metal shows only what it
   *  reflects, so a metallic appearance under plain lights reads black.
   *  Only appearance materials take it — the default CAD look is
   *  unchanged. */
  studioEnvironment() {
    if (this.studio) return this.studio;
    const room = new THREE.Scene();
    const geometry = new THREE.SphereGeometry(10, 32, 16);
    const colors = [];
    const position = geometry.getAttribute("position");
    for (let i = 0; i < position.count; i++) {
      // Bright overhead, warm horizon band, dark floor.
      const t = position.getZ(i) / 10;
      const c = t > 0.35 ? 1.0 : t > -0.1 ? 0.55 + 0.6 * (t + 0.1) : 0.12;
      colors.push(c, c * 0.98, c * 0.95);
    }
    geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
    room.add(new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.BackSide })));
    const pmrem = new THREE.PMREMGenerator(this.renderer);
    this.studio = pmrem.fromScene(room, 0.02).texture;
    pmrem.dispose();
    geometry.dispose();
    return this.studio;
  }

  /** Switch the surface view: cad | parts | shaded | zebra | gaussian | mean. */
  setMode(mode) {
    this.mode = mode;
    this.redraw();
  }

  toggleEdges(on) {
    this.showEdges = on;
    this.redraw();
  }

  /** Rebuild the scene from the last shown model, keeping per-occurrence
   *  visibility and the selection. */
  redraw() {
    if (!this.lastShown) return;
    const hidden = new Set(this.bodies.filter((m) => !m.visible && m.userData.instance).map((m) => m.userData.instance));
    this.show(this.lastShown.bodies, this.lastShown.buffers, this.lastShown.assembly, this.lastShown.ghosts);
    for (const object of [...this.bodies, ...this.edgeLines]) {
      if (hidden.has(object.userData.instance)) object.visible = false;
    }
    if (this.selected != null) this.select(this.selected);
    this.render();
  }

  fitAgain() {
    this.fitted = false;
    const box = this.visibleBounds();
    if (box) {
      this.orbit.fit(box);
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
    const faceHit = this.raycaster.intersectObjects(this.bodies.filter((m) => m.visible), false)[0];
    // World size of the pixel slack at the depth of the hit (or of the
    // orbit target when the ray misses every face).
    const depth = faceHit ? faceHit.distance : this.camera.position.distanceTo(this.orbit.target);
    const perPixel = (2 * depth * Math.tan(THREE.MathUtils.degToRad(this.camera.fov) / 2)) / Math.max(1, rect.height);
    const slack = EDGE_PICK_PX * perPixel;
    this.raycaster.params.Line.threshold = slack;
    const lines = this.edgeLines.filter((l) => l.visible);
    const edgeHit = lines.length ? this.raycaster.intersectObjects(lines, false)[0] : null;
    if (edgeHit && (!faceHit || edgeHit.distance <= faceHit.distance + 2 * slack)) {
      const { edges } = edgeHit.object.userData;
      const edge = edges.segmentEdge[Math.floor(edgeHit.index / 2)];
      const mesh = this.bodies.find((m) => m.userData.lines === edgeHit.object);
      return { kind: "edge", mesh, edge, point: edgeHit.point.toArray() };
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
    else if (hit?.kind === "edge") this.highlight = this.edgeHighlight(hit.mesh, hit.edge);
    if (this.highlight) this.group.add(this.highlight);
    this.canvas.style.cursor = hit ? "pointer" : "";
    this.render();
  }

  /** Glow the faces a history step created and fade them out:
   *  `created` maps a drawn body index to its face indices; every face
   *  of the `occurrences` (assembly paths, with what lies under them)
   *  glows whole. */
  glowFaces(created, occurrences = []) {
    for (const glow of this.glows ?? []) this.group.remove(glow);
    this.glows = [];
    cancelAnimationFrame(this.glowFrame);
    for (const mesh of this.bodies) {
      const path = mesh.userData.instance;
      const whole = path != null && occurrences.some((p) => path === p || path.startsWith(`${p}/`));
      const faces = created.get(mesh.userData.index);
      if (!whole && !faces?.size) continue;
      const { faceIds } = mesh.userData;
      const index = mesh.geometry.index.array;
      const picked = [];
      for (let t = 0; t < faceIds.length; t++) {
        if (whole || faces.has(faceIds[t])) picked.push(index[t * 3], index[t * 3 + 1], index[t * 3 + 2]);
      }
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute("position", mesh.geometry.getAttribute("position"));
      geometry.setAttribute("normal", mesh.geometry.getAttribute("normal"));
      geometry.setIndex(picked);
      const material = new THREE.MeshStandardMaterial({
        color: GLOW,
        emissive: GLOW,
        emissiveIntensity: 0.35,
        transparent: true,
        opacity: 0.85,
        side: mesh.userData.sheet ? THREE.FrontSide : THREE.DoubleSide,
        polygonOffset: true,
        polygonOffsetFactor: HOVER_OFFSET,
        polygonOffsetUnits: HOVER_OFFSET,
      });
      const glow = new THREE.Mesh(geometry, material);
      if (!mesh.matrixAutoUpdate) placeAt(glow, mesh.matrix);
      this.group.add(glow);
      this.glows.push(glow);
    }
    if (!this.glows.length) return this.render();
    const started = performance.now();
    const fade = () => {
      const t = Math.min(1, (performance.now() - started) / GLOW_MS);
      for (const glow of this.glows) glow.material.opacity = 0.85 * (1 - t) + 0.25 * t;
      this.render();
      if (t < 1) this.glowFrame = requestAnimationFrame(fade);
    };
    fade();
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
    const highlight = new THREE.Mesh(geometry, material);
    if (!mesh.matrixAutoUpdate) placeAt(highlight, mesh.matrix);
    return highlight;
  }

  /** The segments of boundary edge `edge` of the body `mesh` draws. */
  edgeHighlight(mesh, edge) {
    const lines = mesh.userData.lines;
    if (!lines) return null;
    const { edges } = lines.userData;
    const points = [];
    for (let s = 0; s < edges.segmentEdge.length; s++) {
      if (edges.segmentEdge[s] !== edge) continue;
      for (let k = 0; k < 6; k++) points.push(edges.positions[s * 6 + k]);
    }
    const highlight = fatLines(THREE, new Float32Array(points), HOVER, HOVER_EDGE_PX * this.renderer.getPixelRatio());
    highlight.renderOrder = 1;
    if (!mesh.matrixAutoUpdate) placeAt(highlight, mesh.matrix);
    return highlight;
  }

  pick(event) {
    const hit = this.hitTest(event.clientX, event.clientY);
    if (!hit) {
      this.onPick(null);
      return;
    }
    const { body, instance } = hit.mesh.userData;
    const name = (face) => (face === null ? "open boundary" : body.face_labels[face] ?? "unnamed");
    if (hit.kind === "edge") {
      const faces = hit.mesh.userData.lines.userData.edges.edges[hit.edge].faces;
      this.onPick({ body, instance, edge: { faces, labels: faces.map(name) }, point: hit.point });
      return;
    }
    this.onPick({ body, instance, face: hit.face, label: name(hit.face), point: hit.point });
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
