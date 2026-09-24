// The three.js viewer: one mesh per body from the engine's buffers,
// per-face colouring and picking through the face-id attribute, a
// small orbit control (no addon), fit-to-bounds on the first render.
import * as THREE from "./vendor/three.module.min.js";

const PALETTE = [0x6c8ebf, 0x8fb56c, 0xd9a25f, 0xb87ab8, 0x66b3b3, 0xc97a7a, 0x9c9c6c, 0x7a8fc9];

function faceColor(index) {
  return new THREE.Color(PALETTE[index % PALETTE.length]);
}

class Orbit {
  constructor(camera, element, onChange) {
    this.camera = camera;
    this.element = element;
    this.onChange = onChange;
    this.target = new THREE.Vector3();
    this.spherical = new THREE.Spherical(5, Math.PI / 3, Math.PI / 4);
    this.dragging = null;
    element.addEventListener("pointerdown", (e) => {
      this.dragging = { x: e.clientX, y: e.clientY, button: e.button, shift: e.shiftKey };
      element.setPointerCapture(e.pointerId);
    });
    element.addEventListener("pointermove", (e) => {
      if (!this.dragging) return;
      const dx = e.clientX - this.dragging.x;
      const dy = e.clientY - this.dragging.y;
      this.dragging.x = e.clientX;
      this.dragging.y = e.clientY;
      if (this.dragging.button === 2 || this.dragging.shift) {
        this.pan(dx, dy);
      } else {
        this.spherical.theta -= dx * 0.008;
        this.spherical.phi = Math.min(Math.PI - 0.05, Math.max(0.05, this.spherical.phi - dy * 0.008));
      }
      this.update();
    });
    const release = () => (this.dragging = null);
    element.addEventListener("pointerup", release);
    element.addEventListener("pointercancel", release);
    element.addEventListener("contextmenu", (e) => e.preventDefault());
    element.addEventListener(
      "wheel",
      (e) => {
        e.preventDefault();
        this.spherical.radius *= Math.exp(e.deltaY * 0.001);
        this.update();
      },
      { passive: false },
    );
    element.addEventListener("keydown", (e) => {
      const step = 0.05;
      if (e.key === "ArrowLeft") this.spherical.theta += step;
      else if (e.key === "ArrowRight") this.spherical.theta -= step;
      else if (e.key === "ArrowUp") this.spherical.phi = Math.max(0.05, this.spherical.phi - step);
      else if (e.key === "ArrowDown") this.spherical.phi = Math.min(Math.PI - 0.05, this.spherical.phi + step);
      else if (e.key === "+" || e.key === "=") this.spherical.radius *= 0.9;
      else if (e.key === "-") this.spherical.radius *= 1.1;
      else return;
      e.preventDefault();
      this.update();
    });
  }

  pan(dx, dy) {
    const scale = this.spherical.radius * 0.0015;
    const right = new THREE.Vector3();
    const up = new THREE.Vector3(0, 0, 1);
    this.camera.getWorldDirection(right);
    right.cross(up).normalize();
    const trueUp = right.clone().cross(new THREE.Vector3().subVectors(this.target, this.camera.position).normalize()).normalize();
    this.target.addScaledVector(right, -dx * scale);
    this.target.addScaledVector(trueUp, dy * scale);
  }

  fit(bounds) {
    const min = new THREE.Vector3(...bounds[0]);
    const max = new THREE.Vector3(...bounds[1]);
    this.target = min.clone().add(max).multiplyScalar(0.5);
    const size = max.clone().sub(min).length() || 1;
    this.spherical.radius = size * 1.4;
    this.update();
  }

  update() {
    // Z is up: spherical about the Z axis.
    const s = this.spherical;
    const x = s.radius * Math.sin(s.phi) * Math.cos(s.theta);
    const y = s.radius * Math.sin(s.phi) * Math.sin(s.theta);
    const z = s.radius * Math.cos(s.phi);
    this.camera.position.set(this.target.x + x, this.target.y + y, this.target.z + z);
    this.camera.up.set(0, 0, 1);
    this.camera.lookAt(this.target);
    this.onChange();
  }
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
    this.mode = "shaded";
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
    this.showEdges = true;
    this.fitted = false;
    this.raycaster = new THREE.Raycaster();
    this.orbit = new Orbit(this.camera, canvas, () => this.render());
    canvas.tabIndex = 0;
    canvas.addEventListener("click", (e) => this.pick(e));
    this.resize();
    window.addEventListener("resize", () => this.resize());
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

  /** Replace the scene with the envelope's bodies and their buffers. */
  show(bodies, buffers) {
    this.lastShown = { bodies, buffers };
    this.group.clear();
    this.bodies = [];
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
        const edges = new THREE.LineSegments(new THREE.EdgesGeometry(geometry, 25), new THREE.LineBasicMaterial({ color: 0x0f1115 }));
        this.group.add(edges);
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
    if (this.mode === "gaussian" || this.mode === "mean") {
      return new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
    }
    return new THREE.MeshStandardMaterial({ vertexColors: true, metalness: 0.1, roughness: 0.65, side: THREE.DoubleSide });
  }

  /** Switch the surface view: shaded | zebra | gaussian | mean. */
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

  pick(event) {
    const rect = this.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((event.clientX - rect.left) / rect.width) * 2 - 1, -((event.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObjects(this.bodies, false);
    if (!hits.length) {
      this.onPick(null);
      return;
    }
    const hit = hits[0];
    const { body, faceIds } = hit.object.userData;
    const face = faceIds[hit.faceIndex];
    this.onPick({ body, face, label: body.face_labels[face] ?? "unnamed", point: hit.point.toArray() });
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
    this.renderer.render(this.scene, this.camera);
  }
}
