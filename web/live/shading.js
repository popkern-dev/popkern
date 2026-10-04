// Surface shading for the viewer: the face palette, the reflection-line
// zebra shaders, and the diverging curvature map with its robust range.
import * as THREE from "./vendor/three.module.min.js";

const PALETTE = [0x6c8ebf, 0x8fb56c, 0xd9a25f, 0xb87ab8, 0x66b3b3, 0xc97a7a, 0x9c9c6c, 0x7a8fc9];

export function faceColor(index) {
  return new THREE.Color(PALETTE[index % PALETTE.length]);
}

// Reflection-line zebra: per-pixel stripes over the reflected view
// direction — the class-A continuity inspection (a G1 defect kinks
// the stripes across a seam; a G2 defect kinks their derivative).
export const ZEBRA_VERTEX = `
  varying vec3 vWorldPos;
  varying vec3 vWorldNormal;
  void main() {
    vec4 world = modelMatrix * vec4(position, 1.0);
    vWorldPos = world.xyz;
    vWorldNormal = normalize(mat3(modelMatrix) * normal);
    gl_Position = projectionMatrix * viewMatrix * world;
  }
`;
export const ZEBRA_FRAGMENT = `
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
export function divergingColor(value, halfRange) {
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
export function robustHalfRange(buffers, key) {
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
