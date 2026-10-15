// The presentation look for films (docs/architecture/presentation-films.md
// §6): flat toon shading in each part's colour, the face edges in black,
// a plain background and a soft contact shadow on a ground plane under
// the model. Applied to the viewer while a film is open, removed when it
// closes; nothing here changes the model.
import * as THREE from "./vendor/three.module.min.js";

/** A three-tone ramp: shadow, mid, lit. */
function toonRamp() {
  const data = new Uint8Array([90, 90, 90, 255, 175, 175, 175, 255, 255, 255, 255, 255]);
  const texture = new THREE.DataTexture(data, 3, 1, THREE.RGBAFormat);
  texture.minFilter = THREE.NearestFilter;
  texture.magFilter = THREE.NearestFilter;
  texture.generateMipmaps = false;
  texture.needsUpdate = true;
  return texture;
}

let ramp = null;

/** A toon material in `color` (a THREE.Color), drawn on both sides with
 *  the viewer's face depth offset so edges stay on top. */
export function toonMaterial(color, offset) {
  ramp ??= toonRamp();
  return new THREE.MeshToonMaterial({ color, gradientMap: ramp, side: THREE.DoubleSide, ...offset });
}

/** The ground shadow and the look's scene state, keyed to one viewer. */
export class FilmLook {
  constructor(viewer) {
    this.viewer = viewer;
    this.ground = null;
    this.saved = null;
  }

  /** Put the look on: `look` ("toon" | "studio" | "cad"), `background`
   *  ("#rrggbb"), `bounds` (the model at its placements). */
  apply({ look, background }, bounds) {
    const viewer = this.viewer;
    if (!this.saved) {
      this.saved = {
        mode: viewer.mode,
        clear: viewer.renderer.getClearColor(new THREE.Color()),
        alpha: viewer.renderer.getClearAlpha(),
        shadows: viewer.renderer.shadowMap.enabled,
      };
    }
    const mode = look === "toon" ? "toon" : "cad";
    if (viewer.mode !== mode) viewer.setMode(mode);
    viewer.renderer.setClearColor(new THREE.Color(background), 1);
    this.background = background;
    this.placeGround(bounds);
  }

  /** A shadow-catching plane just under the model and a key light that
   *  casts onto it, fitted to `bounds`. */
  placeGround(bounds) {
    const viewer = this.viewer;
    if (!bounds) return;
    const [lo, hi] = bounds;
    const c = new THREE.Vector3((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2);
    const r = Math.max(1e-6, Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2);
    if (!this.ground) {
      this.ground = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), new THREE.ShadowMaterial({ opacity: 0.16 }));
      this.ground.receiveShadow = true;
      this.ground.name = "film-ground";
      this.light = new THREE.DirectionalLight(0xffffff, 0.0);
      this.light.castShadow = true;
      this.light.shadow.mapSize.set(2048, 2048);
      this.light.shadow.radius = 6;
      this.light.name = "film-shadow-light";
      viewer.scene.add(this.light, this.light.target);
      viewer.scene.add(this.ground);
    }
    this.ground.scale.set(8 * r, 8 * r, 1);
    this.ground.position.set(c.x, c.y, lo[2] - 1e-3 * r);
    // The shadow falls from above and a little in front, like the key.
    this.light.position.set(c.x + 0.6 * r, c.y - 0.9 * r, c.z + 4 * r);
    this.light.target.position.copy(c);
    const cam = this.light.shadow.camera;
    cam.left = cam.bottom = -2.5 * r;
    cam.right = cam.top = 2.5 * r;
    cam.near = 0.1 * r;
    cam.far = 10 * r;
    cam.updateProjectionMatrix();
    viewer.renderer.shadowMap.enabled = true;
    viewer.renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    this.castShadows();
  }

  /** Every body casts (after a redraw rebuilt the meshes). */
  castShadows() {
    for (const mesh of this.viewer.bodies) mesh.castShadow = true;
  }

  /** Enable shadows on another renderer (the export's). */
  prepare(renderer) {
    renderer.shadowMap.enabled = !!this.ground;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    renderer.setClearColor(new THREE.Color(this.background ?? "#ffffff"), 1);
  }

  /** Take the look off and restore the viewer as it was. */
  remove() {
    const viewer = this.viewer;
    if (this.ground) {
      viewer.scene.remove(this.ground, this.light, this.light.target);
      this.ground.geometry.dispose();
      this.ground.material.dispose();
      this.light.dispose?.();
      this.ground = null;
      this.light = null;
    }
    if (this.saved) {
      viewer.renderer.setClearColor(this.saved.clear, this.saved.alpha);
      viewer.renderer.shadowMap.enabled = this.saved.shadows;
      const mode = this.saved.mode;
      this.saved = null;
      if (viewer.mode !== mode) viewer.setMode(mode);
    }
    viewer.render();
  }
}
