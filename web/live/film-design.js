// A film's generative designs converging (PK-959): for each design an
// evolve track names, a layer of two display meshes — the design at a
// position among the records its solve kept (two neighbouring records
// blended, meshed in the worker), over its design space drawn faint
// with its outline. They are not bodies: no picking, export or
// section; the film player moves, fades, shows and highlights them by
// their selector ("design carve") like any body.
//
// A design's history is in the worker when this page solved it, else
// in a bake's design file, which is fetched and handed over here the
// first time a film needs it.
import * as THREE from "./vendor/three.module.min.js";
import { designsAt } from "./film-state.js";

/** The design space's opacity, and its outline's. */
const SPACE_OPACITY = 0.1;
const OUTLINE_OPACITY = 0.35;
/** Positions are meshed to 1/STEPS of a record: smooth at any film
 *  length (34 records over 8 s at 30 fps is under 1/7 of a record a
 *  frame), and the same frames export after export. */
const STEPS = 32;
/** Meshes kept per design (a scrub back and forth reuses them). */
const KEEP = 64;

function geometryOf(mesh) {
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(mesh.positions, 3));
  geometry.setAttribute("normal", new THREE.BufferAttribute(mesh.normals, 3));
  geometry.setIndex(new THREE.BufferAttribute(mesh.indices, 1));
  return geometry;
}

function boundsOf(positions) {
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += 3) {
    for (let d = 0; d < 3; d++) {
      lo[d] = Math.min(lo[d], positions[i + d]);
      hi[d] = Math.max(hi[d], positions[i + d]);
    }
  }
  return [lo, hi];
}

export class DesignLayers {
  constructor(page) {
    this.page = page;
    this.layers = new Map();
    this.group = null;
  }

  get active() {
    return this.layers.size > 0;
  }

  /** Load what `film`'s evolve tracks need from the worker (seeding the
   *  bakes' design files when it holds no history for a design).
   *  Resolves to null, or the reason a design cannot be shown. */
  async prepare(film, context) {
    this.remove();
    const selectors = new Set(film.tracks.filter((t) => t.channel.kind === "evolve").map((t) => t.channel.target));
    if (!selectors.size) return null;
    const host = this.page.host;
    const group = new THREE.Group();
    group.name = "film-designs";
    const layers = new Map();
    for (const selector of selectors) {
      const design = context?.designs?.find((d) => `design ${d.label}` === selector);
      if (!design) return `${selector} is not in the document`;
      let reply = await host.designHistory(design.design);
      if (!reply.ok) {
        for await (const bytes of this.page.bakes.designFiles()) {
          await host.seedDesignHistory(bytes);
          reply = await host.designHistory(design.design);
          if (reply.ok) break;
        }
      }
      if (!reply.ok) return `the solve history of ${selector} is not here: re-solve the design to record it`;
      const space = await host.designMesh(design.design, 0, true);
      if (!space.ok) return `the design space of ${selector} did not load`;
      const spaceGeometry = geometryOf(space.mesh);
      const bounds = boundsOf(space.mesh.positions);
      const ghost = new THREE.Mesh(
        spaceGeometry,
        new THREE.MeshLambertMaterial({ color: 0x8a94a6, transparent: true, opacity: SPACE_OPACITY, depthWrite: false, side: THREE.DoubleSide }),
      );
      const outline = new THREE.LineSegments(
        new THREE.EdgesGeometry(spaceGeometry, 30),
        new THREE.LineBasicMaterial({ color: 0x4b5563, transparent: true, opacity: OUTLINE_OPACITY, depthWrite: false }),
      );
      // Hidden until its first mesh; styled as the model at its first
      // frame (the film's look is on by then).
      const carved = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshStandardMaterial({ color: 0x9aa4b2 }));
      carved.castShadow = true;
      carved.visible = false;
      for (const [object, scale] of [[carved, 1], [ghost, SPACE_OPACITY], [outline, OUTLINE_OPACITY]]) {
        object.userData = { instance: selector, body: { bounds }, opacityScale: scale, designLayer: true };
        object.matrixAutoUpdate = false;
        group.add(object);
      }
      // Drawn under the model's bodies' order: the faint space last.
      ghost.renderOrder = 2;
      outline.renderOrder = 3;
      layers.set(selector, { selector, design: design.design, frames: reply.history.frames, carved, ghost, outline, meshes: new Map(), shown: null, styled: false });
    }
    this.layers = layers;
    this.group = group;
    this.page.viewer.scene.add(group);
    return null;
  }

  /** The layer objects the player animates, each as a body mesh is:
   *  `{ object, mesh }` (the mesh carries the selector and bounds). */
  objects() {
    const out = [];
    for (const layer of this.layers.values()) {
      for (const object of [layer.carved, layer.ghost, layer.outline]) out.push({ object, mesh: layer.carved });
    }
    return out;
  }

  /** The selectors of the layers. */
  selectors() {
    return [...this.layers.keys()];
  }

  /** Show each design at its position at `t` (meshing it in the worker
   *  when it is not cached). The carved mesh wears the model's own
   *  material, so it takes the film's look. */
  async apply(film, t) {
    if (!this.active) return;
    const positions = designsAt(film, t, (selector) => this.layers.get(selector)?.frames);
    await Promise.all(
      [...this.layers.values()].map(async (layer) => {
        const body = this.page.viewer.bodies[0]?.material;
        if (!layer.styled && body) {
          const material = body.clone();
          material.vertexColors = false;
          material.side = THREE.DoubleSide;
          material.needsUpdate = true;
          layer.carved.material.dispose();
          layer.carved.material = material;
          layer.styled = true;
        }
        const key = Math.round((positions.get(layer.selector) ?? 0) * STEPS) / STEPS;
        layer.wanted = key;
        if (layer.shown === key) return;
        let geometry = layer.meshes.get(key);
        if (!geometry) {
          const reply = await this.page.host.designMesh(layer.design, key);
          if (!reply.ok) return;
          geometry = geometryOf(reply.mesh);
          layer.meshes.set(key, geometry);
          if (layer.meshes.size > KEEP) {
            const [oldest, old] = layer.meshes.entries().next().value;
            layer.meshes.delete(oldest);
            if (old !== layer.carved.geometry) old.dispose();
          }
        }
        // A later frame asked for another position meanwhile.
        if (layer.wanted !== key) return;
        layer.carved.geometry = geometry;
        layer.shown = key;
      }),
    );
  }

  /** Take the layers off the scene and free them. */
  remove() {
    if (this.group) this.page.viewer.scene.remove(this.group);
    for (const layer of this.layers.values()) {
      for (const geometry of layer.meshes.values()) geometry.dispose();
      layer.ghost.geometry.dispose();
      layer.outline.geometry.dispose();
      for (const object of [layer.carved, layer.ghost, layer.outline]) object.material?.dispose();
    }
    this.layers = new Map();
    this.group = null;
  }
}
