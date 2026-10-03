// Pre-built examples (crates/popkern-playground/src/bake.rs): a served,
// fully certified envelope and its display buffers per gallery script,
// built natively at deploy time. A script is looked up by the FNV-1a 64
// key of its text and confirmed by comparing the stored source, so only
// the exact script ever shows a bake. No index beside the page (a dev
// tree that never baked) means every run goes to the kernel.

const MAGIC = "PKBAKE01";
const FORMAT = 1;

/** FNV-1a 64 over the UTF-8 bytes of `text`, as 16 hex digits. */
export function sourceKey(text) {
  const bytes = new TextEncoder().encode(text);
  let hash = 0xcbf29ce484222325n;
  const prime = 0x100000001b3n;
  const mask = 0xffffffffffffffffn;
  for (const byte of bytes) hash = ((hash ^ BigInt(byte)) * prime) & mask;
  return hash.toString(16).padStart(16, "0");
}

export class Bakes {
  constructor() {
    this.entries = new Map();
    this.cache = new Map();
    this.ready = this.#loadIndex();
  }

  async #loadIndex() {
    try {
      const response = await fetch("./bakes/index.json");
      if (!response.ok) return;
      const index = await response.json();
      if (index.format !== FORMAT) return;
      for (const entry of index.entries ?? []) {
        if (entry.bake) this.entries.set(entry.source_key, entry);
      }
    } catch {
      // No bakes: the page runs everything live.
    }
  }

  /** The index entry for a gallery file, or null (for the picker). */
  forFile(file) {
    for (const entry of this.entries.values()) if (entry.file === file) return entry;
    return null;
  }

  /** The bake of exactly `source`: `{ envelope, buffers, native }`,
   *  or null when there is none (or it cannot be read). */
  async find(source) {
    await this.ready;
    const entry = this.entries.get(sourceKey(source));
    if (!entry) return null;
    let bake = this.cache.get(entry.bake);
    if (!bake) {
      bake = await this.#load(entry.bake);
      if (!bake) return null;
      this.cache.set(entry.bake, bake);
    }
    return bake.source === source ? bake : null;
  }

  async #load(name) {
    try {
      const response = await fetch(`./bakes/${name}`);
      if (!response.ok) return null;
      return decode(await response.arrayBuffer());
    } catch {
      return null;
    }
  }
}

/** Decode a bake file (layout in bake.rs). The typed arrays view the
 *  fetched buffer in place; every buffer starts 4-byte aligned. */
export function decode(arrayBuffer) {
  const bytes = new Uint8Array(arrayBuffer);
  if (bytes.length < 12 || new TextDecoder().decode(bytes.subarray(0, 8)) !== MAGIC) return null;
  const length = new DataView(arrayBuffer).getUint32(8, true);
  const header = JSON.parse(new TextDecoder().decode(bytes.subarray(12, 12 + length)));
  if (header.format !== FORMAT) return null;
  let offset = 12 + length;
  const take = (Kind, count) => {
    const view = new Kind(arrayBuffer, offset, count);
    offset += count * 4;
    return view;
  };
  const buffers = header.buffers.map((counts) => ({
    positions: take(Float32Array, counts.positions),
    normals: take(Float32Array, counts.normals),
    indices: take(Uint32Array, counts.indices),
    face_ids: take(Uint32Array, counts.face_ids),
    gaussian: take(Float32Array, counts.gaussian),
    mean: take(Float32Array, counts.mean),
  }));
  if (offset !== arrayBuffer.byteLength) return null;
  return { source: header.source, envelope: header.envelope, buffers, native: header.native };
}
