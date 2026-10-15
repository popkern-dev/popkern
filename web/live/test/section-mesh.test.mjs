// The clipping view's cut outline on kernel display meshes handed in as
// JSON on stdin: `{ meshes: [{positions, indices}], planes: [{normal,
// offset}] }`. Prints one row per (mesh, plane): `{ segments, length,
// open }`, where `open` counts outline endpoints used an odd number of
// times (0 for a closed outline). Run by crates/popkern-live/tests/
// page_section.rs.
import { readFileSync } from "node:fs";
import { cutOutline } from "../section-geometry.js";

const { meshes, planes } = JSON.parse(readFileSync(0, "utf8"));
const rows = [];
for (const m of meshes) {
  const positions = Float32Array.from(m.positions);
  const indices = Uint32Array.from(m.indices);
  for (const plane of planes) {
    const flat = cutOutline(positions, indices, null, [plane], 0);
    const ends = new Map();
    let length = 0;
    for (let i = 0; i < flat.length; i += 6) {
      length += Math.hypot(flat[i + 3] - flat[i], flat[i + 4] - flat[i + 1], flat[i + 5] - flat[i + 2]);
      // Weld at 1e-6 of the model size: the display mesh's shared
      // vertices give bit-identical crossing points on shared edges.
      for (const o of [0, 3]) {
        const k = [flat[i + o], flat[i + o + 1], flat[i + o + 2]].map((c) => Math.round(c * 1e6)).join(",");
        ends.set(k, (ends.get(k) ?? 0) + 1);
      }
    }
    let open = 0;
    for (const n of ends.values()) if (n % 2) open++;
    rows.push({ segments: flat.length / 6, length, open });
  }
}
process.stdout.write(JSON.stringify(rows));
