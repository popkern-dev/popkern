// The viewer's face-boundary edges (edges.js) for meshes handed in as
// JSON on stdin: `[{positions, indices, face_ids}]`. Prints one row per
// mesh, `{edges, segments, sheet}` with every segment's two points, for the
// Rust test (crates/popkern-live/tests/page_edges.rs) to check.
import { readFileSync } from "node:fs";
import { boundaryEdges, enclosesNoVolume } from "../edges.js";

const meshes = JSON.parse(readFileSync(0, "utf8"));
const rows = meshes.map((m) => {
  const b = {
    positions: Float32Array.from(m.positions),
    indices: Uint32Array.from(m.indices),
    face_ids: Uint32Array.from(m.face_ids),
  };
  const out = boundaryEdges(b);
  return { edges: out.edges.length, segments: Array.from(out.positions), sheet: enclosesNoVolume(b) };
});
process.stdout.write(JSON.stringify(rows));
