// PK-410 node smoke: drive the wasm-bindgen (nodejs target) bundle the
// way the page's worker will, and print what the Rust side can compare
// against its native run — outcome classes, feature counts, and FNV-1a
// hashes of every buffer (bit-identity without shipping buffers back).
//
//   node smoke.mjs <pkg-dir> <script.rhai> [max_ops] [timeout_ms] [dump-dir]
//
// Output: one JSON object on stdout. With a dump directory, every
// body's buffers are written raw (`body<i>.positions.f32` …) so the
// Rust side can measure deviations, not just detect them.
import { createRequire } from "node:module";
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { resolve, basename, join } from "node:path";

const [pkgDir, scriptPath, maxOps, timeoutMs, dumpDir] = process.argv.slice(2);
if (!pkgDir || !scriptPath) {
  console.error("usage: node smoke.mjs <pkg-dir> <script.rhai> [max_ops] [timeout_ms]");
  process.exit(2);
}
const require = createRequire(import.meta.url);
const pkg = require(resolve(pkgDir, "popkern_live.js"));

const FNV_OFFSET = 0xcbf29ce484222325n;
const FNV_PRIME = 0x100000001b3n;
const MASK = (1n << 64n) - 1n;
function fnv1a(bytes) {
  let hash = FNV_OFFSET;
  for (const byte of bytes) {
    hash ^= BigInt(byte);
    hash = (hash * FNV_PRIME) & MASK;
  }
  return hash.toString();
}
const bytesOf = (typed) => new Uint8Array(typed.buffer, typed.byteOffset, typed.byteLength);

const clock = pkg.init();
const options = {
  max_operations: maxOps ? Number(maxOps) : null,
  timeout_ms: timeoutMs ? Number(timeoutMs) : null,
};
const session = new pkg.LiveSession(JSON.stringify(options));
const source = readFileSync(scriptPath, "utf8");
const envelope = JSON.parse(session.run(source, basename(scriptPath)));

const bodies = [];
for (let i = 0; i < session.body_count(); i++) {
  const positions = session.body_positions(i);
  const normals = session.body_normals(i);
  const indices = session.body_indices(i);
  const faceIds = session.body_face_ids(i);
  const meta = envelope.bodies[i];
  if (dumpDir) {
    mkdirSync(dumpDir, { recursive: true });
    writeFileSync(join(dumpDir, `body${i}.positions.f32`), bytesOf(positions));
    writeFileSync(join(dumpDir, `body${i}.normals.f32`), bytesOf(normals));
    writeFileSync(join(dumpDir, `body${i}.indices.u32`), bytesOf(indices));
    writeFileSync(join(dumpDir, `body${i}.face_ids.u32`), bytesOf(faceIds));
  }
  bodies.push({
    label: meta.label,
    vertex_count: meta.vertex_count,
    triangle_count: meta.triangle_count,
    face_count: meta.face_count,
    open_edges: meta.open_edges,
    face_labels: meta.face_labels,
    positions_fnv: fnv1a(bytesOf(positions)),
    normals_fnv: fnv1a(bytesOf(normals)),
    indices_fnv: fnv1a(bytesOf(indices)),
    face_ids_fnv: fnv1a(bytesOf(faceIds)),
    positions_length: positions.length,
    typed: positions instanceof Float32Array && indices instanceof Uint32Array,
  });
}

let exports = null;
if (envelope.outcome === "served") {
  const step = session.export_step("smoke");
  const pkn = session.export_pkn();
  const rhai = session.to_rhai();
  if (dumpDir) {
    writeFileSync(join(dumpDir, "smoke.step"), step);
    writeFileSync(join(dumpDir, "smoke.pkn"), pkn);
    writeFileSync(join(dumpDir, "smoke.rhai"), rhai);
  }
  exports = {
    step_fnv: fnv1a(step),
    step_length: step.length,
    pkn_fnv: fnv1a(pkn),
    pkn_length: pkn.length,
    rhai_fnv: fnv1a(new TextEncoder().encode(rhai)),
    rhai_lines: rhai.split("\n").length,
  };
}

console.log(
  JSON.stringify({
    protocol_version: pkg.protocol_version(),
    clock_available: clock,
    outcome: envelope.outcome,
    kind: envelope.kind ?? null,
    step: envelope.step ?? null,
    wall_budget: envelope.wall_budget,
    features: envelope.features,
    elapsed_ms: envelope.elapsed_ms,
    printed: envelope.printed,
    params: envelope.params.map((p) => p.name),
    report_failures: envelope.report ? envelope.report.failures : null,
    payload_budget: envelope.payload ? envelope.payload.budget ?? null : null,
    bodies,
    exports,
  }),
);
