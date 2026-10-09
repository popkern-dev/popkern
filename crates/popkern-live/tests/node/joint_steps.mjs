// PK-725 (six-axis-arm §8, §9.8): an assembly's joint slider in the
// wasm engine, as the page's worker drives it — run the script, then
// step one joint parameter and time each `set_param` reply (JSON
// included); then, optionally, one edit of a parameter that rebuilds a
// part (by `delta` from its value). Prints one JSON object.
//
//   node joint_steps.mjs <pkg-dir> <script.rhai> <param> <from> <steps> [<rebuild-param> <delta>]
import { createRequire } from "node:module";
import { readFileSync } from "node:fs";
import { resolve, basename } from "node:path";

const [pkgDir, scriptPath, param, from, steps, rebuild, delta] = process.argv.slice(2);
const require = createRequire(import.meta.url);
const pkg = require(resolve(pkgDir, "popkern_live.js"));
pkg.init();
const session = new pkg.LiveSession(JSON.stringify({ max_operations: null, timeout_ms: null }));
const t0 = performance.now();
const run = JSON.parse(session.run(readFileSync(scriptPath, "utf8"), basename(scriptPath)));
const runMs = performance.now() - t0;
const bodies = session.body_count();
const times = [];
let unchanged = true;
let served = run.outcome === "served";
for (let k = 1; k <= Number(steps); k++) {
  const t = performance.now();
  const reply = JSON.parse(session.set_param(param, Number(from) + k));
  times.push(performance.now() - t);
  unchanged &&= reply.bodies_unchanged === true;
  served &&= reply.outcome === "served";
}
times.sort((a, b) => a - b);
let rebuilt = null;
if (rebuild) {
  const value = run.params.find((p) => p.name === rebuild)?.value;
  const t = performance.now();
  const reply = JSON.parse(session.set_param(rebuild, value + Number(delta)));
  rebuilt = { ms: performance.now() - t, outcome: reply.outcome, bodies_unchanged: reply.bodies_unchanged };
}
console.log(JSON.stringify({
  outcome: run.outcome,
  run_ms: runMs,
  served,
  bodies_unchanged: unchanged,
  bodies_before: bodies,
  bodies_after: session.body_count(),
  fastest_ms: times[0],
  median_ms: times[Math.floor(times.length / 2)],
  max_ms: times[times.length - 1],
  rebuilt,
}));
