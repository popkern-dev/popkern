// PK-411 headless protocol harness: drives worker.js through
// WorkerHost under node's worker_threads exactly as the page does —
// init, run, set_param, exports, show-as-code, and cancel of a script
// that ignores its budget. Prints one JSON object the Rust test
// compares against the native engine.
//
//   node host.test.mjs <boss_plate.rhai>
import { Worker } from "node:worker_threads";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { performance } from "node:perf_hooks";
import { WorkerHost } from "../worker-host.js";

const [scriptPath] = process.argv.slice(2);
const workerFile = fileURLToPath(new URL("../worker.js", import.meta.url));
const spawn = () => new Worker(workerFile, { type: "module" });

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

const out = {};
const host = new WorkerHost(spawn, { max_operations: null, timeout_ms: null });
const init = await host.start();
out.init = { protocol_version: init.protocol_version, clock_available: init.clock_available };

const source = readFileSync(scriptPath, "utf8");
const run = await host.run(source, "boss_plate.rhai");
out.run = {
  ok: run.ok,
  outcome: run.envelope.outcome,
  features: run.envelope.features,
  bodies: run.envelope.bodies.length,
  transferred: run.buffers.map((b) => b.positions instanceof Float32Array && b.indices instanceof Uint32Array),
  indices_fnv: run.buffers.map((b) => fnv1a(bytesOf(b.indices))),
  face_ids_fnv: run.buffers.map((b) => fnv1a(bytesOf(b.face_ids))),
  face_labels: run.envelope.bodies.map((b) => b.face_labels),
  params: run.envelope.params.length,
};

// Zebra analysis mesh (PK-577): the same bodies at the 6° bound.
const analysis = await host.analysisMesh();
out.analysis = {
  ok: analysis.ok,
  facet_deg: analysis.stats?.facet_deg ?? null,
  triangles: analysis.stats?.triangles ?? null,
  buffer_triangles: (analysis.buffers ?? []).reduce((n, b) => n + b.indices.length / 3, 0),
  display_triangles: run.buffers.reduce((n, b) => n + b.indices.length / 3, 0),
  indices_fnv: (analysis.buffers ?? []).map((b) => fnv1a(bytesOf(b.indices))),
};

const edited = await host.setParam("extrude-1.distance", 0.5);
out.set_param = {
  outcome: edited.envelope.outcome,
  value: edited.envelope.params.find((p) => p.name === "extrude-1.distance")?.value ?? null,
  bodies: edited.envelope.bodies.length,
};
const bad = await host.setParam("nope", 1.0);
out.bad_param = { outcome: bad.envelope.outcome, kind: bad.envelope.kind ?? null };

const step = await host.exportStep("boss");
const pkn = await host.exportPkn();
const code = await host.toRhai();
out.exports = {
  step_ok: step.ok && step.bytes instanceof Uint8Array && step.bytes.length > 0,
  pkn_fnv: pkn.ok ? fnv1a(pkn.bytes) : null,
  rhai_fnv: code.ok ? fnv1a(new TextEncoder().encode(code.code)) : null,
};

// The page's time limit: set_budget changes the budgets of later runs
// without touching the model (an operation budget stops a script loop;
// lifting it lets the same loop finish).
const counting = "let m = model();\nlet n = 0;\nfor i in 0..20000 { n += i; }\nm\n";
const tight = await host.setBudget({ max_operations: 1000, timeout_ms: null });
const stopped = await host.run(counting, "counting.rhai");
const lifted = await host.setBudget({ max_operations: null, timeout_ms: null });
const finished = await host.run(counting, "counting.rhai");
out.set_budget = {
  ok: tight.ok && lifted.ok,
  echoed: tight.budget ?? null,
  stopped_kind: stopped.envelope?.kind ?? null,
  finished_outcome: finished.envelope?.outcome ?? null,
};

// A malformed request has the one error shape.
const bogus = await host.request({ type: "frobnicate" });
out.bogus = bogus.outcome;

// A script that ignores its (absent) budget: the page cancels, the
// worker is replaced, and the next run answers.
const budget = 400;
const started = performance.now();
const runaway = host.run("let m = model();\nloop {}\n", "runaway.rhai");
let runawayResolvedAt = null;
runaway.then(() => {
  runawayResolvedAt = performance.now();
});
await new Promise((resolve) => setTimeout(resolve, budget));
const cancelStarted = performance.now();
const restarted = host.cancel("wall budget defeated");
const runawayReply = await runaway;
const cancelled = await restarted;
const restartedAt = performance.now();
const recovered = await host.run(source, "boss_plate.rhai");
out.cancel = {
  budget_ms: budget,
  cancel_ms: runawayResolvedAt - cancelStarted,
  restart_ms: restartedAt - cancelStarted,
  total_ms: performance.now() - started,
  runaway_outcome: runawayReply.outcome?.outcome ?? null,
  runaway_step: runawayReply.outcome?.step ?? null,
  restarted: cancelled.ok === true && host.restarts === 1,
  recovered: recovered.ok && recovered.envelope.outcome === "served",
};

await host.worker.terminate();
console.log(JSON.stringify(out));
