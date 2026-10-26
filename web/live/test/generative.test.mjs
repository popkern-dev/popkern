// PK-956: a generative design through the worker as the page drives it.
// A plain build defers the solve (the browser never solves a design on
// its own: minutes in wasm), so it declines `NeedsSolve`. Seeded with
// the designs a native run solved (a bake's `generative`), the same
// build serves from the worker's cache: the wasm key matches the native
// one bit for bit. The page's re-solve (`generative: "solve"`) of a
// design nobody solved builds it, reporting each update as a progress
// message, and the next plain build serves it from the cache. Prints
// one JSON object the Rust test asserts on.
//
//   node generative.test.mjs <script.rhai> <designs.json> <variant.rhai>
import { Worker } from "node:worker_threads";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { performance } from "node:perf_hooks";
import { WorkerHost } from "../worker-host.js";

const [scriptPath, designsPath, variantPath] = process.argv.slice(2);
const workerFile = fileURLToPath(new URL("../worker.js", import.meta.url));
const host = new WorkerHost(() => new Worker(workerFile, { type: "module" }), { max_operations: null, timeout_ms: null });
await host.start();
const source = readFileSync(scriptPath, "utf8");
const variant = readFileSync(variantPath, "utf8");
const progress = [];
host.onProgress = (update, id) => progress.push({ ...update, id });

const summary = (reply) => {
  const envelope = reply.envelope ?? {};
  const generative = (envelope.report?.features ?? []).find((f) => f.operation === "generative");
  return {
    ok: reply.ok,
    outcome: envelope.outcome ?? reply.outcome?.outcome,
    generative_kind: envelope.payload?.error?.generative_kind ?? null,
    generative: generative?.generative ?? null,
    bodies: envelope.bodies?.length ?? 0,
  };
};
const timed = async (run) => {
  const started = performance.now();
  const reply = await run();
  return { ...summary(reply), ms: performance.now() - started };
};

const out = {};
out.deferred = await timed(() => host.run(source, "generative.rhai"));
out.seed = await host.seedGenerative(readFileSync(designsPath, "utf8"));
out.seeded = await timed(() => host.run(source, "generative.rhai"));
out.solve = await timed(() => host.run(variant, "variant.rhai", { generative: "solve" }));
out.progress = {
  count: progress.length,
  first: progress[0] ?? null,
  last: progress[progress.length - 1] ?? null,
  ids: [...new Set(progress.map((p) => p.id))],
};
out.cached = await timed(() => host.run(variant, "variant.rhai"));
process.stdout.write(JSON.stringify(out));
process.exit(0);
