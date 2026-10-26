// PK-956: a generative design through the worker as the page drives it.
// A plain build defers the solve (the browser never solves a design on
// its own: minutes in wasm), so it declines `NeedsSolve`. Seeded with
// the designs a native run solved (a bake's `generative`), the same
// build serves from the worker's cache: the wasm key matches the native
// one bit for bit. The page's re-solve (`generative: "solve"`) of a
// design nobody solved builds it, reporting each update as a progress
// message, and the next plain build serves it from the cache. A film
// of the design converging (PK-959) finds its history in the worker
// once the bake's design file is seeded (the re-solved variant's
// without), and meshes it as the native session does. Prints one JSON
// object the Rust test asserts on.
//
//   node generative.test.mjs <script.rhai> <designs.json> <variant.rhai> <designs.bin>
import { Worker } from "node:worker_threads";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { performance } from "node:perf_hooks";
import { WorkerHost } from "../worker-host.js";

const [scriptPath, designsPath, variantPath, historiesPath] = process.argv.slice(2);
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

// The design converging: its history, meshes among its records, and
// a film that evolves it.
const design = out.seeded.generative?.design ?? "";
const meshSummary = (reply) => {
  if (!reply.ok) return { ok: false };
  let sum = 0;
  for (const v of reply.mesh.positions) sum += v;
  return { ok: true, triangles: reply.mesh.indices.length / 3, vertices: reply.mesh.positions.length / 3, sum };
};
out.history_missing = (await host.designHistory(design)).ok;
const bytes = readFileSync(historiesPath);
out.history_seed = await host.seedDesignHistory(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength));
const history = await host.designHistory(design);
out.history = history.ok ? { updates: history.history.updates, frames: history.history.frames.length } : null;
const last = (out.history?.frames ?? 1) - 1;
out.meshes = {
  first: meshSummary(await host.designMesh(design, 0)),
  middle: meshSummary(await host.designMesh(design, last / 2)),
  last: meshSummary(await host.designMesh(design, last)),
  space: meshSummary(await host.designMesh(design, 0, true)),
};
out.variant_history = (await host.designHistory(out.solve.generative?.design ?? "")).ok;
await host.run(source, "generative.rhai");
const film = await host.request({
  type: "film",
  source: 'let film = storyboard();\nfilm.step("Carve", 4.0, |s| { s.show("body 0", false); s.evolve("carve"); });\nfilm',
  name: "films/carve.rhai",
});
out.film = film.ok
  ? { designs: film.context.designs, evolve: film.film.tracks.filter((t) => t.channel.kind === "evolve").map((t) => t.channel) }
  : { outcome: film.outcome };
process.stdout.write(JSON.stringify(out));
process.exit(0);
