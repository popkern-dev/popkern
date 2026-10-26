// The Web Worker entry: one LiveSession per worker, driven by
// PROTOCOL.md messages. Runs unchanged under node's worker_threads for
// the headless harness (the transport is the only difference).
import * as pkg from "./pkg/popkern_live.js";
import { createHandler } from "./worker-core.js";

const inBrowser = typeof self !== "undefined" && typeof self.postMessage === "function" && typeof self.addEventListener === "function";

// Resolves to the module's linear memory (for the `stats` verb).
async function loadModule() {
  if (inBrowser) {
    const exports = await pkg.default();
    return exports?.memory ?? null;
  }
  // node: instantiate from the bytes beside the JS module.
  const { readFileSync } = await import("node:fs");
  const { fileURLToPath } = await import("node:url");
  const wasm = fileURLToPath(new URL("./pkg/popkern_live_bg.wasm", import.meta.url));
  const exports = pkg.initSync({ module: readFileSync(wasm) });
  return exports?.memory ?? null;
}

// node's port, resolved once (progress is posted mid-request).
const nodePort = inBrowser ? null : (await import("node:worker_threads")).parentPort;
// Progress messages go out while a request still runs.
const notify = inBrowser ? (message) => self.postMessage(message) : (message) => nodePort.postMessage(message);

const handle = createHandler(pkg, loadModule, notify);

if (inBrowser) {
  self.addEventListener("message", async (event) => {
    const { reply, transfer } = await handle(event.data);
    self.postMessage(reply, transfer);
  });
} else {
  nodePort.on("message", async (message) => {
    const { reply, transfer } = await handle(message);
    nodePort.postMessage(reply, transfer);
  });
}
