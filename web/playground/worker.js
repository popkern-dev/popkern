// The Web Worker entry: one Playground session per worker, driven by
// PROTOCOL.md messages. Runs unchanged under node's worker_threads for
// the headless harness (the transport is the only difference).
import * as pkg from "./pkg/popkern_playground.js";
import { createHandler } from "./worker-core.js";

const inBrowser = typeof self !== "undefined" && typeof self.postMessage === "function" && typeof self.addEventListener === "function";

async function loadModule() {
  if (inBrowser) {
    await pkg.default();
    return;
  }
  // node: instantiate from the bytes beside the JS module.
  const { readFileSync } = await import("node:fs");
  const { fileURLToPath } = await import("node:url");
  const wasm = fileURLToPath(new URL("./pkg/popkern_playground_bg.wasm", import.meta.url));
  pkg.initSync({ module: readFileSync(wasm) });
}

const handle = createHandler(pkg, loadModule);

if (inBrowser) {
  self.addEventListener("message", async (event) => {
    const { reply, transfer } = await handle(event.data);
    self.postMessage(reply, transfer);
  });
} else {
  const { parentPort } = await import("node:worker_threads");
  parentPort.on("message", async (message) => {
    const { reply, transfer } = await handle(message);
    parentPort.postMessage(reply, transfer);
  });
}
