// PK-411/PK-415 browser harness: serve the playground from a local
// static server, drive headless Chrome over the DevTools protocol
// (node ≥ 22: global WebSocket + fetch, no packages), wait for the
// real kernel run, and record what the launch checklist needs —
// load timings, the served status, the report and parameter panels,
// the error overlay for a broken script and a kernel decline, a
// parameter edit, a face pick, every network request the page made,
// every console message — plus a screenshot.
//
//   node browser.mjs <screenshot.png> [chrome-binary]
//
// Output: one JSON object on stdout. Exit 3 when no browser is found.
import { spawn } from "node:child_process";
import { createServer } from "node:http";
import { readFile, stat, mkdtemp, rm } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, extname, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const [screenshotPath, chromeArg] = process.argv.slice(2);
const root = join(dirname(fileURLToPath(import.meta.url)), "..");
const MIME = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".mjs": "text/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".wasm": "application/wasm",
  ".rhai": "text/plain; charset=utf-8",
  ".png": "image/png",
  ".md": "text/markdown; charset=utf-8",
};

const requests = [];
const server = createServer(async (req, res) => {
  const url = new URL(req.url, "http://127.0.0.1");
  requests.push(url.pathname);
  let path = join(root, decodeURIComponent(url.pathname));
  try {
    if ((await stat(path)).isDirectory()) path = join(path, "index.html");
    const body = await readFile(path);
    res.writeHead(200, { "content-type": MIME[extname(path)] ?? "application/octet-stream", "cache-control": "no-store" });
    res.end(body);
  } catch {
    res.writeHead(404);
    res.end("not found");
  }
});
await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
const origin = `http://127.0.0.1:${server.address().port}`;

const candidates = [chromeArg, process.env.POPKERN_CHROME, "google-chrome", "chromium", "chromium-browser", "google-chrome-stable"].filter(Boolean);
const profile = await mkdtemp(join(tmpdir(), "popkern-playground-chrome-"));
let chrome = null;
let devtools = null;
for (const binary of candidates) {
  try {
    const port = 9222 + Math.floor(Math.random() * 1000);
    const child = spawn(
      binary,
      [
        "--headless=new",
        "--no-sandbox",
        "--disable-gpu",
        "--use-gl=angle",
        "--use-angle=swiftshader",
        "--enable-unsafe-swiftshader",
        "--ignore-gpu-blocklist",
        "--disable-dev-shm-usage",
        "--disable-extensions",
        "--disable-background-networking",
        "--disable-component-update",
        "--no-first-run",
        "--no-default-browser-check",
        `--user-data-dir=${profile}`,
        `--remote-debugging-port=${port}`,
        "--window-size=1440,900",
        "about:blank",
      ],
      { stdio: ["ignore", "ignore", "ignore"] },
    );
    child.on("error", () => {});
    const deadline = Date.now() + 15000;
    let list = null;
    while (Date.now() < deadline && !list) {
      await new Promise((r) => setTimeout(r, 200));
      try {
        list = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json();
      } catch {
        if (child.exitCode !== null) break;
      }
    }
    if (!list) {
      child.kill();
      continue;
    }
    chrome = child;
    devtools = list.find((t) => t.type === "page") ?? list[0];
    break;
  } catch {
    // try the next binary
  }
}
if (!chrome || !devtools) {
  server.close();
  await rm(profile, { recursive: true, force: true });
  console.error("no headless browser found (tried " + candidates.join(", ") + ")");
  process.exit(3);
}

const ws = new WebSocket(devtools.webSocketDebuggerUrl);
await new Promise((resolve, reject) => {
  ws.onopen = resolve;
  ws.onerror = reject;
});
let nextId = 1;
const waiting = new Map();
const consoleLines = [];
const exceptions = [];
const networkUrls = [];
ws.onmessage = (event) => {
  const msg = JSON.parse(event.data);
  if (msg.id && waiting.has(msg.id)) {
    waiting.get(msg.id)(msg);
    waiting.delete(msg.id);
  } else if (msg.method === "Runtime.consoleAPICalled") {
    consoleLines.push(`${msg.params.type}: ${msg.params.args.map((a) => a.value ?? a.description ?? "").join(" ")}`);
  } else if (msg.method === "Runtime.exceptionThrown") {
    exceptions.push(msg.params.exceptionDetails.exception?.description ?? msg.params.exceptionDetails.text);
  } else if (msg.method === "Network.requestWillBeSent") {
    networkUrls.push(msg.params.request.url);
  }
};
function send(method, params = {}) {
  const id = nextId++;
  return new Promise((resolve) => {
    waiting.set(id, resolve);
    ws.send(JSON.stringify({ id, method, params }));
  });
}
async function evaluate(expression) {
  const reply = await send("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true });
  if (reply.result?.exceptionDetails) throw new Error(reply.result.exceptionDetails.exception?.description ?? "evaluate failed");
  return reply.result?.result?.value;
}
async function waitFor(expression, timeoutMs, label) {
  const started = Date.now();
  while (Date.now() - started < timeoutMs) {
    const value = await evaluate(expression);
    if (value) return value;
    await new Promise((r) => setTimeout(r, 100));
  }
  throw new Error(`timed out waiting for ${label}`);
}
const statusText = () => evaluate(`document.getElementById("status").textContent`);

const out = { origin };
try {
  await send("Runtime.enable");
  await send("Page.enable");
  await send("Network.enable");
  const navStarted = Date.now();
  await send("Page.navigate", { url: `${origin}/` });
  out.status_after_first_run = await waitFor(
    `(() => { const s = document.getElementById("status").textContent; return /^(served|script error|kernel decline|gate failure|page error|worker error|kernel failed)/.test(s) ? s : null; })()`,
    120000,
    "the first run",
  );
  out.wall_to_first_status_ms = Date.now() - navStarted;
  out.kernel_ready_ms = await evaluate(`window.playground.kernelReadyMs ?? null`);
  out.first_render_ms = await evaluate(`window.playground.firstRenderMs ?? null`);
  out.report_summary = await evaluate(`document.getElementById("report-summary").textContent`);
  out.report_rows = await evaluate(`document.querySelectorAll("#report tr").length`);
  out.param_controls = await evaluate(`document.querySelectorAll("#params .param").length`);
  out.dialog_hidden = await evaluate(`getComputedStyle(document.getElementById("code-dialog")).display === "none"`);
  out.overlay_hidden = await evaluate(`getComputedStyle(document.getElementById("overlay")).display === "none"`);
  out.webgl = await evaluate(`!!document.getElementById("viewer").getContext("webgl2") || !!document.getElementById("viewer").getContext("webgl")`);
  out.bodies_on_screen = await evaluate(`window.playground.viewer.bodies.length`);

  // Face pick: click the centre of the canvas.
  const box = await evaluate(`(() => { const r = document.getElementById("viewer").getBoundingClientRect(); return { x: r.left + r.width / 2, y: r.top + r.height / 2 }; })()`);
  for (const type of ["mousePressed", "mouseReleased"]) {
    await send("Input.dispatchMouseEvent", { type, x: box.x, y: box.y, button: "left", clickCount: 1 });
  }
  await new Promise((r) => setTimeout(r, 200));
  out.pick_text = await evaluate(`document.getElementById("pick").textContent`);

  // Screenshot of the served state.
  const shot = await send("Page.captureScreenshot", { format: "png" });
  if (screenshotPath && shot.result?.data) {
    const { writeFile } = await import("node:fs/promises");
    await writeFile(screenshotPath, Buffer.from(shot.result.data, "base64"));
    out.screenshot = screenshotPath;
  }

  // A parameter edit through the page's own path.
  const edited = await evaluate(`window.playground.setParam("extrude-1.distance", 0.5).then((r) => ({ outcome: r.envelope?.outcome ?? r.outcome?.outcome, value: r.envelope?.params?.find((p) => p.name === "extrude-1.distance")?.value ?? null }))`);
  out.set_param = edited;

  // A broken script: syntax error → overlay, last good model kept.
  await evaluate(`window.playground.editor.value = "let m = model(\\n"`);
  const broken = await evaluate(`window.playground.run().then((r) => ({ ok: r.ok, outcome: r.envelope?.outcome ?? r.outcome?.outcome, kind: r.envelope?.kind ?? null }))`);
  out.broken = {
    ...broken,
    overlay_visible: await evaluate(`getComputedStyle(document.getElementById("overlay")).display !== "none"`),
    overlay_class: await evaluate(`document.getElementById("overlay-class").textContent`),
    overlay_text: await evaluate(`document.getElementById("overlay-text").textContent`),
    marked_line: await evaluate(`(document.querySelector(".editor-line-error") && [...document.querySelectorAll(".editor-line")].indexOf(document.querySelector(".editor-line-error")) + 1) || null`),
    bodies_still_on_screen: await evaluate(`window.playground.viewer.bodies.length`),
    status: await statusText(),
  };

  // A kernel decline: oversized fillet → overlay with the remediation.
  await evaluate(`window.playground.editor.value = 'let m = model();\\nlet cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));\\ncube.fillet(edges(), 0.9);\\nm\\n'`);
  const decline = await evaluate(`window.playground.run().then((r) => ({ ok: r.ok, outcome: r.envelope?.outcome ?? r.outcome?.outcome, kind: r.envelope?.kind ?? null }))`);
  out.decline = {
    ...decline,
    overlay_class: await evaluate(`document.getElementById("overlay-class").textContent`),
    overlay_text: await evaluate(`document.getElementById("overlay-text").textContent`),
    marked_line: await evaluate(`(document.querySelector(".editor-line-error") && [...document.querySelectorAll(".editor-line")].indexOf(document.querySelector(".editor-line-error")) + 1) || null`),
    status: await statusText(),
  };

  // The report panel states measured continuity on every blend row
  // and flags a downgrade. The real downgrading fixture (a G2 blend on
  // a blend, 27 s native) is beyond the page budget, so the display
  // path is exercised with the row shape `report_value` emits — the
  // data path is pinned natively by popkern-rhai/tests/report_rows.rs.
  out.report_rows_text = await evaluate(`[...document.querySelectorAll("#report tr")].map((r) => r.textContent.trim()).slice(0, 8)`);
  out.downgraded_row = await evaluate(`(() => {
    window.playground.renderReport({ report: { failures: 0, downgraded_blends: 1, is_clean: false, features: [
      { label: "verticals", operation: "fillet", outcome: "Served", elapsed_ms: 12.5,
        blend: { downgraded: true, continuity: { requested: "G2", achieved: "G0" } } },
    ] } });
    const cell = document.querySelector("#report tr td:nth-child(3)");
    return { text: cell.textContent, css: cell.className, summary: document.getElementById("report-summary").textContent };
  })()`);

  // Exports through the page's host (a download click cannot be
  // observed headlessly; the bytes are what the click would save).
  await evaluate(`window.playground.editor.value = ${JSON.stringify(await readFile(join(root, "..", "..", "crates/popkern-rhai/examples/boss_plate.rhai"), "utf8"))}`);
  await evaluate(`window.playground.run().then(() => true)`);
  // The bytes a download click would save, handed back for the Rust
  // side to run through the OCC oracle / reopen / re-run.
  const toB64 = `(bytes) => btoa(Array.from(bytes, (b) => String.fromCharCode(b)).join(""))`;
  out.export_step_b64 = await evaluate(`window.playground.host.exportStep("model").then((r) => (r.ok ? (${toB64})(r.bytes) : null))`);
  out.export_pkn_b64 = await evaluate(`window.playground.host.exportPkn().then((r) => (r.ok ? (${toB64})(r.bytes) : null))`);
  out.show_code = await evaluate(`window.playground.host.toRhai().then((r) => (r.ok ? r.code : null))`);
  out.export_step_bytes = out.export_step_b64 ? Buffer.from(out.export_step_b64, "base64").length : -1;
  out.export_pkn_bytes = out.export_pkn_b64 ? Buffer.from(out.export_pkn_b64, "base64").length : -1;
  out.show_code_lines = out.show_code ? out.show_code.split("\n").length : -1;

  // Runaway script: the page's cancel path stops the worker.
  await evaluate(`window.playground.editor.value = "let m = model();\\nloop {}\\n"`);
  const cancelStarted = Date.now();
  await evaluate(`(async () => { const p = window.playground.run(); await new Promise((r) => setTimeout(r, 300)); await window.playground.cancel("harness cancel"); await p; return true; })()`);
  out.cancel = { elapsed_ms: Date.now() - cancelStarted, status: await statusText(), restarts: await evaluate(`window.playground.host.restarts`) };
  await evaluate(`window.playground.editor.value = ${JSON.stringify(await readFile(join(root, "..", "..", "crates/popkern-rhai/examples/boss_plate.rhai"), "utf8"))}`);
  out.after_cancel = await evaluate(`window.playground.run().then((r) => r.envelope?.outcome ?? r.outcome?.outcome)`);

  out.network = [...new Set(networkUrls)];
  out.foreign_requests = out.network.filter((u) => !u.startsWith(origin));
  out.server_requests = [...new Set(requests)];
  out.console = consoleLines.slice(0, 40);
  out.exceptions = exceptions.slice(0, 20);
} catch (error) {
  out.error = error.message;
  out.status = await statusText().catch(() => null);
  out.console = consoleLines.slice(0, 40);
  out.exceptions = exceptions.slice(0, 20);
}

console.log(JSON.stringify(out));
ws.close();
const exited = new Promise((resolve) => chrome.once("exit", resolve));
chrome.kill();
await Promise.race([exited, new Promise((r) => setTimeout(r, 5000))]);
server.close();
for (let attempt = 0; attempt < 5; attempt++) {
  try {
    await rm(profile, { recursive: true, force: true });
    break;
  } catch {
    await new Promise((r) => setTimeout(r, 300));
  }
}
process.exit(out.error ? 1 : 0);
