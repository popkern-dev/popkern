// PK-411/PK-415 browser harness: serve popkern live from a local
// static server, drive headless Chrome over the DevTools protocol
// (node ≥ 22: global WebSocket + fetch, no packages), wait for the
// real kernel run, and record what the launch checklist needs —
// load timings, the served status, the report and parameter panels,
// the error overlay for a broken script and a kernel decline, a
// parameter edit, a face pick, every network request the page made,
// every console message — plus a screenshot.
//
//   node browser.mjs <screenshot.png> [chrome-binary]
//   POPKERN_LAUNCH=1 node browser.mjs <screenshot.png>   # + the launch checklist rows:
//     throttled cold load (Fast 3G), the three runaway classes under
//     small budgets, a 200-feature document's memory, keyboard reach
//   POPKERN_PAGE_ORIGIN=<url> node browser.mjs <screenshot.png>
//     drive an already-served copy of the page (e.g. the sandbox
//     /routes/ deploy) instead of this harness's own static server;
//     foreign-request accounting is relative to that origin.
//
// Output: one JSON object on stdout. Exit 3 when no browser is found.
import { spawn } from "node:child_process";
import { createServer } from "node:http";
import { readFile, stat, mkdtemp, rm } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, extname, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const [screenshotPath, chromeArg] = process.argv.slice(2);
// The assembly rows' document: four pins and a two-pin sub-assembly,
// one part definition with a radius parameter (fast to build live).
const ASSEMBLY_SCRIPT = `let a = assembly("pins");
let pin = a.part("pin", |m| {
    m.param("radius", 0.1, "m");
    m.profile(xy()).circle(p("radius")).sketch("pin").extrude(blind(1.0));
});
for i in 0..4 { a.place(pin, "pin_" + i).translate(i * 0.5, 0.0, 0.0); }
let pair = a.subassembly("pair", |s| {
    s.place(pin, "low");
    s.place(pin, "high").translate(0.0, 0.0, 1.5);
});
a.place(pair).translate(0.0, 2.0, 0.0);
a
`;
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
// An externally served copy (the sandbox /routes/ deploy) takes over as
// the page under test; the local server stays up but unvisited, so
// server_requests being empty distinguishes the two modes in the JSON.
const external = (process.env.POPKERN_PAGE_ORIGIN ?? "").replace(/\/+$/, "");
const origin = external || `http://127.0.0.1:${server.address().port}`;

const candidates = [chromeArg, process.env.POPKERN_CHROME, "google-chrome", "chromium", "chromium-browser", "google-chrome-stable"].filter(Boolean);
const profile = await mkdtemp(join(tmpdir(), "popkern-live-chrome-"));
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
// `stallMs` is how long the page may show the same status before the
// wait fails (PK-828): suite load slows a run, it does not stop the
// status advancing, so no fixed wall-clock cap is load-proof.
async function waitFor(expression, stallMs, label) {
  return waitForProgress(expression, `document.getElementById("status")?.textContent ?? ""`, stallMs, 900000, label);
}
// A wait that tolerates load but not a stall (PK-784): it keeps
// waiting while `progress` (an expression read each poll) keeps
// changing, fails when it has not changed for `stallMs`, and caps the
// whole wait at `capMs`; the error names the last progress seen. A
// fixed wall-clock deadline failed the certification wait under full
// test-suite load while the page was still advancing seam by seam.
async function waitForProgress(expression, progress, stallMs, capMs, label) {
  const started = Date.now();
  let last = null;
  let changed = Date.now();
  while (Date.now() - started < capMs) {
    const value = await evaluate(expression);
    if (value) return value;
    const now = await evaluate(progress);
    if (now !== last) {
      last = now;
      changed = Date.now();
    } else if (Date.now() - changed > stallMs) {
      throw new Error(`stalled waiting for ${label} at ${JSON.stringify(last)}`);
    }
    await new Promise((r) => setTimeout(r, 100));
  }
  throw new Error(`timed out waiting for ${label} at ${JSON.stringify(last)}`);
}
// Polls read the status null-safely: under load the first polls can
// land while the navigation still shows about:blank (no #status yet).
const statusText = () => evaluate(`document.getElementById("status")?.textContent ?? ""`);

const out = { origin };
try {
  await send("Runtime.enable");
  await send("Page.enable");
  await send("Network.enable");
  const navStarted = Date.now();
  await send("Page.navigate", { url: `${origin}/` });
  out.status_after_first_run = await waitFor(
    `(() => { const s = document.getElementById("status")?.textContent ?? ""; return /^(served|script error|kernel decline|gate failure|page error|worker error|kernel failed)/.test(s) ? s : null; })()`,
    120000,
    "the first run",
  );
  out.wall_to_first_status_ms = Date.now() - navStarted;
  out.kernel_ready_ms = await evaluate(`window.popkernLive.kernelReadyMs ?? null`);
  out.first_render_ms = await evaluate(`window.popkernLive.firstRenderMs ?? null`);
  out.report_summary = await evaluate(`document.getElementById("report-summary").textContent`);
  out.report_rows = await evaluate(`document.querySelectorAll("#report tr").length`);
  out.param_controls = await evaluate(`document.querySelectorAll("#params .param").length`);
  out.dialog_hidden = await evaluate(`getComputedStyle(document.getElementById("code-dialog")).display === "none"`);
  out.overlay_hidden = await evaluate(`getComputedStyle(document.getElementById("overlay")).display === "none"`);
  out.webgl = await evaluate(`!!document.getElementById("viewer").getContext("webgl2") || !!document.getElementById("viewer").getContext("webgl")`);
  out.bodies_on_screen = await evaluate(`window.popkernLive.viewer.bodies.length`);
  out.first_run_live = await evaluate(`window.popkernLive.live`);
  out.host = await evaluate(`document.body.dataset.host ?? null`);
  // The default script is pre-built where the bundle was baked; build
  // it live too, so the kernel run and its deferred certification are
  // what the rows below measure.
  if (!out.first_run_live) {
    await waitFor(`window.popkernLive.host.info ? true : null`, 120000, "the kernel");
    const live = await evaluate(`window.popkernLive.run({ live: true }).then((r) => r?.envelope?.outcome ?? null)`);
    out.live_run = { outcome: live, status: await statusText() };
  }
  // Deferred continuity certification (PK-561): the page certifies the
  // served fillet's seams after render and swaps in the certified rows.
  const certifyStarted = Date.now();
  out.status_after_certify = await waitForProgress(
    `(() => { const s = document.getElementById("status")?.textContent ?? ""; return /continuity (certified|certification stopped)/.test(s) ? s : null; })()`,
    `document.getElementById("status")?.textContent ?? ""`,
    120000,
    900000,
    "continuity certification",
  );
  out.certify_wall_ms = Date.now() - certifyStarted;
  out.certified_rows_text = await evaluate(`[...document.querySelectorAll("#report tr")].map((r) => r.textContent.trim()).slice(0, 8)`);

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

  // History scrubber (PK-663): first, middle and last steps of the live
  // boss plate, a visited step redrawn from cache, play to the end, a
  // deep link; the camera never refits while scrubbing.
  out.scrubber = await evaluate(`(async () => {
    const page = window.popkernLive;
    const s = page.scrubber;
    // Load-proof polls (PK-828): wait while the watched value advances,
    // fail only when it stalls.
    const until = async (done, progress, stallMs = 60000) => {
      let last = progress(), changed = performance.now();
      while (!done()) {
        await new Promise((r) => setTimeout(r, 50));
        const now = progress();
        if (now !== last) { last = now; changed = performance.now(); }
        else if (performance.now() - changed > stallMs) throw new Error("scrubber stalled at " + JSON.stringify(now));
      }
    };
    await until(() => s.steps.length > 0, () => s.steps.length);
    const camera = () => page.viewer.camera.position.toArray().map((v) => v.toFixed(6)).join(",");
    const before = camera();
    const ticks = document.querySelectorAll("#scrubber-ticks .tick").length;
    const visible = !document.getElementById("scrubber").hidden;
    await s.go(1);
    const first = { bodies: page.viewer.bodies.length, label: document.getElementById("scrubber-label").textContent, line: page.editor.stepLine };
    await s.go(3);
    const middle = { bodies: page.viewer.bodies.length, glows: page.viewer.glows.length, row: document.querySelector("#report tr.step-current")?.dataset.feature ?? null, hash: location.hash };
    await s.go(5);
    // A visited step redraws from the cache: no kernel request (the
    // redraw time is reported, not gated — load can stretch it).
    const request = page.host.request;
    let revisit_requests = 0;
    page.host.request = (...args) => { revisit_requests += 1; return request.apply(page.host, args); };
    const revisit = await s.go(3);
    page.host.request = request;
    await s.go(s.steps.length);
    const last = { bodies: page.viewer.bodies.length, scrubbed: s.scrubbed, hash: location.hash };
    s.play();
    await until(() => !s.playing, () => s.at);
    const played = { at: s.at, playing: Boolean(s.playing) };
    location.hash = "step=2";
    await s.load(page.lastGood.envelope);
    const linked = s.at;
    await s.go(s.steps.length);
    return { visible, ticks, first, middle, revisit_ms: revisit, revisit_requests, last, played, linked, camera_kept: camera() === before };
  })()`);

  // Edit-to-refresh latency (PK-827): a parameter of the first feature
  // replays the whole tree; timed from the call to the frame after the
  // reply is drawn, the kernel's own share from the envelope. Fastest
  // of six (suite load only adds).
  out.edit_latency = await evaluate(`(async () => {
    const page = window.popkernLive;
    const frame = () => new Promise((r) => requestAnimationFrame(() => r()));
    const samples = [];
    for (const value of [0.45, 0.5, 0.45, 0.5, 0.45, 0.5]) {
      const t0 = performance.now();
      const reply = await page.setParam("extrude-1.distance", value);
      await frame();
      const total = performance.now() - t0;
      samples.push({ total, kernel: reply.envelope?.elapsed_ms ?? null, outcome: reply.envelope?.outcome ?? null });
    }
    const best = samples.reduce((a, b) => (b.total < a.total ? b : a));
    return { fastest_ms: best.total, kernel_ms: best.kernel, rest_ms: best.kernel == null ? null : best.total - best.kernel, outcomes: samples.map((x) => x.outcome), samples: samples.map((x) => Math.round(x.total)) };
  })()`);

  // A parameter edit through the page's own path.
  const edited = await evaluate(`window.popkernLive.setParam("extrude-1.distance", 0.5).then((r) => ({ outcome: r.envelope?.outcome ?? r.outcome?.outcome, value: r.envelope?.params?.find((p) => p.name === "extrude-1.distance")?.value ?? null }))`);
  out.set_param = edited;

  // A broken script: syntax error → overlay, last good model kept.
  await evaluate(`window.popkernLive.editor.value = "let m = model(\\n"`);
  const broken = await evaluate(`window.popkernLive.run().then((r) => ({ ok: r.ok, outcome: r.envelope?.outcome ?? r.outcome?.outcome, kind: r.envelope?.kind ?? null }))`);
  out.broken = {
    ...broken,
    overlay_visible: await evaluate(`getComputedStyle(document.getElementById("overlay")).display !== "none"`),
    overlay_class: await evaluate(`document.getElementById("overlay-class").textContent`),
    overlay_text: await evaluate(`document.getElementById("overlay-text").textContent`),
    marked_line: await evaluate(`(document.querySelector(".editor-line-error") && [...document.querySelectorAll(".editor-line")].indexOf(document.querySelector(".editor-line-error")) + 1) || null`),
    bodies_still_on_screen: await evaluate(`window.popkernLive.viewer.bodies.length`),
    status: await statusText(),
  };

  // A kernel decline: oversized fillet → overlay with the remediation.
  await evaluate(`window.popkernLive.editor.value = 'let m = model();\\nlet cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));\\ncube.fillet(edges(), 0.9);\\nm\\n'`);
  const decline = await evaluate(`window.popkernLive.run().then((r) => ({ ok: r.ok, outcome: r.envelope?.outcome ?? r.outcome?.outcome, kind: r.envelope?.kind ?? null }))`);
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
    window.popkernLive.renderReport({ report: { failures: 0, downgraded_blends: 1, is_clean: false, features: [
      { label: "verticals", operation: "fillet", outcome: "Served", elapsed_ms: 12.5,
        blend: { downgraded: true, continuity: { requested: "G2", achieved: "G0" } } },
    ] } });
    const cell = document.querySelector("#report tr td:nth-child(3)");
    return { text: cell.textContent, css: cell.className, summary: document.getElementById("report-summary").textContent };
  })()`);

  // Exports through the page's host (a download click cannot be
  // observed headlessly; the bytes are what the click would save).
  await evaluate(`window.popkernLive.editor.value = ${JSON.stringify(await readFile(join(root, "..", "..", "crates/popkern-rhai/examples/boss_plate.rhai"), "utf8"))}`);
  // live: the boss plate is also a pre-built example, and the exports
  // read the worker's model.
  await evaluate(`window.popkernLive.run({ live: true }).then(() => true)`);
  // The bytes a download click would save, handed back for the Rust
  // side to run through the OCC oracle / reopen / re-run.
  const toB64 = `(bytes) => btoa(Array.from(bytes, (b) => String.fromCharCode(b)).join(""))`;
  out.export_step_b64 = await evaluate(`window.popkernLive.host.exportStep("model").then((r) => (r.ok ? (${toB64})(r.bytes) : null))`);
  out.export_pkn_b64 = await evaluate(`window.popkernLive.host.exportPkn().then((r) => (r.ok ? (${toB64})(r.bytes) : null))`);
  out.show_code = await evaluate(`window.popkernLive.host.toRhai().then((r) => (r.ok ? r.code : null))`);
  out.export_step_bytes = out.export_step_b64 ? Buffer.from(out.export_step_b64, "base64").length : -1;
  out.export_pkn_bytes = out.export_pkn_b64 ? Buffer.from(out.export_pkn_b64, "base64").length : -1;
  out.show_code_lines = out.show_code ? out.show_code.split("\n").length : -1;

  // Runaway script: the page's cancel path stops the worker.
  await evaluate(`window.popkernLive.editor.value = "let m = model();\\nloop {}\\n"`);
  const cancelStarted = Date.now();
  await evaluate(`(async () => { const p = window.popkernLive.run(); await new Promise((r) => setTimeout(r, 300)); await window.popkernLive.cancel("harness cancel"); await p; return true; })()`);
  out.cancel = { elapsed_ms: Date.now() - cancelStarted, status: await statusText(), restarts: await evaluate(`window.popkernLive.host.restarts`) };
  await evaluate(`window.popkernLive.editor.value = ${JSON.stringify(await readFile(join(root, "..", "..", "crates/popkern-rhai/examples/boss_plate.rhai"), "utf8"))}`);
  out.after_cancel = await evaluate(`window.popkernLive.run({ live: true }).then((r) => r.envelope?.outcome ?? r.outcome?.outcome)`);

  // Pre-built examples (bakes): choosing one shows its certified build
  // without the kernel; anything needing the live model rebuilds it.
  out.bakes_present = await evaluate(`fetch("./bakes/index.json").then((r) => r.ok)`);
  if (out.bakes_present) {
    const choose = async (file) => {
      const before = await statusText();
      const started = Date.now();
      await evaluate(`window.popkernLive.loadExample("${file}")`);
      const status = await waitFor(`(() => { const s = document.getElementById("status")?.textContent ?? ""; return s.startsWith("served (pre-built)") && s !== ${JSON.stringify(before)} ? s : null; })()`, 30000, `the ${file} example`);
      return { status, wall_ms: Date.now() - started, live: await evaluate(`window.popkernLive.live`) };
    };
    out.prebuilt = await choose("pyramid_apex.rhai");
    out.prebuilt.report_text = await evaluate(`[...document.querySelectorAll("#report tr")].map((r) => r.textContent.trim()).join(" | ")`);
    out.prebuilt.bodies_on_screen = await evaluate(`window.popkernLive.viewer.bodies.length`);
    out.prebuilt.read_only = await evaluate(`window.popkernLive.editor.readOnly && !document.getElementById("readonly-bar").hidden`);
    out.prebuilt_cheap = await choose("revolved_knob.rhai");
    out.prebuilt_cheap.made_live = await evaluate(`window.popkernLive.ensureLive()`);
    out.prebuilt_cheap.status_after_live = await statusText();

    // Projects (PK-662): copy the read-only example to edit it; the
    // copy saves in this browser and the example stays as it was.
    out.copy_to_edit = await evaluate(`(async () => {
      const page = window.popkernLive;
      const example = page.editor.value;
      await page.copyToEdit();
      const copy = page.projects.current;
      page.editor.textarea.value = "// my edit\\n" + example;
      page.editor.textarea.dispatchEvent(new Event("input"));
      await page.projects.flush();
      const stored = await page.projects.store.get(copy.id);
      const path = copy.path;
      await page.loadExample("revolved_knob.rhai");
      return {
        copy_editable: !copy.readOnly,
        copy_saved: stored.files[path].startsWith("// my edit"),
        example_unchanged: page.editor.value === example && page.editor.readOnly,
        projects: (await page.projects.store.list()).length,
      };
    })()`);
  }

  // A project example (three migrations) builds through run_project.
  out.project_example = await evaluate(`(async () => {
    const page = window.popkernLive;
    await page.projects.openExample("boss_plate_in_steps");
    await page.runDone;
    const reply = await page.run({ live: true });
    return {
      outcome: reply?.envelope?.outcome ?? null,
      features: reply?.envelope?.features ?? null,
      versions: (reply?.envelope?.project?.versions ?? []).map((v) => v.name),
      read_only: page.editor.readOnly,
      tab: document.getElementById("doc-name").textContent,
    };
  })()`);
  // Later rows build one script: back to Scratch.
  await evaluate(`window.popkernLive.projects.openStored("scratch").then(() => window.popkernLive.runDone)`);

  // The time limit: presets up to none, remembered, applied to the
  // worker's budgets.
  out.time_limit = await evaluate(`(async () => {
    const select = document.getElementById("time-limit");
    const options = [...select.options].map((o) => o.textContent);
    const initial = select.value;
    select.value = "none";
    select.dispatchEvent(new Event("change"));
    await new Promise((r) => setTimeout(r, 200));
    const none = { host_timeout: window.popkernLive.host.options.timeout_ms, host_ops: window.popkernLive.host.options.max_operations, stored: localStorage.getItem("popkern-live-time-limit"), status: document.getElementById("status").textContent };
    select.value = initial;
    select.dispatchEvent(new Event("change"));
    await new Promise((r) => setTimeout(r, 200));
    return { options, initial, none, restored_timeout: window.popkernLive.host.options.timeout_ms };
  })()`);

  // An assembly (PK-629): the tree, instanced drawing (one geometry per
  // part definition), hide, open-part scoping, a qualified parameter
  // edit moving every instance, instance-path picks, product-structure
  // STEP and visibility saved into the .pkn.
  out.assembly = await evaluate(`(async () => {
    const page = window.popkernLive;
    page.editor.value = ${JSON.stringify(ASSEMBLY_SCRIPT)};
    const run = await page.run({ live: true });
    const rows = () => [...document.querySelectorAll("#assembly-tree li")];
    const drawn = () => page.viewer.bodies.filter((m) => m.visible).length;
    const result = {
      outcome: run?.envelope?.outcome ?? run?.outcome?.outcome ?? null,
      status: document.getElementById("status").textContent,
      section_visible: !document.getElementById("assembly-section").hidden,
      summary: document.getElementById("assembly-summary").textContent,
      rows: rows().map((r) => r.dataset.path),
      meshes: page.viewer.bodies.length,
      geometries: new Set(page.viewer.bodies.map((m) => m.geometry)).size,
    };
    rows().find((r) => r.dataset.path === "pair_1").querySelector(".tree-eye").click();
    result.drawn_after_hide = drawn();
    result.hidden_rows = rows().filter((r) => r.classList.contains("tree-hidden")).map((r) => r.dataset.path);
    rows().find((r) => r.dataset.path === "pin_2").querySelector(".tree-open").click();
    result.open_name = document.getElementById("assembly-open-name").textContent;
    result.open_params = [...document.querySelectorAll("#params .param-name")].map((e) => e.textContent);
    result.open_report_rows = document.querySelectorAll("#report tr").length;
    // Pick the centre of pin_2 through the real hit test.
    const mesh = page.viewer.bodies.find((m) => m.userData.instance === "pin_2");
    const box = mesh.geometry.boundingBox ?? (mesh.geometry.computeBoundingBox(), mesh.geometry.boundingBox);
    const centre = box.getCenter(mesh.position.clone()).applyMatrix4(mesh.matrixWorld).project(page.viewer.camera);
    const rect = page.viewer.canvas.getBoundingClientRect();
    page.viewer.pick({ clientX: rect.left + (centre.x + 1) / 2 * rect.width, clientY: rect.top + (1 - centre.y) / 2 * rect.height });
    result.pick = document.getElementById("pick").textContent;
    const before = page.viewer.bodies.map((m) => m.userData.body.bounds[1][0]);
    const edited = await page.setParam("pin:radius", 0.2);
    result.param_outcome = edited?.envelope?.outcome ?? edited?.outcome?.outcome ?? null;
    result.radius_grew_everywhere = page.viewer.bodies.every((m, i) => m.userData.body.bounds[1][0] > before[i] - 1e-12) && page.viewer.bodies[0].userData.body.bounds[1][0] > before[0];
    result.drawn_after_edit = drawn();
    const step = await page.host.exportStep("pins");
    const text = step.ok ? new TextDecoder().decode(step.bytes) : "";
    result.step_occurrences = text.split("=NEXT_ASSEMBLY_USAGE_OCCURRENCE(").length - 1;
    result.step_solids = text.split("=MANIFOLD_SOLID_BREP(").length - 1;
    // Clearance (PK-630): after the edit the pins are radius 0.2 at a
    // 0.5 pitch — three neighbouring pairs 0.1 apart.
    document.getElementById("clearance").value = "0.2";
    const check = await page.checkClearance();
    result.clearance = check?.ok ? check.report.pairs.map((p) => [p.a, p.b, p.status, Math.round(p.distance * 1e9) / 1e9]) : check;
    result.clearance_rows = document.querySelectorAll("#clearance-list li").length;
    // Exploded view (PK-630): transforms only — the per-frame cost of
    // moving every occurrence and drawing, over a sweep of the slider.
    const slider = document.getElementById("explode");
    result.explode_visible = !document.getElementById("explode-row").hidden;
    const frames = [];
    for (let k = 0; k <= 60; k++) {
      const t0 = performance.now();
      slider.value = String(k / 60);
      slider.dispatchEvent(new Event("input"));
      frames.push(performance.now() - t0);
    }
    frames.sort((x, y) => x - y);
    result.explode_frame_ms = { median: frames[30], p95: frames[57], max: frames[60] };
    const moved = page.viewer.bodies.find((m) => m.userData.instance === "pin_0");
    result.explode_moved = moved.matrix.elements[12] !== moved.userData.base.elements[12] || moved.matrix.elements[13] !== moved.userData.base.elements[13] || moved.matrix.elements[14] !== moved.userData.base.elements[14];
    slider.value = "0";
    slider.dispatchEvent(new Event("input"));
    await page.host.setVisible("pair_1", false);
    const pkn = await page.host.exportPkn();
    const doc = pkn.ok ? JSON.parse(new TextDecoder().decode(pkn.bytes)) : null;
    result.pkn_version = doc?.version ?? null;
    result.pkn_hidden = doc?.assembly?.root?.filter((i) => !i.visible).map((i) => i.name) ?? null;
    document.getElementById("assembly-close").click();
    return result;
  })()`);

  // Assembly history (PK-666): the construction steps scrub — the first
  // step draws one occurrence glowing, the last is the live assembly.
  out.assembly_history = await evaluate(`(async () => {
    const page = window.popkernLive;
    const s = page.scrubber;
    await page.run({ live: true });
    await s.loading;
    const visible = !document.getElementById("scrubber").hidden;
    const kinds = s.steps.map((x) => x.kind);
    await s.go(1);
    const first = {
      meshes: page.viewer.bodies.length,
      glows: page.viewer.glows.length,
      label: document.getElementById("scrubber-label").textContent,
      insert_hidden: document.getElementById("scrubber-insert").hidden,
    };
    await s.go(s.steps.length);
    const last = { meshes: page.viewer.bodies.length, scrubbed: s.scrubbed };
    return { visible, kinds, first, last, recorded: s.assembly?.recorded ?? null };
  })()`);

  // Open STEP (PK-631): a third-party writer's assembly (OpenCASCADE
  // XCAF, tools/occ-oracle/occ_write_assembly.py) opens with its tree,
  // names and shared parts.
  const stepFixture = (await readFile(join(root, "..", "..", "crates/popkern-assembly/tests/fixtures/occ_bench_assembly.step"))).toString("base64");
  // Part appearance (six-axis-arm §4.5): a part's colour and finish
  // are its CAD-view material (with the studio environment so a metal
  // reads as metal); a part without one keeps the default CAD face.
  out.appearance = await evaluate(`(async () => {
    const page = window.popkernLive;
    page.editor.value = ${JSON.stringify(`let a = assembly("colours");
let blue = a.part("housing", |m| { m.profile(xy()).circle(0.3).sketch("s").extrude(blind(0.5)); m.appearance(#{ color: "#2b5cab", metal: 0.6, rough: 0.35 }); });
let plain = a.part("plain", |m| { m.profile(xy()).rect(0.2, 0.2).sketch("s").extrude(blind(0.2)); });
a.place(blue, "housing_1");
a.place(plain, "plain_1").translate(1.0, 0.0, 0.0);
a`)};
    const run = await page.run({ live: true });
    const material = (instance) => {
      const m = page.viewer.bodies.find((b) => b.userData.instance === instance).material;
      return { color: m.color.getHexString(), metal: m.metalness, rough: m.roughness, env: !!m.envMap };
    };
    return {
      outcome: run?.envelope?.outcome ?? run?.outcome?.outcome ?? null,
      housing: material("housing_1"),
      plain: material("plain_1"),
    };
  })()`);

  // The section view (six-axis-arm §4.6): a clipping plane on every
  // material, one stencil-capped cut per visible body — the caps draw
  // real pixels — moving the cut is a transform, and "off" restores.
  out.section = await evaluate(`(async () => {
    const page = window.popkernLive;
    const viewer = page.viewer;
    const select = document.getElementById("section-axis");
    // 0.2 of the bounds puts the cut through the housing (the other
    // part lies wholly beyond it and is removed).
    document.getElementById("section-at").value = "0.2";
    select.value = "x";
    select.dispatchEvent(new Event("change"));
    const gl = viewer.renderer.getContext();
    const pixels = () => {
      viewer.renderer.render(viewer.scene, viewer.camera);
      const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight;
      const buf = new Uint8Array(w * h * 4);
      gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, buf);
      return buf;
    };
    const visible = viewer.bodies.filter((m) => m.visible).length;
    const crossed = (viewer, plane) => viewer.bodies.filter((m) => {
      if (!m.visible) return false;
      const [lo, hi] = m.userData.body.bounds;
      const box = new viewer.camera.position.constructor(0, 0, 0);
      let min = Infinity, max = -Infinity;
      for (let k = 0; k < 8; k++) {
        box.set(k & 1 ? hi[0] : lo[0], k & 2 ? hi[1] : lo[1], k & 4 ? hi[2] : lo[2]).applyMatrix4(m.matrixWorld);
        const d = box.x * plane.normal[0] + box.y * plane.normal[1] + box.z * plane.normal[2];
        min = Math.min(min, d); max = Math.max(max, d);
      }
      return min < plane.offset && plane.offset < max;
    }).length;
    const result = {
      slider_shown: !document.getElementById("section-at").hidden,
      caps: viewer.section.caps.children.length,
      visible,
      crossed: crossed(viewer, viewer.section.planes[0]),
      clipped: viewer.bodies.every((m) => (m.material.clippingPlanes ?? []).length === 1),
    };
    // Look at the cut face: from the removed (+x) side.
    const [lo, hi] = viewer.visibleBounds();
    const centre = lo.map((v, i) => (v + hi[i]) / 2);
    const span = Math.hypot(...hi.map((v, i) => v - lo[i]));
    const saved = viewer.camera.position.clone();
    const savedUp = viewer.camera.up.clone();
    viewer.camera.up.set(0, 0, 1);
    viewer.camera.position.set(centre[0] + 2 * span, centre[1] + 0.3 * span, centre[2] + 0.3 * span);
    viewer.camera.lookAt(...centre);
    viewer.camera.updateMatrixWorld();
    const withCaps = pixels();
    viewer.section.caps.visible = false;
    const without = pixels();
    viewer.section.caps.visible = true;
    viewer.camera.position.copy(saved);
    viewer.camera.up.copy(savedUp);
    viewer.orbit.update();
    let changed = 0;
    for (let i = 0; i < withCaps.length; i += 4) if (withCaps[i] !== without[i] || withCaps[i + 1] !== without[i + 1]) changed++;
    result.cap_pixels = changed;
    // The fastest move is its cost: suite load only adds to a sample
    // (a 20-move mean read 16.7 ms beside the whole suite, PK-874).
    const moves = [];
    for (let k = 0; k < 20; k++) {
      const t0 = performance.now();
      viewer.section.move(k / 20);
      moves.push(performance.now() - t0);
    }
    result.move_ms = Math.min(...moves);
    select.value = "";
    select.dispatchEvent(new Event("change"));
    result.cleared = viewer.section.caps.children.length === 0 && viewer.bodies.every((m) => !(m.material.clippingPlanes ?? []).length);
    return result;
  })()`);

  // The clipping view (docs/architecture/clipping-view.md §9): a tilted
  // plane clips and caps, flip changes what is drawn, two planes in
  // "all" remove only their corner, an excluded occurrence stays whole,
  // the gizmo's arrow drag slides the plane without orbiting, the link
  // restores the cut, and clear restores every material.
  out.clipping = await evaluate(`(async () => {
    const page = window.popkernLive;
    const viewer = page.viewer;
    const section = viewer.section;
    const gl = viewer.renderer.getContext();
    const pixels = () => {
      viewer.renderer.render(viewer.scene, viewer.camera);
      const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight;
      const buf = new Uint8Array(w * h * 4);
      gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, buf);
      return buf;
    };
    const differ = (a, b) => { let n = 0; for (let i = 0; i < a.length; i += 4) if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 3] !== b[i + 3]) n++; return n; };
    const [lo, hi] = viewer.visibleBounds();
    const centre = lo.map((v, i) => (v + hi[i]) / 2);
    const span = Math.hypot(...hi.map((v, i) => v - lo[i]));
    const housing = viewer.bodies.find((m) => m.userData.instance?.startsWith("housing"));
    const plain = viewer.bodies.find((m) => m.userData.instance?.startsWith("plain"));
    const [hlo, hhi] = housing.userData.body.bounds;
    const hc = hlo.map((v, i) => (v + hhi[i]) / 2);
    const look = (dir) => {
      viewer.camera.up.set(0, 0, 1);
      viewer.camera.position.set(hc[0] + 3 * span * dir[0], hc[1] + 3 * span * dir[1], hc[2] + 3 * span * dir[2]);
      viewer.camera.lookAt(...hc);
      viewer.camera.updateMatrixWorld();
    };
    const saved = viewer.camera.position.clone();
    const savedUp = viewer.camera.up.clone();
    const result = {};
    // A tilted plane through the housing's centre.
    const n = [0.6, 0, 0.8];
    section.setPlanes([{ normal: n, offset: n[0] * hc[0] + n[1] * hc[1] + n[2] * hc[2] }]);
    const visible = viewer.bodies.filter((m) => m.visible).length;
    const crossed = (viewer, plane) => viewer.bodies.filter((m) => {
      if (!m.visible) return false;
      const [lo, hi] = m.userData.body.bounds;
      const box = new viewer.camera.position.constructor(0, 0, 0);
      let min = Infinity, max = -Infinity;
      for (let k = 0; k < 8; k++) {
        box.set(k & 1 ? hi[0] : lo[0], k & 2 ? hi[1] : lo[1], k & 4 ? hi[2] : lo[2]).applyMatrix4(m.matrixWorld);
        const d = box.x * plane.normal[0] + box.y * plane.normal[1] + box.z * plane.normal[2];
        min = Math.min(min, d); max = Math.max(max, d);
      }
      return min < plane.offset && plane.offset < max;
    }).length;
    result.tilted = {
      clipped: viewer.bodies.every((m) => (m.material.clippingPlanes ?? []).length === 1),
      caps: section.caps.children.length,
      visible,
      crossed: crossed(viewer, section.planes[0]),
      outline_segments: section.describe().outline_segments,
    };
    look([0.6, -0.3, 0.8]);
    const before = pixels();
    section.flip();
    const after = pixels();
    result.flip_pixels = differ(before, after);
    // Two planes through the housing's axis: x > cx and y > cy.
    section.setPlanes([
      { normal: [1, 0, 0], offset: hc[0] },
      { normal: [0, 1, 0], offset: hc[1] },
    ], { combine: "any" });
    const r = (hhi[0] - hlo[0]) / 2;
    // A surface point in the x > cx, y < cy quadrant: removed by "any",
    // kept by "all".
    const probe = [hc[0] + r * Math.SQRT1_2, hc[1] - r * Math.SQRT1_2, hc[2]];
    look([0.7, -0.7, 0.15]);
    const screen = (p) => {
      const v = new viewer.camera.position.constructor(...p).project(viewer.camera);
      const rect = viewer.canvas.getBoundingClientRect();
      return [rect.left + ((v.x + 1) / 2) * rect.width, rect.top + ((1 - v.y) / 2) * rect.height];
    };
    const hitNear = () => {
      const [x, y] = screen(probe);
      const hit = viewer.hitTest(x, y);
      if (!hit) return false;
      return Math.hypot(...hit.point.map((c, i) => c - probe[i])) < 0.05 * span;
    };
    result.any_keeps_probe = hitNear();
    section.setCombine("all");
    result.all_keeps_probe = hitNear();
    result.all_caps = section.caps.children.length;
    // Exclude the plain part: no clipping, no caps on it.
    // Through the plain part's own centre, so the cut crosses it.
    const [plo, phi] = plain.userData.body.bounds;
    const pc = new viewer.camera.position.constructor((plo[0] + phi[0]) / 2, (plo[1] + phi[1]) / 2, (plo[2] + phi[2]) / 2).applyMatrix4(plain.matrixWorld);
    section.setPlanes([{ normal: [1, 0, 0], offset: pc.x }], { combine: "any" });
    const capsBefore = section.caps.children.length;
    const crossedBefore = crossed(viewer, section.planes[0]);
    section.setExcluded(plain.userData.instance, true);
    result.excluded = {
      plain_unclipped: !(plain.material.clippingPlanes ?? []).length,
      housing_clipped: (housing.material.clippingPlanes ?? []).length === 1,
      caps_before: capsBefore,
      crossed_before: crossedBefore,
      caps_after: section.caps.children.length,
      link: section.describe().link,
      toggle_shown: document.querySelectorAll(".tree-cut:not([hidden])").length,
    };
    section.setExcluded(plain.userData.instance, false);
    // The gizmo: drag the arrow tip along the normal on screen.
    look([0.3, -1, 0.4]);
    viewer.render();
    const g = section.gizmo;
    g.group.updateMatrixWorld(true);
    const tipWorld = g.arrow.children[1].getWorldPosition(new viewer.camera.position.constructor()).toArray();
    const base = g.group.position.toArray();
    const [tx, ty] = screen(tipWorld);
    const [bx, by] = screen(base);
    const len = Math.hypot(tx - bx, ty - by);
    const ux = (tx - bx) / len, uy = (ty - by) / len;
    const offset0 = section.planes[0].offset;
    const camera0 = viewer.camera.position.toArray();
    const fire = (type, x, y) => viewer.canvas.dispatchEvent(new PointerEvent(type, { pointerId: 7, pointerType: "mouse", button: 0, buttons: type === "pointerup" ? 0 : 1, clientX: x, clientY: y, bubbles: true }));
    fire("pointerdown", tx, ty);
    const grabbed = !!g.dragging;
    fire("pointermove", tx + ux * 20, ty + uy * 20);
    fire("pointermove", tx + ux * 40, ty + uy * 40);
    fire("pointerup", tx + ux * 40, ty + uy * 40);
    // 40 px along the arrow is 40 / (pixels per unit) world units.
    const tipDistance = Math.hypot(...tipWorld.map((c, i) => c - base[i]));
    const perUnit = len / tipDistance;
    result.gizmo = {
      grabbed,
      moved: section.planes[0].offset - offset0,
      expected: 40 / perUnit,
      orbited: viewer.camera.position.toArray().some((c, i) => Math.abs(c - camera0[i]) > 1e-9),
    };
    // A plane facing the viewer: the arrow is a dot, so the square
    // slides it — dragging up pushes it away from the camera.
    section.setPlanes([]);
    section.addPlane("view");
    viewer.render();
    g.group.updateMatrixWorld(true);
    const sq = g.group.position.toArray();
    const [qx, qy] = screen(sq);
    const camDistance = () => {
      const p = section.planes[0];
      const c = viewer.camera.position;
      return Math.abs(p.normal[0] * c.x + p.normal[1] * c.y + p.normal[2] * c.z - p.offset);
    };
    const d0 = camDistance();
    const fireAt = (type, x, y) => viewer.canvas.dispatchEvent(new PointerEvent(type, { pointerId: 8, pointerType: "mouse", button: 0, buttons: type === "pointerup" ? 0 : 1, clientX: x, clientY: y, bubbles: true }));
    // Off-centre on the square, clear of the arrow and the rings' hit tubes.
    fireAt("pointerdown", qx + 15, qy + 15);
    const squareGrabbed = !!g.dragging;
    fireAt("pointermove", qx + 15, qy - 15);
    fireAt("pointerup", qx + 15, qy - 15);
    result.facing = { grabbed: squareGrabbed, away: camDistance() - d0 };

    // A plane exactly on a face (the housing's top) removes nothing: the
    // face draws the same pixels as with no cut (it shimmered when its
    // fragments sat at a distance of ±0 from the clip plane).
    look([0.4, -0.5, 0.8]);
    section.setPlanes([]);
    const uncut = pixels();
    section.setPlanes([{ normal: [0, 0, 1], offset: hhi[2] }]);
    section.gizmo.setEnabled(false);
    result.on_face_pixels = differ(uncut, pixels());
    section.gizmo.setEnabled(true);

    // The viewport's section button: on (a cut facing the viewer, with
    // handles), then off.
    section.setPlanes([]);
    document.getElementById("section-toggle").click();
    result.toggle = { on: section.active, pressed: document.getElementById("section-toggle").getAttribute("aria-pressed"), handles: section.gizmo.group.visible };
    document.getElementById("section-toggle").click();
    result.toggle.off = !section.active;

    // The link restores the cut.
    section.setPlanes([{ normal: [0, 0.6, 0.8], offset: 0.1 }, { normal: [1, 0, 0], offset: hc[0] }], { combine: "all" });
    const link = new URLSearchParams(location.hash.slice(1)).get("section");
    section.clear();
    const clearedLink = new URLSearchParams(location.hash.slice(1)).get("section");
    section.pendingLink = link;
    section.apply();
    result.link = {
      written: link,
      cleared_link: clearedLink,
      restored: section.planes.length,
      combine: section.combine,
      normal: section.planes[0]?.normal,
    };
    // Outline cost on this scene, one plane.
    section.setPlanes([{ normal: n, offset: n[0] * hc[0] + n[2] * hc[2] }]);
    const outlines = [];
    for (let k = 0; k < 10; k++) {
      const t0 = performance.now();
      section.outlineNow();
      outlines.push(performance.now() - t0);
    }
    result.outline_ms = Math.min(...outlines);
    section.clear();
    result.cleared = section.caps.children.length === 0 && section.outlines.children.length === 0 &&
      viewer.bodies.every((m) => !(m.material.clippingPlanes ?? []).length) && !section.gizmo.group.visible;
    viewer.camera.position.copy(saved);
    viewer.camera.up.copy(savedUp);
    viewer.orbit.update();
    return result;
  })()`);

  // Presentation films (docs/architecture/presentation-films.md §10): a
  // film script runs against the document on screen (no rebuild), its
  // state moves an occurrence, the look and caption apply, and the
  // export writes MP4, WebM and GIF with the film's frame count.
  out.film = await evaluate(`(async () => {
    const page = window.popkernLive;
    const viewer = page.viewer;
    const film = page.film;
    const source = ${JSON.stringify(`let doc = document();
let film = storyboard(#{ size: [320, 180], fps: 10 });
film.title("Colours", "a housing and a plain part", 0.5, |s| { s.orbit(20.0); });
film.step("Lift", "the plain part rises", 1.0, |s| {
    s.move("plain_1", #{ to: [0.0, 0.0, 0.6] });
    s.fade("housing_1", 1.0, 0.4);
});
film`)};
    const result = {};
    const t0 = performance.now();
    const reply = await film.run(source, "films/lift.rhai");
    result.run_ms = performance.now() - t0;
    result.ok = reply.ok;
    result.active = film.active;
    result.duration = film.film?.duration ?? null;
    result.bar_shown = !document.getElementById("film-bar").hidden;
    result.mode = viewer.mode;
    const plain = viewer.bodies.find((m) => m.userData.instance?.startsWith("plain"));
    const housing = viewer.bodies.find((m) => m.userData.instance?.startsWith("housing"));
    const z = () => new viewer.camera.position.constructor().setFromMatrixPosition(plain.matrixWorld).z;
    await film.seek(0);
    const z0 = z();
    await film.seek(1.5);
    result.lifted = z() - z0;
    result.housing_opacity = housing.material.opacity;
    await film.seek(0.25);
    result.caption = film.caption ? { title: film.caption.title, opacity: film.caption.opacity } : null;
    const { exportFilm } = await import("./film-export.js");
    const header = (bytes, n) => Array.from(bytes.slice(0, n)).map((b) => String.fromCharCode(b)).join("");
    for (const key of ["mp4-720p", "webm", "gif"]) {
      try {
        const e0 = performance.now();
        const o = await exportFilm(page, key);
        result[key] = { bytes: o.bytes.length, frames: o.frames, head: header(o.bytes, 12), ms: performance.now() - e0 };
        if (key === "gif") {
          let n = 0;
          for (let i = 0; i + 2 < o.bytes.length; i++) if (o.bytes[i] === 0x21 && o.bytes[i + 1] === 0xf9 && o.bytes[i + 2] === 0x04) n++;
          result[key].gif_frames = n;
        }
      } catch (e) {
        result[key] = { error: String(e) };
      }
    }
    // A broken film reports at its line and leaves the open one alone.
    const broken = await film.run(${JSON.stringify('let film = storyboard();\nfilm.step("a", 1.0, |s| { s.move("lid", #{ to: [0.0, 0.0, 1.0] }); });\nfilm')}, "films/broken.rhai");
    result.broken = { ok: broken.ok, kind: broken.outcome?.kind ?? null, rendered: broken.outcome?.rendered ?? null, still_open: film.active };
    // The viewport's own close chip (on a phone the bar's sheet may be
    // hidden): shown while a film is open, and it closes the film.
    result.exit_shown = !document.getElementById("film-exit").hidden;
    document.getElementById("film-exit").click();
    result.exit_closed = !film.active;
    await film.run(source, "films/lift.rhai");
    // Escape closes the film (the page's own key binding).
    document.body.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true }));
    await new Promise((r) => setTimeout(r, 50));
    result.escape_closed = !film.active;
    film.close();
    result.closed = { exit_hidden: document.getElementById("film-exit").hidden, active: film.active, mode: viewer.mode, bar_hidden: document.getElementById("film-bar").hidden, lifted_back: Math.abs(z() - z0) < 1e-9 };
    return result;
  })()`);

  // A joint slider (six-axis-arm §9.8, PK-725): an assembly parameter
  // that drives no part re-places the occurrences — the page keeps its
  // geometries (the same objects), the link and its follower move, the
  // ground stays — and the round trip is timed over a drag.
  out.joint = await evaluate(`(async () => {
    const page = window.popkernLive;
    page.editor.value = ${JSON.stringify(`let a = assembly("arm");
let block = a.part("block", |m| { m.profile(xy()).rect(4.0, 1.0).sketch("s").extrude(blind(1.0)); });
let base = a.place(block, "base");
let link = a.place(block, "link").translate(1.0, 2.0, 0.0);
let rotor = a.place(block, "rotor").translate(-3.0, 0.0, 0.0);
a.ground(base);
a.revolute("j1", base.axis([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]), link.datum("z"),
    #{ zero: [base.datum("x"), link.datum("x")], angle: 20.0, limits: [-170.0, 170.0], offset: 1.0 });
a.revolute("j1_rotor", base.datum("z"), rotor.datum("z"),
    #{ zero: [base.datum("x"), rotor.datum("x")], follow: "j1", ratio: -21.0 });
a`)};
    const runStarted = performance.now();
    const run = await page.run({ live: true });
    const runMs = performance.now() - runStarted;
    const at = (instance) => page.viewer.bodies.find((b) => b.userData.instance === instance).matrix.elements.slice();
    const geometries = page.viewer.bodies.map((m) => m.geometry);
    const sliderRow = [...document.querySelectorAll("#params .param")].find((r) => r.querySelector(".param-name").textContent === "j1");
    const before = { base: at("base"), link: at("link"), rotor: at("rotor") };
    // The step's own cost (PK-874): the worker's reply and the
    // placement, without SwiftShader's frame, whose backpressure grows
    // with machine load (the fastest full step read 37.8 ms beside the
    // whole suite, 20 ms alone).
    // The worker's own rebuild of the same script, fastest of 5 under
    // the same load: the step must cost a fraction of it (a re-mesh or
    // a full describe inside the step reads as one). The page's live
    // run is no yardstick — on this three-block arm it is mostly
    // page and frame overhead, and an injected rebuild still read under
    // half of it.
    const rebuilds = [];
    for (let k = 0; k < 5; k++) {
      const t0 = performance.now();
      await page.host.run(page.editor.value, "arm.rhai");
      rebuilds.push(performance.now() - t0);
    }
    const rebuildMs = Math.min(...rebuilds);
    const steps = [];
    for (let k = 1; k <= 20; k++) {
      const t0 = performance.now();
      const reply = await page.host.setParam("j1", 60 + k);
      page.viewer.place(reply.envelope.assembly, { render: false });
      steps.push(performance.now() - t0);
    }
    const stepMs = Math.min(...steps);
    const times = [];
    let last = null;
    for (let k = 1; k <= 20; k++) {
      const t0 = performance.now();
      last = await page.setParam("j1", 20 + k);
      times.push(performance.now() - t0);
    }
    times.sort((x, y) => x - y);
    const moved = (name) => at(name).some((v, i) => Math.abs(v - before[name][i]) > 1e-6);
    return {
      outcome: run?.envelope?.outcome ?? run?.outcome?.outcome ?? null,
      placed: last?.envelope?.bodies_unchanged ?? null,
      sent_buffers: last?.buffers?.length ?? null,
      same_geometries: page.viewer.bodies.every((m, i) => m.geometry === geometries[i]),
      base_moved: moved("base"),
      link_moved: moved("link"),
      rotor_moved: moved("rotor"),
      slider_value: page.lastGood.envelope.params.find((p) => p.name === "j1")?.value ?? null,
      // The panel row is the same element, showing the new value.
      row_kept: sliderRow != null && sliderRow.isConnected,
      row_value: sliderRow ? Number(sliderRow.querySelector('input[type="number"]').value) : null,
      bodies_kept: page.lastGood.envelope.bodies.length,
      report_kept: page.lastGood.envelope.report != null,
      frame_ms: { fastest: times[0], median: times[10], max: times[19] },
      step_ms: stepMs,
      run_ms: runMs,
      rebuild_ms: rebuildMs,
    };
  })()`);

  // A large model (millimetres, two metres across) fits and draws:
  // the depth range follows the model's size.
  out.large_fit = await evaluate(`(async () => {
    const page = window.popkernLive;
    page.editor.value = ${JSON.stringify(`let m = model();
m.profile(xy()).rect(2000.0, 600.0).sketch("s").extrude(blind(400.0));
m`)};
    const run = await page.run({ live: true });
    page.viewer.fitAgain();
    const viewer = page.viewer;
    viewer.renderer.render(viewer.scene, viewer.camera);
    const gl = viewer.renderer.getContext();
    const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight;
    const buf = new Uint8Array(w * h * 4);
    gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, buf);
    const [r0, g0, b0] = [buf[0], buf[1], buf[2]];
    let drawn = 0;
    for (let i = 0; i < buf.length; i += 4) if (Math.abs(buf[i] - r0) + Math.abs(buf[i + 1] - g0) + Math.abs(buf[i + 2] - b0) > 30) drawn++;
    // The feature tree shows a structured warning by its message, and
    // a script's label as text (never markup).
    page.renderReport({ assembly: null, report: { failures: 0, downgraded_blends: 0, is_clean: false, features: [
      { label: "<b>fins</b>", operation: "fillet", outcome: "Served", warnings: [{ message: "sliver retained", payload: { edge: 7 } }] },
    ] } });
    const warnRow = document.querySelector("#report tr.row-warn td");
    const labelCell = document.querySelector("#report tr td");
    const report = { warning: warnRow?.textContent ?? null, label: labelCell?.textContent ?? null, markup: document.querySelectorAll("#report b").length };
    page.renderReport(page.lastGood.envelope);
    return { outcome: run?.envelope?.outcome ?? run?.outcome?.outcome ?? null, drawn_fraction: drawn / (w * h), far: viewer.camera.far, report };
  })()`);

  out.open_step = await evaluate(`(async () => {
    const page = window.popkernLive;
    const bytes = Uint8Array.from(atob(${JSON.stringify(stepFixture)}), (c) => c.charCodeAt(0)).buffer;
    const reply = await page.openStep("occ_bench_assembly.step", bytes);
    return {
      outcome: reply?.envelope?.outcome ?? reply?.outcome?.outcome ?? null,
      printed: reply?.envelope?.printed ?? null,
      rows: [...document.querySelectorAll("#assembly-tree li")].map((r) => r.dataset.path + "=" + r.querySelector(".tree-kind").textContent),
      geometries: new Set(page.viewer.bodies.map((m) => m.geometry)).size,
      meshes: page.viewer.bodies.length,
      tab: document.getElementById("doc-name").textContent,
      build_kept: await page.run({ live: true }).then((r) => r === null),
    };
  })()`);

  // A returning visitor with an edited script saved: the page still
  // opens pre-built (nothing builds on load) and offers the edit back.
  if (out.bakes_present) {
    // The draft is the default script edited (the editor may hold the
    // view-only placeholder of the STEP opened above). It is stored
    // under the key the page used before the popkern-live rename, so
    // this also proves the one-time key adoption.
    await evaluate(`window.popkernLive.defaultScript().then((d) => { localStorage.removeItem("popkern-live-draft"); localStorage.setItem("popkern-playground-draft", "// edited last time\\n" + d); })`);
    await send("Page.navigate", { url: `${origin}/` });
    const reopened = await waitFor(`(() => { const s = document.getElementById("status")?.textContent ?? ""; return s.startsWith("served") ? s : null; })()`, 60000, "the reopened page");
    out.returning = {
      status: reopened,
      live: await evaluate(`window.popkernLive.live`),
      notice_visible: await evaluate(`!document.getElementById("draft-notice").hidden`),
      playground_draft_key_cleared: await evaluate(`localStorage.getItem("popkern-playground-draft") === null`),
      editor_is_default: await evaluate(`window.popkernLive.defaultScript().then((d) => window.popkernLive.editor.value === d)`),
    };
    await evaluate(`document.getElementById("draft-restore").click()`);
    out.returning.restored = await waitFor(`(() => { const s = document.getElementById("status")?.textContent ?? ""; return /^served: /.test(s) ? s : null; })()`, 60000, "the restored draft");
    out.returning.editor_restored = await evaluate(`window.popkernLive.editor.value.startsWith("// edited last time")`);
    out.returning.notice_after = await evaluate(`!document.getElementById("draft-notice").hidden`);
    out.returning.scratch_open = await evaluate(`window.popkernLive.projects.current?.id === "scratch"`);

    // Another tab saves Scratch first: this tab's next save is refused
    // and the conflict bar offers both versions; nothing is overwritten.
    out.conflict = await evaluate(`(async () => {
      const page = window.popkernLive;
      const p = page.projects.current;
      await page.projects.store.save(p.id, { ...p.files, [p.path]: "// the other tab\\n" + p.files[p.path] }, p.revision);
      page.editor.textarea.value = "// this tab\\n" + page.editor.value;
      page.editor.textarea.dispatchEvent(new Event("input"));
      await page.projects.flush();
      const bar = !document.getElementById("conflict-bar").hidden;
      const stored = await page.projects.store.get(p.id);
      const theirsKept = stored.files[p.path].startsWith("// the other tab");
      await page.projects.loadTheirs();
      return { bar, theirs_kept: theirsKept, loaded_theirs: page.editor.value.startsWith("// the other tab"), bar_after: !document.getElementById("conflict-bar").hidden };
    })()`);
  }

  // Insert a feature here (PK-666): on a scrubbed-back step of an
  // editable project the scrubber offers an insert; it adds a migration
  // anchored on that step, which — once its closure records a feature —
  // builds with the new feature right after the anchor.
  out.insert_here = await evaluate(`(async () => {
    const page = window.popkernLive;
    const wait = () => new Promise((r) => setTimeout(r, 50));
    const block = 'let m = model();\\nlet block = m.profile(xy()).rect(2.0, 2.0).sketch("block_sketch").extrude(blind(1.0));\\nlet bore = m.tool(|m| m.profile(xy()).circle(0.4).sketch("bore_sketch").extrude(blind(1.0)));\\nblock.subtract(bore);\\nm.active().fillet(edges().circular().at_height(1.0), 0.05);\\nm\\n';
    await page.projects.create("Insert test", block);
    await page.runDone;
    // This project's own build (a baked page may still be finishing
    // another), and its timeline.
    await page.run({ live: true });
    const s = page.scrubber;
    await s.loading;
    const hiddenAtEnd = document.getElementById("scrubber-insert").hidden;
    const anchor = s.steps.find((x) => x.operation === "boolean");
    await s.go(anchor.index - 1);
    await s.go(anchor.index);
    const shown = !document.getElementById("scrubber-insert").hidden;
    document.getElementById("scrubber-insert").click();
    await s.inserting;
    const statusAfterInsert = document.getElementById("status").textContent;
    const path = page.projects.path;
    const template = page.editor.value;
    page.editor.textarea.value = template.replace("    body\\n", "    body.chamfer(edges().straight().at_height(1.0), 0.05)\\n");
    page.editor.textarea.dispatchEvent(new Event("input"));
    await page.projects.flush();
    clearTimeout(page.idleTimer);
    await page.runDone;
    const built = await page.run({ live: true });
    await s.loading;
    const labels = s.steps.map((x) => x.label);
    const result = {
      hidden_at_end: hiddenAtEnd,
      shown,
      path,
      anchored: template.includes('m.feature("' + anchor.label + '")'),
      statuses: (built?.envelope?.project?.versions ?? []).map((v) => v.status),
      after_anchor: labels[anchor.index] ?? null,
      status_after_insert: statusAfterInsert,
      anchor: anchor.label,
    };
    await page.projects.openStored("scratch");
    await page.runDone;
    return result;
  })()`);

  // Folders (PK-668): a project folder opened through the File System
  // Access API — an origin-private folder stands in for the user's pick
  // (headless Chrome shows no picker) — edited and saved in place; a
  // change on disk reloads a clean editor and conflicts with a dirty
  // one; forgetting keeps the files; a browser without the API says so.
  out.folder = await evaluate(`(async () => {
    const page = window.popkernLive;
    const root = await navigator.storage.getDirectory();
    try { await root.removeEntry("bracket", { recursive: true }); } catch {}
    const dir = await root.getDirectoryHandle("bracket", { create: true });
    const at = async (path, create) => {
      const parts = path.split("/");
      let d = dir;
      for (const p of parts.slice(0, -1)) d = await d.getDirectoryHandle(p, { create });
      return d.getFileHandle(parts.at(-1), { create });
    };
    const put = async (path, text) => { const w = await (await at(path, true)).createWritable(); await w.write(text); await w.close(); };
    const read = async (path) => (await (await at(path, false)).getFile()).text();
    const base = "migrations/001_base.rhai";
    await put("popkern.toml", '[project]\\nname = "Bracket"\\n');
    await put(base, 'let m = model();\\nm.profile(xy()).rect(1.0, 0.5).sketch("plate").extrude(blind(0.1));\\nm\\n');
    window.showDirectoryPicker = async () => dir;
    const opened = await page.projects.openFolder();
    await page.runDone;
    const result = {
      opened: opened?.name ?? null,
      kind: page.projects.current?.kind ?? null,
      outcome: page.lastGood?.envelope?.outcome ?? null,
      listed: (await page.projects.store.list()).some((p) => p.id === opened?.id && p.kind === "folder"),
    };
    const edit = (from, to) => {
      page.editor.textarea.value = page.editor.value.replace(from, to);
      page.editor.textarea.dispatchEvent(new Event("input"));
    };
    edit("blind(0.1)", "blind(0.2)");
    await page.projects.flush();
    result.written = (await read(base)).includes("blind(0.2)");
    await put(base, (await read(base)).replace("blind(0.2)", "blind(0.3)"));
    await page.projects.checkFolder();
    await page.runDone;
    result.reloaded = page.editor.value.includes("blind(0.3)");
    edit("blind(0.3)", "blind(0.4)");
    await put(base, (await read(base)).replace("blind(0.3)", "blind(0.5)"));
    await page.projects.flush();
    result.conflict_shown = !document.getElementById("conflict-bar").hidden;
    result.conflict_text = document.getElementById("conflict-text").textContent;
    result.disk_untouched = (await read(base)).includes("blind(0.5)");
    await page.projects.loadTheirs();
    await page.runDone;
    result.theirs = page.editor.value.includes("blind(0.5)");
    await page.projects.remove(opened.id);
    result.forgotten = !(await page.projects.store.list()).some((p) => p.id === opened.id);
    result.files_kept = (await read("popkern.toml")).includes("Bracket");
    window.showDirectoryPicker = undefined;
    page.projects.store.folders = null;
    result.unsupported_opened = await page.projects.openFolder();
    result.unsupported_status = document.getElementById("status").textContent;
    delete window.showDirectoryPicker;
    await page.projects.openStored("scratch");
    await page.runDone;
    return result;
  })()`);

  // Insert into an assembly part (PK-678): a part opened in the
  // assembly panel scrubs its own features (its bodies alone); insert
  // there adds an edit_part migration that builds into every instance.
  out.part_insert = await evaluate(`(async () => {
    const page = window.popkernLive;
    const s = page.scrubber;
    const wait = () => new Promise((r) => setTimeout(r, 50));
    const script = 'let a = assembly("pins");\\nlet pin = a.part("pin", |m| {\\n    m.profile(xy()).circle(0.1).sketch("pin").extrude(blind(1.0));\\n    m.active().fillet(edges().circular().at_height(0.0), 0.02);\\n});\\na.place(pin, "left");\\na.place(pin, "right").translate(0.5, 0.0, 0.0);\\na\\n';
    await page.projects.create("Part insert", script);
    await page.runDone;
    // This project's own build (a baked page may still be finishing
    // another), and its timeline.
    await page.run({ live: true });
    await s.loading;
    page.assembly.open("pin");
    await s.loading;
    const anchor = s.steps.find((x) => x.operation === "extrude");
    await s.go(anchor.index);
    const result = {
      scoped: s.part,
      meshes: page.viewer.bodies.length,
      shown: !document.getElementById("scrubber-insert").hidden,
    };
    document.getElementById("scrubber-insert").click();
    await s.inserting;
    result.path = page.projects.path;
    const template = page.editor.value;
    result.edits_part = template.includes('a.edit_part("pin"');
    page.editor.textarea.value = template.replace("        body\\n", "        body.chamfer(edges().circular().at_height(1.0), 0.02)\\n");
    page.editor.textarea.dispatchEvent(new Event("input"));
    await page.projects.flush();
    clearTimeout(page.idleTimer);
    await page.runDone;
    const built = await page.run({ live: true });
    result.statuses = (built?.envelope?.project?.versions ?? []).map((v) => v.status);
    result.occurrences = page.viewer.bodies.length;
    await s.loading;
    page.assembly.open(null);
    await s.loading;
    page.assembly.open("pin");
    await s.loading;
    result.labels = s.steps.map((x) => x.label);
    page.assembly.open(null);
    await page.projects.openStored("scratch");
    await page.runDone;
    return result;
  })()`);

  // Migrations (PK-664): add, rename and reorder files through the
  // rail; an upstream edit that drops a label breaks the downstream
  // migration, which the overlay names with its line and the upstream
  // change, opening that file.
  out.migrations = await evaluate(`(async () => {
    const page = window.popkernLive;
    window.prompt = () => "rim_round";
    window.confirm = () => true;
    const block = 'let m = model();\\nlet block = m.profile(xy()).rect(2.0, 2.0).sketch("block_sketch").extrude(blind(1.0));\\nlet bore = m.tool(|m| m.profile(xy()).circle(0.4).sketch("bore_sketch").extrude(blind(1.0)));\\nblock.subtract(bore);\\nm.label_edge("bore_rim", edges().circular().at_height(1.0));\\nm\\n';
    await page.projects.create("Rail test", block);
    await page.runDone;
    const single = { rail_hidden: document.getElementById("migrations-rail").hidden, add_shown: !document.getElementById("migration-add").hidden };
    await page.rail.add();
    const added = page.projects.path;
    page.editor.textarea.value = 'let m = previous();\\nm.active().fillet(edges().labeled("bore_rim"), 0.05, #{ label: "rim_round" });\\nm\\n';
    page.editor.textarea.dispatchEvent(new Event("input"));
    await page.projects.flush();
    clearTimeout(page.idleTimer);
    await page.runDone;
    const built = await page.run({ live: true });
    const statuses = (built?.envelope?.project?.versions ?? []).map((v) => [v.name, v.status]);
    window.prompt = () => "fillet the rim";
    await page.rail.rename(added);
    const renamed = page.projects.path;
    // Drop the label upstream: the fillet migration fails on its line.
    page.rail.open("migrations/001_base.rhai");
    page.editor.textarea.value = page.editor.value.split("\\n").filter((l) => !l.includes("label_edge")).join("\\n");
    page.editor.textarea.dispatchEvent(new Event("input"));
    await page.projects.flush();
    clearTimeout(page.idleTimer);
    await page.runDone;
    const broken = await page.run({ live: true });
    const overlay = document.getElementById("overlay-text").textContent;
    const failedBadge = document.querySelector("#migrations-list li[data-path='" + renamed + "'] .badge")?.textContent ?? null;
    const result = {
      single,
      added,
      statuses,
      renamed,
      broken_outcome: broken?.envelope?.outcome ?? null,
      overlay_names_file: overlay.includes("in 010_fillet_the_rim.rhai, line 2"),
      overlay_names_upstream: overlay.includes("changed upstream since it last built: 001_base"),
      editor_on_failing_file: page.projects.path === renamed,
      marked_line: [...document.querySelectorAll(".editor-line-error")].length === 1,
      failed_badge: failedBadge,
    };
    // Reorder: the fillet migration moves first (numbers swap).
    await page.rail.move(renamed, -1);
    result.reordered = Object.keys(page.projects.current.files).filter((p) => p.startsWith("migrations/")).sort();
    await page.projects.openStored("scratch");
    await page.runDone;
    return result;
  })()`);

  // The standard views (view-orient.js): each turns the eye to look
  // along its axis, Z up, and refits.
  out.views = await evaluate(`(async () => {
    const page = window.popkernLive;
    await page.projects.create("Views", 'let m = model();\\nm.profile(xy()).rect(3.0, 1.0).sketch("s").extrude(blind(0.5));\\nm\\n');
    await page.runDone;
    const { orientView } = await import("./view-orient.js");
    const look = {};
    for (const name of ["front", "back", "left", "right", "top", "bottom", "iso"]) {
      orientView(page.viewer, name, 0);
      const d = page.viewer.camera.getWorldDirection(new page.viewer.camera.position.constructor());
      look[name] = [d.x, d.y, d.z].map((v) => Math.round(v * 1000) / 1000);
    }
    return { look, pad_buttons: document.querySelectorAll("#view-orient [data-view]").length };
  })()`);

  // A phone (PK-867, PK-868): a project's films are listed on the
  // History sheet; a film opens as a stage over the model view (the
  // sheet stays where it was) with exactly one Close and one Play, a
  // timeline of one chapter per step, and while it is prepared a card
  // naming what it waits on — also when a build is still running as it
  // opens (it waits, then plays). The film frame's shade stays in the
  // model view, the status chip clears the corner's tools, and the view
  // modes keep their clustered card.
  await send("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 1, mobile: true });
  await new Promise((r) => setTimeout(r, 400));
  out.phone = await evaluate(`(async () => {
    const page = window.popkernLive;
    const rect = (el) => el.getBoundingClientRect();
    const shown = (el) => el && getComputedStyle(el).display !== "none" && !el.hidden && el.getClientRects().length > 0;
    const stages = [];
    const loading = document.getElementById("film-loading");
    new MutationObserver(() => { if (!loading.hidden) stages.push(document.getElementById("film-loading-stage").textContent.replace(/[.…].*$/, "")); })
      .observe(loading, { attributes: true, childList: true, subtree: true, characterData: true });
    const waitActive = async () => {
      const t0 = performance.now();
      while (!page.film.active && performance.now() - t0 < 120000) await new Promise((r) => setTimeout(r, 100));
      await new Promise((r) => setTimeout(r, 250));
    };
    window.prompt = () => "spin";
    document.querySelector(".sheet-tabs [data-sheet=code]").click();
    await page.rail.addFilm();
    await waitActive();
    const pane = rect(document.querySelector(".pane-viewer"));
    const hole = rect(document.querySelector("#film-frame .film-hole"));
    const frame = rect(document.getElementById("film-frame"));
    const bar = rect(document.getElementById("film-bar"));
    const visible = (selector) => [...document.querySelectorAll(selector)].filter(shown).length;
    const open = {
      active: page.film.active,
      sheet: document.body.dataset.sheet,
      bar_shown: shown(document.getElementById("film-bar")),
      bar_in_pane: bar.left >= pane.left - 1 && bar.right <= pane.right + 1 && bar.bottom <= pane.bottom + 1,
      frame_above_bar: frame.bottom <= bar.top + 1,
      closes: [...document.querySelectorAll("button")].filter((b) => shown(b) && /close film/i.test(b.textContent + (b.getAttribute("aria-label") ?? ""))).length,
      plays: [...document.querySelectorAll("button")].filter((b) => shown(b) && /play the film|pause the film/i.test(b.getAttribute("aria-label") ?? "")).length,
      chapters: document.querySelectorAll("#film-chapters .film-chapter").length,
      steps: page.film.film.steps.length,
      corner_hidden: visible(".view-corner") === 0,
      frame_in_pane: frame.left >= pane.left - 1 && frame.right <= pane.right + 1 && frame.top >= pane.top - 1 && frame.bottom <= pane.bottom + 1,
      hole_in_pane: hole.left >= pane.left - 1 && hole.right <= pane.right + 1 && hole.top >= pane.top - 1 && hole.bottom <= pane.bottom + 1,
      clipped: getComputedStyle(document.getElementById("film-frame")).overflow === "hidden",
      stages: [...new Set(stages)],
    };
    document.getElementById("film-play").click();
    await new Promise((r) => setTimeout(r, 300));
    open.played = page.film.playing;
    page.film.pause();
    // Seek on the timeline: a click three quarters along.
    const line = rect(document.getElementById("film-timeline"));
    const at = line.left + line.width * 0.75;
    document.getElementById("film-timeline").dispatchEvent(new PointerEvent("pointerdown", { clientX: at, clientY: line.top + line.height / 2, pointerId: 7, bubbles: true }));
    document.getElementById("film-timeline").dispatchEvent(new PointerEvent("pointerup", { clientX: at, pointerId: 7, bubbles: true }));
    await new Promise((r) => setTimeout(r, 100));
    open.seek_fraction = page.film.time / page.film.film.duration;
    document.getElementById("film-exit").click();
    await new Promise((r) => setTimeout(r, 200));
    const closed = { active: page.film.active, stage_gone: !document.body.classList.contains("film-stage") };
    // Opening a film while a build runs: the stage waits for the build
    // (saying so), then plays.
    stages.length = 0;
    const building = page.run({ live: true });
    page.projects.select("films/spin.rhai");
    await building;
    await waitActive();
    const during_build = { active: page.film.active, stages: [...new Set(stages)] };
    document.getElementById("film-exit").click();
    await new Promise((r) => setTimeout(r, 200));
    document.querySelector(".sheet-tabs [data-sheet=history]").click();
    // The settings menu opens over everything, the sheet tabs too.
    document.getElementById("settings-menu-button").click();
    await new Promise((r) => setTimeout(r, 250));
    const menuBox = rect(document.getElementById("settings-menu"));
    const tabsBox = rect(document.querySelector(".sheet-tabs"));
    const probeY = Math.max(menuBox.top + 8, Math.min(menuBox.bottom - 8, tabsBox.top + tabsBox.height / 2));
    const topmost = document.elementFromPoint(menuBox.left + menuBox.width / 2, probeY);
    const settings_on_top = { open: !document.getElementById("settings-menu").hidden, topmost: Boolean(topmost?.closest("#settings-menu")), scrim: document.body.classList.contains("menu-open") };
    document.getElementById("settings-done").click();
    settings_on_top.closed = document.getElementById("settings-menu").hidden && !document.body.classList.contains("menu-open");
    const list = document.getElementById("film-list");
    const status = rect(document.querySelector(".dock-status"));
    const tools = rect(document.querySelector(".view-tools"));
    const view = document.querySelector(".bar-view");
    const mode = document.querySelector(".bar-view .mode");
    return {
      open,
      closed: { ...closed, list_shown: shown(list), list_items: list.querySelectorAll(".film-list-item").length },
      during_build,
      settings_on_top,
      status_clear_of_tools: status.right <= tools.left + 0.5 || status.bottom <= tools.top,
      modes_clustered: getComputedStyle(view).backgroundColor !== "rgba(0, 0, 0, 0)" && parseFloat(getComputedStyle(mode).borderTopLeftRadius) < 20,
      orient_shown: shown(document.getElementById("view-orient")),
    };
  })()`);
  // A mid desktop width (PK-871): the film tab keeps its name in view
  // beside the add buttons (at 1280 px it hid behind "+ migration
  // + film" until they collapsed to glyphs in a narrow card).
  await send("Emulation.setDeviceMetricsOverride", { width: 1280, height: 800, deviceScaleFactor: 1, mobile: false });
  await new Promise((r) => setTimeout(r, 400));
  out.rail_desktop = await evaluate(`(async () => {
    if (window.popkernLive.film.active) document.getElementById("film-exit")?.click();
    document.querySelector(".migration:not(.migration-film) .migration-name")?.click();
    await new Promise((r) => setTimeout(r, 300));
    const tabs = document.querySelector(".file-tabs").getBoundingClientRect();
    const film = document.querySelector(".migration-film .migration-name");
    const r = film?.getBoundingClientRect();
    return {
      film_text: film?.textContent ?? null,
      visible: Boolean(r) && r.width > 0 && r.left >= tabs.left - 0.5 && r.right <= tabs.right + 0.5,
      add_labels_hidden: [...document.querySelectorAll(".rail-add .rail-add-label")].every((l) => getComputedStyle(l).display === "none"),
    };
  })()`);
  await send("Emulation.clearDeviceMetricsOverride");

  // Themes, Scratch's film, deleting projects, and a film opened the
  // moment its example opens (PK-869): it must play on that example's
  // model, never the previous one's.
  out.projects_ui = await evaluate(`(async () => {
    const page = window.popkernLive;
    const theme = await import("./theme.js");
    const bg = () => getComputedStyle(document.documentElement).getPropertyValue("--bg").trim();
    const before = bg();
    theme.setTheme("light");
    const themes = { default: before, light: bg(), stored: localStorage.getItem("popkern.theme"), options: document.querySelectorAll("#theme option").length };
    theme.setTheme("cad");
    themes.back = bg() === before;
    const scratch = await page.projects.store.get("scratch");
    const scratch_film = Boolean(scratch?.files?.["films/tour.rhai"]);
    // Delete the open Scratch: another project opens; Scratch is gone.
    await page.projects.openStored("scratch");
    await page.runDone;
    await page.projects.removeAndReopen("scratch", await page.defaultScript(), await page.defaultFilm());
    const deleted = { scratch_gone: (await page.projects.store.get("scratch")) == null, reopened: page.projects.current != null && page.projects.current.id !== "scratch" };
    // A fresh Scratch, with the default film, comes back when asked for.
    const fresh = await page.projects.scratch(await page.defaultScript(), await page.defaultFilm());
    deleted.fresh_film = Boolean(fresh.files["films/tour.rhai"]);
    // A film selected the moment its example opens, twice in a row.
    const race = [];
    for (const [id, film] of [["example:boss_plate_in_steps", "films/section_sweep.rhai"], ["example:door_hinge.rhai", "films/swing.rhai"]]) {
      await page.projects.openExample(id);
      page.projects.select(film);
      const t0 = performance.now();
      while (!page.film.active && performance.now() - t0 < 120000) await new Promise((r) => setTimeout(r, 100));
      // Every film chosen starts at zero (the last film was left at 4 s).
      const start = page.film.time;
      // A seek moves the film's bodies on screen: three.js rebuilt an
      // auto-updating mesh's matrix at every render, dropping them.
      // (from the start: the player keeps its time across films)
      await page.film.seek(0);
      page.viewer.render();
      const before = page.viewer.bodies.map((m) => m.matrixWorld.elements.join());
      await page.film.seek(4.0);
      page.viewer.render();
      const moved = page.viewer.bodies.some((m, i) => m.matrixWorld.elements.join() !== before[i]);
      race.push({ id, active: page.film.active, moved, start, overlay: document.getElementById("overlay").hidden ? null : document.getElementById("overlay-text").textContent.slice(0, 200), bodies: page.film.context?.bodies ?? null, kind: page.film.context?.kind ?? null });
      document.getElementById("film-exit").click();
    }
    // Settings: one menu in the top bar; the section's options are a
    // panel that shows with the section.
    const settings = document.getElementById("settings-menu");
    const settings_row = {
      has: ["theme", "time-limit", "edges"].every((id) => settings.contains(document.getElementById(id))),
      in_topbar: Boolean(settings.closest(".topbar")),
      corner_buttons: document.querySelectorAll(".view-tools button").length,
    };
    const panel = document.getElementById("section-panel");
    settings_row.panel_hidden_off = panel.hidden;
    document.getElementById("section-toggle").click();
    await new Promise((r) => setTimeout(r, 100));
    settings_row.panel_shown_on = !panel.hidden;
    document.getElementById("section-toggle").click();
    return { themes, scratch_film, deleted, race, settings: settings_row };
  })()`);

  if (process.env.POPKERN_LAUNCH === "1") {
    out.launch = await launchRows();
  }
  out.network = [...new Set(networkUrls)];
  // Foreign means another scheme://host:port — CSP 'self' semantics.
  // (Under a subpath deploy the page origin is host + prefix, but a
  // same-host request outside the prefix still never leaves the box.)
  const pageHost = new URL(`${origin}/`).origin;
  out.foreign_requests = out.network.filter((u) => !u.startsWith(`${pageHost}/`));
  out.server_requests = [...new Set(requests)];
  out.console = consoleLines.slice(0, 40);
  out.exceptions = exceptions.slice(0, 20);
} catch (error) {
  out.error = error.message;
  out.status = await statusText().catch(() => null);
  out.console = consoleLines.slice(0, 40);
  out.exceptions = exceptions.slice(0, 20);
}

/// The launch-checklist rows (PK-415).
async function launchRows() {
  const rows = {};
  // Keyboard reach and reduced motion, by inspection of the live page.
  rows.keyboard = await evaluate(`({
    canvas_focusable: document.getElementById("viewer").tabIndex >= 0,
    run_button: !!document.getElementById("run"),
    editor_is_textarea: !!document.querySelector("#editor textarea"),
    reduced_motion_rule: [...document.styleSheets].some((s) => { try { return [...s.cssRules].some((r) => r.conditionText && r.conditionText.includes("prefers-reduced-motion")); } catch { return false; } }),
    focus_visible_rule: [...document.styleSheets].some((s) => { try { return [...s.cssRules].some((r) => r.selectorText && r.selectorText.includes("focus-visible")); } catch { return false; } }),
  })`);
  // The three runaway classes under small budgets on a fresh host.
  const spawnHost = `(() => { const Host = window.popkernLive.host.constructor; const h = new Host(() => new Worker("./worker.js", { type: "module" }), { max_operations: 200000, timeout_ms: 400 }); return h.start().then(() => { window.__launchHost = h; return h.info.protocol_version; }); })()`;
  await evaluate(spawnHost);
  const timed = async (label, script, grace) => {
    const started = Date.now();
    const reply = await evaluate(`(async () => { const h = window.__launchHost; const p = h.run(${JSON.stringify(script)}, "runaway.rhai"); const cancel = new Promise((r) => setTimeout(() => r("timeout"), ${grace})); const first = await Promise.race([p, cancel]); if (first === "timeout") { await h.cancel("launch harness"); const r = await p; return { cancelled: true, outcome: r.outcome?.outcome ?? r.envelope?.outcome, restarts: h.restarts }; } return { cancelled: false, outcome: first.envelope?.outcome ?? first.outcome?.outcome, kind: first.envelope?.kind ?? null, budget: first.envelope?.payload?.budget ?? null }; })()`);
    return { label, elapsed_ms: Date.now() - started, ...reply };
  };
  rows.runaway_ops = await timed("op-bound loop", "let m = model();\nlet i = 0;\nloop { i += 1; }\n", 15000);
  rows.runaway_wall = await timed("wall-bound feature loop", 'let m = model();\nlet i = 0;\nloop { m.profile(xy_at(i * 1.0)).rect(1.0, 1.0).sketch("s" + i).extrude(blind(0.5)); i += 1; }\n', 15000);
  rows.runaway_defeats_budgets = await timed("one long kernel call (G2 cube, no progress ticks)", 'let m = model();\nlet cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));\ncube.fillet(edges(), 0.1, #{ continuity: "g2" });\nm\n', 2000);
  // Memory: a 200-feature document on a generous host.
  await evaluate(`(() => { const Host = window.popkernLive.host.constructor; const h = new Host(() => new Worker("./worker.js", { type: "module" }), { max_operations: null, timeout_ms: null }); return h.start().then(() => { window.__bigHost = h; return true; }); })()`);
  const big = `let m = model();\nlet i = 0;\nwhile i < 100 { m.profile(xy_at(i * 1.0)).rect(1.0, 1.0).sketch("s" + i).extrude(blind(0.5)); i += 1; }\nm\n`;
  const memStart = Date.now();
  rows.big_document = await evaluate(`(async () => { const h = window.__bigHost; const before = await h.stats(); const r = await h.run(${JSON.stringify(big)}, "big.rhai"); const after = await h.stats(); return { outcome: r.envelope?.outcome ?? r.outcome?.outcome, features: r.envelope?.features ?? null, bodies: r.envelope?.bodies?.length ?? null, wasm_memory_before: before.wasm_memory_bytes, wasm_memory_after: after.wasm_memory_bytes, js_heap: performance.memory ? performance.memory.usedJSHeapSize : null, elapsed_ms: r.envelope?.elapsed_ms ?? null }; })()`);
  rows.big_document.wall_ms = Date.now() - memStart;
  // Throttled cold load: Fast 3G, a fresh page load.
  await send("Network.emulateNetworkConditions", { offline: false, latency: 562.5, downloadThroughput: (1.6 * 1024 * 1024) / 8, uploadThroughput: (750 * 1024) / 8 });
  await send("Network.setCacheDisabled", { cacheDisabled: true });
  const t0 = Date.now();
  await send("Page.navigate", { url: `${origin}/?cold` });
  const status = await waitFor(
    `(() => { const s = document.getElementById("status"); if (!s) return null; const t = s.textContent; return /^(served|script error|kernel decline|gate failure|page error|worker error|kernel failed)/.test(t) ? t : null; })()`,
    240000,
    "the throttled first run",
  );
  rows.cold_load_fast3g = {
    status,
    wall_to_first_status_ms: Date.now() - t0,
    kernel_ready_ms: await evaluate(`window.popkernLive.kernelReadyMs ?? null`),
    first_render_ms: await evaluate(`window.popkernLive.firstRenderMs ?? null`),
    bundle_gzip_bytes_note: "the wasm is served uncompressed by this harness; a static host with gzip/brotli halves the transfer",
  };
  await send("Network.emulateNetworkConditions", { offline: false, latency: 0, downloadThroughput: -1, uploadThroughput: -1 });
  return rows;
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
