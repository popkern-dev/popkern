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
  out.status_after_certify = await waitFor(
    `(() => { const s = document.getElementById("status").textContent; return /continuity (certified|certification stopped)/.test(s) ? s : null; })()`,
    120000,
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
    for (let i = 0; i < 100 && s.steps.length === 0; i++) await new Promise((r) => setTimeout(r, 50));
    const camera = () => page.viewer.camera.position.toArray().map((v) => v.toFixed(6)).join(",");
    const before = camera();
    const ticks = document.querySelectorAll("#scrubber-ticks .tick").length;
    const visible = !document.getElementById("scrubber").hidden;
    await s.go(1);
    const first = { bodies: page.viewer.bodies.length, label: document.getElementById("scrubber-label").textContent, line: page.editor.stepLine };
    await s.go(3);
    const middle = { bodies: page.viewer.bodies.length, glows: page.viewer.glows.length, row: document.querySelector("#report tr.step-current")?.dataset.feature ?? null, hash: location.hash };
    await s.go(5);
    const revisit = await s.go(3);
    await s.go(s.steps.length);
    const last = { bodies: page.viewer.bodies.length, scrubbed: s.scrubbed, hash: location.hash };
    s.play();
    for (let i = 0; i < 200 && s.playing; i++) await new Promise((r) => setTimeout(r, 50));
    const played = { at: s.at, playing: Boolean(s.playing) };
    location.hash = "step=2";
    await s.load(page.lastGood.envelope);
    const linked = s.at;
    await s.go(s.steps.length);
    return { visible, ticks, first, middle, revisit_ms: revisit, last, played, linked, camera_kept: camera() === before };
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
      const status = await waitFor(`(() => { const s = document.getElementById("status").textContent; return s.startsWith("served (pre-built)") && s !== ${JSON.stringify(before)} ? s : null; })()`, 30000, `the ${file} example`);
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

  // Open STEP (PK-631): a third-party writer's assembly (OpenCASCADE
  // XCAF, tools/occ-oracle/occ_write_assembly.py) opens with its tree,
  // names and shared parts.
  const stepFixture = (await readFile(join(root, "..", "..", "crates/popkern-assembly/tests/fixtures/occ_bench_assembly.step"))).toString("base64");
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
    out.returning.restored = await waitFor(`(() => { const s = document.getElementById("status").textContent; return /^served: /.test(s) ? s : null; })()`, 60000, "the restored draft");
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
