// Render the lander's model images in popkern live itself: each gallery
// example opened pre-built, every card hidden, the transparent canvas
// filling the window, screenshotted at 2x by headless Chrome (CDP).
// tools/cutout.py then crops and downscales each one into assets/img/.
//
//   tools/build-live.sh --web && tools/bake-live.sh      (repo root)
//   node web/lander/tools/render-gallery.mjs <outdir> [chrome]
//   for each JOBS entry: python3 web/lander/tools/cutout.py \
//       <outdir>/<name>.png web/lander/assets/img/<out>.png <width>
//
// A job: { name, out, width, file, shade?, color?, theta?, phi?, zoom?,
// edges?, w?, h?, live?, settle? } — shade is a Shade button (cad,
// zebra, gaussian, mean), color a Color button (parts, faces), theta
// and phi the orbit angles (radians, Z up), live builds the model in
// the page first (zebra needs the live analysis mesh).
import { createServer } from "node:http";
import { readFile, stat, mkdtemp, writeFile } from "node:fs/promises";
import { join, extname } from "node:path";
import { tmpdir } from "node:os";
import { spawn } from "node:child_process";

import { fileURLToPath } from "node:url";
const [outdir, chromeArg] = process.argv.slice(2);
const chromeBin = chromeArg ?? process.env.POPKERN_CHROME ?? "google-chrome";
const root = fileURLToPath(new URL("../../live", import.meta.url));
const JOBS = JSON.parse(await readFile(new URL("./gallery-jobs.json", import.meta.url), "utf8"));
const jobs = JOBS;
const MIME = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".wasm": "application/wasm", ".rhai": "text/plain", ".png": "image/png", ".svg": "image/svg+xml", ".bin": "application/octet-stream" };
const server = createServer(async (req, res) => {
  const url = new URL(req.url, "http://x");
  let path = join(root, decodeURIComponent(url.pathname));
  try {
    if ((await stat(path)).isDirectory()) path = join(path, "index.html");
    res.writeHead(200, { "content-type": MIME[extname(path)] ?? "application/octet-stream", "cache-control": "no-store" });
    res.end(await readFile(path));
  } catch { res.writeHead(404); res.end(); }
});
await new Promise((r) => server.listen(0, "127.0.0.1", r));
const origin = `http://127.0.0.1:${server.address().port}`;
const port = 9300 + Math.floor(Math.random() * 500);
const profile = await mkdtemp(join(tmpdir(), "render-"));
const chrome = spawn(chromeBin, ["--headless=new", "--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist", "--disable-dev-shm-usage", `--user-data-dir=${profile}`, `--remote-debugging-port=${port}`, "about:blank"], { stdio: "ignore" });
let list = null;
for (let i = 0; i < 75 && !list; i++) {
  await new Promise((r) => setTimeout(r, 200));
  try { list = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); } catch {}
}
const ws = new WebSocket(list.find((t) => t.type === "page").webSocketDebuggerUrl);
await new Promise((r) => (ws.onopen = r));
let id = 1;
const waiting = new Map();
ws.onmessage = (e) => { const m = JSON.parse(e.data); if (m.id && waiting.has(m.id)) { waiting.get(m.id)(m); waiting.delete(m.id); } };
const send = (method, params = {}) => new Promise((r) => { const i = id++; waiting.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
const evaluate = async (expression) => {
  const r = await send("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true });
  if (r.result?.exceptionDetails) throw new Error(r.result.exceptionDetails.exception?.description ?? "eval failed");
  return r.result?.result?.value;
};
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
await send("Page.enable");
await send("Emulation.setDefaultBackgroundColorOverride", { color: { r: 0, g: 0, b: 0, a: 0 } });
let size = null;
let navigated = false;
for (const job of jobs) {
  const w = job.w ?? 1400, h = job.h ?? 1050;
  if (size !== `${w}x${h}`) {
    await send("Emulation.setDeviceMetricsOverride", { width: w, height: h, deviceScaleFactor: 2, mobile: false });
    size = `${w}x${h}`;
  }
  if (!navigated) {
    navigated = true;
    await send("Page.navigate", { url: `${origin}/index.html` });
    await sleep(2500);
    await evaluate(`window.popkernLive.setTimeLimit(null)`);
  }
  const t0 = Date.now();
  const result = await evaluate(`(async () => {
    const p = window.popkernLive;
    location.hash = "";
    await p.loadExample(${JSON.stringify(job.file)});
    await p.runDone;
    ${job.live ? "await p.ensureLive();" : ""}
    const shade = ${JSON.stringify(job.shade ?? "cad")};
    document.querySelector('[data-shade="' + shade + '"]').click();
    const color = ${JSON.stringify(job.color ?? null)};
    for (const b of document.querySelectorAll("[data-color]")) if (b.getAttribute("aria-pressed") === "true") b.click();
    if (color) document.querySelector('[data-color="' + color + '"]').click();
    if (shade === "zebra") await p.ensureAnalysisMesh();
    p.viewer.toggleEdges(${job.edges ?? true});
    // Only the model: every card hidden, the viewport fills the window.
    for (const el of document.querySelectorAll("body > *")) if (!el.matches("main")) el.style.display = "none";
    for (const el of document.querySelectorAll("main > *")) if (!el.matches(".pane-viewer")) el.style.display = "none";
    for (const el of document.querySelectorAll(".pane-viewer > *")) if (!el.matches("canvas")) el.style.display = "none";
    const main = document.querySelector("main");
    main.style.display = "block"; main.style.padding = "0";
    const pane = document.querySelector(".pane-viewer");
    pane.style.position = "fixed"; pane.style.inset = "0";
    document.body.style.background = "transparent";
    document.body.style.backgroundImage = "none";
    await new Promise((r) => setTimeout(r, 300));
    p.viewer.resize();
    p.viewer.fitAgain();
    const o = p.viewer.orbit;
    ${job.theta != null ? `o.spherical.theta = ${job.theta};` : ""}
    ${job.phi != null ? `o.spherical.phi = ${job.phi};` : ""}
    o.spherical.radius *= ${job.zoom ?? 1};
    o.update();
    return document.getElementById("status").textContent;
  })()`);
  await sleep(job.settle ?? 1500);
  const shot = await send("Page.captureScreenshot", { format: "png" });
  await writeFile(join(outdir, `${job.name}.png`), Buffer.from(shot.result.data, "base64"));
  console.log(job.name, `${((Date.now() - t0) / 1000).toFixed(1)} s`, result.slice(0, 110));
}
chrome.kill();
server.close();
process.exit(0);
