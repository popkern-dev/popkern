// PK-660 project smoke: build a project folder through the wasm-bindgen
// (nodejs target) bundle's run_project, the way the page hands a
// project over, and print what the Rust side compares against the
// native runner: the outcome, feature count, each migration's row, and
// every body's mesh volume computed from its buffers.
//
//   node project_smoke.mjs <pkg-dir> <project-dir>
import { createRequire } from "node:module";
import { readFileSync, readdirSync, statSync } from "node:fs";
import { resolve, join, relative } from "node:path";

const [pkgDir, projectDir] = process.argv.slice(2);
if (!pkgDir || !projectDir) {
  console.error("usage: node project_smoke.mjs <pkg-dir> <project-dir>");
  process.exit(2);
}
const require = createRequire(import.meta.url);
const pkg = require(resolve(pkgDir, "popkern_live.js"));

const files = {};
const base64 = {};
function walk(dir) {
  for (const name of readdirSync(dir)) {
    const path = join(dir, name);
    if (statSync(path).isDirectory()) {
      walk(path);
      continue;
    }
    const key = relative(projectDir, path).split("\\").join("/");
    if (/\.(rhai|toml)$/.test(name)) files[key] = readFileSync(path, "utf8");
    else base64[key] = readFileSync(path).toString("base64");
  }
}
walk(projectDir);

pkg.init();
const session = new pkg.LiveSession("{}");
const envelope = JSON.parse(session.run_project(JSON.stringify({ files, base64 })));

const volumes = [];
for (let i = 0; i < session.body_count(); i++) {
  const p = session.body_positions(i);
  const t = session.body_indices(i);
  let v = 0;
  for (let k = 0; k < t.length; k += 3) {
    const a = 3 * t[k], b = 3 * t[k + 1], c = 3 * t[k + 2];
    v += p[a] * (p[b + 1] * p[c + 2] - p[b + 2] * p[c + 1])
       - p[a + 1] * (p[b] * p[c + 2] - p[b + 2] * p[c])
       + p[a + 2] * (p[b] * p[c + 1] - p[b + 1] * p[c]);
  }
  volumes.push(v / 6);
}
console.log(JSON.stringify({
  outcome: envelope.outcome,
  features: envelope.features,
  project: envelope.project,
  volumes,
}));
