// PK-662: the page's project zip (web/live/zip.js) under node.
//   node zip_roundtrip.mjs write <project-dir> <out.zip>   — zip a folder
//   node zip_roundtrip.mjs read <in.zip>                   — print its files
import { readFileSync, readdirSync, statSync, writeFileSync } from "node:fs";
import { join, relative } from "node:path";
import { unzipFiles, zipFiles } from "../../../../web/live/zip.js";

const [mode, a, b] = process.argv.slice(2);
if (mode === "write") {
  const files = {};
  const walk = (dir) => {
    for (const name of readdirSync(dir)) {
      const path = join(dir, name);
      if (statSync(path).isDirectory()) walk(path);
      else files[relative(a, path).split("\\").join("/")] = new Uint8Array(readFileSync(path));
    }
  };
  walk(a);
  writeFileSync(b, zipFiles(files));
} else if (mode === "read") {
  try {
    const files = unzipFiles(new Uint8Array(readFileSync(a)));
    const out = {};
    for (const [path, data] of Object.entries(files)) out[path] = Buffer.from(data).toString("base64");
    console.log(JSON.stringify({ ok: true, files: out }));
  } catch (error) {
    console.log(JSON.stringify({ ok: false, error: error.message }));
  }
} else {
  console.error("usage: zip_roundtrip.mjs write <dir> <out.zip> | read <in.zip>");
  process.exit(2);
}
