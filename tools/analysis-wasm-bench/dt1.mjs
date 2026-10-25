// DT1, wasm half: solves the exported cases with the wasm build in node
// and compares each result line with the native one.
import { readFileSync } from "node:fs";
import init, { run_case } from "../../target/wasm-bench/pkg/popkern_analysis_wasm_bench.js";

const dir = new URL("../../target/wasm-bench/", import.meta.url);
await init({ module_or_path: readFileSync(new URL("pkg/popkern_analysis_wasm_bench_bg.wasm", dir)) });
const native = readFileSync(new URL("native.txt", dir), "utf8").trim().split("\n");
let failed = 0;
for (const line of native) {
  const [name, ...expect] = line.split(" ");
  const bytes = readFileSync(new URL(`${name}.bin`, dir));
  const got = run_case(bytes, 4000);
  const ok = got === expect.join(" ");
  if (!ok) failed += 1;
  console.log(`${ok ? "same" : "DIFFERENT"} ${name}: wasm ${got} native ${expect.join(" ")}`);
}
process.exit(failed ? 1 : 0);
