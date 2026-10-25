"""Writes <out>/e3-<budget>.html for each budget: the wasm build and the
bench case embedded as base64, solved synchronously while the page
parses, so `chrome --dump-dom` returns the time on the real clock. One
page per budget: wasm memory never shrinks, so solves share no page."""
import base64
import pathlib
import sys

out = pathlib.Path(sys.argv[1])
glue = (out / "pkg-classic" / "popkern_analysis_wasm_bench.js").read_text()
wasm = base64.b64encode((out / "pkg-classic" / "popkern_analysis_wasm_bench_bg.wasm").read_bytes()).decode()
case = base64.b64encode((out / "block_two_bores.bin").read_bytes()).decode()
page = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>E3 solve timing</title></head>
<body><pre id="out">running</pre>
<script>{glue}</script>
<script>
const bytes = s => Uint8Array.from(atob(s), c => c.charCodeAt(0));
const lines = [];
try {{
  const wasm = wasm_bindgen.initSync({{ module: bytes("{wasm}") }});
  const input = bytes("{case}");
  for (const cells of [BUDGET]) {{
    const t0 = performance.now();
    const r = wasm_bindgen.run_case(input, cells);
    const t = (performance.now() - t0) / 1000;
    const mb = wasm.memory.buffer.byteLength / 1048576;
    lines.push(`budget ${{cells}}: ${{t.toFixed(2)}} s, wasm memory ${{mb.toFixed(0)}} MB, result ${{r}}`);
  }}
}} catch (e) {{
  lines.push(`error: ${{e}}`);
}}
document.getElementById("out").textContent = lines.join(" | ");
</script></body></html>
"""
for budget in sys.argv[2:]:
    (out / f"e3-{budget}.html").write_text(page.replace("BUDGET", budget))
