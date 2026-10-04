# popkern live (static page)

popkern live, the browser workbench (crate `popkern-live`; the working
name until 2026-10-04 was "playground"): a script editor, the kernel compiled to
WebAssembly running in a Web Worker, a three.js viewport, the feature
tree, parameter sliders, STEP/`.pkn` exports and show-as-code.
Everything runs in the tab — no server-side component, no network
use beyond loading these files.

## Build and serve

```
tools/build-live.sh --web        # → web/live/pkg/ (gitignored)
tools/bake-live.sh               # → web/live/bakes/ (optional: pre-built gallery)
cd web/live && python3 -m http.server 8000
```

Any static file server works; the page needs `application/wasm` for
`.wasm` (every common server has it) and must be served over HTTP —
`file://` cannot load module workers. There is no bundler and no
build step for the JS: the files here are what the browser loads.
Every URL in the page is relative and the CSP lives in meta tags, so
the page is subpath-agnostic — it serves identically from `/` and
from a mounted prefix.

## Sandbox deploy (`/routes/`)

```
tools/deploy-live-routes.sh              # build, bake, copy, verify
tools/deploy-live-routes.sh --no-build   # redeploy the built bundle
```

One command builds the bundle, copies the runtime assets (not the
harnesses or npm metadata) to `~/code/.routes/static/popkern-live/`,
and verifies the served copy: content types by `curl`, then the
`test/browser.mjs` harness pointed at the deployed URL via
`POPKERN_PAGE_ORIGIN` — a real kernel run, request accounting, and
recorded load numbers in `docs/reports/live/routes-deploy.json`.
The sandbox's Caddy serves the tree at `/routes/popkern-live/`
(locally on port 7681 and on the tunnel hostname); the script's
header documents the `conf.d` reverse-proxy variant a live
`popkern-daemon` session will use instead. The former route,
`/routes/popkern-live/`, is rewritten on every deploy to
redirect pages pointing here.

## Files

| File | Role |
|---|---|
| `index.html`, `style.css` | The page and its CSP (`'self'` for everything; `'wasm-unsafe-eval'` for the module) |
| `app.js` | The page controller (`window.popkernLive`): owns the build lifecycle — run, bake shortcut, cancel, parameter replay — and the state the modules below share |
| `controls.js` | Wires the top bar, viewer bar, draft notice and dialogs to the controller |
| `local-store.js` | Everything kept in `localStorage` behind one guarded accessor: the draft, the held-aside draft, the time limit, and the one-time `popkern-playground-*` key adoption |
| `drafts.js` | Which document a visit opens (pre-built first, an edited script held aside), restore, reset, opening an example |
| `document-tab.js` | The open document: its title, file stem and edited mark; an imported STEP placeholder |
| `gallery.js` | The examples dropdown from `examples.json`, example and default-script fetches |
| `time-limit.js` | The build time-limit presets, the worker budgets they map to, and applying a change |
| `status-line.js` | The status line, the pick readout, served-model wording, page-error reporting |
| `outcome-overlay.js` | The viewport overlay for a build that did not serve, with the line jump |
| `report-table.js` | The feature tree rows (blend continuity measured vs requested) and its summary |
| `params-panel.js` | Parameter sliders and number boxes, scoped to an opened assembly part |
| `model-display.js` | A served model or pre-built bake landing on screen (viewer, overlay, status) |
| `continuity-certification.js` | Deferred seam certification in slices after a served build |
| `analysis-mesh.js` | The zebra view's finer analysis mesh, fetched once per served model |
| `exports.js` | STEP / `.pkn` downloads (tree visibility saved into the `.pkn`) and show-as-code |
| `clearance-panel.js` | Interference and clearance check and its pair list |
| `step-import.js` | Open STEP as a view-only assembly |
| `assembly-panel.js` | The assembly tree: visibility, selection, open part |
| `instances.js` | Assembly occurrences in the viewer: per-occurrence drawing, placement, explode directions |
| `shading.js` | Viewer shading: face palette, zebra shaders, curvature map |
| `editor.js` | Dependency-free editor (textarea + highlighted mirror + gutter) |
| `viewer.js` | three.js scene: the view modes (CAD grey with hover pre-highlight, face colours, zebra, curvature), face and edge picking |
| `edges.js` | Face-boundary edges from the mesh buffers, grouped by the face pair they separate |
| `orbit.js` | Orbit control: drag, pinch, wheel, arrow keys |
| `worker.js`, `worker-core.js` | The Web Worker owning one `LiveSession` |
| `worker-host.js` | The page's side of the protocol: request matching, cancel = terminate + restart |
| `PROTOCOL.md` | The worker message shapes (versioned; the daemon's socket schema grows from it) |
| `vendor/` | Pinned third-party code with hashes in `vendor/LICENSES.md` |
| `test/host.test.mjs` | Headless protocol harness under node `worker_threads` |
| `test/browser.mjs` | Headless Chrome harness over the DevTools protocol (serves the page itself) |
| `examples.json`, `examples/` | The gallery manifest and scripts (phase 4c) |
| `default.rhai` | The script a fresh page and reset show (pre-built with the gallery) |
| `bakes.js`, `bakes/` | Pre-built examples: the page shows the default or a gallery script's certified build without the kernel (`bakes/` is generated by `tools/bake-live.sh`, gitignored) |

The Rust side runs both harnesses from
`crates/popkern-live/tests/page.rs` and skips, with a recorded
reason, when the bundle, node or a browser is absent.

## Budgets and runaway scripts

The worker runs scripts under an operation budget and a wall-clock
budget (`performance.now` is registered as the kernel's clock). A
script that defeats both is stopped by the page: the worker is
terminated and a fresh one started, the last good model stays on
screen.
