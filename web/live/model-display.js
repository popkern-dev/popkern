// How a served model lands on screen: drawn in the viewer with the
// assembly tree's visibility and selection re-applied, the overlay
// cleared, and the status line saying what served — for a kernel run
// (then certified in the background) and for a pre-built bake (already
// certified natively).
import { clearOutcome } from "./outcome-overlay.js";
import { meshClosure, modelContents } from "./status-line.js";

const $ = (id) => document.getElementById(id);

/** Draw a served model, re-applying the tree's visibility and
 *  selection (a rebuild keeps what the user hid). */
export function drawModel(page, envelope, buffers) {
  // The explode slider belongs to assemblies; a new document starts
  // assembled.
  const exploding = envelope.assembly != null;
  if (!exploding || $("explode-row").hidden) {
    $("explode").value = "0";
    page.viewer.explode = 0;
  }
  $("explode-row").hidden = !exploding;
  page.viewer.show(envelope.bodies, buffers, envelope.assembly ?? null);
  for (const path of page.assembly.hidden) page.viewer.setVisible(path, false);
  if (page.assembly.selected != null) page.viewer.select(page.assembly.selected);
}

/** " — first render N ms after page start" the first time a model is
 *  drawn (recording when on `page.firstRenderMs`), else "". */
function firstRenderNote(page) {
  return page.firstRenderMs == null ? ` — first render ${(page.firstRenderMs = performance.now()).toFixed(0)} ms after page start` : "";
}

/** Show a kernel run's served model and start certifying its deferred
 *  seams. The zebra view's analysis mesh is fetched first when zebra
 *  is the current mode. */
export async function showServed(page, envelope, buffers) {
  page.lastGood = { envelope, buffers };
  page.viewer.setAnalysis(null);
  if (page.viewer.mode === "zebra") await page.ensureAnalysisMesh();
  drawModel(page, envelope, buffers);
  clearOutcome(page.editor);
  const ms = envelope.elapsed_ms == null ? "" : ` in ${envelope.elapsed_ms.toFixed(0)} ms`;
  page.servedText = `served: ${envelope.features} feature(s), ${modelContents(envelope)}${ms}${meshClosure(envelope)}${firstRenderNote(page)}`;
  page.status(page.servedText);
  page.certify(envelope);
  // The history of the model now live (scrubber.js).
  page.scrubber.load(envelope);
}

/** Show a pre-built example: the certified report, parameters and
 *  display meshes exactly as the native build left them (the worker's
 *  model is not the one on screen until something rebuilds it). */
export function showBake(page, bake, source) {
  page.live = false;
  page.scrubber.reset();
  page.certifyGeneration += 1;
  const { envelope, buffers } = bake;
  page.lastGood = { envelope, buffers, bakedSource: source };
  page.assembly.render(envelope.assembly ?? null);
  page.renderInspector();
  page.viewer.setAnalysis(null);
  drawModel(page, envelope, buffers);
  clearOutcome(page.editor);
  const native = ((bake.native?.run_ms ?? 0) + (bake.native?.certify_ms ?? 0)) / 1000;
  page.servedText = `served (pre-built): ${envelope.features} feature(s), ${modelContents(envelope)}${meshClosure(envelope)}${firstRenderNote(page)}`;
  page.status(`${page.servedText} — continuity certified; builds natively in ${native.toFixed(1)} s, and edits, sliders and exports rebuild it live`);
  if (page.viewer.mode === "zebra") page.ensureAnalysisMesh();
}
