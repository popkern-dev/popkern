// The zebra view renders from a finer ANALYSIS mesh (6°, PK-577) than
// the display mesh: fetched from the worker once per served model and
// kept until the next run. On a failure the zebra view falls back to
// the display mesh and the status line says why.

/** Fetch the analysis mesh for `page`'s model on screen (building the
 *  live model first if a pre-built one is shown) and hand it to the
 *  viewer. A newer model on screen drops the reply. */
export async function ensureAnalysisMesh(page) {
  if (!page.lastGood || page.lastGood.analysis) return;
  if (!(await page.ensureLive())) return;
  const good = page.lastGood;
  if (good.analysis) return;
  const reply = await page.host.analysisMesh();
  if (page.lastGood !== good) return;
  if (!reply?.ok) {
    const why = reply?.outcome?.detail ?? "unavailable";
    page.status(`zebra: ${why}`, reply?.outcome?.step !== "analysis_mesh");
    return;
  }
  good.analysis = reply.buffers;
  page.viewer.setAnalysis(reply.buffers);
  const { triangles, mesh_ms, facet_deg } = reply.stats;
  const ms = mesh_ms == null ? "" : ` in ${mesh_ms.toFixed(0)} ms`;
  page.status(`zebra: ${triangles} triangles at ${facet_deg}°${ms}`);
}
