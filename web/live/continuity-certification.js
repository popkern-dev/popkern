// Deferred continuity certification (PROTOCOL.md, "Deferred continuity
// certification"): fillets serve on sampled continuity evidence, and
// the page certifies their seams afterwards in short slices, so a G2
// part is on screen long before its continuity is proven.
import { deferredSeams } from "./report-table.js";

/** The worker time one certification slice may take. */
const CERTIFY_SLICE_MS = 250;

/** Certify `envelope`'s deferred seams slice by slice on `page`'s
 *  worker, swapping in each upgraded report. Abandoned when a newer
 *  run, parameter edit or stop bumps `page.certifyGeneration`. */
export async function certifyContinuity(page, envelope) {
  const total = deferredSeams(envelope.report);
  if (!total) return;
  const generation = ++page.certifyGeneration;
  const started = performance.now();
  let remaining = total;
  while (remaining > 0) {
    page.status(`${page.servedText} — certifying continuity: ${total - remaining}/${total} seams`);
    const reply = await page.host.certify(CERTIFY_SLICE_MS);
    if (generation !== page.certifyGeneration) return;
    const step = reply.ok ? reply.step : { failure: reply.outcome };
    if (step.failure) {
      page.status(`${page.servedText} — continuity certification stopped: ${step.failure.detail ?? step.failure.kind ?? "failed"}`, true);
      return;
    }
    remaining = step.remaining;
    if (step.report) {
      envelope.report = step.report;
      page.renderInspector(envelope);
    }
  }
  page.status(`${page.servedText} — continuity certified (${total} seams, ${((performance.now() - started) / 1000).toFixed(1)} s)`);
}
