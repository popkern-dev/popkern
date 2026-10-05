// The page's controls wired to the page controller: the top bar
// (build, time limit, stop, open STEP, downloads, show-as-code; the
// projects controls are wired in projects-dialog.js), the viewer bar (fit, view mode, edges, explode,
// clearance), the draft notice and the dialogs.
import { dismissDraft, restoreDraft } from "./drafts.js";
import { downloadModel, showModelCode } from "./exports.js";

const $ = (id) => document.getElementById(id);
const DIALOGS = ["about-dialog", "code-dialog", "projects-dialog"];

/** Wire every control to `page` (the Page in app.js). */
export function wireControls(page) {
  $("run").addEventListener("click", () => page.run({ live: true }));
  $("time-limit").addEventListener("change", (e) => page.setTimeLimit(e.target.value === "none" ? null : Number(e.target.value)));
  $("cancel").addEventListener("click", () => page.cancel("cancelled by the user"));
  $("fit").addEventListener("click", () => page.viewer.fitAgain());
  $("view-mode").addEventListener("change", async (e) => {
    const mode = e.target.value;
    if (mode === "zebra") await page.ensureAnalysisMesh();
    page.viewer.setMode(mode);
  });
  $("edges").addEventListener("change", (e) => page.viewer.toggleEdges(e.target.checked));
  $("explode").addEventListener("input", (e) => page.viewer.setExplode(Number(e.target.value)));
  $("clearance-check").addEventListener("click", () => page.checkClearance());
  $("open-step").addEventListener("click", () => $("open-step-file").click());
  $("open-step-file").addEventListener("change", async (e) => {
    const file = e.target.files?.[0];
    if (file) await page.openStep(file.name, await file.arrayBuffer());
    e.target.value = "";
  });
  $("download-step").addEventListener("click", () => downloadModel(page, "step"));
  $("download-pkn").addEventListener("click", () => downloadModel(page, "pkn"));
  $("show-code").addEventListener("click", () => showModelCode(page));
  $("code-close").addEventListener("click", () => ($("code-dialog").hidden = true));
  $("code-use").addEventListener("click", () => {
    if (page.editor.readOnly) return;
    page.editor.value = $("code-text").textContent;
    page.tab.show();
    page.projects.edited(page.editor.value);
    $("code-dialog").hidden = true;
    page.run();
  });
  $("draft-restore").addEventListener("click", () => restoreDraft(page));
  $("about").addEventListener("click", () => ($("about-dialog").hidden = false));
  $("about-close").addEventListener("click", () => ($("about-dialog").hidden = true));
  for (const dialog of DIALOGS) {
    // A click on the backdrop (not the box) or Escape closes it.
    $(dialog).addEventListener("click", (e) => {
      if (e.target === e.currentTarget) e.currentTarget.hidden = true;
    });
  }
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") for (const dialog of DIALOGS) $(dialog).hidden = true;
  });
  $("draft-dismiss").addEventListener("click", () => dismissDraft());
}
