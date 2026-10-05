// Which project the page opens (live-projects-and-history.md §12). A
// visitor types into the Scratch project, saved in this browser. A
// fresh visit opens Scratch when it is still pre-built (the default, or
// an example it was set to) so nothing waits on the kernel; a Scratch
// with edits from last time stays in the store while the page opens on
// the pre-built default, with a notice offering it back. The draft
// kept by pages before projects existed seeds Scratch once.
import { DEFAULT_TITLE } from "./document-tab.js";
import { clearLegacyDrafts, readLegacyDraft } from "./local-store.js";
import { migrationPaths, singleScriptProject } from "./project-store.js";
import { SCRATCH } from "./projects.js";

const $ = (id) => document.getElementById(id);

/** The default script as a read-only project. */
function defaultProject(source) {
  return {
    id: "example:default",
    name: DEFAULT_TITLE,
    files: singleScriptProject(DEFAULT_TITLE, "boss_plate", source),
    revision: 0,
    readOnly: true,
  };
}

/** Open the first project of a visit. */
export async function openFirstScript(page) {
  const projects = await page.projectsReady;
  const fallback = await page.defaultScript();
  const legacy = readLegacyDraft();
  let scratch = await projects.scratch(legacy ?? fallback);
  if (legacy != null) {
    // The pre-project draft fills an unedited Scratch, else becomes a
    // project of its own — never silently dropped.
    const path = migrationPaths(scratch.files).at(-1);
    if (scratch.files[path] === fallback && legacy !== fallback) {
      scratch = await projects.store.save(SCRATCH, { ...scratch.files, [path]: legacy }, scratch.revision);
    } else if (scratch.files[path] !== legacy) {
      await projects.store.create("Recovered draft", singleScriptProject("Recovered draft", "draft", legacy));
    }
    // Forget the pre-project keys (unless storage forgets everything
    // anyway; then they seed the next visit too).
    if (projects.store.persistent) clearLegacyDrafts();
  }
  const source = scratch.files[migrationPaths(scratch.files).at(-1)] ?? "";
  if (source === fallback || (await page.bakes.find(source))) {
    $("draft-notice").hidden = true;
    await projects.open(scratch);
    return;
  }
  $("draft-notice").hidden = false;
  await projects.open(defaultProject(fallback));
}

/** Show the notice while Scratch waits to be reopened. */
export function showDraftNotice() {
  return !$("draft-notice").hidden;
}

/** Open Scratch with its edits and build it. */
export async function restoreDraft(page) {
  $("draft-notice").hidden = true;
  await page.projects.openStored(SCRATCH);
}

/** Hide the notice (Scratch is kept; the projects dialog opens it). */
export function dismissDraft() {
  $("draft-notice").hidden = true;
}

/** Open a gallery example (`file` in examples/) read-only. */
export async function openExample(page, file) {
  if (!file) return;
  await page.projects.openExample(file);
}
