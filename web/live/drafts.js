// Which document the page opens, and the draft behind it. A refresh
// reopens the working draft; a visit opens on something pre-built so it
// never waits on the kernel, holding an edited script from last time
// aside with a notice offering to restore it. Choosing an example or
// resetting replaces the draft.
import { DEFAULT_TITLE, UNTITLED } from "./document-tab.js";
import { fetchExample } from "./gallery.js";
import { clearSavedDraft, readDraft, readSavedDraft, writeDraft, writeSavedDraft } from "./local-store.js";

const $ = (id) => document.getElementById(id);

/** The page opens on something pre-built: the last script when it has
 *  a bake (the default or an example), else the default — and an
 *  edited script from last time is held aside. */
export async function openFirstScript(page) {
  const draft = readDraft();
  const fallback = await page.defaultScript();
  const prebuilt = draft ? await page.bakes.find(draft) : null;
  if (draft && !prebuilt && draft !== fallback) writeSavedDraft(draft);
  const name = prebuilt ? await page.bakes.nameOf(draft) : null;
  const titles = await page.examplesReady;
  page.tab.open(titles?.get(name) ?? DEFAULT_TITLE, prebuilt ? draft : fallback);
  showDraftNotice();
  page.run();
}

/** Show the restore notice while a held-aside script exists. */
export function showDraftNotice() {
  $("draft-notice").hidden = !readSavedDraft();
}

/** Put the held-aside script back as the working draft and build it. */
export function restoreDraft(page) {
  const saved = readSavedDraft();
  clearSavedDraft();
  showDraftNotice();
  if (!saved) return;
  page.tab.open(UNTITLED, saved);
  writeDraft(saved);
  page.viewer.fitted = false;
  page.run();
}

/** Forget the held-aside script. */
export function dismissDraft() {
  clearSavedDraft();
  showDraftNotice();
}

/** Back to the default script, remembered as the draft. */
export async function resetToDefault(page) {
  page.tab.open(DEFAULT_TITLE, await page.defaultScript());
  writeDraft(page.editor.value);
  page.run();
}

/** Open a gallery example (`file` in examples/) under its title. */
export async function openExample(page, file) {
  if (!file) return;
  const text = await fetchExample(file);
  if (text == null) return;
  const titles = await page.examplesReady;
  page.tab.open(titles?.get(file) ?? file.replace(/\.rhai$/, ""), text);
  // A refresh reopens the chosen example (pre-built, so instant).
  writeDraft(page.editor.value);
  page.viewer.fitted = false;
  page.run();
}
