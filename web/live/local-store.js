// What the page remembers in this browser: the working draft, an
// edited script held aside from an earlier visit, and the build time
// limit. Everything lives in `localStorage` on this origin (design §7:
// nothing leaves the tab), and every access is guarded — a browser
// that blocks storage (private mode, a sandboxed frame) gets a page
// that simply forgets.

const DRAFT_KEY = "popkern-live-draft";
// An edited script from an earlier visit, held aside while the page
// opens on the pre-built default (restore puts it back).
const SAVED_DRAFT_KEY = "popkern-live-saved-draft";
const TIME_LIMIT_KEY = "popkern-live-time-limit";

/** Run `action` against `localStorage`; null when storage is blocked
 *  or the action throws. */
function storage(action) {
  try {
    return action(window.localStorage);
  } catch {
    return null;
  }
}

/** The page stored its keys as "popkern-playground-*" until the crate
 *  became popkern-live (2026-10-04). Both routes share one origin, so a
 *  visitor's draft and time limit are moved to the current keys once,
 *  on the first load that finds them. */
export function adoptPlaygroundKeys() {
  storage((s) => {
    for (const key of [DRAFT_KEY, SAVED_DRAFT_KEY, TIME_LIMIT_KEY]) {
      const earlier = key.replace(/^popkern-live-/, "popkern-playground-");
      const value = s.getItem(earlier);
      if (value != null && s.getItem(key) == null) s.setItem(key, value);
      s.removeItem(earlier);
    }
  });
}

/** The working draft (the script a refresh reopens), or null. */
export function readDraft() {
  return storage((s) => s.getItem(DRAFT_KEY));
}

/** Remember `text` as the working draft. */
export function writeDraft(text) {
  storage((s) => s.setItem(DRAFT_KEY, text));
}

/** The edited script held aside from an earlier visit, or null. */
export function readSavedDraft() {
  return storage((s) => s.getItem(SAVED_DRAFT_KEY));
}

/** Hold `text` aside as the edited script a restore puts back. */
export function writeSavedDraft(text) {
  storage((s) => s.setItem(SAVED_DRAFT_KEY, text));
}

/** Forget the held-aside script. */
export function clearSavedDraft() {
  storage((s) => s.removeItem(SAVED_DRAFT_KEY));
}

/** The stored time-limit choice: "none", a number of seconds as text,
 *  or null when nothing (or nothing readable) is stored. */
export function readTimeLimit() {
  return storage((s) => s.getItem(TIME_LIMIT_KEY));
}

/** Remember the time limit in seconds (null = no limit). */
export function writeTimeLimit(limit) {
  storage((s) => s.setItem(TIME_LIMIT_KEY, limit == null ? "none" : String(limit)));
}
