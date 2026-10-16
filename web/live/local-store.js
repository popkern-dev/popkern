// What the page keeps in localStorage: the build time limit, the
// colour theme (theme.js; read before the first paint by
// theme-boot.js, the one other reader), and the
// draft keys pages used before projects existed (read once to seed the
// Scratch project; projects themselves live in IndexedDB,
// project-store.js). Everything lives in `localStorage` on this origin (design §7:
// nothing leaves the tab), and every access is guarded — a browser
// that blocks storage (private mode, a sandboxed frame) gets a page
// that simply forgets.

const DRAFT_KEY = "popkern-live-draft";
// An edited script from an earlier visit, held aside while the page
// opens on the pre-built default (restore puts it back).
const SAVED_DRAFT_KEY = "popkern-live-saved-draft";
const TIME_LIMIT_KEY = "popkern-live-time-limit";
/** The colour theme's key (theme-boot.js reads the same one). */
export const THEME_KEY = "popkern.theme";

/** The remembered colour theme id, or null. */
export function readTheme() {
  return storage((s) => s.getItem(THEME_KEY));
}

/** Remember colour theme `id` (forgotten when storage is blocked). */
export function writeTheme(id) {
  storage((s) => s.setItem(THEME_KEY, id));
}

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

/** The script pages kept before projects existed (the held-aside
 *  edited script, else the working draft), or null. It seeds the
 *  Scratch project once. */
export function readLegacyDraft() {
  return storage((s) => s.getItem(SAVED_DRAFT_KEY) ?? s.getItem(DRAFT_KEY));
}

/** Forget the pre-project draft keys (adopted into Scratch). */
export function clearLegacyDrafts() {
  storage((s) => {
    s.removeItem(DRAFT_KEY);
    s.removeItem(SAVED_DRAFT_KEY);
  });
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

const PARAM_FILTER_KEY = "popkern-live-param-filter";

/** Which parameter kinds the panel shows: `{ named, implicit }`, both
 *  shown unless the visitor turned one off. */
export function readParamFilter() {
  const stored = storage((s) => JSON.parse(s.getItem(PARAM_FILTER_KEY) ?? "null"));
  return { named: stored?.named !== false, implicit: stored?.implicit !== false };
}

export function writeParamFilter(filter) {
  storage((s) => s.setItem(PARAM_FILTER_KEY, JSON.stringify(filter)));
}
