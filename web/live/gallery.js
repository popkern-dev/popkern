// The default script (`default.rhai`): what a fresh Scratch project
// holds. (The examples are read-only projects: example-store.js.)

/** The script a fresh Scratch project holds. It lives in default.rhai
 *  so tools/bake-live.sh pre-builds it with the gallery; a failed fetch
 *  resolves to a placeholder script that says why. */
/** The default model's film (`films/tour.rhai` of a fresh Scratch), or
 *  null when it cannot be loaded (Scratch then starts without one). */
export function fetchDefaultFilm() {
  return fetch("./default.film.rhai")
    .then((r) => (r.ok ? r.text() : null))
    .catch(() => null);
}

export function fetchDefaultScript() {
  return fetch("./default.rhai")
    .then((r) => (r.ok ? r.text() : Promise.reject(new Error(`HTTP ${r.status}`))))
    .catch((error) => `// The default script could not be loaded (${error.message}).\n// Pick an example above, or write a script here.\nlet m = model();\nm\n`);
}
