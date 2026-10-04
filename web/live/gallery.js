// The example gallery beside the page (`examples.json` + `examples/`)
// and the default script (`default.rhai`): the scripts a visitor can
// open without writing one. All same-origin, relative fetches.

const $ = (id) => document.getElementById(id);

/** Fill the examples dropdown from the gallery manifest. Resolves to a
 *  map of example file -> gallery title, or undefined/null when there
 *  is no gallery (the picker then stays hidden). */
export async function loadGallery() {
  try {
    const response = await fetch("./examples.json");
    if (!response.ok) return;
    const manifest = await response.json();
    const select = $("examples");
    const titles = new Map();
    for (const entry of manifest.examples ?? []) {
      titles.set(entry.file, entry.title);
      const option = document.createElement("option");
      option.value = entry.file;
      option.textContent = entry.title;
      select.append(option);
    }
    $("examples-row").hidden = select.options.length <= 1;
    return titles;
  } catch {
    // No gallery beside the page: the example picker stays hidden.
    return null;
  }
}

/** An example script's text, or null when it cannot be fetched. */
export async function fetchExample(file) {
  const response = await fetch(`./examples/${file}`);
  if (!response.ok) return null;
  return response.text();
}

/** The script a fresh page (and reset) shows. It lives in default.rhai
 *  so tools/bake-live.sh pre-builds it with the gallery; a failed fetch
 *  resolves to a placeholder script that says why. */
export function fetchDefaultScript() {
  return fetch("./default.rhai")
    .then((r) => (r.ok ? r.text() : Promise.reject(new Error(`HTTP ${r.status}`))))
    .catch((error) => `// The default script could not be loaded (${error.message}).\n// Pick an example above, or write a script here.\nlet m = model();\nm\n`);
}
