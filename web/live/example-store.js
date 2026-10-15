// The examples as read-only projects (live-projects-and-history.md
// §12). `examples.json` lists them: an entry with `file` is one script
// (served as a one-migration project), an entry with `project` is a
// folder beside the page whose files the entry lists. Opening one never
// changes it; "make a copy to edit" copies it into the browser store.
import { singleScriptProject } from "./project-store.js";

/** The text of a response, or an Error naming what failed. */
async function text(response, what) {
  if (!response.ok) throw new Error(`${what}: HTTP ${response.status}`);
  return response.text();
}

/** The stem of a script file name (`boss_plate.rhai` -> `boss_plate`). */
function stemOf(file) {
  return file.replace(/\.rhai$/, "").replace(/[^A-Za-z0-9_-]/g, "_");
}

export class ExampleStore {
  constructor(entries) {
    this.entries = entries;
    this.kind = "example";
  }

  /** The gallery (an empty store when there is none beside the page). */
  static async open() {
    try {
      const response = await fetch("./examples.json");
      if (!response.ok) return new ExampleStore([]);
      const manifest = await response.json();
      return new ExampleStore(manifest.examples ?? []);
    } catch {
      return new ExampleStore([]);
    }
  }

  /** Rows for the projects dialog. */
  list() {
    return this.entries.map((e) => ({
      id: `example:${e.project ?? e.file}`,
      name: e.title,
      blurb: e.blurb ?? "",
      demonstrates: e.demonstrates ?? "",
      migrations: e.project ? (e.files ?? []).filter((f) => f.startsWith("migrations/")).length : 1,
    }));
  }

  /** The entry with gallery file `file` (single-script examples). */
  byFile(file) {
    return this.entries.find((e) => e.file === file) ?? null;
  }

  /** The read-only project `id` (`example:<file or folder>`). */
  async get(id) {
    const key = id.replace(/^example:/, "");
    const entry = this.entries.find((e) => (e.project ?? e.file) === key);
    if (!entry) return null;
    let files;
    if (entry.project) {
      files = {};
      for (const path of entry.files ?? []) files[path] = await text(await fetch(`./examples/${entry.project}/${path}`), path);
    } else {
      files = singleScriptProject(entry.title, stemOf(entry.file), await text(await fetch(`./examples/${entry.file}`), entry.file));
    }
    // A single script's films (`stem.<name>.film.rhai` beside it) open
    // as the project's `films/<name>.rhai`.
    for (const film of entry.films ?? []) {
      files[`films/${film.name}.rhai`] = await text(await fetch(`./examples/${film.file}`), film.file);
    }
    return { id, name: entry.title, files, revision: 0, readOnly: true, file: entry.file ?? null };
  }
}
