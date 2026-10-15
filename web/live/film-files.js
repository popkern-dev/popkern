// Film files in a project (docs/architecture/presentation-films.md §2):
// `films/<name>.rhai`, any number, beside `migrations/`; a single
// script's `stem.film.rhai` / `stem.<name>.film.rhai` siblings open as
// `films/<name>.rhai` of its one-migration project.

/** Whether `path` names a film of a project. */
export function isFilmPath(path) {
  return typeof path === "string" && /^films\/[A-Za-z0-9_-]+\.rhai$/.test(path);
}

/** The film's name from its path. */
export function filmName(path) {
  return path.replace(/^films\//, "").replace(/\.rhai$/, "");
}

/** The path of the film called `name` (letters, digits, _ -). */
export function filmPath(name) {
  const clean = String(name).trim().toLowerCase().replace(/\s+/g, "_").replace(/[^a-z0-9_-]/g, "") || "film";
  return `films/${clean}.rhai`;
}

/** The starting text of a new film. */
export const NEW_FILM = `// A film of the model (docs: the guide's Films section). It reads the
// built document and never changes it: edit and press Build to replay.
let doc = document();
let film = storyboard(#{ look: "toon" });
film.title(doc.name, "", 2.0, |s| { s.orbit(40.0); });
film.step("Turn", "", 4.0, |s| { s.orbit(120.0); });
film
`;
