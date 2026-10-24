// The migrations rail (live-projects-and-history.md §12): a project's
// migration files in order as the script card's tabs, each with the
// status the last build gave it (built, cached, failed, stale; the
// stylesheet shows only failed and stale), the open one highlighted. A migration opens in the editor with a click; an
// editable project adds, renames, reorders and deletes migrations
// (renaming and reordering rename files: the number orders them). The
// rail stays hidden for a one-migration project until a second
// migration is added, so a single script looks as it always did.
import { NEXT_MIGRATION, insertAfterMigration, insertAfterName, insertAfterPartMigration } from "./migration-templates.js";
import { migrationPaths } from "./project-store.js";
import { NEW_FILM, filmName, filmPath } from "./film-files.js";
import { iconLabel, iconSvg } from "./icons.js";

const $ = (id) => document.getElementById(id);

/** `migrations/012_name.rhai` -> { number: 12, name: "name" }. */
export function parseMigration(path) {
  const m = /^migrations\/(\d+)_([A-Za-z0-9_-]+)\.rhai$/.exec(path);
  return m ? { number: Number(m[1]), name: m[2] } : null;
}

/** The path of migration `number` called `name`. */
export function migrationPath(number, name) {
  return `migrations/${String(number).padStart(3, "0")}_${name}.rhai`;
}

/** A migration name from what a person typed (letters, digits, _ -). */
export function cleanName(text) {
  return text.trim().toLowerCase().replace(/\s+/g, "_").replace(/[^a-z0-9_-]/g, "") || "step";
}


export class MigrationsRail {
  constructor(page) {
    this.page = page;
    this.versions = new Map();
    $("migration-add").addEventListener("click", () => this.add());
    $("film-add").addEventListener("click", () => this.addFilm());
  }

  /** Statuses from a project build's envelope (`project.versions`). */
  setVersions(envelope) {
    this.versions = new Map((envelope?.project?.versions ?? []).map((v) => [v.name, v.status]));
    this.render();
  }

  /** The phone's History sheet lists the project's films (a phone
   *  shows the rail only on the Code sheet): one button each, the
   *  open film marked. */
  renderFilms() {
    const films = this.page.projects.current ? this.page.projects.filmPaths() : [];
    const box = $("film-list");
    box.hidden = films.length === 0;
    const openPath = this.page.projects.filmOpen ? this.page.projects.path : null;
    for (const b of box.querySelectorAll("button")) b.classList.toggle("open", b.dataset.path === openPath);
    if (box.dataset.films === films.join("\n")) return;
    box.dataset.films = films.join("\n");
    box.querySelectorAll("button").forEach((b) => b.remove());
    for (const path of films) {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "film-list-item";
      b.dataset.path = path;
      b.classList.toggle("open", path === openPath);
      b.replaceChildren(...iconLabel("play", filmName(path)));
      b.title = `play the film ${filmName(path)}`;
      b.addEventListener("click", () => this.open(path));
      box.append(b);
    }
  }

  render() {
    this.renderFilms();
    const p = this.page.projects.current;
    const list = $("migrations-list");
    list.replaceChildren();
    const paths = p ? migrationPaths(p.files) : [];
    const films = p ? this.page.projects.filmPaths() : [];
    // The rail shows for several migrations, or for any film (then the
    // one migration shows too, to get back to the model).
    const shown = p && (paths.length >= 2 || films.length > 0);
    $("migrations-rail").hidden = !shown;
    $("migration-add").hidden = !p || p.readOnly;
    $("film-add").hidden = !p || p.readOnly;
    if (!shown) return;
    paths.forEach((path, i) => {
      const parsed = parseMigration(path);
      const stem = path.split("/").pop().replace(/\.rhai$/, "");
      const li = document.createElement("li");
      li.className = "migration";
      li.dataset.path = path;
      if (path === p.path) li.classList.add("open");
      const status = this.versions.get(stem);
      const open = document.createElement("button");
      open.type = "button";
      open.className = "migration-name";
      open.textContent = `${String(parsed.number).padStart(3, "0")} ${parsed.name}`;
      open.title = `open ${stem}.rhai`;
      open.addEventListener("click", () => this.open(path));
      li.append(open);
      if (status) {
        const badge = document.createElement("span");
        badge.className = `badge badge-${status}`;
        badge.textContent = status;
        li.append(badge);
      }
      if (!p.readOnly) {
        for (const [icon, title, act, hidden] of [
          ["up", "move earlier", () => this.move(path, -1), i === 0],
          ["down", "move later", () => this.move(path, 1), i === paths.length - 1],
          ["rename", "rename", () => this.rename(path), false],
          ["close", "delete this migration", () => this.remove(path), false],
        ]) {
          if (hidden) continue;
          const b = document.createElement("button");
          b.type = "button";
          b.className = "migration-tool";
          b.innerHTML = iconSvg(icon);
          b.title = title;
          b.setAttribute("aria-label", title);
          b.addEventListener("click", act);
          li.append(b);
        }
      }
      list.append(li);
    });
    for (const path of films) {
      const li = document.createElement("li");
      li.className = "migration migration-film";
      li.dataset.path = path;
      if (path === p.path) li.classList.add("open");
      const open = document.createElement("button");
      open.type = "button";
      open.className = "migration-name";
      open.replaceChildren(...iconLabel("play", filmName(path)));
      open.title = `play the film ${filmName(path)} (films/${filmName(path)}.rhai)`;
      open.addEventListener("click", () => this.open(path));
      li.append(open);
      if (!p.readOnly) {
        for (const [icon, title, act] of [
          ["rename", "rename", () => this.renameFilm(path)],
          ["close", "delete this film", () => this.removeFilm(path)],
        ]) {
          const b = document.createElement("button");
          b.type = "button";
          b.className = "migration-tool";
          b.innerHTML = iconSvg(icon);
          b.title = title;
          b.setAttribute("aria-label", title);
          b.addEventListener("click", act);
          li.append(b);
        }
      }
      list.append(li);
    }
  }

  /** Add a film (films/<name>.rhai from a starter) and open it. */
  async addFilm() {
    const p = this.page.projects.current;
    if (!p || p.readOnly) return;
    const name = window.prompt("Name of the new film", "turntable");
    if (name == null) return;
    const path = filmPath(name);
    if (path in p.files) return this.page.status(`a film named ${path} exists`, true);
    await this.page.projects.replaceFiles({ ...p.files, [path]: NEW_FILM }, path);
    this.render();
  }

  async renameFilm(path) {
    const p = this.page.projects.current;
    const name = window.prompt("New name of the film", filmName(path));
    if (name == null) return;
    const target = filmPath(name);
    if (target === path) return;
    if (target in p.files) return this.page.status(`a film named ${target} exists`, true);
    const files = { ...p.files, [target]: p.files[path] };
    delete files[path];
    await this.page.projects.replaceFiles(files, p.path === path ? target : p.path);
    this.render();
  }

  async removeFilm(path) {
    const p = this.page.projects.current;
    if (!window.confirm(`Delete the film ${filmName(path)}?`)) return;
    const files = { ...p.files };
    delete files[path];
    if (p.path === path) this.page.film.close();
    await this.page.projects.replaceFiles(files, p.path === path ? migrationPaths(files).at(-1) : p.path);
    this.render();
  }

  /** Show migration `path` in the editor (no build: the project as a
   *  whole is what builds). */
  open(path) {
    const p = this.page.projects.current;
    if (!p || !(path in p.files)) return;
    this.page.projects.select(path);
    this.render();
  }

  /** Add the next migration (`previous()` template) and open it. */
  async add() {
    const p = this.page.projects.current;
    if (!p || p.readOnly) return;
    const name = window.prompt("Name of the new migration (what it does)", "next_step");
    if (name == null) return;
    const last = migrationPaths(p.files).map(parseMigration).at(-1);
    const number = last ? Math.floor(last.number / 10) * 10 + 10 : 1;
    const path = migrationPath(number, cleanName(name));
    await this.page.projects.replaceFiles({ ...p.files, [path]: NEXT_MIGRATION }, path);
    this.render();
  }

  /** "Insert a feature here" (design §12): a new last migration that
   *  inserts one feature right after history step `step` of the final
   *  tree `steps` (of assembly part `part`, when given), opened in the
   *  editor. */
  async insertAt(step, steps, part = null) {
    const p = this.page.projects.current;
    if (!p || p.readOnly) return this.page.status("make a copy to edit before inserting a feature", true);
    const unique = steps.filter((s) => s.label === step.label).length === 1;
    const last = migrationPaths(p.files).map(parseMigration).at(-1);
    const number = last ? Math.floor(last.number / 10) * 10 + 10 : 1;
    const path = migrationPath(number, insertAfterName(part ? `${part} ${step.label}` : step.label));
    if (path in p.files) return this.page.status(`a migration named ${path} exists`, true);
    const source = part ? insertAfterPartMigration(part, step.label, step.index, unique) : insertAfterMigration(step.label, step.index, unique);
    await this.page.projects.replaceFiles({ ...p.files, [path]: source }, path);
    this.render();
  }

  async rename(path) {
    const p = this.page.projects.current;
    const parsed = parseMigration(path);
    const name = window.prompt("New name of this migration", parsed.name);
    if (name == null) return;
    const target = migrationPath(parsed.number, cleanName(name));
    if (target === path) return;
    if (target in p.files) return this.page.status(`a migration named ${target} exists`, true);
    const files = { ...p.files, [target]: p.files[path] };
    delete files[path];
    await this.page.projects.replaceFiles(files, p.path === path ? target : p.path);
    this.render();
  }

  /** Swap `path` with its neighbour `delta` (-1 earlier, +1 later) by
   *  exchanging their numbers. */
  async move(path, delta) {
    const p = this.page.projects.current;
    const paths = migrationPaths(p.files);
    const i = paths.indexOf(path);
    const j = i + delta;
    if (j < 0 || j >= paths.length) return;
    const a = parseMigration(paths[i]);
    const b = parseMigration(paths[j]);
    const files = { ...p.files };
    const aTo = migrationPath(b.number, a.name);
    const bTo = migrationPath(a.number, b.name);
    delete files[paths[i]];
    delete files[paths[j]];
    files[aTo] = p.files[paths[i]];
    files[bTo] = p.files[paths[j]];
    const openPath = p.path === paths[i] ? aTo : p.path === paths[j] ? bTo : p.path;
    await this.page.projects.replaceFiles(files, openPath);
    this.render();
    this.page.run();
  }

  async remove(path) {
    const p = this.page.projects.current;
    if (migrationPaths(p.files).length < 2) return;
    if (!window.confirm(`Delete ${path.split("/").pop()}? (Export the project first to keep a copy.)`)) return;
    const files = { ...p.files };
    delete files[path];
    const open = p.path === path ? migrationPaths(files).at(-1) : p.path;
    await this.page.projects.replaceFiles(files, open);
    this.render();
    this.page.run();
  }
}
