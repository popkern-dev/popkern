// The open project (live-projects-and-history.md §5.7, §12): which
// project the editor shows, autosaving a browser project as it is
// edited (a write names the revision it edited; a stale one shows the
// conflict bar instead of overwriting), examples opened read-only with
// "make a copy to edit", the Scratch project a visitor types into, and
// project export/import as a .zip of the folder.
import { ExampleStore } from "./example-store.js";
import { isFilmPath } from "./film-files.js";
import { FolderStore, RoutedStore } from "./folder-store.js";
import { ServerStore } from "./host-mode.js";
import { BrowserStore, ConflictError, manifestName, manifestText, migrationPaths, singleScriptProject } from "./project-store.js";
import { saveBytes } from "./exports.js";
import { fileStem } from "./document-tab.js";
import { unzipFiles, zipFiles } from "./zip.js";

const $ = (id) => document.getElementById(id);

/** The fixed id of the project a visitor types into. */
export const SCRATCH = "scratch";
/** Where the default model's film lives in Scratch. */
export const DEFAULT_FILM_PATH = "films/tour.rhai";
// An edit is saved once typing pauses this long.
const SAVE_MS = 400;
// How often the open server project's revision is checked.
const DISK_POLL_MS = 2000;

export class Projects {
  constructor(page) {
    this.page = page;
    this.current = null;
    this.saveTimer = null;
    this.saving = Promise.resolve();
  }

  /** Open both stores; the page shows a "not saved" mark when the
   *  browser blocks persistent storage. */
  async init() {
    const mode = await this.page.modeReady;
    const store = mode?.storage === "server" ? Promise.resolve(new ServerStore()) : BrowserStore.open();
    const [primary, folders, examples] = await Promise.all([store, FolderStore.open(), ExampleStore.open()]);
    this.store = new RoutedStore(primary, folders);
    this.examples = examples;
    $("storage-note").hidden = this.store.persistent;
    // A folder edited elsewhere (an editor, `git pull`, an agent over
    // MCP) is noticed when the page comes back to the front; a project
    // on the server is also polled while the page is visible, so an
    // agent's edit shows up as it lands.
    window.addEventListener("focus", () => this.checkFolder());
    if (this.store.kind === "server") {
      setInterval(() => {
        if (document.visibilityState === "visible") this.checkFolder();
      }, DISK_POLL_MS);
    }
    return this;
  }

  /** The Scratch project, created from `seed` (the legacy draft or the
   *  default script) the first time, with the default model's `film`
   *  (`films/tour.rhai`) when one is given. */
  async scratch(seed, film = null) {
    const files = { ...singleScriptProject("Scratch", "scratch", seed), ...(film ? { [DEFAULT_FILM_PATH]: film } : {}) };
    return (await this.store.get(SCRATCH)) ?? this.store.create("Scratch", files, SCRATCH);
  }

  /** The open migration's path (the last one: where a project grows). */
  get path() {
    return this.current?.path ?? null;
  }

  /** Whether the open file is a film (`films/<name>.rhai`). */
  get filmOpen() {
    return isFilmPath(this.current?.path);
  }

  /** The open project's films' paths, by name. */
  filmPaths() {
    return Object.keys(this.current?.files ?? {}).filter(isFilmPath).sort();
  }

  /** The geometry a film plays: the migrations' text (the editor may be
   *  showing a film), as the key a film rerun checks before rebuilding. */
  geometrySource() {
    const p = this.current;
    const migrations = migrationPaths(p?.files ?? {});
    if (!migrations.length) return "";
    if (migrations.length === 1) return p.files[migrations[0]];
    return JSON.stringify(migrations.map((m) => [m, p.files[m]]));
  }

  /** What to build: one script, or the whole project with the open
   *  migration as `source` (the editor's text is always what builds);
   *  with a film open, the project's own migrations. */
  buildRequest(source) {
    const p = this.current;
    if (!p?.path) return null;
    const migrations = migrationPaths(p.files);
    if (this.filmOpen) {
      if (migrations.length > 1) return { kind: "project", files: { ...p.files, [p.path]: source } };
      return migrations.length ? { kind: "script", source: p.files[migrations[0]], name: migrations[0].split("/").pop() } : null;
    }
    if (migrations.length > 1) return { kind: "project", files: { ...p.files, [p.path]: source } };
    return { kind: "script", source, name: p.path.split("/").pop() };
  }

  /** Show `project` (from either store) in the editor and build it.
   *  Pending edits of the project being left are saved first. */
  async open(project, { build = true, path: keep = null } = {}) {
    await this.flush();
    // A film plays on the project it came from.
    this.page.film?.close();
    const migrations = migrationPaths(project.files);
    const path = keep && keep in project.files ? keep : migrations.at(-1) ?? null;
    this.current = { ...project, path, readOnly: Boolean(project.readOnly), unsaved: false };
    this.page.tab.importedSource = null;
    this.page.tab.open(project.name, path ? project.files[path] : "");
    this.page.editor.readOnly = this.current.readOnly;
    $("readonly-bar").hidden = !this.current.readOnly;
    $("conflict-bar").hidden = true;
    this.showName();
    this.page.rail?.setVersions(null);
    this.page.viewer.fitted = false;
    // The project's saved simulation setup, if any (PK-948).
    this.page.simulate?.restore(project.files);
    // The model on screen is the previous project's until this one builds.
    this.page.geometryKey = null;
    if (build) this.page.run();
  }

  /** Show migration `path` of the open project in the editor (the
   *  project builds as a whole, so nothing rebuilds). */
  select(path) {
    const p = this.current;
    if (!p || !(path in p.files)) return;
    const wasFilm = this.filmOpen;
    p.path = path;
    this.page.tab.open(p.name, p.files[path]);
    this.showName();
    this.page.rail?.render();
    // Opening a film plays it on the model; leaving one closes it.
    if (this.filmOpen) this.page.run();
    else if (wasFilm) this.page.film.close();
  }

  /** Replace the open project's files (migrations added, renamed,
   *  reordered or deleted), save at once, and show `openPath`. */
  async replaceFiles(files, openPath) {
    const p = this.current;
    if (!p || p.readOnly) return;
    clearTimeout(this.saveTimer);
    this.saveTimer = null;
    await this.saving;
    p.files = files;
    p.unsaved = true;
    this.select(openPath in files ? openPath : migrationPaths(files).at(-1));
    await this.save();
  }

  /** Open a browser project by id. */
  async openStored(id) {
    const project = await this.store.get(id);
    if (project) await this.open(project);
  }

  /** Pick a project folder on this computer and open it (Chromium);
   *  elsewhere say why not and what to do instead. */
  async openFolder() {
    const folders = await this.store.ensureFolders();
    if (!folders) {
      this.page.status("this browser cannot open folders (Chrome and Edge can); projects stay in this browser — export one as a .zip to take it elsewhere", true);
      return null;
    }
    let project;
    try {
      project = await folders.pick();
    } catch (error) {
      this.page.status(`cannot open the folder: ${error.message}`, true);
      return null;
    }
    if (project) await this.open(project);
    return project;
  }

  /** The open folder project changed on disk: reload it when nothing
   *  is pending here, show the conflict bar when edits are. */
  async checkFolder() {
    const p = this.current;
    const source = p?.kind === "folder" ? this.store.folders : p?.kind === "server" ? this.store.primary : null;
    if (!source || this.checking || this.conflict) return;
    this.checking = true;
    try {
      await this.saving;
      // While typing, the next save's revision check decides.
      if (this.current !== p || this.saveTimer) return;
      const saves = this.saves;
      const revision = await source.revision(p.id).catch(() => null);
      // A save that started meanwhile moved the revision itself.
      if (saves !== this.saves || this.current !== p) return;
      if (revision == null || revision === p.revision) return;
      const theirs = await source.get(p.id, { ask: false });
      if (!theirs) return;
      if (!p.unsaved) {
        await this.open(theirs, { path: p.path });
        this.page.status(`reloaded ${p.name}: its folder changed on disk`);
        return;
      }
      this.showConflict(theirs, p);
    } finally {
      this.checking = false;
    }
  }

  showConflict(theirs, p) {
    this.conflict = theirs;
    $("conflict-text").textContent = `"${p.name}" changed on disk (another editor, git, or an agent) while you were editing it here`;
    $("conflict-bar").hidden = false;
  }

  /** Open example `id` (`example:<file>`) read-only. */
  async openExample(id) {
    const project = await this.examples.get(id.startsWith("example:") ? id : `example:${id}`);
    if (project) await this.open(project);
  }

  /** The editor changed: keep the text and save it shortly. */
  edited(text) {
    const p = this.current;
    if (!p || p.readOnly || !p.path) return;
    p.files = { ...p.files, [p.path]: text };
    p.unsaved = true;
    this.showName();
    clearTimeout(this.saveTimer);
    this.saveTimer = setTimeout(() => this.save(), SAVE_MS);
  }

  /** Save pending edits now (before switching projects). */
  async flush() {
    if (this.saveTimer) {
      clearTimeout(this.saveTimer);
      this.saveTimer = null;
      await this.save();
    }
    await this.saving;
  }

  /** Write the open project; a stale revision shows the conflict bar. */
  save() {
    const p = this.current;
    this.saveTimer = null;
    if (!p || p.readOnly || !p.unsaved) return this.saving;
    this.saves = (this.saves ?? 0) + 1;
    this.saving = this.store
      .save(p.id, p.files, p.revision)
      .then((saved) => {
        if (this.current?.id !== saved.id) return;
        this.current.revision = saved.revision;
        this.current.unsaved = false;
        this.showName();
      })
      .catch((error) => {
        if (error instanceof ConflictError) {
          this.conflict = error.project;
          $("conflict-text").textContent = error.message;
          $("conflict-bar").hidden = false;
        } else {
          this.page.status(`not saved: ${error.message}`, true);
        }
      });
    return this.saving;
  }

  /** Conflict resolution: take the other tab's version… */
  async loadTheirs() {
    if (this.conflict) await this.open(this.conflict);
    this.conflict = null;
  }

  /** …or overwrite it with this tab's (an explicit choice). */
  async keepMine() {
    if (!this.conflict || !this.current) return;
    this.current.revision = this.conflict.revision;
    this.conflict = null;
    $("conflict-bar").hidden = true;
    this.current.unsaved = true;
    await this.save();
  }

  /** Copy the open (read-only) project into the browser store and open
   *  the copy; the example itself is never changed. */
  async copyToEdit() {
    const p = this.current;
    if (!p) return;
    const name = p.readOnly ? p.name : `${p.name} (copy)`;
    const files = { ...p.files, "popkern.toml": manifestText(name) };
    await this.open(await this.store.create(name, files));
  }

  /** A new one-migration project. */
  async create(name, source) {
    await this.open(await this.store.create(name, singleScriptProject(name, "base", source)));
  }

  async remove(id) {
    if (this.current?.id === id) this.current = null;
    await this.store.remove(id);
  }

  /** Delete project `id`; when it was the open one, open the most
   *  recent remaining project, or a fresh Scratch from `seed` (the
   *  default script, with its `film`) when none is left. A deleted Scratch is never
   *  recreated behind the person's back otherwise: the next visit
   *  starts a fresh one (openFirstScript). */
  async removeAndReopen(id, seed, film = null) {
    const wasOpen = this.current?.id === id;
    await this.remove(id);
    if (!wasOpen) return;
    const rest = (await this.store.list()).filter((p) => p.kind !== "folder");
    if (rest.length) await this.openStored(rest[0].id);
    else await this.open(await this.scratch(seed, film));
  }

  /** Download the open project as `<name>.zip`. */
  exportZip() {
    const p = this.current;
    if (!p) return;
    saveBytes(zipFiles(p.files), `${fileStem(p.name)}.zip`);
  }

  /** Import a project .zip as a new browser project and open it. */
  async importZip(fileName, bytes) {
    let files;
    try {
      const raw = unzipFiles(new Uint8Array(bytes));
      const decoder = new TextDecoder("utf-8", { fatal: true });
      files = Object.fromEntries(Object.entries(raw).map(([path, data]) => [path, decoder.decode(data)]));
    } catch (error) {
      this.page.status(`cannot import ${fileName}: ${error.message}`, true);
      return null;
    }
    if (!files["popkern.toml"] || migrationPaths(files).length === 0) {
      this.page.status(`cannot import ${fileName}: a project has popkern.toml and migrations/NNN_name.rhai`, true);
      return null;
    }
    const name = manifestName(files["popkern.toml"]) ?? fileName.replace(/\.zip$/i, "");
    const project = await this.store.create(name, files);
    await this.open(project);
    return project;
  }

  /** The project in the top bar (its name, where it is kept) and the
   *  editor tab: the open file, and the save state beside it. */
  showName() {
    const p = this.current;
    if (!p) return;
    const file = p.path ? p.path.split("/").pop() : "";
    $("doc-name").textContent = file;
    const where = p.kind === "folder" ? `saved in the folder ${p.folder}` : p.kind === "server" ? `saved in ${p.path ?? "the server's projects folder"}` : "saved in this browser";
    $("doc-name").title = p.readOnly ? "an example: make a copy to edit" : p.unsaved ? "saving…" : where;
    $("doc-name").closest(".editor-head").dataset.state = p.readOnly ? "readonly" : p.unsaved ? "unsaved" : "saved";
    $("project-name").textContent = p.name;
    const chip = $("project-where");
    chip.hidden = false;
    chip.textContent = p.readOnly ? "example" : p.kind === "folder" ? "folder" : p.kind === "server" ? "server" : this.store.persistent ? "browser" : "this tab";
    chip.title = p.readOnly ? "read only: make a copy to edit" : where;
  }
}
