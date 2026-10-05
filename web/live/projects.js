// The open project (live-projects-and-history.md §5.7, §12): which
// project the editor shows, autosaving a browser project as it is
// edited (a write names the revision it edited; a stale one shows the
// conflict bar instead of overwriting), examples opened read-only with
// "make a copy to edit", the Scratch project a visitor types into, and
// project export/import as a .zip of the folder.
import { ExampleStore } from "./example-store.js";
import { BrowserStore, ConflictError, manifestName, manifestText, migrationPaths, singleScriptProject } from "./project-store.js";
import { saveBytes } from "./exports.js";
import { fileStem } from "./document-tab.js";
import { unzipFiles, zipFiles } from "./zip.js";

const $ = (id) => document.getElementById(id);

/** The fixed id of the project a visitor types into. */
export const SCRATCH = "scratch";
// An edit is saved once typing pauses this long.
const SAVE_MS = 400;

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
    [this.store, this.examples] = await Promise.all([BrowserStore.open(), ExampleStore.open()]);
    $("storage-note").hidden = this.store.persistent;
    return this;
  }

  /** The Scratch project, created from `seed` (the legacy draft or the
   *  default script) the first time. */
  async scratch(seed) {
    return (await this.store.get(SCRATCH)) ?? this.store.create("Scratch", singleScriptProject("Scratch", "scratch", seed), SCRATCH);
  }

  /** The open migration's path (the last one: where a project grows). */
  get path() {
    return this.current?.path ?? null;
  }

  /** What to build: one script, or the whole project with the open
   *  migration as `source` (the editor's text is always what builds). */
  buildRequest(source) {
    const p = this.current;
    if (!p?.path) return null;
    const migrations = migrationPaths(p.files);
    if (migrations.length > 1) return { kind: "project", files: { ...p.files, [p.path]: source } };
    return { kind: "script", source, name: p.path.split("/").pop() };
  }

  /** Show `project` (from either store) in the editor and build it.
   *  Pending edits of the project being left are saved first. */
  async open(project, { build = true } = {}) {
    await this.flush();
    const migrations = migrationPaths(project.files);
    const path = migrations.at(-1) ?? null;
    this.current = { ...project, path, readOnly: Boolean(project.readOnly), unsaved: false };
    this.page.tab.importedSource = null;
    this.page.tab.open(project.name, path ? project.files[path] : "");
    this.page.editor.readOnly = this.current.readOnly;
    $("readonly-bar").hidden = !this.current.readOnly;
    $("conflict-bar").hidden = true;
    this.showName();
    this.page.viewer.fitted = false;
    if (build) this.page.run();
  }

  /** Open a browser project by id. */
  async openStored(id) {
    const project = await this.store.get(id);
    if (project) await this.open(project);
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

  /** The editor tab: project, migration file, and a mark while unsaved. */
  showName() {
    const p = this.current;
    if (!p) return;
    const file = p.path ? p.path.split("/").pop() : "";
    const mark = p.readOnly ? " (example, read only)" : p.unsaved ? " ●" : "";
    $("doc-name").textContent = `${p.name} · ${file}${mark}`;
    $("doc-name").title = p.readOnly ? "make a copy to edit" : p.unsaved ? "saving…" : "saved in this browser";
  }
}
