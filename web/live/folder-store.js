// Folder-backed projects (live-projects-and-history.md §5.7, E-fsapi,
// PK-668): a project folder on this computer — `popkern.toml` and
// `migrations/NNN_name.rhai` — opened through the File System Access
// API (Chromium) and edited in place, no server. The page remembers
// opened folders (their handles, in IndexedDB) and asks again for
// permission after a reload, on the click that reopens one.
//
// A revision is the SHA-256 of the project's files as they are on
// disk, so a save is checked against the folder itself: when the files
// changed underneath (a `git pull`, another editor), the save is a
// ConflictError carrying what is on disk, never a silent overwrite.
// Each file is written through a swap file the browser commits on
// close; a save of several files is not atomic as a whole, and a save
// interrupted half-way is seen as a conflict by the next one.
//
// "Forget" drops the remembered folder; it never deletes files.
import { ConflictError, manifestName, manifestText } from "./project-store.js";
import { firstMigration } from "./migration-templates.js";

const DB_NAME = "popkern-live-folders";
const STORE = "folders";
/** Ids of folder projects start with this. */
export const FOLDER_PREFIX = "folder:";
const MIGRATION = /^\d+_[A-Za-z0-9_-]+\.rhai$/;

function promise(request) {
  return new Promise((resolve, reject) => {
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
  });
}

/** Whether this browser can open folders. */
export function foldersSupported() {
  return typeof globalThis.showDirectoryPicker === "function";
}

/** The revision of `files`: SHA-256 over the sorted paths and texts. */
export async function revisionOf(files) {
  const text = Object.keys(files)
    .sort()
    .map((path) => `${path}\u0000${files[path]}\u0000`)
    .join("");
  const digest = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text));
  return [...new Uint8Array(digest)].map((b) => b.toString(16).padStart(2, "0")).join("");
}

/** The project files in `dir`: the manifest and the migrations. */
async function readFolder(dir) {
  const files = {};
  for await (const [name, entry] of dir.entries()) {
    if (entry.kind === "file" && name === "popkern.toml") {
      files[name] = await (await entry.getFile()).text();
    } else if (entry.kind === "directory" && name === "migrations") {
      for await (const [file, handle] of entry.entries()) {
        if (handle.kind === "file" && MIGRATION.test(file)) files[`migrations/${file}`] = await (await handle.getFile()).text();
      }
    }
  }
  return files;
}

async function writeText(dir, path, text) {
  const parts = path.split("/");
  let parent = dir;
  for (const part of parts.slice(0, -1)) parent = await parent.getDirectoryHandle(part, { create: true });
  const file = await parent.getFileHandle(parts.at(-1), { create: true });
  const out = await file.createWritable();
  await out.write(text);
  await out.close();
}

async function removeFile(dir, path) {
  const parts = path.split("/");
  let parent = dir;
  for (const part of parts.slice(0, -1)) parent = await parent.getDirectoryHandle(part);
  await parent.removeEntry(parts.at(-1));
}

/** Read and write access to `handle`, asking when `ask` (a click). */
async function access(handle, ask) {
  const options = { mode: "readwrite" };
  if ((await handle.queryPermission?.(options)) === "granted" || !handle.queryPermission) return true;
  return ask && (await handle.requestPermission(options)) === "granted";
}

/** The remembered project folders. */
export class FolderStore {
  constructor(db) {
    this.db = db;
    this.persistent = true;
    this.kind = "folder";
  }

  /** The store, or null where folders cannot be opened or remembered. */
  static async open() {
    if (!foldersSupported()) return null;
    try {
      const request = indexedDB.open(DB_NAME, 1);
      request.onupgradeneeded = () => request.result.createObjectStore(STORE, { keyPath: "id" });
      return new FolderStore(await promise(request));
    } catch {
      return null;
    }
  }

  async entry(id) {
    return (await promise(this.db.transaction(STORE).objectStore(STORE).get(id))) ?? null;
  }

  async list() {
    const all = await promise(this.db.transaction(STORE).objectStore(STORE).getAll());
    return all.map((e) => ({ id: e.id, name: e.name, updated: e.opened, files: null, kind: "folder" }));
  }

  /** The project in a remembered folder; `ask` permits a permission
   *  prompt (the call comes from a click). */
  async get(id, { ask = true } = {}) {
    const entry = await this.entry(id);
    if (!entry) return null;
    if (!(await access(entry.handle, ask))) throw new Error(`popkern live may not read the folder "${entry.handle.name}" (permission not granted)`);
    const files = await readFolder(entry.handle);
    return { id, name: manifestName(files["popkern.toml"]) ?? entry.name, files, revision: await revisionOf(files), kind: "folder", folder: entry.handle.name };
  }

  /** Pick a folder and open it as a project: a folder without a
   *  manifest becomes one (a manifest and a first migration are
   *  written). Resolves to the project, or null when the picker is
   *  dismissed. */
  async pick() {
    let handle;
    try {
      handle = await globalThis.showDirectoryPicker({ mode: "readwrite", id: "popkern-project" });
    } catch (error) {
      if (error?.name === "AbortError") return null;
      throw error;
    }
    let files = await readFolder(handle);
    if (!files["popkern.toml"]) {
      const name = handle.name;
      const seed = { "popkern.toml": manifestText(name) };
      if (!Object.keys(files).length) seed["migrations/001_base.rhai"] = firstMigration(name);
      for (const [path, text] of Object.entries(seed)) await writeText(handle, path, text);
      files = await readFolder(handle);
    }
    // The same folder picked again reopens its project.
    for (const e of await promise(this.db.transaction(STORE).objectStore(STORE).getAll())) {
      if (await e.handle.isSameEntry?.(handle)) return this.get(e.id);
    }
    const id = `${FOLDER_PREFIX}${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
    const name = manifestName(files["popkern.toml"]) ?? handle.name;
    await promise(this.db.transaction(STORE, "readwrite").objectStore(STORE).put({ id, name, handle, opened: new Date().toISOString() }));
    return this.get(id);
  }

  /** Write `files` to the folder, provided its files are still those of
   *  `revision`; migrations no longer in `files` are removed (a rename
   *  or a delete in the rail). */
  async save(id, files, revision) {
    const entry = await this.entry(id);
    if (!entry) throw new Error("the folder is no longer remembered here");
    if (!(await access(entry.handle, false))) throw new Error(`popkern live may no longer write to the folder "${entry.handle.name}"; reopen it from Projects`);
    const disk = await readFolder(entry.handle);
    const current = await revisionOf(disk);
    if (current !== revision) {
      const theirs = { id, name: manifestName(disk["popkern.toml"]) ?? entry.name, files: disk, revision: current, kind: "folder", folder: entry.handle.name };
      throw new ConflictError(theirs, "its folder changed on disk: another editor, or git");
    }
    for (const [path, text] of Object.entries(files)) {
      if (disk[path] !== text) await writeText(entry.handle, path, text);
    }
    for (const path of Object.keys(disk)) {
      if (!(path in files)) await removeFile(entry.handle, path);
    }
    const name = manifestName(files["popkern.toml"]) ?? entry.name;
    return { id, name, files: { ...files }, revision: await revisionOf(files), kind: "folder", folder: entry.handle.name };
  }

  /** The folder's revision on disk now (an external change shows up as
   *  a different one), or null when it cannot be read. */
  async revision(id) {
    const entry = await this.entry(id);
    if (!entry || !(await access(entry.handle, false))) return null;
    return revisionOf(await readFolder(entry.handle));
  }

  /** Forget the folder (its files stay where they are). */
  async remove(id) {
    await promise(this.db.transaction(STORE, "readwrite").objectStore(STORE).delete(id));
  }
}

/** One store over the page's own (browser or server) and the folders:
 *  `folder:` ids go to the folders, everything else (and every create)
 *  to the primary. */
export class RoutedStore {
  constructor(primary, folders) {
    this.primary = primary;
    this.folders = folders;
    this.persistent = primary.persistent;
    this.kind = primary.kind;
  }

  /** The folder store, opened on first use (null where folders cannot
   *  be opened). */
  async ensureFolders() {
    if (!this.folders) this.folders = await FolderStore.open();
    return this.folders;
  }

  route(id) {
    return this.folders && String(id).startsWith(FOLDER_PREFIX) ? this.folders : this.primary;
  }

  async list() {
    const own = await this.primary.list();
    return this.folders ? [...(await this.folders.list()), ...own] : own;
  }

  get(id, options) {
    return this.route(id).get(id, options);
  }

  create(name, files, id) {
    return this.primary.create(name, files, id);
  }

  save(id, files, revision, name) {
    return this.route(id).save(id, files, revision, name);
  }

  remove(id) {
    return this.route(id).remove(id);
  }
}
