// Where projects live (live-projects-and-history.md §5.7). A project is
// a folder — `popkern.toml` and `migrations/NNN_name.rhai`, text files
// by project-relative path — plus a revision every write must name:
// a write against a stale revision (another tab saved first) is refused
// with a ConflictError instead of silently overwriting.
//
// BrowserStore keeps projects in this origin's IndexedDB; a browser
// that blocks it (private mode, a sandboxed frame) gets MemoryStore,
// which forgets on reload and says so (`persistent: false`).

/** A write named a revision the store has moved past. */
export class ConflictError extends Error {
  constructor(project) {
    super(`"${project.name}" changed elsewhere (another tab saved it) since it was opened here`);
    this.project = project;
  }
}

const DB_NAME = "popkern-live";
const STORE = "projects";

function promise(request) {
  return new Promise((resolve, reject) => {
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error);
  });
}

function newId() {
  return `p-${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
}

/** A summary row for a project list. */
function summary(project) {
  return { id: project.id, name: project.name, updated: project.updated, files: Object.keys(project.files).length };
}

/** Projects in this origin's IndexedDB. */
export class BrowserStore {
  constructor(db) {
    this.db = db;
    this.persistent = true;
    this.kind = "browser";
  }

  /** The store, or a MemoryStore when IndexedDB is unavailable. */
  static async open() {
    try {
      const request = indexedDB.open(DB_NAME, 1);
      request.onupgradeneeded = () => request.result.createObjectStore(STORE, { keyPath: "id" });
      return new BrowserStore(await promise(request));
    } catch {
      return new MemoryStore();
    }
  }

  async list() {
    const all = await promise(this.db.transaction(STORE).objectStore(STORE).getAll());
    return all.map(summary).sort((a, b) => (b.updated ?? "").localeCompare(a.updated ?? ""));
  }

  async get(id) {
    return (await promise(this.db.transaction(STORE).objectStore(STORE).get(id))) ?? null;
  }

  /** A new project called `name` holding `files`; resolves to it. `id`
   *  may be given (the scratch project has a fixed one). */
  async create(name, files, id = newId()) {
    const project = { id, name, files: { ...files }, revision: 1, updated: new Date().toISOString() };
    await promise(this.db.transaction(STORE, "readwrite").objectStore(STORE).add(project));
    return project;
  }

  /** Replace `id`'s files, provided the stored revision is still
   *  `revision`; resolves to the updated project or rejects with a
   *  ConflictError carrying the stored one. One transaction: the check
   *  and the write cannot interleave with another tab's. */
  save(id, files, revision, name = null) {
    return new Promise((resolve, reject) => {
      const tx = this.db.transaction(STORE, "readwrite");
      const store = tx.objectStore(STORE);
      let result = null;
      const read = store.get(id);
      read.onsuccess = () => {
        const current = read.result;
        if (!current) {
          reject(new Error("the project no longer exists in this browser"));
          tx.abort();
          return;
        }
        if (current.revision !== revision) {
          reject(new ConflictError(current));
          tx.abort();
          return;
        }
        result = { ...current, files: { ...files }, name: name ?? current.name, revision: revision + 1, updated: new Date().toISOString() };
        store.put(result);
      };
      tx.oncomplete = () => resolve(result);
      tx.onerror = () => reject(tx.error);
    });
  }

  async remove(id) {
    await promise(this.db.transaction(STORE, "readwrite").objectStore(STORE).delete(id));
  }
}

/** The same interface in memory: nothing survives a reload. */
export class MemoryStore {
  constructor() {
    this.projects = new Map();
    this.persistent = false;
    this.kind = "browser";
  }

  static async open() {
    return new MemoryStore();
  }

  async list() {
    return [...this.projects.values()].map(summary).sort((a, b) => (b.updated ?? "").localeCompare(a.updated ?? ""));
  }

  async get(id) {
    const p = this.projects.get(id);
    return p ? structuredClone(p) : null;
  }

  async create(name, files, id = newId()) {
    if (this.projects.has(id)) throw new Error(`project ${id} exists`);
    const project = { id, name, files: { ...files }, revision: 1, updated: new Date().toISOString() };
    this.projects.set(id, project);
    return structuredClone(project);
  }

  async save(id, files, revision, name = null) {
    const current = this.projects.get(id);
    if (!current) throw new Error("the project no longer exists");
    if (current.revision !== revision) throw new ConflictError(structuredClone(current));
    const next = { ...current, files: { ...files }, name: name ?? current.name, revision: revision + 1, updated: new Date().toISOString() };
    this.projects.set(id, next);
    return structuredClone(next);
  }

  async remove(id) {
    this.projects.delete(id);
  }
}

/** The manifest of a project called `name`. */
export function manifestText(name) {
  return `[project]\nname = ${JSON.stringify(name)}\n`;
}

/** The project name a manifest declares, or null. */
export function manifestName(text) {
  const m = /^\s*name\s*=\s*"((?:[^"\\]|\\.)*)"/m.exec(text ?? "");
  return m ? JSON.parse(`"${m[1]}"`) : null;
}

/** The migration paths of `files`, in migration order. */
export function migrationPaths(files) {
  return Object.keys(files)
    .filter((p) => /^migrations\/\d+_[A-Za-z0-9_-]+\.rhai$/.test(p))
    .sort((a, b) => Number(/\/(\d+)_/.exec(a)[1]) - Number(/\/(\d+)_/.exec(b)[1]));
}

/** A one-migration project holding `source`. */
export function singleScriptProject(name, stem, source) {
  return { "popkern.toml": manifestText(name), [`migrations/001_${stem}.rhai`]: source };
}
