// Which deployment the page is in (live-projects-and-history.md §5.7):
// the public page has no `host.json` (projects in this browser, kernel
// in a Web Worker); `popkern serve` answers it with where projects live
// and where the kernel runs. One static page serves both.
import { ConflictError } from "./project-store.js";

/** `{ storage: "browser" | "server", compute: "browser" | "native",
 *  version, root, mcp }`: under `popkern serve`, `root` is the projects
 *  folder on the serving machine and `mcp` its MCP endpoint path. */
export async function hostMode() {
  try {
    const response = await fetch("./host.json", { cache: "no-store" });
    if (!response.ok) return { storage: "browser", compute: "browser" };
    const host = await response.json();
    return {
      storage: host.storage === "server" ? "server" : "browser",
      compute: host.compute === "native" ? "native" : "browser",
      version: host.version ?? null,
      root: typeof host.root === "string" ? host.root : null,
      mcp: typeof host.mcp === "string" ? host.mcp : null,
    };
  } catch {
    return { storage: "browser", compute: "browser" };
  }
}

async function answer(response) {
  if (response.status === 409) throw new ConflictError((await response.json()).project);
  if (!response.ok) throw new Error((await response.json().catch(() => ({}))).error ?? `HTTP ${response.status}`);
  return response.status === 204 ? null : response.json();
}

const JSON_HEADERS = { "Content-Type": "application/json" };

/** A project as the server answers it, marked as a folder on the
 *  serving machine (`path` is that folder). */
function onServer(project) {
  return { ...project, kind: "server" };
}

/** Projects in the server's folder (`./api/projects`): the same
 *  interface as BrowserStore; revisions are content hashes the server
 *  checks, so a stale save (another tab, `git pull`, an editor) is a
 *  ConflictError carrying the stored project. */
export class ServerStore {
  constructor() {
    this.persistent = true;
    this.kind = "server";
  }

  async list() {
    return (await answer(await fetch("./api/projects"))).map(onServer);
  }

  async get(id) {
    const response = await fetch(`./api/projects/${encodeURIComponent(id)}`);
    return response.status === 404 ? null : onServer(await answer(response));
  }

  /** `id` is ignored: the server names folders after the project. The
   *  scratch project is the one called "Scratch". */
  async create(name, files) {
    return onServer(await answer(await fetch("./api/projects", { method: "POST", headers: JSON_HEADERS, body: JSON.stringify({ name, files }) })));
  }

  async save(id, files, revision) {
    return onServer(await answer(await fetch(`./api/projects/${encodeURIComponent(id)}`, { method: "PUT", headers: JSON_HEADERS, body: JSON.stringify({ files, revision }) })));
  }

  /** The stored revision of `id` (null when it is gone): what the page
   *  polls to notice an edit made on disk (an editor, git, an agent
   *  through MCP). */
  async revision(id) {
    return (await this.get(id))?.revision ?? null;
  }

  async remove(id) {
    await answer(await fetch(`./api/projects/${encodeURIComponent(id)}`, { method: "DELETE" }));
  }
}
