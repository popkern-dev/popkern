// The projects dialog: the visitor's projects in this browser (open,
// delete, new, import a .zip) and the examples (opened read-only).
// Replaces the example dropdown (live-projects-and-history.md §12).
import { SCRATCH } from "./projects.js";

const $ = (id) => document.getElementById(id);

function row(title, detail, actions) {
  const li = document.createElement("li");
  li.className = "project-row";
  const text = document.createElement("div");
  text.className = "project-text";
  const name = document.createElement("strong");
  name.textContent = title;
  text.append(name);
  if (detail) {
    const span = document.createElement("span");
    span.className = "muted";
    span.textContent = detail;
    text.append(span);
  }
  li.append(text);
  for (const [label, act, extra] of actions) {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = label;
    if (extra) button.className = extra;
    button.addEventListener("click", act);
    li.append(button);
  }
  return li;
}

function when(iso) {
  if (!iso) return "";
  const date = new Date(iso);
  return Number.isNaN(date.getTime()) ? "" : `saved ${date.toLocaleString()}`;
}

/** Fill and show the dialog. */
export async function showProjects(page) {
  const projects = page.projects;
  const close = () => ($("projects-dialog").hidden = true);
  const mine = $("projects-list");
  mine.replaceChildren();
  for (const p of await projects.store.list()) {
    const actions = [
      [
        "open",
        async () => {
          close();
          await projects.openStored(p.id);
        },
      ],
    ];
    if (p.kind === "folder") {
      actions.push([
        "forget",
        async () => {
          await projects.remove(p.id);
          await showProjects(page);
        },
      ]);
    } else if (p.id !== SCRATCH) {
      actions.push([
        "delete",
        async () => {
          if (!window.confirm(`Delete "${p.name}" from this browser? This cannot be undone (export it first to keep a copy).`)) return;
          await projects.remove(p.id);
          await showProjects(page);
        },
        "danger",
      ]);
    }
    const detail = p.kind === "folder" ? "folder on this computer" : when(p.updated);
    mine.append(row(p.name, detail, actions));
  }
  const examples = $("examples-list");
  examples.replaceChildren();
  for (const e of projects.examples.list()) {
    const detail = e.migrations > 1 ? `${e.migrations} migrations — ${e.blurb}` : e.blurb;
    examples.append(
      row(e.name, detail, [
        [
          "open",
          async () => {
            close();
            await projects.openExample(e.id);
          },
        ],
      ]),
    );
  }
  $("projects-dialog").hidden = false;
}

/** Wire the dialog's own buttons. */
export function wireProjects(page) {
  $("projects-open").addEventListener("click", () => showProjects(page));
  $("projects-close").addEventListener("click", () => ($("projects-dialog").hidden = true));
  $("project-new").addEventListener("click", async () => {
    const name = window.prompt("Name of the new project", "My part");
    if (!name) return;
    $("projects-dialog").hidden = true;
    await page.projects.create(name.trim() || "My part", await page.defaultScript());
  });
  $("project-import").addEventListener("click", () => $("project-import-file").click());
  $("project-folder").addEventListener("click", async () => {
    $("projects-dialog").hidden = true;
    await page.projects.openFolder();
  });
  $("project-import-file").addEventListener("change", async (e) => {
    const file = e.target.files?.[0];
    e.target.value = "";
    if (!file) return;
    $("projects-dialog").hidden = true;
    await page.projects.importZip(file.name, await file.arrayBuffer());
  });
  $("export-zip").addEventListener("click", () => page.projects.exportZip());
  $("copy-to-edit").addEventListener("click", () => page.projects.copyToEdit());
  $("conflict-theirs").addEventListener("click", () => page.projects.loadTheirs());
  $("conflict-mine").addEventListener("click", () => page.projects.keepMine());
  // Typing into a read-only example points at the copy button.
  page.editor.textarea.addEventListener("keydown", (e) => {
    if (!page.editor.readOnly || e.ctrlKey || e.metaKey || e.key.length !== 1) return;
    const bar = $("readonly-bar");
    bar.classList.remove("flash");
    void bar.offsetWidth;
    bar.classList.add("flash");
  });
}
