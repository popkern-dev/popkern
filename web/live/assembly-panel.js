// The inspector's assembly tree (docs/architecture/assemblies.md §6.5):
// every occurrence indented by depth, a visibility toggle, a click that
// selects (emphasises) the occurrence in the viewport, and "open" on a
// part occurrence that scopes the feature tree and the parameters to
// that part's definition (a slider then rebuilds the part once and every
// instance of it moves with it).

import { iconSvg } from "./icons.js";
import { paletteCss } from "./shading.js";

const $ = (id) => document.getElementById(id);

export class AssemblyPanel {
  /** `on.select(path)`, `on.visible(path, visible)`, `on.open(partName | null)`. */
  constructor(on) {
    this.on = on;
    this.view = null;
    this.hidden = new Set();
    this.selected = null;
    this.openPart = null;
    $("assembly-close").addEventListener("click", () => this.open(null));
  }

  /** Show `view` (the envelope's assembly view, or null for a part
   *  document). Visibility the user set survives a rebuild of the same
   *  tree; a different tree starts fresh. */
  render(view) {
    const section = $("assembly-section");
    const same = this.view && view && sameTree(this.view, view);
    this.view = view;
    if (!view) {
      section.hidden = true;
      this.hidden.clear();
      this.selected = null;
      this.setOpen(null);
      return;
    }
    if (!same) {
      // A different tree: earlier clearance results no longer apply.
      $("clearance-list").hidden = true;
      $("clearance-summary").hidden = true;
      this.hidden = new Set(view.instances.filter((i) => !i.visible).map((i) => i.path));
      this.selected = null;
      this.setOpen(null);
    }
    section.hidden = false;
    const parts = view.definitions.filter((d) => d.kind === "part").length;
    const occurrences = view.instances.filter((i) => i.is_part).length;
    $("assembly-summary").textContent = `${view.name}: ${parts} part(s), ${occurrences} occurrence(s)`;
    const list = $("assembly-tree");
    list.innerHTML = "";
    for (const instance of view.instances) {
      const definition = view.definitions[instance.definition];
      const row = document.createElement("li");
      row.className = "tree-row";
      row.dataset.path = instance.path;
      row.style.paddingLeft = `${0.3 + instance.depth * 1.1}em`;
      const eye = document.createElement("button");
      eye.type = "button";
      eye.className = "tree-eye";
      eye.title = "Show or hide";
      eye.addEventListener("click", (e) => {
        e.stopPropagation();
        this.toggle(instance.path);
      });
      // The part's colour in the Parts colouring (one per definition).
      const dot = document.createElement("span");
      dot.className = "tree-dot";
      dot.style.background = paletteCss(instance.definition);
      const name = document.createElement("span");
      name.className = "tree-name";
      name.textContent = instance.name;
      const kind = document.createElement("span");
      kind.className = "tree-kind";
      kind.textContent = definition ? definition.name : "?";
      row.append(eye, dot, name, kind);
      if (definition?.kind === "part") {
        const open = document.createElement("button");
        open.type = "button";
        open.className = "tree-open";
        open.textContent = "open";
        open.title = `Show the features and parameters of ${definition.name}`;
        open.addEventListener("click", (e) => {
          e.stopPropagation();
          this.open(definition.name);
        });
        row.append(open);
      }
      row.addEventListener("click", () => this.select(this.selected === instance.path ? null : instance.path));
      list.append(row);
    }
    this.paint();
  }

  /** Visibility is inherited: a row is shown only when it and every
   *  ancestor are. */
  shown(path) {
    const names = path.split("/");
    for (let i = 1; i <= names.length; i++) {
      if (this.hidden.has(names.slice(0, i).join("/"))) return false;
    }
    return true;
  }

  toggle(path) {
    if (this.hidden.has(path)) this.hidden.delete(path);
    else this.hidden.add(path);
    this.on.visible(path, !this.hidden.has(path));
    this.paint();
  }

  select(path) {
    this.selected = path;
    this.on.select(path);
    this.paint();
  }

  open(partName) {
    this.setOpen(partName);
    this.on.open(partName);
  }

  setOpen(partName) {
    this.openPart = partName;
    $("assembly-open").hidden = partName == null;
    $("assembly-open-name").textContent = partName ?? "";
  }

  paint() {
    for (const row of $("assembly-tree").children) {
      const path = row.dataset.path;
      const own = !this.hidden.has(path);
      row.classList.toggle("tree-hidden", !this.shown(path));
      row.classList.toggle("tree-selected", path === this.selected);
      const eye = row.querySelector(".tree-eye");
      eye.innerHTML = iconSvg(own ? "eye" : "eye-off");
      eye.setAttribute("aria-pressed", String(own));
      eye.setAttribute("aria-label", `${own ? "Hide" : "Show"} ${path}`);
    }
  }
}

function sameTree(a, b) {
  return a.instances.length === b.instances.length && a.instances.every((x, i) => x.path === b.instances[i].path && x.definition === b.instances[i].definition);
}
