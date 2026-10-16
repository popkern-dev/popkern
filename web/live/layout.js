// The page's layout controls: popover menus (export, settings, the
// film's export), the view-mode segments (a shade, and on CAD shading a
// colouring), the inspector's tabs, the phone layout's bottom-sheet
// tabs and size, and hiding the script card. All of it is view state:
// nothing here builds, saves or changes the document.

const $ = (id) => document.getElementById(id);

// Popover menus: button id -> menu id.
const MENUS = [
  ["export-menu-button", "export-menu"],
  ["settings-menu-button", "settings-menu"],
  ["film-export-menu-button", "film-export-menu"],
];

/** Close every open menu (but `keep`). */
function closeMenus(keep = null) {
  for (const [button, menu] of MENUS) {
    if (menu === keep) continue;
    $(menu).hidden = true;
    $(button).setAttribute("aria-expanded", "false");
  }
  // On a phone an open menu is a sheet over a scrim (mobile.css).
  document.body.classList.toggle("menu-open", keep != null && !$(keep).hidden);
}

function wireMenus() {
  for (const [buttonId, menuId] of MENUS) {
    const button = $(buttonId);
    const menu = $(menuId);
    button.addEventListener("click", (e) => {
      e.stopPropagation();
      const open = menu.hidden;
      closeMenus(menuId);
      menu.hidden = !open;
      button.setAttribute("aria-expanded", String(open));
      document.body.classList.toggle("menu-open", open);
    });
    // A menu item that acts closes its menu; a field (a select, a
    // slider, a checkbox) keeps it open while it is adjusted.
    menu.addEventListener("click", (e) => {
      e.stopPropagation();
      if (e.target.closest("button[role='menuitem']")) closeMenus();
    });
  }
  document.addEventListener("click", () => closeMenus());
  $("settings-done").addEventListener("click", (e) => {
    e.stopPropagation();
    closeMenus();
  });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") closeMenus();
  });
}

/** The viewer mode for a shade and a colouring: a colouring applies
 *  to CAD shading only. */
export function viewMode(shade, color) {
  if (shade !== "cad") return shade;
  return color === "parts" ? "parts" : color === "faces" ? "shaded" : "cad";
}

function wireViewModes(page) {
  const state = { shade: "cad", color: null };
  const paint = () => {
    for (const b of document.querySelectorAll("[data-shade]")) b.setAttribute("aria-pressed", String(b.dataset.shade === state.shade));
    for (const b of document.querySelectorAll("[data-color]")) {
      b.setAttribute("aria-pressed", String(state.shade === "cad" && b.dataset.color === state.color));
      b.disabled = state.shade !== "cad";
    }
  };
  const apply = async () => {
    paint();
    const mode = viewMode(state.shade, state.color);
    if (mode === "zebra") await page.ensureAnalysisMesh();
    page.viewer.setMode(mode);
  };
  for (const b of document.querySelectorAll("[data-shade]")) {
    b.addEventListener("click", () => {
      state.shade = b.dataset.shade;
      apply();
    });
  }
  for (const b of document.querySelectorAll("[data-color]")) {
    // A second press returns to the parts' own appearance.
    b.addEventListener("click", () => {
      state.color = state.color === b.dataset.color ? null : b.dataset.color;
      apply();
    });
  }
  paint();
}

/** Show inspector tab `name` (parts | features | params | checks). */
export function showInspectorTab(name) {
  document.querySelector(".inspector").dataset.inspect = name;
  for (const t of document.querySelectorAll("[data-inspect].tab")) t.setAttribute("aria-selected", String(t.dataset.inspect === name));
}

/** Show bottom-sheet tab `name` (history | code | model) on the phone
 *  layout; on a wide screen every panel is already in view. */
export function showSheet(name) {
  document.body.dataset.sheet = name;
  for (const t of document.querySelectorAll(".sheet-tabs [data-sheet]")) t.setAttribute("aria-selected", String(t.dataset.sheet === name));
}

/** Whether the phone layout (the bottom sheet) is in use. */
export function sheetLayout() {
  return window.matchMedia("(max-width: 60rem)").matches;
}

function setEditorCollapsed(collapsed) {
  document.body.classList.toggle("editor-collapsed", collapsed);
  $("editor-expand").hidden = !collapsed;
  $("editor-collapse").setAttribute("aria-expanded", String(!collapsed));
}

/** Wire the layout controls once, at page start. */
export function wireLayout(page) {
  wireMenus();
  wireViewModes(page);
  for (const t of document.querySelectorAll("[data-inspect].tab")) t.addEventListener("click", () => showInspectorTab(t.dataset.inspect));
  for (const t of document.querySelectorAll(".sheet-tabs [data-sheet]")) t.addEventListener("click", () => showSheet(t.dataset.sheet));
  $("sheet-size").addEventListener("click", () => document.body.classList.toggle("sheet-tall"));
  $("editor-collapse").addEventListener("click", () => setEditorCollapsed(true));
  $("editor-expand").addEventListener("click", () => {
    setEditorCollapsed(false);
    page.editor.textarea.focus();
  });
}
