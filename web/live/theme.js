// The colour theme (themes.css): CAD (the default), Midnight or Light,
// chosen in the view menu and remembered in this browser. View state
// only: nothing here builds or changes the document.

import { writeTheme } from "./local-store.js";

const $ = (id) => document.getElementById(id);

/** The themes, in menu order: `[id, label]`. */
export const THEMES = [
  ["cad", "CAD"],
  ["midnight", "Midnight"],
  ["light", "Light"],
];

/** The theme on screen. */
export function currentTheme() {
  const theme = document.documentElement.dataset.theme;
  return THEMES.some(([id]) => id === theme) ? theme : "cad";
}

/** Put theme `id` on screen and remember it. */
export function setTheme(id) {
  if (!THEMES.some(([t]) => t === id)) return;
  document.documentElement.dataset.theme = id;
  // The browser chrome follows the page's ground.
  const meta = document.querySelector('meta[name="theme-color"]');
  if (meta) meta.content = getComputedStyle(document.documentElement).getPropertyValue("--bg").trim() || meta.content;
  writeTheme(id);
  const select = $("theme");
  if (select) select.value = id;
}

/** Wire the Theme field once, at page start. */
export function wireTheme() {
  const select = $("theme");
  for (const [id, label] of THEMES) {
    const option = document.createElement("option");
    option.value = id;
    option.textContent = label;
    select.append(option);
  }
  select.value = currentTheme();
  select.addEventListener("change", () => setTheme(select.value));
  setTheme(currentTheme());
}
