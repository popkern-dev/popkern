// Resizable panels on a wide screen: a splitter in the gap between the
// script card and the model, and one between the model and the
// inspector. Drag one (or focus it and use the arrow keys) to move the
// panel's edge; double-click resets it. Widths are view state,
// remembered in this browser (local-store.js) and clamped so the model
// always keeps room. The phone layout has no columns, and no splitters.

import { readPanelSizes, writePanelSizes } from "./local-store.js";

const $ = (id) => document.getElementById(id);

// Each splitter: the panel it sizes, the CSS variable that carries the
// width (style.css reads it in the workbench's columns), which way a
// drag grows it, and the panel's narrowest width in rem.
const SPLITS = [
  { id: "split-editor", panel: "editor", pane: ".pane-editor", variable: "--editor-w", grow: 1, minRem: 16 },
  { id: "split-side", panel: "side", pane: ".pane-side", variable: "--side-w", grow: -1, minRem: 15 },
];
// The model keeps at least this much width, and neither panel takes
// more than this share of the workbench.
const MODEL_MIN_REM = 20;
const PANEL_MAX_SHARE = 0.5;
const KEY_STEP = 16;

function rem() {
  return parseFloat(getComputedStyle(document.documentElement).fontSize) || 16;
}

/** Wire the splitters once, at page start. */
export function wirePanelResize() {
  const bench = document.querySelector(".workbench");
  if (!bench) return;
  const sizes = readPanelSizes();

  // The widest `split`'s panel may be now: the workbench less the
  // model's minimum, the other panel and the two gaps.
  const widest = (split) => {
    const gap = parseFloat(getComputedStyle(bench).columnGap) || 0;
    const inner = bench.clientWidth - parseFloat(getComputedStyle(bench).paddingLeft) - parseFloat(getComputedStyle(bench).paddingRight);
    const other = SPLITS.find((s) => s !== split);
    const otherPane = document.querySelector(other.pane);
    const otherWidth = otherPane && otherPane.offsetParent ? otherPane.getBoundingClientRect().width : 2.75 * rem();
    return Math.max(split.minRem * rem(), Math.min(inner * PANEL_MAX_SHARE, inner - otherWidth - 2 * gap - MODEL_MIN_REM * rem()));
  };
  const apply = (split, width) => {
    if (width == null) {
      bench.style.removeProperty(split.variable);
      return;
    }
    const clamped = Math.round(Math.min(widest(split), Math.max(split.minRem * rem(), width)));
    bench.style.setProperty(split.variable, `${clamped}px`);
    const handle = $(split.id);
    handle.setAttribute("aria-valuenow", String(clamped));
  };
  const save = () => writePanelSizes(sizes);
  const current = (split) => document.querySelector(split.pane)?.getBoundingClientRect().width ?? 0;

  for (const split of SPLITS) {
    const handle = $(split.id);
    if (!handle) continue;
    handle.setAttribute("aria-valuemin", String(Math.round(split.minRem * rem())));
    if (sizes[split.panel] != null) apply(split, sizes[split.panel]);

    let start = null;
    handle.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      e.preventDefault();
      handle.setPointerCapture(e.pointerId);
      start = { x: e.clientX, width: current(split) };
      document.body.classList.add("resizing");
      handle.classList.add("dragging");
    });
    handle.addEventListener("pointermove", (e) => {
      if (!start) return;
      sizes[split.panel] = start.width + split.grow * (e.clientX - start.x);
      apply(split, sizes[split.panel]);
    });
    const end = () => {
      if (!start) return;
      start = null;
      document.body.classList.remove("resizing");
      handle.classList.remove("dragging");
      sizes[split.panel] = current(split);
      save();
    };
    handle.addEventListener("pointerup", end);
    handle.addEventListener("pointercancel", end);
    handle.addEventListener("dblclick", () => {
      delete sizes[split.panel];
      apply(split, null);
      save();
    });
    handle.addEventListener("keydown", (e) => {
      const step = e.shiftKey ? KEY_STEP * 4 : KEY_STEP;
      const by = e.key === "ArrowLeft" ? -step : e.key === "ArrowRight" ? step : null;
      if (by == null) return;
      e.preventDefault();
      sizes[split.panel] = current(split) + split.grow * by;
      apply(split, sizes[split.panel]);
      save();
    });
  }

  // A smaller window re-clamps the remembered widths (the model keeps
  // its minimum); a larger one gives them back up to what was chosen.
  if (typeof ResizeObserver === "function") {
    new ResizeObserver(() => {
      for (const split of SPLITS) if (sizes[split.panel] != null) apply(split, sizes[split.panel]);
    }).observe(bench);
  }
}
