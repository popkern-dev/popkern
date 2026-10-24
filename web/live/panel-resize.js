// Resizable panels. On a wide screen: a splitter in the gap between the
// script card and the model, one between the model and the inspector,
// and one between the model and the history dock. Drag one (or focus it
// and use the arrow keys) to move the panel's edge; double-click resets
// it. On a phone the bottom sheet's grabber drags the sheet taller or
// shorter (a tap still toggles it). Sizes are view state, remembered in
// this browser (local-store.js) and clamped so the model always keeps
// room.

import { readPanelSizes, writePanelSizes } from "./local-store.js";

const $ = (id) => document.getElementById(id);

// Each splitter: the panel it sizes, the CSS variable that carries the
// size (style.css reads it in the workbench's tracks), the axis, which
// way a drag grows it, and the panel's smallest size in rem.
const SPLITS = [
  { id: "split-editor", panel: "editor", pane: ".pane-editor", variable: "--editor-w", axis: "x", grow: 1, minRem: 16 },
  { id: "split-side", panel: "side", pane: ".pane-side", variable: "--side-w", axis: "x", grow: -1, minRem: 15 },
  // The dock down to its title row (the step cards scroll in what is left).
  { id: "split-dock", panel: "dock", pane: ".dock", variable: "--dock-h", axis: "y", grow: -1, minRem: 3.25 },
];
// A dock lower than this shows its title row only (the step cards
// would be a clipped sliver).
const DOCK_SHORT_REM = 6;
// The model keeps at least this much room, and a side panel takes at
// most this share of the workbench's width.
const MODEL_MIN_REM = { x: 20, y: 12 };
const PANEL_MAX_SHARE = 0.5;
const KEY_STEP = 16;
// The phone sheet's range: its smallest height in rem, and the share
// of the screen it may take (the model keeps its 12rem row besides).
const SHEET_MIN_REM = 4;
const SHEET_MAX_SHARE = 0.85;

function rem() {
  return parseFloat(getComputedStyle(document.documentElement).fontSize) || 16;
}

/** Wire the splitters and the phone sheet's drag once, at page start. */
export function wirePanelResize() {
  const bench = document.querySelector(".workbench");
  if (!bench) return;
  const sizes = readPanelSizes();
  const save = () => writePanelSizes(sizes);
  wireSplitters(bench, sizes, save);
  wireSheetDrag(sizes, save);
}

function wireSplitters(bench, sizes, save) {
  const box = (el) => el.getBoundingClientRect();
  const size = (split) => {
    const pane = document.querySelector(split.pane);
    if (!pane) return 0;
    return split.axis === "x" ? box(pane).width : box(pane).height;
  };
  // The largest `split`'s panel may be now: the workbench less the
  // model's minimum (and on the width, the other side panel), the gaps
  // and padding.
  const largest = (split) => {
    const style = getComputedStyle(bench);
    if (split.axis === "y") {
      const inner = bench.clientHeight - parseFloat(style.paddingBottom) - parseFloat(style.paddingTop);
      return Math.max(split.minRem * rem(), inner - (parseFloat(style.rowGap) || 0) - MODEL_MIN_REM.y * rem());
    }
    const gap = parseFloat(style.columnGap) || 0;
    const inner = bench.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight);
    const other = SPLITS.find((s) => s.axis === "x" && s !== split);
    const otherPane = document.querySelector(other.pane);
    const otherWidth = otherPane && otherPane.offsetParent ? box(otherPane).width : 2.75 * rem();
    return Math.max(split.minRem * rem(), Math.min(inner * PANEL_MAX_SHARE, inner - otherWidth - 2 * gap - MODEL_MIN_REM.x * rem()));
  };
  const apply = (split, value) => {
    if (value == null) {
      bench.style.removeProperty(split.variable);
      if (split.panel === "dock") document.body.classList.remove("dock-sized", "dock-short");
      return;
    }
    const clamped = Math.round(Math.min(largest(split), Math.max(split.minRem * rem(), value)));
    bench.style.setProperty(split.variable, `${clamped}px`);
    if (split.panel === "dock") {
      document.body.classList.add("dock-sized");
      document.body.classList.toggle("dock-short", clamped < DOCK_SHORT_REM * rem());
    }
    $(split.id).setAttribute("aria-valuenow", String(clamped));
  };

  for (const split of SPLITS) {
    const handle = $(split.id);
    if (!handle) continue;
    handle.setAttribute("aria-valuemin", String(Math.round(split.minRem * rem())));
    if (sizes[split.panel] != null) apply(split, sizes[split.panel]);
    const along = (e) => (split.axis === "x" ? e.clientX : e.clientY);

    let start = null;
    handle.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      e.preventDefault();
      handle.setPointerCapture(e.pointerId);
      start = { at: along(e), size: size(split) };
      document.body.classList.add("resizing", `resizing-${split.axis}`);
      handle.classList.add("dragging");
    });
    handle.addEventListener("pointermove", (e) => {
      if (!start) return;
      sizes[split.panel] = start.size + split.grow * (along(e) - start.at);
      apply(split, sizes[split.panel]);
    });
    const end = () => {
      if (!start) return;
      start = null;
      document.body.classList.remove("resizing", `resizing-${split.axis}`);
      handle.classList.remove("dragging");
      sizes[split.panel] = size(split);
      save();
    };
    handle.addEventListener("pointerup", end);
    handle.addEventListener("pointercancel", end);
    handle.addEventListener("dblclick", () => {
      delete sizes[split.panel];
      apply(split, null);
      save();
    });
    const keys = split.axis === "x" ? ["ArrowLeft", "ArrowRight"] : ["ArrowUp", "ArrowDown"];
    handle.addEventListener("keydown", (e) => {
      const step = e.shiftKey ? KEY_STEP * 4 : KEY_STEP;
      // Right and down move the edge right and down.
      const by = e.key === keys[0] ? -step : e.key === keys[1] ? step : null;
      if (by == null) return;
      e.preventDefault();
      sizes[split.panel] = size(split) + split.grow * by;
      apply(split, sizes[split.panel]);
      save();
    });
  }

  // A smaller window re-clamps the remembered sizes (the model keeps
  // its minimum); a larger one gives them back up to what was chosen.
  if (typeof ResizeObserver === "function") {
    new ResizeObserver(() => {
      for (const split of SPLITS) if (sizes[split.panel] != null) apply(split, sizes[split.panel]);
    }).observe(bench);
  }
}

/** The phone sheet: dragging its grabber sets the sheet's height (one
 *  height for all three sheets); a tap without a drag keeps the
 *  grabber's toggle (taller or shorter), and clears a dragged height. */
function wireSheetDrag(sizes, save) {
  const grip = $("sheet-size");
  const tabs = document.querySelector(".sheet-tabs");
  if (!grip || !tabs) return;
  const range = () => [SHEET_MIN_REM * rem(), window.innerHeight * SHEET_MAX_SHARE];
  const apply = (height) => {
    if (height == null) {
      document.body.style.removeProperty("--sheet-h");
      document.body.classList.remove("sheet-sized");
      return;
    }
    const [lo, hi] = range();
    document.body.style.setProperty("--sheet-h", `${Math.round(Math.min(hi, Math.max(lo, height)))}px`);
    document.body.classList.add("sheet-sized");
    document.body.classList.remove("sheet-tall");
  };
  if (sizes.sheet != null) apply(sizes.sheet);
  // The sheet's height is what lies under the tabs: the screen's
  // bottom less the tabs' bottom edge.
  const sheetHeight = () => window.innerHeight - tabs.getBoundingClientRect().bottom;

  let start = null;
  let dragged = false;
  grip.addEventListener("pointerdown", (e) => {
    if (e.button !== 0) return;
    grip.setPointerCapture(e.pointerId);
    start = { y: e.clientY, height: sheetHeight() };
    dragged = false;
  });
  grip.addEventListener("pointermove", (e) => {
    if (!start) return;
    const dy = e.clientY - start.y;
    if (!dragged && Math.abs(dy) < 6) return;
    if (!dragged) document.body.classList.add("resizing", "resizing-y");
    dragged = true;
    sizes.sheet = start.height - dy;
    apply(sizes.sheet);
  });
  const end = () => {
    if (!start) return;
    start = null;
    document.body.classList.remove("resizing", "resizing-y");
    if (dragged) {
      sizes.sheet = sheetHeight();
      save();
    }
  };
  grip.addEventListener("pointerup", end);
  grip.addEventListener("pointercancel", end);
  // A tap (no drag): a dragged height goes back to the layout's own;
  // otherwise the toggle (layout.js) runs as before.
  grip.addEventListener(
    "click",
    (e) => {
      if (dragged) {
        e.stopImmediatePropagation();
        dragged = false;
        return;
      }
      if (sizes.sheet != null) {
        e.stopImmediatePropagation();
        delete sizes.sheet;
        apply(null);
        save();
      }
    },
    true,
  );
  window.addEventListener("resize", () => sizes.sheet != null && apply(sizes.sheet));
}
