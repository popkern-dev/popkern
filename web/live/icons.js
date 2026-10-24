// Line icons for history steps (feature operations and assembly steps)
// and the few interface marks the modules draw. Each is a 16x16 stroke
// path set drawn with `currentColor`; an operation without its own
// icon gets the generic mark, so a new feature kind never breaks the
// history strip.

const PATHS = {
  extrude: "M2.5 12.5h11M4.5 12.5V7h7v5.5M8 1.8v3.6M6.4 3.4 8 1.8l1.6 1.6",
  revolve: "M8 1.5v13M3 9.5a5 2.5 0 0 0 10 0M3 9.5V6M13 9.5V6M3 6a5 2.5 0 0 1 10 0",
  fillet: "M3 13.5V8a5 5 0 0 1 5-5h5.5M3 13.5h10.5",
  chamfer: "M3 13.5V7.5L7.5 3h6M3 13.5h10.5",
  shell: "M2.5 3.5v10h11v-10M5 3.5v7.5h6V3.5",
  boolean: "M9 8a3.5 3.5 0 1 1-7 0 3.5 3.5 0 0 1 7 0zM14 8a3.5 3.5 0 1 1-7 0 3.5 3.5 0 0 1 7 0z",
  linear_pattern: "M2 6h3v4H2zM6.5 6h3v4h-3zM11 6h3v4h-3z",
  circular_pattern: "M8 2.5v1.5M8 12v1.5M2.5 8H4M12 8h1.5M4.1 4.1l1 1M10.9 10.9l1 1M4.1 11.9l1-1M10.9 5.1l1-1M9.2 8a1.2 1.2 0 1 1-2.4 0 1.2 1.2 0 0 1 2.4 0z",
  sketch: "M3 13l.9-3.4L10.6 3l2.4 2.4-6.6 6.7zM9.3 4.3l2.4 2.4",
  reference: "M1.8 11.5 4.8 4.5h9.4l-3 7z",
  hole: "M12.5 8a4.5 4.5 0 1 1-9 0 4.5 4.5 0 0 1 9 0zM8 1.5v2.5M8 12v2.5M1.5 8H4M12 8h2.5",
  sweep: "M2 12.5c3 0 3-9 6-9s3 9 6 9",
  helix: "M3 3.5c5 0 5 3 10 3M3 6.5c5 0 5 3 10 3M3 9.5c5 0 5 3 10 3",
  loft: "M2.5 13.5h11M5 2.5h6M2.5 13.5 5 2.5M13.5 13.5 11 2.5",
  draft: "M4 13.5 5.6 2.5h4.8L12 13.5zM2.5 13.5h11",
  transform: "M2 8h12M11 5l3 3-3 3M5 5 2 8l3 3",
  mirror: "M8 1.5v13M6 4 2 12h4zM10 4l4 8h-4z",
  thicken: "M2 10c3-4 9-4 12 0M2 13c3-4 9-4 12 0",
  fastener: "M5 2.5h6l3 5.5-3 5.5H5L2 8z",
  hinge: "M2.5 3h4.5v10H2.5zM9 3h4.5v10H9zM7 6h2M7 10h2",
  body: "M8 1.8l5.5 3v6.4L8 14.2l-5.5-3V4.8zM8 8l5.5-3.2M8 8v6.2M8 8 2.5 4.8",
  place: "M8 1.8l5.5 3v6.4L8 14.2l-5.5-3V4.8zM8 8l5.5-3.2M8 8v6.2M8 8 2.5 4.8",
  move: "M8 1.5v13M1.5 8h13M6 3.5 8 1.5l2 2M6 12.5l2 2 2-2M3.5 6l-2 2 2 2M12.5 6l2 2-2 2",
  mate: "M6.6 9.4l2.8-2.8M5.2 7.2 3.8 8.6a2.5 2.5 0 0 0 3.6 3.6l1.4-1.4M10.8 8.8l1.4-1.4a2.5 2.5 0 0 0-3.6-3.6L7.2 5.2",
  step: "M10.5 8a2.5 2.5 0 1 1-5 0 2.5 2.5 0 0 1 5 0z",
  eye: "M1.5 8s2.4-4.5 6.5-4.5S14.5 8 14.5 8 12.1 12.5 8 12.5 1.5 8 1.5 8zM10 8a2 2 0 1 1-4 0 2 2 0 0 1 4 0z",
  section: "M2.5 2.5h11v11h-11zM2.5 13.5 13.5 2.5M5.5 13.5l8-8M2.5 10.5l8-8",
  "section-off": "M2.5 2.5h11v11h-11z",
  "eye-off": "M2 2l12 12M6.6 3.7A6.6 6.6 0 0 1 8 3.5c4.1 0 6.5 4.5 6.5 4.5a11 11 0 0 1-1.9 2.4M4.2 5.1A11 11 0 0 0 1.5 8s2.4 4.5 6.5 4.5c1 0 1.9-.3 2.7-.7",
  // Interface marks: drawn, not typed, so they sit centred in their
  // buttons whatever font the device has (a typed ✎ × ▶ › rode the
  // font's own metrics, and ▶ can render as a colour emoji).
  rename: "M10.5 2.5l3 3-8 8h-3v-3z",
  close: "M4 4l8 8M12 4l-8 8",
  up: "M8 13V3M4.5 6.5 8 3l3.5 3.5",
  down: "M8 3v10M4.5 9.5 8 13l3.5-3.5",
  plus: "M8 3v10M3 8h10",
  next: "m6 3.5 4.5 4.5L6 12.5",
};

// Filled marks (drawn with `fill`, no stroke).
const SOLID = {
  play: "M5 3v10l8-5z",
};

// Operation names (the history's `operation_label`, an assembly
// step's kind) to icon names.
const ALIASES = {
  boolean: "boolean",
  pattern: "linear_pattern",
  tapped_hole: "hole",
  clearance_hole: "hole",
  holes: "hole",
  pipe: "sweep",
  helical_sweep: "helix",
  offset_solid: "thicken",
  raw_body: "body",
  hex_bolt: "fastener",
  socket_cap_screw: "fastener",
  countersunk_screw: "fastener",
  threaded_rod: "fastener",
  hex_nut: "fastener",
  living_hinge: "hinge",
};

/** The icon name drawn for history operation `operation`. */
export function iconName(operation) {
  if (!operation) return "step";
  if (PATHS[operation]) return operation;
  if (ALIASES[operation]) return ALIASES[operation];
  if (operation.includes("hinge")) return "hinge";
  return "step";
}

/** An icon's SVG markup (constant text, safe to assign as HTML). */
export function iconSvg(name) {
  if (SOLID[name]) return `<svg class="icon solid" viewBox="0 0 16 16" aria-hidden="true"><path d="${SOLID[name]}"/></svg>`;
  const d = PATHS[name] ?? PATHS.step;
  return `<svg class="icon" viewBox="0 0 16 16" aria-hidden="true"><path d="${d}"/></svg>`;
}

/** A label after an icon: `<icon><span>text</span>` as nodes, the
 *  text set as text (never parsed as HTML). */
export function iconLabel(name, text) {
  const box = document.createElement("template");
  box.innerHTML = iconSvg(name);
  const label = document.createElement("span");
  label.textContent = text;
  return [box.content.firstChild, label];
}

/** The SVG markup for history operation `operation`. */
export function operationIcon(operation) {
  return iconSvg(iconName(operation));
}
