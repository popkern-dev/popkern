// The status line: one line saying what the page is doing or what the
// last build served, and the viewer's pick readout (a face's or edge's
// stable name). Uncaught page errors land on the status line too.

const $ = (id) => document.getElementById(id);

/** Show `text` on the status line, styled as a failure when `failed`. */
export function setStatus(text, failed = false) {
  const el = $("status");
  el.textContent = text;
  // One line in the dock: the whole text is the tooltip.
  el.title = text;
  el.classList.toggle("status-failed", failed);
  document.body.classList.toggle("status-failed", failed);
}

/** "3 body(ies)" for a part document, "4 part(s), 9 occurrence(s)" for an assembly. */
export function modelContents(envelope) {
  const view = envelope.assembly;
  if (!view) return `${envelope.bodies.length} body(ies)`;
  const parts = view.definitions.filter((d) => d.kind === "part").length;
  return `${parts} part(s), ${view.instances.filter((i) => i.is_part).length} occurrence(s)`;
}

/** ", watertight" or ", N open mesh edges" for a served envelope. */
export function meshClosure(envelope) {
  const open = envelope.bodies.reduce((n, b) => n + b.open_edges, 0);
  return open ? `, ${open} open mesh edges` : ", watertight";
}

/** The pick readout's text for a face or edge pick (null: the hint). */
export function pickText(pick) {
  if (!pick) return "click a face or edge to see its stable name";
  if (pick.instance) {
    // An assembly pick is an entity path: the occurrence, then the
    // part's own stable name (design §4): "planet_1:planet/…".
    const at = (label) => `${pick.instance}:${label}`;
    if (pick.edge) {
      const [la, lb] = pick.edge.labels;
      return pick.edge.faces[1] === null ? `open edge of ${at(la)}` : `edge of ${at(la)} | ${at(lb)}`;
    }
    return `face ${at(pick.label)}`;
  }
  const body = pick.body.label ? `${pick.body.label} · ` : "";
  if (pick.edge) {
    const [a, b] = pick.edge.faces;
    const [la, lb] = pick.edge.labels;
    return b === null ? `${body}open edge of face ${a}: ${la}` : `${body}edge of faces ${a} | ${b}: ${la} | ${lb}`;
  }
  return `${body}face ${pick.face}: ${pick.label}`;
}

/** Show a viewer pick in the viewport's corner. */
export function showPick(pick) {
  const el = $("pick");
  el.textContent = pickText(pick);
  // One line (the full name is the tooltip), so a long pick never
  // reflows the viewer under the pointer.
  el.title = el.textContent;
}

/** Report uncaught errors and rejections on the status line, so a page
 *  fault is never silent. */
export function reportPageErrors() {
  window.addEventListener("error", (event) => {
    const status = $("status");
    if (status) status.textContent = `page error: ${event.message}`;
  });
  window.addEventListener("unhandledrejection", (event) => {
    const status = $("status");
    if (status) status.textContent = `page error: ${event.reason?.message ?? event.reason}`;
  });
}
