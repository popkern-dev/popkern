// The parameters panel: a slider and a number box per editable model
// parameter. Read-only and expression-driven parameters are not
// offered; an opened assembly part shows only its own ("part:name").
// Two toggles filter by kind: parameters the script declares by name
// (`param("wall", …)`) and feature arguments the history exposes
// implicitly (`fillet-3.radius`); the choice is kept in this browser.
import { readParamFilter, writeParamFilter } from "./local-store.js";

const $ = (id) => document.getElementById(id);

// The last render's arguments, so a toggle re-renders without a build.
let last = null;
// The last render's slider and number box per parameter name.
const controls = new Map();

/** Whether `p` was declared by name. Envelopes from before the
 *  `explicit` field (older bakes) fall back to the naming rule: a
 *  feature's parameters are qualified `feature.role`. */
export function isNamed(p) {
  if (typeof p.explicit === "boolean") return p.explicit;
  return !p.name.split(":").pop().includes(".");
}

/** Wire the two kind toggles (once, at page start). */
export function wireParamFilter() {
  const filter = readParamFilter();
  for (const kind of ["named", "implicit"]) {
    const box = $(`params-show-${kind}`);
    box.checked = filter[kind];
    box.addEventListener("change", () => {
      writeParamFilter({ named: $("params-show-named").checked, implicit: $("params-show-implicit").checked });
      if (last) renderParams(...last);
    });
  }
}

/** Render `params` (an envelope's parameter list) scoped to `part`
 *  (null: every part); a committed edit calls `setParam(name, value)`. */
export function renderParams(params, part, setParam) {
  last = [params, part, setParam];
  const show = { named: $("params-show-named").checked, implicit: $("params-show-implicit").checked };
  const list = $("params");
  list.innerHTML = "";
  controls.clear();
  const hidden = { named: 0, implicit: 0 };
  for (const p of params) {
    if (p.read_only || p.expression) continue;
    // An assembly's parameters are qualified "part:name"; an opened
    // part shows only its own.
    if (part != null && !p.name.startsWith(`${part}:`)) continue;
    const kind = isNamed(p) ? "named" : "implicit";
    if (!show[kind]) {
      hidden[kind] += 1;
      continue;
    }
    const row = document.createElement("label");
    row.className = `param param-${kind}`;
    const bounds = p.bounds ?? [p.value > 0 ? p.value / 4 : p.value * 2 - 1, p.value > 0 ? p.value * 3 : p.value / 2 + 1];
    const lo = Math.min(bounds[0], p.value);
    const hi = Math.max(bounds[1], p.value);
    row.innerHTML = `<span class="param-name"></span>
        <input type="range" min="${lo}" max="${hi}" step="${(hi - lo) / 200}" value="${p.value}">
        <input type="number" value="${p.value}" step="any">`;
    const name = row.querySelector(".param-name");
    name.textContent = p.name;
    name.title = kind === "named" ? `${p.name}: declared by name in the script` : `${p.name}: a feature argument, exposed implicitly`;
    const [slider, number] = row.querySelectorAll("input");
    controls.set(p.name, { slider, number });
    slider.setAttribute("aria-label", `${p.name} slider`);
    number.setAttribute("aria-label", `${p.name} value`);
    const apply = (value) => {
      if (!Number.isFinite(value)) return;
      number.value = value;
      slider.value = value;
      setParam(p.name, value);
    };
    slider.addEventListener("change", () => apply(Number(slider.value)));
    number.addEventListener("change", () => apply(Number(number.value)));
    list.append(row);
  }
  const off = hidden.named + hidden.implicit;
  $("params-empty").hidden = list.children.length > 0 || off > 0;
  const note = $("params-hidden");
  note.hidden = off === 0;
  const kinds = [hidden.named && `${hidden.named} named`, hidden.implicit && `${hidden.implicit} implicit`].filter(Boolean);
  note.textContent = `${off} hidden (${kinds.join(", ")})`;
}

/** Set the values of the rendered rows named in `changed` in place —
 *  a joint slider (a placement-only change) re-renders nothing, so a
 *  dragged slider keeps its focus and ten thousand part rows are not
 *  rebuilt per step. `params` is the full list now current. */
export function refreshParamValues(changed, params) {
  if (last) last = [params, last[1], last[2]];
  for (const p of changed) {
    const control = controls.get(p.name);
    if (!control) continue;
    const { slider, number } = control;
    if (p.value < Number(slider.min)) slider.min = String(p.value);
    if (p.value > Number(slider.max)) slider.max = String(p.value);
    slider.value = String(p.value);
    number.value = String(p.value);
  }
}
