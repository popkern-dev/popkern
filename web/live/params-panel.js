// The parameters panel: a slider and a number box per editable model
// parameter. Read-only and expression-driven parameters are not
// offered; an opened assembly part shows only its own ("part:name").

const $ = (id) => document.getElementById(id);

/** Render `params` (an envelope's parameter list) scoped to `part`
 *  (null: every part); a committed edit calls `setParam(name, value)`. */
export function renderParams(params, part, setParam) {
  const list = $("params");
  list.innerHTML = "";
  for (const p of params) {
    if (p.read_only || p.expression) continue;
    // An assembly's parameters are qualified "part:name"; an opened
    // part shows only its own.
    if (part != null && !p.name.startsWith(`${part}:`)) continue;
    const row = document.createElement("label");
    row.className = "param";
    const bounds = p.bounds ?? [p.value > 0 ? p.value / 4 : p.value * 2 - 1, p.value > 0 ? p.value * 3 : p.value / 2 + 1];
    const lo = Math.min(bounds[0], p.value);
    const hi = Math.max(bounds[1], p.value);
    row.innerHTML = `<span class="param-name">${p.name}</span>
        <input type="range" min="${lo}" max="${hi}" step="${(hi - lo) / 200}" value="${p.value}" aria-label="${p.name} slider">
        <input type="number" value="${p.value}" step="any" aria-label="${p.name} value">`;
    const [slider, number] = row.querySelectorAll("input");
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
  $("params-empty").hidden = list.children.length > 0;
}
