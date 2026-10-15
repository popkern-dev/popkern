// The view menu's section row (docs/architecture/clipping-view.md §10):
// a preset select for the current plane, its position slider, plane
// chips (up to three), flip / add / remove, the combine mode, and the
// in-viewport handles toggle. The controls follow the section's state
// through `section.onChange`, so a link, a film or a gizmo drag shows.
import { MAX_PLANES } from "./section-geometry.js";

const $ = (id) => document.getElementById(id);

/** The preset an orientation matches: "x" | "y" | "z" (either sense),
 *  else "custom". */
function presetOf(plane) {
  const [x, y, z] = plane.normal.map(Math.abs);
  if (x > 1 - 1e-9) return "x";
  if (y > 1 - 1e-9) return "y";
  if (z > 1 - 1e-9) return "z";
  return "custom";
}

export function wireSection(page) {
  const section = page.viewer.section;
  const select = $("section-axis");
  const slider = $("section-at");

  const pickFace = () => {
    $("section-hint").hidden = false;
    $("section-tools").hidden = false;
    page.viewer.pickHook = (hit) => {
      $("section-hint").hidden = true;
      if (!hit.normal) return;
      // On the face, then a quarter of the model's depth along it into
      // the material, so the first click already shows a cut.
      const bounds = page.viewer.visibleBounds();
      const n = hit.normal;
      const depth = bounds ? Math.abs((bounds[1][0] - bounds[0][0]) * n[0]) + Math.abs((bounds[1][1] - bounds[0][1]) * n[1]) + Math.abs((bounds[1][2] - bounds[0][2]) * n[2]) : 0;
      const point = hit.point.map((c, i) => c - 0.25 * depth * n[i]);
      section.addPlane({ point, normal: n }, { replace: true });
    };
  };

  // The viewport's section button: a cut facing the viewer with its
  // handles, or none.
  $("section-toggle").addEventListener("click", () => {
    if (section.active) {
      section.clear();
      return;
    }
    section.gizmo.setEnabled(true);
    $("section-handles").checked = true;
    section.addPlane("view");
  });

  select.addEventListener("change", () => {
    const value = select.value;
    page.viewer.pickHook = null;
    $("section-hint").hidden = true;
    if (!value) return section.clear();
    if (value === "face") return pickFace();
    if (value === "custom") return;
    if (!section.active && "xyz".includes(value)) {
      // The axis preset at the slider's position (the original row).
      section.set(value, Number(slider.value));
      return;
    }
    section.addPlane(value, { replace: true });
  });
  slider.addEventListener("input", () => section.move(Number(slider.value)));
  $("section-flip").addEventListener("click", () => section.flip());
  $("section-add").addEventListener("click", () => section.addPlane("view"));
  $("section-remove").addEventListener("click", () => section.removePlane());
  $("section-combine").addEventListener("change", (e) => section.setCombine(e.target.value));
  $("section-handles").addEventListener("change", (e) => section.gizmo.setEnabled(e.target.checked));

  const paint = () => {
    const planes = section.planes;
    const active = planes.length > 0;
    $("section-toggle").setAttribute("aria-pressed", String(active));
    slider.hidden = !active;
    $("section-tools").hidden = !active && $("section-hint").hidden;
    if (!active) {
      if (select.value !== "face") select.value = "";
      page.assembly?.paintSection?.(false);
      return;
    }
    const preset = presetOf(planes[section.current]);
    $("section-axis").querySelector('option[value="custom"]').hidden = preset !== "custom";
    if (select.value !== "face" || $("section-hint").hidden) select.value = preset;
    if (document.activeElement !== slider) slider.value = String(section.sliderValue());
    const chips = $("section-planes");
    chips.innerHTML = "";
    if (planes.length > 1) {
      planes.forEach((_, i) => {
        const chip = document.createElement("button");
        chip.type = "button";
        chip.textContent = `plane ${i + 1}`;
        chip.setAttribute("aria-pressed", String(i === section.current));
        chip.addEventListener("click", () => {
          section.current = i;
          section.apply();
        });
        chips.append(chip);
      });
    }
    $("section-add").disabled = planes.length >= MAX_PLANES;
    $("section-combine-row").hidden = planes.length < 2;
    $("section-combine").value = section.combine;
    page.assembly?.paintSection?.(true);
  };
  section.onChange(paint);
  paint();
}
