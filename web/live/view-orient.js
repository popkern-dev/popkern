// The standard views (the pad in the model view's top-right corner):
// isometric, top, bottom, front, back, left and right, each refitted
// to the visible model. Z is up and the front looks along +Y, so the
// front view has +X to the right and +Z up; top has +X right, +Y up.
// View state only: nothing here builds or changes the document.

const $ = (id) => document.getElementById(id);

// Polar angles a hair off the poles: looking straight down Z leaves
// the camera's up (+Z) parallel to the view and lookAt undefined.
const POLE = 1e-4;

/** Azimuth `theta` (from +X towards +Y) and polar angle `phi` (from
 *  +Z) of the eye for each standard view. */
export const VIEWS = {
  iso: [-Math.PI / 4, Math.acos(1 / Math.sqrt(3))],
  front: [-Math.PI / 2, Math.PI / 2],
  back: [Math.PI / 2, Math.PI / 2],
  right: [0, Math.PI / 2],
  left: [Math.PI, Math.PI / 2],
  top: [-Math.PI / 2, POLE],
  bottom: [-Math.PI / 2, Math.PI - POLE],
};

/** Turn the viewer to standard view `name`. */
export function orientView(viewer, name, ms) {
  const view = VIEWS[name];
  const box = viewer.visibleBounds();
  if (!view || !box) return false;
  viewer.orbit.orient(view[0], view[1], box, ms);
  viewer.fitted = true;
  return true;
}

/** Wire the standard-view pad once, at page start. */
export function wireViewOrient(page) {
  for (const button of $("view-orient").querySelectorAll("[data-view]")) {
    button.addEventListener("click", () => orientView(page.viewer, button.dataset.view));
  }
}
