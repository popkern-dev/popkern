// Golden film states from web/live/film-state.js for the Rust port's
// parity test (crates/popkern-cli/tests/film_state_parity.rs;
// docs/architecture/cli-render-export.md §7).
//
// Usage: node tools/film-state-fixtures.mjs <fixture.json>...
// Each fixture holds `{ film, paths, bounds, times }` (written by the
// test under POPKERN_FILM_FIXTURES=1); this adds `states`, one per time:
// the camera, the items that differ from an untouched one, joints,
// explode, section, its plane and the caption, with Maps and objects as
// ordered [key, value] pairs.
import { readFileSync, writeFileSync } from "node:fs";
import { cameraAt, sectionAt, stateAt } from "../web/live/film-state.js";

/** An item a film has touched (an untouched one is left out). */
const moved = (item) =>
  item.offset.some((c) => c !== 0) || item.turns.length > 0 || item.opacity !== 1 || !item.visible || item.highlight !== 0;

for (const file of process.argv.slice(2)) {
  const fixture = JSON.parse(readFileSync(file, "utf8"));
  const { film, paths, bounds, times } = fixture;
  fixture.states = times.map((t) => {
    const camera = cameraAt(film, t);
    const state = stateAt(film, t, paths);
    return {
      t,
      camera: { azimuth: camera.azimuth, elevation: camera.elevation, zoom: camera.zoom, targets: Object.entries(camera.targets) },
      items: [...state.items.entries()].filter(([, item]) => moved(item)),
      joints: Object.entries(state.joints),
      explode: state.explode,
      section: state.section,
      plane: state.section ? sectionAt(state.section, bounds) : null,
      caption: state.caption,
    };
  });
  writeFileSync(file, `${JSON.stringify(fixture)}\n`);
  console.log(`${file}: ${times.length} states`);
}
