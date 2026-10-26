// The film state rules (film-state.js; docs/architecture/
// presentation-films.md §4, §10), run by crates/popkern-live/tests/
// page_film.rs. Exits non-zero with the failing check's message;
// prints the check count on success.
import assert from "node:assert/strict";
import { CAPTION_FADE, cameraAt, captionAt, designsAt, ease, frameTimes, recordPosition, sectionAt, selects, stateAt } from "../film-state.js";

let checks = 0;
const check = (name, fn) => {
  try {
    fn();
    checks++;
  } catch (e) {
    console.error(`FAIL ${name}: ${e.message}`);
    process.exit(1);
  }
};
const near = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps, `${a} vs ${b}`);

const track = (start, duration, channel, easeKind = "linear") => ({ start, duration, ease: easeKind, line: null, channel });
const film = (tracks, steps = [{ title: "A", subtitle: "", start: 0, duration: 10, kind: "step", line: 1 }]) => ({
  name: "t",
  size: [1920, 1080],
  fps: 30,
  look: "toon",
  background: "#ffffff",
  duration: steps.reduce((d, s) => Math.max(d, s.start + s.duration), 0),
  camera: { azimuth: -35, elevation: 25, zoom: 1, target: "all" },
  steps,
  tracks,
});

check("easing curves hit their ends and shapes", () => {
  for (const k of ["linear", "in", "out", "inout"]) {
    near(ease(k, 0), 0);
    near(ease(k, 1), 1);
  }
  near(ease("inout", 0.5), 0.5);
  assert.ok(ease("in", 0.5) < 0.5 && ease("out", 0.5) > 0.5);
});

check("a move holds its from before, its to after, and eases between", () => {
  const f = film([track(2, 2, { kind: "move", target: "top", from: [0, 0, 1], to: [0, 0, 0] })]);
  near(stateAt(f, 0, ["top"]).items.get("top").offset[2], 1);
  near(stateAt(f, 3, ["top"]).items.get("top").offset[2], 0.5);
  near(stateAt(f, 9, ["top"]).items.get("top").offset[2], 0);
});

check("moves on one target add; a prefix moves what lies under it", () => {
  const f = film([
    track(0, 1, { kind: "move", target: "pair", from: [1, 0, 0], to: [1, 0, 0] }),
    track(0, 1, { kind: "move", target: "pair/high", from: [0, 2, 0], to: [0, 2, 0] }),
  ]);
  const s = stateAt(f, 0.5, ["pair/low", "pair/high", "base"]);
  assert.deepEqual(s.items.get("pair/high").offset, [1, 2, 0]);
  assert.deepEqual(s.items.get("pair/low").offset, [1, 0, 0]);
  assert.deepEqual(s.items.get("base").offset, [0, 0, 0]);
  assert.ok(selects("all", "x") && selects("a", "a/b") && !selects("a", "ab"));
});

check("fades hold the latest started value; shows flip at their start", () => {
  const f = film([
    track(1, 1, { kind: "fade", target: "a", from: 1, to: 0.2 }),
    track(5, 1, { kind: "fade", target: "a", from: 0.2, to: 1 }),
    track(3, 0, { kind: "show", target: "b", visible: true }),
  ]);
  near(stateAt(f, 0, ["a"]).items.get("a").opacity, 1);
  near(stateAt(f, 4, ["a"]).items.get("a").opacity, 0.2);
  near(stateAt(f, 5.5, ["a"]).items.get("a").opacity, 0.6);
  assert.equal(stateAt(f, 2, ["b"]).items.get("b").visible, false, "hidden before it is shown");
  assert.equal(stateAt(f, 3, ["b"]).items.get("b").visible, true);
});

check("joints, explode and highlight", () => {
  const f = film([
    track(0, 2, { kind: "joint", param: "j1", from: 0, to: 90 }),
    track(4, 2, { kind: "joint", param: "j1", from: 90, to: 45 }),
    track(1, 2, { kind: "explode", from: 0, to: 1 }),
    track(0, 4, { kind: "highlight", target: "a" }),
  ]);
  near(stateAt(f, 1, ["a"]).joints.j1, 45);
  near(stateAt(f, 3, ["a"]).joints.j1, 90);
  near(stateAt(f, 5, ["a"]).joints.j1, 67.5);
  near(stateAt(f, 2, ["a"]).explode, 0.5);
  near(stateAt(f, 2, ["a"]).items.get("a").highlight, 1);
  near(stateAt(f, 0, ["a"]).items.get("a").highlight, 0);
  near(stateAt(f, 5, ["a"]).items.get("a").highlight, 0);
});

check("camera: orbits add, a camera track eases to its pose absolutely, zoom geometric, targets blend", () => {
  const f = film([
    track(0, 2, { kind: "orbit", degrees: 40 }),
    track(2, 2, { kind: "camera", to: { azimuth: 90, elevation: null, zoom: 4, target: "top" } }, "linear"),
    track(4, 2, { kind: "orbit", degrees: 10 }),
  ]);
  near(cameraAt(f, 1).azimuth, -15);
  const mid = cameraAt(f, 3);
  near(mid.azimuth, (5 + 90) / 2);
  near(mid.elevation, 25);
  near(mid.zoom, 2, 1e-12);
  near(mid.targets.all, 0.5);
  near(mid.targets.top, 0.5);
  const end = cameraAt(f, 9);
  near(end.azimuth, 100);
  assert.deepEqual(Object.keys(end.targets), ["top"]);
});

check("captions fade in and out, carry over identical steps, and holds have none", () => {
  const steps = [
    { title: "One", subtitle: "a", start: 0, duration: 2, kind: "title", line: 1 },
    { title: "One", subtitle: "a", start: 2, duration: 2, kind: "title", line: 2 },
    { title: "", subtitle: "", start: 4, duration: 1, kind: "hold", line: 3 },
    { title: "Two", subtitle: "", start: 5, duration: 2, kind: "step", line: 4 },
  ];
  const f = film([], steps);
  near(captionAt(f, 0).opacity, 0);
  near(captionAt(f, CAPTION_FADE / 2).opacity, 0.5);
  near(captionAt(f, 1.99).opacity, 1, 1e-9);
  near(captionAt(f, 2.01).opacity, 1, 1e-9);
  assert.equal(captionAt(f, 4.5), null);
  assert.equal(captionAt(f, 6.9).title, "Two");
  near(captionAt(f, 7).opacity, 1, 1e-9);
});

check("section keys blend: an axis preset across bounds, off cuts nothing", () => {
  const bounds = [[0, 0, 0], [2, 4, 6]];
  const p = sectionAt({ from: { kind: "axis", axis: "y", at: 1 }, to: { kind: "axis", axis: "y", at: 0.5 }, k: 0.5, combine: "any" }, bounds);
  assert.deepEqual(p.normal, [0, 1, 0]);
  near(p.offset, 3);
  const off = sectionAt({ from: { kind: "off" }, to: { kind: "axis", axis: "z", at: 0.5 }, k: 0, combine: "any" }, bounds);
  assert.ok(off.offset > 6, "off sits past the bounds");
  const tilt = sectionAt({ from: { kind: "plane", normal: [1, 0, 0], offset: 1 }, to: { kind: "plane", normal: [0, 1, 0], offset: 1 }, k: 0.5, combine: "any" }, bounds);
  near(Math.hypot(...tilt.normal), 1);
  // Halfway: 0.5x + 0.5y = 1, i.e. (x + y)/sqrt2 = sqrt2.
  near(tilt.offset, Math.SQRT2, 1e-12);
  assert.equal(sectionAt(null, bounds), null);
});

check("frame times cover the film at its fps", () => {
  const f = film([], [{ title: "A", subtitle: "", start: 0, duration: 2.5, kind: "step", line: 1 }]);
  const times = frameTimes(f);
  assert.equal(times.length, 75);
  near(times[74], 74 / 30);
  assert.equal(frameTimes(f, 15).length, 38);
});

// A design converging (PK-959): records at updates 1, 2, 4, 10.
const records = [1, 2, 4, 10].map((update) => ({ update }));
check("an update's position among the records is piecewise linear and clamped", () => {
  near(recordPosition(records, 0), 0);
  near(recordPosition(records, 1), 0);
  near(recordPosition(records, 3), 1.5);
  near(recordPosition(records, 7), 2.5);
  near(recordPosition(records, 10), 3);
  near(recordPosition(records, 50), 3);
  near(recordPosition([{ update: 5 }], 9), 0);
  near(recordPosition([], 9), 0);
});

check("evolve spreads its time over the records, holding before and after", () => {
  const f = film([
    track(2, 4, { kind: "evolve", target: "design carve", from: 0, to: 10 }),
    track(8, 1, { kind: "evolve", target: "design carve", from: 10, to: 4 }),
  ]);
  const at = (t) => designsAt(f, t, (selector) => (selector === "design carve" ? records : null)).get("design carve");
  near(at(0), 0);
  near(at(4), 1.5);
  near(at(6), 3);
  near(at(7), 3);
  near(at(8.5), 2.5);
  near(at(10), 2);
  // A design whose records are not loaded has no position.
  assert.equal(designsAt(f, 4, () => null).size, 0);
  // Its layer is a selector like any body's: fades reach it.
  const state = stateAt(film([track(0, 2, { kind: "fade", target: "design carve", from: 1, to: 0 })]), 1, ["design carve", "body 0"]);
  near(state.items.get("design carve").opacity, 0.5);
  near(state.items.get("body 0").opacity, 1);
});

console.log(JSON.stringify({ checks }));
