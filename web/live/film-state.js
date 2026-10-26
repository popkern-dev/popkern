// A film's state at a time (docs/architecture/presentation-films.md §4,
// §6): pure — no three.js, no DOM — so the node test
// (test/film.test.mjs) checks the hold, ease, compose and caption rules
// the player and the export both use.
//
// The film is the worker's `film` reply's `film`: steps and tracks at
// absolute times, lengths SI. Before a track starts it holds its first
// value; after it ends, its last.

/** Ease a 0..1 progress. */
export function ease(kind, p) {
  const t = Math.min(1, Math.max(0, p));
  switch (kind) {
    case "linear":
      return t;
    case "in":
      return t * t;
    case "out":
      return 1 - (1 - t) * (1 - t);
    default:
      return t * t * (3 - 2 * t);
  }
}

/** A track's eased progress at `t` (0 before its start, 1 after its end). */
export function progress(track, t) {
  if (t < track.start) return 0;
  if (track.duration <= 0 || t >= track.start + track.duration) return 1;
  return ease(track.ease, (t - track.start) / track.duration);
}

const lerp = (a, b, k) => a + (b - a) * k;
const kind = (track) => track.channel.kind;

/** Whether a track's target selects occurrence `path` (a prefix
 *  selects what lies under it; "all" everything). */
export function selects(target, path) {
  return target === "all" || path === target || (path != null && path.startsWith(`${target}/`));
}

/** The value of an absolute channel (fade, explode, joint, section) for
 *  the tracks `list`: the latest track that has started, at its
 *  progress; before any has, the earliest one's `from`. Returns
 *  `{ track, k }` (k: that track's progress) or null. */
function latest(list, t) {
  let best = null;
  for (const track of list) if (track.start <= t && (!best || track.start >= best.start)) best = track;
  if (best) return { track: best, k: progress(best, t) };
  let first = null;
  for (const track of list) if (!first || track.start < first.start) first = track;
  return first ? { track: first, k: 0 } : null;
}

/** The camera at `t`: `{ azimuth, elevation, zoom, targets }` where
 *  `targets` maps a target name to its weight (a camera move between
 *  two targets blends them). Orbits add to the azimuth; a camera track
 *  eases from wherever the camera is to its pose, absolutely. */
export function cameraAt(film, t) {
  const pose = {
    azimuth: film.camera.azimuth,
    elevation: film.camera.elevation,
    zoom: film.camera.zoom,
    targets: { [film.camera.target]: 1 },
  };
  const tracks = film.tracks
    .filter((tr) => kind(tr) === "orbit" || kind(tr) === "camera")
    .sort((a, b) => a.start - b.start);
  for (const track of tracks) {
    if (track.start > t) break;
    const k = progress(track, t);
    if (kind(track) === "orbit") {
      pose.azimuth += track.channel.degrees * k;
      continue;
    }
    const to = track.channel.to;
    if (to.azimuth != null) pose.azimuth = lerp(pose.azimuth, to.azimuth, k);
    if (to.elevation != null) pose.elevation = lerp(pose.elevation, to.elevation, k);
    if (to.zoom != null) pose.zoom = Math.exp(lerp(Math.log(pose.zoom), Math.log(to.zoom), k));
    if (to.target) {
      const next = {};
      for (const [name, w] of Object.entries(pose.targets)) next[name] = (next[name] ?? 0) + w * (1 - k);
      next[to.target] = (next[to.target] ?? 0) + k;
      for (const name of Object.keys(next)) if (next[name] < 1e-9) delete next[name];
      pose.targets = next;
    }
  }
  return pose;
}

/** Everything a frame at `t` needs, except the camera (`cameraAt`).
 *  `paths`: the occurrence paths (or "body N" names) the scene draws.
 *  Returns per-path `{ offset: [x,y,z], turns: [{origin, direction,
 *  degrees}], opacity, visible, highlight }`, plus `joints` (name →
 *  value), `explode`, `section` (`{ from, to, k, combine }` or null)
 *  and `caption`. */
export function stateAt(film, t, paths) {
  const items = new Map();
  for (const path of paths) items.set(path, { offset: [0, 0, 0], turns: [], opacity: 1, visible: true, highlight: 0 });
  const byKind = (name) => film.tracks.filter((tr) => kind(tr) === name);

  for (const track of byKind("move")) {
    const k = progress(track, t);
    const { target, from, to } = track.channel;
    const d = [lerp(from[0], to[0], k), lerp(from[1], to[1], k), lerp(from[2], to[2], k)];
    for (const [path, item] of items) {
      if (!selects(target, path)) continue;
      item.offset[0] += d[0];
      item.offset[1] += d[1];
      item.offset[2] += d[2];
    }
  }
  for (const track of byKind("turn")) {
    const k = progress(track, t);
    const { target, origin, direction, from, to } = track.channel;
    const degrees = lerp(from, to, k);
    for (const [path, item] of items) if (selects(target, path)) item.turns.push({ origin, direction, degrees });
  }
  // Fades: per path, the latest started fade among those selecting it.
  const fades = byKind("fade");
  const shows = byKind("show").sort((a, b) => a.start - b.start);
  for (const [path, item] of items) {
    const mine = fades.filter((tr) => selects(tr.channel.target, path));
    const at = latest(mine, t);
    if (at) item.opacity = lerp(at.track.channel.from, at.track.channel.to, at.k);
    const myShows = shows.filter((tr) => selects(tr.channel.target, path));
    if (myShows.length) {
      // Before the first show, the opposite of what it sets.
      let visible = !myShows[0].channel.visible;
      for (const tr of myShows) if (tr.start <= t) visible = tr.channel.visible;
      item.visible = visible;
    }
    for (const tr of byKind("highlight")) {
      if (!selects(tr.channel.target, path) || t < tr.start || t > tr.start + tr.duration) continue;
      const u = tr.duration > 0 ? (t - tr.start) / tr.duration : 1;
      item.highlight = Math.max(item.highlight, ease("inout", Math.min(1, u / 0.25, (1 - u) / 0.25)));
    }
  }
  const joints = {};
  const jointTracks = byKind("joint");
  for (const name of new Set(jointTracks.map((tr) => tr.channel.param))) {
    const at = latest(jointTracks.filter((tr) => tr.channel.param === name), t);
    joints[name] = lerp(at.track.channel.from, at.track.channel.to, at.k);
  }
  const explodeAt = latest(byKind("explode"), t);
  const explode = explodeAt ? lerp(explodeAt.track.channel.from, explodeAt.track.channel.to, explodeAt.k) : 0;
  const sectionAt = latest(byKind("section"), t);
  const section = sectionAt
    ? { from: sectionAt.track.channel.from, to: sectionAt.track.channel.to, k: sectionAt.k, combine: sectionAt.track.channel.combine }
    : null;
  return { items, joints, explode, section, caption: captionAt(film, t) };
}

/** The position among a design's records (`frames`: `[{ update }]`,
 *  ascending) of update `update`: piecewise linear between the records
 *  either side, clamped to the first and the last (PK-959). */
export function recordPosition(frames, update) {
  const last = frames.length - 1;
  if (last <= 0 || !(update > frames[0].update)) return 0;
  if (update >= frames[last].update) return last;
  let i = 0;
  while (frames[i + 1].update < update) i += 1;
  const a = frames[i].update;
  const b = frames[i + 1].update;
  return i + (update - a) / (b - a);
}

/** Each evolving design's position among its records at `t`: a Map of
 *  selector ("design carve") to position. `framesOf(selector)` gives a
 *  design's records. As with fades, the latest evolve track that has
 *  started holds; before any has, the earliest one's start. The time
 *  spreads evenly over the records, which are denser early, where the
 *  design changes most. */
export function designsAt(film, t, framesOf) {
  const out = new Map();
  const tracks = film.tracks.filter((tr) => kind(tr) === "evolve");
  for (const target of new Set(tracks.map((tr) => tr.channel.target))) {
    const frames = framesOf(target);
    if (!frames?.length) continue;
    const at = latest(tracks.filter((tr) => tr.channel.target === target), t);
    const { from, to } = at.track.channel;
    out.set(target, lerp(recordPosition(frames, from), recordPosition(frames, to), at.k));
  }
  return out;
}

/** Seconds a caption takes to fade in and out. */
export const CAPTION_FADE = 0.35;

/** The caption at `t`: `{ title, subtitle, kind, opacity, index }` or
 *  null. A caption fades in at its step's start and out at its end,
 *  except where the next step carries the same caption on. */
export function captionAt(film, t) {
  const steps = film.steps;
  const index = steps.findIndex((s) => t >= s.start && t < s.start + s.duration);
  const i = index < 0 ? (t >= film.duration && steps.length ? steps.length - 1 : -1) : index;
  if (i < 0) return null;
  const step = steps[i];
  if (step.kind === "hold" || !step.title) return null;
  const same = (a, b) => a && b && a.title === b.title && a.subtitle === b.subtitle && b.kind !== "hold";
  const fadeIn = same(steps[i - 1], step) ? 1 : Math.min(1, (t - step.start) / CAPTION_FADE);
  const end = step.start + step.duration;
  const last = i === steps.length - 1;
  const fadeOut = same(step, steps[i + 1]) || last ? 1 : Math.min(1, (end - t) / CAPTION_FADE);
  return { title: step.title, subtitle: step.subtitle, kind: step.kind, index: i, opacity: Math.max(0, Math.min(fadeIn, fadeOut)) };
}

/** A section key as a plane `{ normal, offset }` in world space, the
 *  axis presets resolved across `bounds` (removing the high side). */
export function sectionPlane(key, bounds) {
  if (!key || key.kind === "off") return null;
  if (key.kind === "plane") return { normal: key.normal, offset: key.offset };
  const i = "xyz".indexOf(key.axis);
  const normal = [0, 0, 0];
  normal[i] = 1;
  return { normal, offset: bounds[0][i] + (bounds[1][i] - bounds[0][i]) * key.at };
}

/** A section state's plane at its progress: the two keys' planes
 *  blended (normals renormalised), or null when both are off. A key
 *  that is off blends from or to a plane past the bounds (nothing cut). */
export function sectionAt(section, bounds) {
  if (!section) return null;
  let a = sectionPlane(section.from, bounds);
  let b = sectionPlane(section.to, bounds);
  if (!a && !b) return null;
  const beyond = (plane) => {
    // The same plane moved past the bounds' far corner: removes nothing.
    let max = -Infinity;
    for (let c = 0; c < 8; c++) {
      const p = [c & 1 ? bounds[1][0] : bounds[0][0], c & 2 ? bounds[1][1] : bounds[0][1], c & 4 ? bounds[1][2] : bounds[0][2]];
      max = Math.max(max, plane.normal[0] * p[0] + plane.normal[1] * p[1] + plane.normal[2] * p[2]);
    }
    return { normal: plane.normal, offset: max + 1e-6 * (1 + Math.abs(max)) };
  };
  a ??= beyond(b);
  b ??= beyond(a);
  const k = section.k;
  const n = a.normal.map((c, i) => lerp(c, b.normal[i], k));
  const len = Math.hypot(...n) || 1;
  return { normal: n.map((c) => c / len), offset: lerp(a.offset, b.offset, k) / (len || 1) };
}

/** The frame times of an export: `count` frames at `fps`. */
export function frameTimes(film, fps = film.fps) {
  const count = Math.max(1, Math.round(film.duration * fps));
  return Array.from({ length: count }, (_, k) => k / fps);
}
