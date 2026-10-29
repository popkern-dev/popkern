//! A film's state at a time: the Rust twin of `web/live/film-state.js`
//! (docs/architecture/presentation-films.md §4, §6;
//! cli-render-export.md §5), so a native renderer plays a film exactly
//! as the page does. Pure: no scene, no clock. Before a track starts it
//! holds its first value; after it ends, its last. Every rule below is
//! the JavaScript's, in the same order (a stable sort by start, later
//! ties winning, first appearance ordering joints and targets), and the
//! CLI's parity test checks them against the JavaScript on golden
//! fixtures.

use popkern_core::geometry::PortableMath;
use serde::Serialize;

use super::model::{Channel, Ease, Film, SectionKey, Track};

/// Seconds a caption takes to fade in and out (`CAPTION_FADE`).
pub const CAPTION_FADE: f64 = 0.35;

/// Ease a 0..1 progress.
pub fn ease(kind: Ease, p: f64) -> f64 {
    let t = p.clamp(0.0, 1.0);
    match kind {
        Ease::Linear => t,
        Ease::In => t * t,
        Ease::Out => 1.0 - (1.0 - t) * (1.0 - t),
        Ease::InOut => t * t * (3.0 - 2.0 * t),
    }
}

/// A track's eased progress at `t` (0 before its start, 1 after its end).
pub fn progress(track: &Track, t: f64) -> f64 {
    if t < track.start {
        return 0.0;
    }
    if track.duration <= 0.0 || t >= track.start + track.duration {
        return 1.0;
    }
    ease(track.ease, (t - track.start) / track.duration)
}

fn lerp(a: f64, b: f64, k: f64) -> f64 {
    a + (b - a) * k
}

/// Whether `target` selects `path` (a prefix selects what lies under
/// it; `"all"` everything).
pub fn selects(target: &str, path: &str) -> bool {
    target == "all"
        || path == target
        || path
            .strip_prefix(target)
            .is_some_and(|rest| rest.starts_with('/'))
}

/// The latest track of `list` that has started, at its progress; before
/// any has, the earliest one at progress 0.
fn latest<'a>(list: &[&'a Track], t: f64) -> Option<(&'a Track, f64)> {
    let mut best: Option<&Track> = None;
    for &track in list {
        if track.start <= t && best.is_none_or(|b| track.start >= b.start) {
            best = Some(track);
        }
    }
    if let Some(best) = best {
        return Some((best, progress(best, t)));
    }
    let mut first: Option<&Track> = None;
    for &track in list {
        if first.is_none_or(|f| track.start < f.start) {
            first = Some(track);
        }
    }
    first.map(|f| (f, 0.0))
}

/// The camera at a time: orbit angles, zoom, and the targets it frames
/// with their weights (a camera move between two targets blends them).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CameraPose {
    /// Degrees about +Z from +X.
    pub azimuth: f64,
    /// Degrees above the XY plane.
    pub elevation: f64,
    /// 1 fits the target.
    pub zoom: f64,
    /// `(target, weight)` in first-appearance order.
    pub targets: Vec<(String, f64)>,
}

/// The camera at `t`: orbits add to the azimuth; a camera track eases
/// from wherever the camera is to its pose (zoom geometrically).
pub fn camera_at(film: &Film, t: f64) -> CameraPose {
    let mut pose = CameraPose {
        azimuth: film.camera.azimuth,
        elevation: film.camera.elevation,
        zoom: film.camera.zoom,
        targets: vec![(film.camera.target.clone(), 1.0)],
    };
    let mut tracks: Vec<&Track> = film
        .tracks
        .iter()
        .filter(|tr| matches!(tr.channel, Channel::Orbit { .. } | Channel::Camera { .. }))
        .collect();
    tracks.sort_by(|a, b| a.start.total_cmp(&b.start));
    for track in tracks {
        if track.start > t {
            break;
        }
        let k = progress(track, t);
        match &track.channel {
            Channel::Orbit { degrees } => pose.azimuth += degrees * k,
            Channel::Camera { to } => {
                if let Some(a) = to.azimuth {
                    pose.azimuth = lerp(pose.azimuth, a, k);
                }
                if let Some(e) = to.elevation {
                    pose.elevation = lerp(pose.elevation, e, k);
                }
                if let Some(z) = to.zoom {
                    pose.zoom = lerp(pose.zoom.pln(), z.pln(), k).pexp();
                }
                if let Some(target) = to.target.as_deref().filter(|s| !s.is_empty()) {
                    let mut next: Vec<(String, f64)> = Vec::new();
                    let mut add =
                        |name: &str, w: f64| match next.iter_mut().find(|(n, _)| n == name) {
                            Some((_, v)) => *v += w,
                            None => next.push((name.to_string(), w)),
                        };
                    for (name, w) in &pose.targets {
                        add(name, w * (1.0 - k));
                    }
                    add(target, k);
                    next.retain(|(_, w)| *w >= 1e-9);
                    pose.targets = next;
                }
            }
            _ => {}
        }
    }
    pose
}

/// One turn applied to a selection.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TurnState {
    /// A point on the axis.
    pub origin: [f64; 3],
    /// The axis direction.
    pub direction: [f64; 3],
    /// Degrees about it.
    pub degrees: f64,
}

/// What a frame does to one occurrence path or body selector.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ItemState {
    /// Translation added to its placement.
    pub offset: [f64; 3],
    /// Turns, applied in order before the offset.
    pub turns: Vec<TurnState>,
    /// 0..1.
    pub opacity: f64,
    /// Shown.
    pub visible: bool,
    /// 0..1 emphasis.
    pub highlight: f64,
}

/// The section a frame shows: two keys and the progress between them.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SectionState {
    /// The plane at progress 0.
    pub from: SectionKey,
    /// The plane at progress 1.
    pub to: SectionKey,
    /// Progress.
    pub k: f64,
    /// `"any"` or `"all"`.
    pub combine: String,
}

/// The caption a frame draws.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Caption {
    /// Title line.
    pub title: String,
    /// Subtitle line.
    pub subtitle: String,
    /// `"title"` (drawn large) or `"step"`.
    pub kind: String,
    /// Its step's index.
    pub index: usize,
    /// 0..1.
    pub opacity: f64,
}

/// Everything a frame at a time needs except the camera.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FilmState {
    /// Per path, in the order given.
    pub items: Vec<(String, ItemState)>,
    /// Joint values, in first-appearance order.
    pub joints: Vec<(String, f64)>,
    /// The explode factor.
    pub explode: f64,
    /// The section, if any.
    pub section: Option<SectionState>,
    /// The caption, if any.
    pub caption: Option<Caption>,
}

fn of_kind(film: &Film, pick: fn(&Channel) -> bool) -> Vec<&Track> {
    film.tracks.iter().filter(|tr| pick(&tr.channel)).collect()
}

/// The state of `film` at `t` for the scene's `paths` (occurrence paths
/// or `"body N"` selectors).
pub fn state_at(film: &Film, t: f64, paths: &[String]) -> FilmState {
    let mut items: Vec<(String, ItemState)> = paths
        .iter()
        .map(|p| {
            (
                p.clone(),
                ItemState {
                    offset: [0.0; 3],
                    turns: Vec::new(),
                    opacity: 1.0,
                    visible: true,
                    highlight: 0.0,
                },
            )
        })
        .collect();
    for track in &film.tracks {
        if let Channel::Move { target, from, to } = &track.channel {
            let k = progress(track, t);
            let d = [0, 1, 2].map(|i| lerp(from[i], to[i], k));
            for (_, item) in items.iter_mut().filter(|(p, _)| selects(target, p)) {
                for (o, d) in item.offset.iter_mut().zip(d) {
                    *o += d;
                }
            }
        }
    }
    for track in &film.tracks {
        if let Channel::Turn {
            target,
            origin,
            direction,
            from,
            to,
        } = &track.channel
        {
            let degrees = lerp(*from, *to, progress(track, t));
            for (_, item) in items.iter_mut().filter(|(p, _)| selects(target, p)) {
                item.turns.push(TurnState {
                    origin: *origin,
                    direction: *direction,
                    degrees,
                });
            }
        }
    }
    let fades = of_kind(film, |c| matches!(c, Channel::Fade { .. }));
    let mut shows = of_kind(film, |c| matches!(c, Channel::Show { .. }));
    shows.sort_by(|a, b| a.start.total_cmp(&b.start));
    let highlights = of_kind(film, |c| matches!(c, Channel::Highlight { .. }));
    for (path, item) in &mut items {
        let mine: Vec<&Track> = fades
            .iter()
            .copied()
            .filter(
                |tr| matches!(&tr.channel, Channel::Fade { target, .. } if selects(target, path)),
            )
            .collect();
        if let Some((track, k)) = latest(&mine, t)
            && let Channel::Fade { from, to, .. } = &track.channel
        {
            item.opacity = lerp(*from, *to, k);
        }
        let my_shows: Vec<(f64, bool)> = shows
            .iter()
            .filter_map(|tr| match &tr.channel {
                Channel::Show { target, visible } if selects(target, path) => {
                    Some((tr.start, *visible))
                }
                _ => None,
            })
            .collect();
        if let Some(&(_, first)) = my_shows.first() {
            // Before the first show, the opposite of what it sets.
            let mut visible = !first;
            for &(start, v) in &my_shows {
                if start <= t {
                    visible = v;
                }
            }
            item.visible = visible;
        }
        for tr in &highlights {
            let Channel::Highlight { target } = &tr.channel else {
                continue;
            };
            if !selects(target, path) || t < tr.start || t > tr.start + tr.duration {
                continue;
            }
            let u = if tr.duration > 0.0 {
                (t - tr.start) / tr.duration
            } else {
                1.0
            };
            let h = ease(Ease::InOut, 1f64.min(u / 0.25).min((1.0 - u) / 0.25));
            item.highlight = item.highlight.max(h);
        }
    }
    let joint_tracks = of_kind(film, |c| matches!(c, Channel::Joint { .. }));
    let mut joints: Vec<(String, f64)> = Vec::new();
    for tr in &joint_tracks {
        let Channel::Joint { param, .. } = &tr.channel else {
            continue;
        };
        if joints.iter().any(|(n, _)| n == param) {
            continue;
        }
        let mine: Vec<&Track> = joint_tracks
            .iter()
            .copied()
            .filter(|t| matches!(&t.channel, Channel::Joint { param: p, .. } if p == param))
            .collect();
        if let Some((track, k)) = latest(&mine, t)
            && let Channel::Joint { from, to, .. } = &track.channel
        {
            joints.push((param.clone(), lerp(*from, *to, k)));
        }
    }
    let explode = latest(&of_kind(film, |c| matches!(c, Channel::Explode { .. })), t)
        .and_then(|(track, k)| match &track.channel {
            Channel::Explode { from, to } => Some(lerp(*from, *to, k)),
            _ => None,
        })
        .unwrap_or(0.0);
    let section = latest(&of_kind(film, |c| matches!(c, Channel::Section { .. })), t).and_then(
        |(track, k)| match &track.channel {
            Channel::Section { from, to, combine } => Some(SectionState {
                from: from.clone(),
                to: to.clone(),
                k,
                combine: combine.clone(),
            }),
            _ => None,
        },
    );
    FilmState {
        items,
        joints,
        explode,
        section,
        caption: caption_at(film, t),
    }
}

/// The caption at `t`: it fades in at its step's start and out at its
/// end, except where the next step carries the same caption on.
pub fn caption_at(film: &Film, t: f64) -> Option<Caption> {
    let steps = &film.steps;
    let index = steps
        .iter()
        .position(|s| t >= s.start && t < s.start + s.duration);
    let i = match index {
        Some(i) => i,
        None if t >= film.duration && !steps.is_empty() => steps.len() - 1,
        None => return None,
    };
    let step = &steps[i];
    if step.kind == "hold" || step.title.is_empty() {
        return None;
    }
    let same = |a: Option<&super::model::FilmStep>, b: Option<&super::model::FilmStep>| {
        matches!((a, b), (Some(a), Some(b))
            if a.title == b.title && a.subtitle == b.subtitle && b.kind != "hold")
    };
    let fade_in = if same(i.checked_sub(1).and_then(|j| steps.get(j)), Some(step)) {
        1.0
    } else {
        1f64.min((t - step.start) / CAPTION_FADE)
    };
    let end = step.start + step.duration;
    let last = i == steps.len() - 1;
    let fade_out = if same(Some(step), steps.get(i + 1)) || last {
        1.0
    } else {
        1f64.min((end - t) / CAPTION_FADE)
    };
    Some(Caption {
        title: step.title.clone(),
        subtitle: step.subtitle.clone(),
        kind: step.kind.to_string(),
        index: i,
        opacity: fade_in.min(fade_out).clamp(0.0, 1.0),
    })
}

/// A plane `normal · x = offset`; the side the normal points to is
/// removed.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Plane {
    /// Unit normal.
    pub normal: [f64; 3],
    /// Offset along it.
    pub offset: f64,
}

/// A section key as a plane in world space, the axis presets resolved
/// across `bounds`.
pub fn section_plane(key: &SectionKey, bounds: &[[f64; 3]; 2]) -> Option<Plane> {
    match key {
        SectionKey::Off => None,
        SectionKey::Plane { normal, offset } => Some(Plane {
            normal: *normal,
            offset: *offset,
        }),
        SectionKey::Axis { axis, at } => {
            let i = "xyz".find(axis.as_str())?;
            let mut normal = [0.0; 3];
            normal[i] = 1.0;
            Some(Plane {
                normal,
                offset: bounds[0][i] + (bounds[1][i] - bounds[0][i]) * at,
            })
        }
    }
}

/// A section state's plane at its progress: the two keys blended
/// (normals renormalised); a key that is off blends from or to the same
/// plane moved past the bounds, which removes nothing.
pub fn section_at(section: &SectionState, bounds: &[[f64; 3]; 2]) -> Option<Plane> {
    let a = section_plane(&section.from, bounds);
    let b = section_plane(&section.to, bounds);
    let beyond = |plane: Plane| {
        let mut max = f64::NEG_INFINITY;
        for c in 0..8 {
            let p = [
                bounds[usize::from(c & 1 != 0)][0],
                bounds[usize::from(c & 2 != 0)][1],
                bounds[usize::from(c & 4 != 0)][2],
            ];
            max = max.max(plane.normal[0] * p[0] + plane.normal[1] * p[1] + plane.normal[2] * p[2]);
        }
        Plane {
            normal: plane.normal,
            offset: max + 1e-6 * (1.0 + max.abs()),
        }
    };
    let (a, b) = match (a, b) {
        (None, None) => return None,
        (Some(a), Some(b)) => (a, b),
        (Some(a), None) => (a, beyond(a)),
        (None, Some(b)) => (beyond(b), b),
    };
    let k = section.k;
    let n = [0, 1, 2].map(|i| lerp(a.normal[i], b.normal[i], k));
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    let len = if len > 0.0 { len } else { 1.0 };
    Some(Plane {
        normal: n.map(|c| c / len),
        offset: lerp(a.offset, b.offset, k) / len,
    })
}

/// The frame times of an export: `round(duration · fps)` frames (at
/// least one) at `k / fps`.
pub fn frame_times(film: &Film, fps: f64) -> Vec<f64> {
    let count = ((film.duration * fps).round() as usize).max(1);
    (0..count).map(|k| k as f64 / fps).collect()
}
