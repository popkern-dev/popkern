//! Film overlap warnings (PK-855): whether a film's move and turn
//! tracks drive one part through another. The tracks are evaluated as
//! the page plays them (`web/live/film-state.js` `stateAt`: every track
//! holds its first value before its start and its last after its end;
//! a path's offsets add, its turns apply in track order before the
//! offset), and at each frame where a motion track runs the moved
//! occurrences are checked against the assembly with
//! [`Assembly::overlaps_under`] — a proven overlap of the base meshes,
//! the document unchanged. A warning, not an error: a film may show a
//! part passing through a cut-away housing on purpose.

use std::collections::BTreeMap;

use popkern_assembly::{Assembly, AssemblyError, Placement};
use popkern_core::geometry::{Point3, Vec3};
use serde::Serialize;

use super::model::{Channel, Ease, Film, Track};

/// The first frame at which a film drives two occurrences into each
/// other.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FilmOverlap {
    /// Seconds into the film.
    pub time: f64,
    /// The frame at the film's rate.
    pub frame: u32,
    /// The two occurrence paths.
    pub pair: [String; 2],
    /// 1-based script line of a motion track moving one of them then,
    /// when known.
    pub line: Option<usize>,
}

/// Every occurrence pair `film` drives into overlap, at the first frame
/// it does (ordered by time). `unit_scale` is metres per document unit
/// (the film's lengths are SI).
///
/// # Errors
///
/// The interference check's errors (a part that cannot be built or
/// tessellated).
pub fn film_overlaps(
    assembly: &mut Assembly,
    film: &Film,
    unit_scale: f64,
) -> Result<Vec<FilmOverlap>, AssemblyError> {
    let motions: Vec<&Track> = film
        .tracks
        .iter()
        .filter(|t| match &t.channel {
            // Moving everything together changes no pair.
            Channel::Move { target, .. } | Channel::Turn { target, .. } => target != "all",
            _ => false,
        })
        .collect();
    if motions.is_empty() {
        return Ok(Vec::new());
    }
    let to_doc = if unit_scale.is_finite() && unit_scale > 0.0 {
        1.0 / unit_scale
    } else {
        1.0
    };
    let fps = f64::from(film.fps.max(1));
    let last = (film.duration * fps).ceil().max(0.0) as u32;
    // Frame 0 and every frame a motion track runs over: between them the
    // scene holds the pose an earlier frame already checked.
    let frames: Vec<u32> = (0..=last)
        .filter(|&f| {
            let t = f64::from(f) / fps;
            f == 0
                || motions
                    .iter()
                    .any(|m| t >= m.start && t <= m.start + m.duration + 0.5 / fps)
        })
        .collect();
    let leaves: Vec<String> = assembly
        .leaves()
        .into_iter()
        .map(|row| row.path.to_string())
        .collect();
    let mut first: BTreeMap<[String; 2], FilmOverlap> = BTreeMap::new();
    for frame in frames {
        let t = f64::from(frame) / fps;
        let pose = pose_at(&motions, t, &leaves, to_doc)?;
        if pose.is_empty() {
            continue;
        }
        let pairs =
            assembly.overlaps_under(&|p| pose.get(p).cloned(), &|p| pose.contains_key(p))?;
        for (a, b) in pairs {
            let pair = if a <= b { [a, b] } else { [b, a] };
            if first.contains_key(&pair) {
                continue;
            }
            let line = motions
                .iter()
                .filter(|m| m.start <= t && pair.iter().any(|p| selects(target_of(m), p)))
                .filter_map(|m| m.line)
                .next_back();
            first.insert(
                pair.clone(),
                FilmOverlap {
                    time: t,
                    frame,
                    pair,
                    line,
                },
            );
        }
    }
    let mut out: Vec<FilmOverlap> = first.into_values().collect();
    out.sort_by(|x, y| x.time.total_cmp(&y.time).then_with(|| x.pair.cmp(&y.pair)));
    Ok(out)
}

/// The world motion `film` gives occurrence `path` at `t` (seconds),
/// in document units (`unit_scale` metres each): its tracks' offsets
/// and turns as the page composes them. `None` where no track moves it.
///
/// # Errors
///
/// [`AssemblyError`] for a turn about a degenerate axis.
pub fn occurrence_motion(
    film: &Film,
    path: &str,
    t: f64,
    unit_scale: f64,
) -> Result<Option<Placement>, AssemblyError> {
    let motions: Vec<&Track> = film
        .tracks
        .iter()
        .filter(|t| matches!(t.channel, Channel::Move { .. } | Channel::Turn { .. }))
        .collect();
    let to_doc = if unit_scale.is_finite() && unit_scale > 0.0 {
        1.0 / unit_scale
    } else {
        1.0
    };
    Ok(pose_at(&motions, t, &[path.to_string()], to_doc)?.remove(path))
}

/// Each moved leaf occurrence's world motion at `t` (absent: in place).
fn pose_at(
    motions: &[&Track],
    t: f64,
    leaves: &[String],
    to_doc: f64,
) -> Result<BTreeMap<String, Placement>, AssemblyError> {
    let mut out = BTreeMap::new();
    for path in leaves {
        let mine: Vec<&&Track> = motions
            .iter()
            .filter(|m| selects(target_of(m), path))
            .collect();
        if mine.is_empty() {
            continue;
        }
        let mut motion = Placement::identity();
        let mut offset = [0.0_f64; 3];
        for track in mine {
            let k = progress(track, t);
            match &track.channel {
                Channel::Move { from, to, .. } => {
                    for (axis, value) in offset.iter_mut().enumerate() {
                        *value += (from[axis] + (to[axis] - from[axis]) * k) * to_doc;
                    }
                }
                Channel::Turn {
                    origin,
                    direction,
                    from,
                    to,
                    ..
                } => {
                    let degrees = from + (to - from) * k;
                    let turn = Placement::rotation(
                        Point3::new(origin[0] * to_doc, origin[1] * to_doc, origin[2] * to_doc),
                        Vec3::new(direction[0], direction[1], direction[2]),
                        degrees,
                    )?;
                    motion = motion.then(&turn);
                }
                _ => {}
            }
        }
        let motion = motion.then(&Placement::translation(offset[0], offset[1], offset[2]));
        if !motion.is_identity() {
            out.insert(path.clone(), motion);
        }
    }
    Ok(out)
}

fn target_of(track: &Track) -> &str {
    match &track.channel {
        Channel::Move { target, .. } | Channel::Turn { target, .. } => target,
        _ => "",
    }
}

/// `film-state.js` `selects`: a prefix selects what lies under it.
fn selects(target: &str, path: &str) -> bool {
    target == "all"
        || path == target
        || path
            .strip_prefix(target)
            .is_some_and(|rest| rest.starts_with('/'))
}

/// `film-state.js` `progress`: 0 before the start, 1 after the end.
fn progress(track: &Track, t: f64) -> f64 {
    if t < track.start {
        return 0.0;
    }
    if track.duration <= 0.0 || t >= track.start + track.duration {
        return 1.0;
    }
    let p = ((t - track.start) / track.duration).clamp(0.0, 1.0);
    match track.ease {
        Ease::Linear => p,
        Ease::In => p * p,
        Ease::Out => 1.0 - (1.0 - p) * (1.0 - p),
        Ease::InOut => p * p * (3.0 - 2.0 * p),
    }
}
