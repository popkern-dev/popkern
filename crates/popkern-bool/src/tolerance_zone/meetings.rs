//! Where intersection records meet one another on a face: the face
//! curve arrangement's vertices (PK-1127).
//!
//! Every record lies on one face of each operand. Two records on one
//! face can cross or touch inside both: Villarceau's two circles cross
//! twice on the cutting plane, the six circles a cube cuts from its
//! midsphere touch at the edge midpoints, Viviani's figure eight
//! crosses itself. Such a meeting is a vertex of the result, yet no
//! record end and no edge-face hit lies there, so no candidate did,
//! and the records ran through it unsplit (measured: the complex
//! reported one edge's radial order changing between two stations).
//!
//! This stage finds them before the vertex table resolves identity:
//! each meeting becomes a candidate with an interior incidence on both
//! records (two on one record that meets itself), which resolution
//! decides with every other candidate and turns into paves. A meeting
//! is computed once per pair of records, so its identity is shared by
//! construction.
//!
//! **Finding them.** For records `C₁`, `C₂` the meetings are minima of
//! `d(s) = dist(C₁(s), C₂)` at most the tolerance. Searching over one
//! parameter, not two, keeps a tangential touch well posed: there the
//! two-parameter distance is stationary to fourth order along the
//! common tangent, while `d(s)` has an ordinary quadratic minimum. Each
//! sampled local minimum of the distance between the two curves'
//! polylines is refined by golden section on `d`.
//!
//! **Their extent.** A crossing at angle θ lies within the radius `r`
//! of both curves over `r / sin θ` along them, a touch with relative
//! curvature κ over about `√(2r / κ)`: either way the zone is measured,
//! the stretch of `C₁` within `r` of `C₂`, and carried by the candidate
//! so resolution clusters by it (two computations of one touch land
//! anywhere in that stretch).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use popkern_core::topology::{BRep, FaceHandle};

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use crate::intersect::IntersectionRecord;

/// Polyline stations per record.
const SAMPLES: usize = 128;
/// Golden-section iterations per refinement.
pub(super) const REFINE_STEPS: usize = 60;
/// A meeting whose zone covers more than this fraction of the shorter
/// curve is an overlap of the two, not a meeting.
pub(super) const OVERLAP_FRACTION: f64 = 0.05;
/// Stations of the other record per pruning group.
const GROUP: usize = 8;
/// Records leaving a shared end at more than this sine of an angle
/// (about 6 degrees) cannot meet again within a station of it.
const END_ANGLE_SIN: f64 = 0.1;

/// Two records meeting inside at least one of them.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecordMeeting {
    /// The records (equal for a record meeting itself).
    pub records: [usize; 2],
    /// The parameter of the meeting on each.
    pub params: [f64; 2],
    /// Where they meet: the midpoint of their two points.
    pub point: Point3,
    /// The distance between the two points.
    pub gap: f64,
    /// The stretch of the first record within the radius of the second
    /// about the meeting (its larger side).
    pub zone: f64,
    /// Sine of the angle between the curves there.
    pub sin_angle: f64,
}

/// The faces a record lies on, as `(operand, face bits)`.
fn faces_of(r: &IntersectionRecord) -> [(usize, u64); 2] {
    use slotmap::Key;
    [(0, r.face_a.data().as_ffi()), (1, r.face_b.data().as_ffi())]
}

/// Every meeting of `records` that share a face, within `radius`, that
/// lies inside at least one of the two (meetings at both records' ends
/// are their end candidates already).
pub(crate) fn record_meetings(records: &[IntersectionRecord], radius: f64) -> Vec<RecordMeeting> {
    // Sampled only for a pair the closed form leaves (PK-1161).
    let samples: Vec<std::cell::OnceCell<Vec<(f64, Point3)>>> =
        records.iter().map(|_| std::cell::OnceCell::new()).collect();
    let sampled = |k: usize| samples[k].get_or_init(|| sample(&records[k].curve.curve_3d));
    let mut out = Vec::new();
    for i in 0..records.len() {
        if records[i].residue_owner_a.is_some() {
            continue;
        }
        for j in i..records.len() {
            if records[j].residue_owner_a.is_some() {
                continue;
            }
            let shares = faces_of(&records[i])
                .iter()
                .any(|f| faces_of(&records[j]).contains(f));
            if !shares {
                continue;
            }
            let (ci, cj): (&NurbsCurve, &NurbsCurve) =
                (&records[i].curve.curve_3d, &records[j].curve.curve_3d);
            // A segment never meets itself.
            if i == j && super::meetings_line::straight(ci).is_some() {
                continue;
            }
            let found = match (i != j)
                .then(|| super::meetings_line::line_meeting(ci, cj, radius))
                .flatten()
            {
                Some(line) => line.into_iter().collect(),
                None if i != j && !boxes_meet(sampled(i), sampled(j), radius) => Vec::new(),
                None => meetings_between(ci, sampled(i), cj, sampled(j), i == j, radius),
            };
            for m in found {
                let (si, sj, gap, zone, sin_angle) = m;
                // At an end along the curve, not merely where its end is:
                // a closed record's end lies on it, and a crossing there
                // with another branch of the same record is interior to
                // that branch.
                let at_end = |c: &NurbsCurve, s: f64| {
                    let (a, b) = c.parameter_range();
                    [a, b].iter().any(|&e| {
                        super::meetings_line::length_between(c, e.min(s), e.max(s))
                            <= zone.max(radius) + radius
                    })
                };
                if at_end(ci, si) && at_end(cj, sj) {
                    continue;
                }
                let (pi, pj) = (ci.point_at(si), cj.point_at(sj));
                out.push(RecordMeeting {
                    records: [i, j],
                    params: [si, sj],
                    point: Point3::new(
                        0.5 * (pi.x + pj.x),
                        0.5 * (pi.y + pj.y),
                        0.5 * (pi.z + pj.z),
                    ),
                    gap,
                    zone,
                    sin_angle,
                });
            }
        }
    }
    out
}

fn sample(c: &dyn CurveEval) -> Vec<(f64, Point3)> {
    let (a, b) = c.parameter_range();
    sample_between(c, a, b)
}

/// `SAMPLES` stations of `c` over `[a, b]`.
fn sample_between(c: &dyn CurveEval, a: f64, b: f64) -> Vec<(f64, Point3)> {
    (0..=SAMPLES)
        .map(|k| {
            let s = a + (b - a) * k as f64 / SAMPLES as f64;
            (s, c.point_at(s))
        })
        .collect()
}

/// Whether the two polylines' boxes, grown by their chords and the
/// radius, overlap.
fn boxes_meet(p: &[(f64, Point3)], q: &[(f64, Point3)], radius: f64) -> bool {
    let bounds = |pts: &[(f64, Point3)]| {
        let chord = pts
            .windows(2)
            .map(|w| w[0].1.distance_to(&w[1].1))
            .fold(0.0_f64, f64::max);
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (_, x) in pts {
            for (k, v) in [x.x, x.y, x.z].into_iter().enumerate() {
                lo[k] = lo[k].min(v - chord);
                hi[k] = hi[k].max(v + chord);
            }
        }
        (lo, hi)
    };
    let ((pl, ph), (ql, qh)) = (bounds(p), bounds(q));
    (0..3).all(|k| pl[k] <= qh[k] + radius && ql[k] <= ph[k] + radius)
}

/// Meetings of `ci` with `cj` (with itself when `same`): sampled local
/// minima of the stations' distance to the other polyline, refined.
/// Each is `(s on ci, s on cj, gap, zone, sin angle)`.
fn meetings_between(
    ci: &dyn CurveEval,
    pi: &[(f64, Point3)],
    cj: &dyn CurveEval,
    pj: &[(f64, Point3)],
    same: bool,
    radius: f64,
) -> Vec<(f64, f64, f64, f64, f64)> {
    let chord_of = |pts: &[(f64, Point3)]| {
        pts.windows(2)
            .map(|w| w[0].1.distance_to(&w[1].1))
            .fold(0.0_f64, f64::max)
    };
    let chord = chord_of(pi).max(chord_of(pj));
    // A record meets itself only away from its own neighbourhood: more
    // than a tenth of its parameter range apart, and not across its
    // closing seam.
    let (a, b) = (pi[0].0, pi[pi.len() - 1].0);
    let far = |k: usize, m: usize| {
        if !same {
            return true;
        }
        let gap = k.abs_diff(m).min(SAMPLES - k.abs_diff(m).min(SAMPLES));
        gap > SAMPLES / 10
    };
    // Per station of ci, its nearest station of cj, among those within
    // `chord + radius` (only those can start a meeting): cj's stations
    // in groups whose boxes lie farther are never measured, and a
    // station with none near reads infinitely far, which decides the
    // local-minimum tests below as its true distance would (both lie
    // past the bound every candidate is within).
    let reach = chord + radius;
    let groups: Vec<(usize, usize, [f64; 3], [f64; 3])> = (0..pj.len())
        .step_by(GROUP)
        .map(|g0| {
            let g1 = (g0 + GROUP).min(pj.len());
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for (_, q) in &pj[g0..g1] {
                for (k, v) in [q.x, q.y, q.z].into_iter().enumerate() {
                    lo[k] = lo[k].min(v);
                    hi[k] = hi[k].max(v);
                }
            }
            (g0, g1, lo, hi)
        })
        .collect();
    let nearest: Vec<(usize, f64)> = pi
        .iter()
        .enumerate()
        .map(|(k, (_, p))| {
            let mut best = (0, f64::INFINITY);
            for &(g0, g1, lo, hi) in &groups {
                let gap = [p.x, p.y, p.z]
                    .into_iter()
                    .enumerate()
                    .map(|(a, v)| (lo[a] - v).max(v - hi[a]).max(0.0))
                    .fold(0.0_f64, |acc, x| acc + x * x)
                    .sqrt();
                if gap > reach {
                    continue;
                }
                for (m, (_, q)) in pj.iter().enumerate().take(g1).skip(g0) {
                    if !far(k, m) {
                        continue;
                    }
                    let d = p.distance_to(q);
                    if d < best.1 {
                        best = (m, d);
                    }
                }
            }
            best
        })
        .collect();
    let mut out: Vec<(f64, f64, f64, f64, f64)> = Vec::new();
    for k in 0..pi.len() {
        let d = nearest[k].1;
        if d > chord + radius {
            continue;
        }
        let before = if k > 0 {
            nearest[k - 1].1
        } else {
            f64::INFINITY
        };
        let after = nearest.get(k + 1).map_or(f64::INFINITY, |n| n.1);
        if d > before || d > after {
            continue;
        }
        let m = nearest[k].0;
        // Two records leaving one shared end at an angle meet there
        // and nowhere else within a station of it: their end meeting is
        // dropped below (at an end of both), so it is not refined
        // (records of a box meet so at every corner; refining each cost
        // a box union most of its time).
        let last = |pts: &[(f64, Point3)]| pts.len() - 1;
        if !same && (k == 0 || k == last(pi)) && (m == 0 || m == last(pj)) && d <= radius {
            let (ti, tj) = (ci.tangent_at(pi[k].0), cj.tangent_at(pj[m].0));
            let sin = ti.cross(&tj).length() / (ti.length() * tj.length()).max(1e-300);
            if sin > END_ANGLE_SIN {
                continue;
            }
        }
        let lo = pi[k.saturating_sub(1)].0;
        let hi = pi[(k + 1).min(pi.len() - 1)].0;
        let (wj0, wj1) = (pj[m.saturating_sub(2)].0, pj[(m + 2).min(pj.len() - 1)].0);
        let seed = std::cell::Cell::new(pj[m].0);
        let dist = |s: f64| {
            let p = ci.point_at(s);
            let t = super::meetings_foot::warm_foot(cj, &p, (wj0, wj1), &seed);
            (p.distance_to(&cj.point_at(t)), t)
        };
        let s = super::meetings_foot::golden_min(lo, hi, |s| dist(s).0);
        let (gap, t) = dist(s);
        // A record meeting itself is found from both of its sides: keep
        // the one with the first parameter lower.
        if gap > radius || (same && s > t) {
            continue;
        }
        let zone = super::meetings_foot::zone_along(ci, s, (a, b), |x| dist(x).0, radius);
        // A stretch where the curves coincide (a record lying along a
        // boundary edge, cubes sharing an edge) is an overlap, a common
        // block of the two, not a vertex: every sampled minimum on it
        // would otherwise split the edge (measured: two cubes touching
        // along an edge grew from 16 to 92 vertices). Probed against
        // the whole of the other curve, a twentieth of this one's range
        // either side: a touch has left the radius long before.
        if !same {
            let (c0, c1) = (pj[0].0, pj[pj.len() - 1].0);
            let step = OVERLAP_FRACTION * (b - a);
            let still_on = |x: f64| {
                let x = x.clamp(a.min(b), a.max(b));
                let p = ci.point_at(x);
                p.distance_to(&cj.point_at(closest_point_in_window(cj, &p, c0, c1))) <= radius
            };
            if still_on(s - step) && still_on(s + step) {
                continue;
            }
        }
        // One meeting found from neighbouring stations (inside a touch's
        // flat stretch the refinement lands anywhere): one is kept.
        let here = ci.point_at(s);
        if out
            .iter()
            .any(|o| ci.point_at(o.0).distance_to(&here) <= o.3.max(zone) + radius)
        {
            continue;
        }
        let (ti, tj) = (ci.tangent_at(s), cj.tangent_at(t));
        let sin_angle = ti.cross(&tj).length() / (ti.length() * tj.length()).max(1e-300);
        out.push((s, t, gap, zone, sin_angle));
    }
    out
}

/// A face's tolerance, the model's linear where it has none.
fn face_tol(brep: &BRep, face: FaceHandle, linear: f64) -> f64 {
    brep.faces
        .get(face)
        .and_then(|f| f.tolerance)
        .unwrap_or(linear)
}

/// The candidates of the records' meetings: each lies inside its
/// records (an interior incidence on each; on both branches of a record
/// meeting itself) and on their faces, with the radius of the faces'
/// tolerances and the zone the meeting measured.
pub(super) fn meeting_candidates(
    (a, b): (&BRep, &BRep),
    records: &[IntersectionRecord],
    linear: f64,
) -> Vec<CandidateVertex> {
    let mut out = Vec::new();
    let radius_of =
        |r: &IntersectionRecord| face_tol(a, r.face_a, linear) + face_tol(b, r.face_b, linear);
    let widest = records.iter().map(radius_of).fold(2.0 * linear, f64::max);
    let found = record_meetings(records, widest);
    if popkern_core::trace::trace_enabled("arrangement") {
        for (i, r) in records.iter().enumerate() {
            let c = &r.curve.curve_3d;
            let (lo, hi) = c.parameter_range();
            eprintln!(
                "  [arrangement] record {i}: degree {}, {} control points, closed {}, from {:?} to {:?}",
                c.degree,
                c.control_points.len(),
                r.curve.is_closed,
                c.point_at(lo),
                c.point_at(hi)
            );
        }
        for m in &found {
            eprintln!(
                "  [arrangement] records {:?} at {:?}: gap {:.2e}, zone {:.2e}, sin {:.2e}",
                m.records, m.point, m.gap, m.zone, m.sin_angle
            );
        }
    }
    out.extend(found.into_iter().filter_map(|m| {
        let [i, j] = m.records;
        let radius = radius_of(&records[i]).max(radius_of(&records[j]));
        // Two records meet where they are one place: within the linear
        // reach of each other, plus each record's own fit deviation
        // there (its point off its two surfaces), as an edge meets a
        // face within the reach (PK-1180).
        let deviation = |k: usize, s: f64| {
            let r = &records[k];
            super::candidates::record_fit_deviation(a, b, r, false, &r.curve.curve_3d.point_at(s))
        };
        let reach = popkern_core::tolerance::Tolerance::default()
            .with_linear(linear)
            .linear_reach();
        if m.gap > radius.min(reach + deviation(i, m.params[0]) + deviation(j, m.params[1])) {
            return None;
        }
        // A meeting at one record's end (another record passing
        // through the end of this one) is that end, whose candidate
        // exists already: the meeting adds only the other record's
        // pave, and joins the end by proximity. Carrying the end too
        // put it a whole circle along a closed record from that
        // record's other end, and the two stayed apart.
        let end_of = |r: &IntersectionRecord, s: f64| {
            let c = &r.curve.curve_3d;
            let (lo, hi) = c.parameter_range();
            let reach = m.zone.max(radius) + radius;
            if super::meetings_line::length_between(c, lo, s).abs() <= reach {
                Some(false)
            } else if super::meetings_line::length_between(c, s, hi).abs() <= reach {
                Some(true)
            } else {
                None
            }
        };
        let mut incidence: Vec<Incidence> = [(i, m.params[0]), (j, m.params[1])]
            .into_iter()
            .filter(|&(r, s)| end_of(&records[r], s).is_none())
            .map(|(record, s)| Incidence::Record {
                record,
                s,
                end: None,
            })
            .collect();
        for r in [&records[i], &records[j]] {
            for (operand, face) in [(0, r.face_a), (1, r.face_b)] {
                let inc = Incidence::Face {
                    operand,
                    face,
                    uv: None,
                };
                if !incidence.contains(&inc) {
                    incidence.push(inc);
                }
            }
        }
        Some(CandidateVertex {
            position: m.point,
            radius,
            zone: m.zone.max(radius),
            source: SourceKind::RecordMeeting,
            incidence,
            key: (
                SourceKind::RecordMeeting,
                i,
                j as u64,
                m.params[0].to_bits(),
                m.params[1].to_bits(),
            ),
        })
    }));
    out
}

#[cfg(test)]
#[path = "meetings_tests.rs"]
mod tests;
