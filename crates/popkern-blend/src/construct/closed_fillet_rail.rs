//! Rails of a closed fillet whose contacts are not circles (PK-888):
//! the rim of a hole drilled at an angle, a plane x oblique cylinder
//! ellipse.
//!
//! The walker's contacts are exact at its stations; the trim the
//! geometry pipeline fits through them is a degree-1 polyline, whose
//! chords sag off a curved support between stations (measured: 5.0e-3
//! off an r5 hole wall at 77 stations, 6.8e-2 at 20). A cubic through
//! the same contacts lies within 2e-6..4e-6 of the wall there. The rail
//! also has to START at the surgery's anchor: on the side whose loop
//! carries the support's seam, the anchor is where the rail crosses the
//! seam's generator line, so the shortened seam stays a generator —
//! the contact nearest the old seam vertex is off the generator on a
//! tilted rim, and the seam line left the wall by 6.4e-4 at 8 deg and
//! 2.4e-2 at 25 deg.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::construct::EdgeGeometry;

/// Samples of the rail when locating its seam-line crossing.
const CROSSING_SAMPLES: usize = 2048;

/// The median spacing of consecutive `pts`.
fn median_step(pts: &[Point3]) -> f64 {
    let mut steps: Vec<f64> = pts.windows(2).map(|w| w[0].distance_to(&w[1])).collect();
    steps.sort_by(f64::total_cmp);
    steps.get(steps.len() / 2).copied().unwrap_or(0.0)
}

/// The Newell axis of the closed polygon `pts` (unnormalised): the
/// polygon runs counter-clockwise about it.
fn newell_axis(pts: &[Point3]) -> Vec3 {
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    for (i, a) in pts.iter().enumerate() {
        let b = pts[(i + 1) % pts.len()];
        n = n + Vec3::new(
            (a.y - b.y) * (a.z + b.z),
            (a.z - b.z) * (a.x + b.x),
            (a.x - b.x) * (a.y + b.y),
        );
    }
    n
}

/// Points closer than this fraction of the station spacing are one
/// station: a near-duplicate becomes a near-zero chord, and the cubic
/// rings through it (measured: the walk's closing contact repeated
/// its first ~1e-7 away and the rail rang 9e-3 off the wall).
const SAME_STATION: f64 = 1e-3;

/// The exact contacts behind a closed trim (its control points), the
/// repeated closing contact dropped.
fn closed_contacts(trim: &NurbsCurve) -> Vec<Point3> {
    let mut pts = trim.control_points.clone();
    let same = median_step(&pts) * SAME_STATION;
    while pts.len() > 2
        && let (Some(first), Some(last)) = (pts.first(), pts.last())
        && first.distance_to(last) <= same
    {
        pts.pop();
    }
    pts
}

/// A cubic through `pts` (closed back to its first point) by chord
/// length; `None` when the points cannot carry one.
fn closed_cubic(pts: &[Point3]) -> Option<NurbsCurve> {
    let mut ring: Vec<Point3> = pts.to_vec();
    ring.push(*pts.first()?);
    let mut params = vec![0.0];
    for w in ring.windows(2) {
        let step = w[0].distance_to(&w[1]);
        if step <= 0.0 {
            return None;
        }
        params.push(params.last()? + step);
    }
    popkern_core::geometry::fitting::interpolate_curve(&ring, &params, 3.min(ring.len() - 1)).ok()
}

/// The rail of a non-circular closed trim as an edge curve starting
/// and ending at `start` (a point on the rail): a cubic through the
/// contacts from `start` round. Falls back to the fitted trim when the
/// rail cannot be rebuilt.
///
/// The rail runs counter-clockwise about `winding_axis` (the old rim
/// half-edge's traversal axis on this support) so the support's loop
/// keeps its direction; `None` keeps the stations' order.
pub(super) fn rail_on_edge(
    trim: &NurbsCurve,
    start: Point3,
    winding_axis: Option<Vec3>,
) -> CurveOnEdge {
    let fallback = || CurveOnEdge {
        geometry: CurveGeometry::Nurbs(trim.clone()),
        t_start: 0.0,
        t_end: 1.0,
    };
    let mut pts = closed_contacts(trim);
    let n = pts.len();
    if n < 4 {
        return fallback();
    }
    if winding_axis.is_some_and(|axis| newell_axis(&pts).dot(&axis) < 0.0) {
        pts.reverse();
    }
    // The chord the start lies on: its contacts follow it round.
    let chord_distance = |k: usize| {
        let (a, b) = (pts[k], pts[(k + 1) % n]);
        let ab = b - a;
        let t = ((start - a).dot(&ab) / ab.dot(&ab).max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
        (a + ab * t).distance_to(&start)
    };
    let Some(k) = (0..n).min_by(|&i, &j| chord_distance(i).total_cmp(&chord_distance(j))) else {
        return fallback();
    };
    let gap = median_step(&pts) * SAME_STATION;
    let mut rotated = vec![start];
    rotated.extend(
        (1..=n)
            .map(|i| pts[(k + i) % n])
            .filter(|p| p.distance_to(&start) > gap),
    );
    match closed_cubic(&rotated) {
        Some(curve) => {
            let (t_start, t_end) = curve.parameter_range();
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve),
                t_start,
                t_end,
            }
        }
        None => fallback(),
    }
}

/// The anchor on the seam side of a non-circular closed rim: where the
/// rail on `seam_face`'s support crosses `seam_edge`'s generator line.
/// `None` when that rail is an exact circle (the circular anchor owns
/// it) or the seam is not a line.
pub(super) fn seam_generator_anchor(
    brep: &BRep,
    eg: &EdgeGeometry,
    seam_edge: Option<EdgeHandle>,
    seam_face: FaceHandle,
) -> Option<Point3> {
    let (trim, circle) = if seam_face == eg.analysis.face_1 {
        (&eg.trim_3d_1, eg.trim_circle_1)
    } else {
        (&eg.trim_3d_2, eg.trim_circle_2)
    };
    if circle.is_some() {
        return None;
    }
    let CurveGeometry::Line { origin, direction } =
        brep.edges.get(seam_edge?)?.curve.as_ref()?.geometry
    else {
        return None;
    };
    let axis: Vec3 = direction.normalized()?;
    let rail = closed_cubic(&closed_contacts(trim))?;
    let off_line = |p: Point3| {
        let d = p - origin;
        (d - axis * d.dot(&axis)).length()
    };
    let (lo, hi) = rail.parameter_range();
    let at = |i: usize| lo + (hi - lo) * i as f64 / CROSSING_SAMPLES as f64;
    let best = (0..CROSSING_SAMPLES).min_by(|&i, &j| {
        off_line(rail.point_at(at(i))).total_cmp(&off_line(rail.point_at(at(j))))
    })?;
    // Golden-section refinement on the bracketing samples (the closed
    // rail wraps, so the bracket is taken modulo its range).
    let span = (hi - lo) / CROSSING_SAMPLES as f64;
    let wrap = |t: f64| lo + (t - lo).rem_euclid(hi - lo);
    let (mut a, mut b) = (at(best) - span, at(best) + span);
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..80 {
        let (c, d) = (b - g * (b - a), a + g * (b - a));
        if off_line(rail.point_at(wrap(c))) < off_line(rail.point_at(wrap(d))) {
            b = d;
        } else {
            a = c;
        }
    }
    Some(rail.point_at(wrap(0.5 * (a + b))))
}

/// The anchors of a non-circular closed rim, both ON the strip at the
/// section through the seam side's generator crossing
/// ([`seam_generator_anchor`]). The strip's
/// seam joins them along that section's iso line, so they share its
/// station (the contact nearest the old vertex sat a station-fraction
/// away and the seam ran across the strip, its pcurve 0.13 off), and
/// they lie on the strip itself (the fitted strip parts from the
/// support between stations by a sub-precision residual — 4.6e-6 at a
/// 25 deg rim — and an iso seam ending that far from its vertex fell
/// back to a polyline sagging 1.8e-3). The residual moves to the
/// support side, where the declared-accuracy pass records it. `None`
/// keeps the anchors.
pub(super) fn closed_rim_anchors(
    brep: &BRep,
    eg: &EdgeGeometry,
    seam_edge: Option<EdgeHandle>,
    seam_face: FaceHandle,
) -> Option<(Point3, Point3)> {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let anchor = seam_generator_anchor(brep, eg, seam_edge, seam_face)?;
    let SurfaceGeometry::Nurbs(strip) = &eg.surface else {
        return None;
    };
    let (u, v) = strip.closest_point(&anchor)?;
    let (v0, v1) = (*strip.knots_v.first()?, *strip.knots_v.last()?);
    let (own, other) = if (v - v0).abs() <= (v - v1).abs() {
        (v0, v1)
    } else {
        (v1, v0)
    };
    Some((strip.point_at(u, own), strip.point_at(u, other)))
}
