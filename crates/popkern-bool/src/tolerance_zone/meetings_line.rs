//! Meetings of two straight records in closed form (PK-1161).
//!
//! Most records of planar booleans are segments, degree-1 curves of two
//! control points. Two segments meet at most once, where they pass
//! closest, which the segment-segment closest-point solution gives
//! exactly; the sampled search (129 stations each, a nearest-station
//! pass and a golden-section refinement per local minimum) cost a box
//! union 0.17 ms of its 1.3, and a 50-face prism 0.26 of 2.4. Parallel
//! segments, which meet along a stretch or not at all, are left to the
//! sampled search and its overlap rule. A meeting's zone, the stretch
//! within the radius of the other segment, and a straight record's
//! length between two parameters are closed form too: measured by
//! doubling and bisection and by quadrature they still cost the prism
//! 0.05 ms for its five T-junctions.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

/// Below this sine of the angle between two segments they are parallel
/// here and left to the sampled search.
const PARALLEL_SIN: f64 = 1e-9;

/// A straight record: its parameter range and end points, when it is a
/// degree-1 curve of two equally weighted control points on a clamped
/// knot vector (its point then moves linearly with its parameter).
pub(super) fn straight(c: &NurbsCurve) -> Option<((f64, f64), Point3, Point3)> {
    let clamped = c.knots.len() == 4 && c.knots[0] == c.knots[1] && c.knots[2] == c.knots[3];
    let even = c.weights.len() == 2 && (c.weights[0] - c.weights[1]).abs() <= 1e-12 * c.weights[0];
    if c.degree != 1 || c.control_points.len() != 2 || !clamped || !even {
        return None;
    }
    Some((
        (c.knots[0], c.knots[3]),
        c.control_points[0],
        c.control_points[1],
    ))
}

/// A meeting as the sampled search reports it: `(s on ci, s on cj,
/// gap, zone, sin angle)`.
type Meeting = (f64, f64, f64, f64, f64);

/// The meeting of two straight records within `radius`, as the sampled
/// search reports it: `(s on ci, s on cj, gap, zone, sin angle)`.
/// `None` when the two are not both straight or are parallel within the
/// radius (the sampled search decides those); `Some(None)` when they do
/// not meet.
pub(super) fn line_meeting(
    ci: &NurbsCurve,
    cj: &NurbsCurve,
    radius: f64,
) -> Option<Option<Meeting>> {
    let ((a, b), p0, p1) = straight(ci)?;
    let ((c, d), q0, q1) = straight(cj)?;
    let (d1, d2) = (p1 - p0, q1 - q0);
    let (l1, l2) = (d1.dot(&d1), d2.dot(&d2));
    if l1 <= 0.0 || l2 <= 0.0 {
        return None;
    }
    let sin_angle = d1.cross(&d2).length() / (l1 * l2).sqrt();
    if sin_angle < PARALLEL_SIN {
        // Parallel lines farther apart than the radius never meet;
        // nearer, they share a stretch or touch end to end, which the
        // sampled search's overlap and end rules decide.
        let r = q0 - p0;
        let across = (r - d1 * (r.dot(&d1) / l1)).length();
        return (across > radius).then_some(None);
    }
    let (u, v) = closest_fractions(p0, d1, q0, d2);
    // Closest at an end of both: records leaving a shared end (every
    // corner of a box), whose meeting there is their end candidates'
    // and is dropped by the caller; not measured.
    let at_end = |w: f64| w == 0.0 || w == 1.0;
    if at_end(u) && at_end(v) {
        return Some(None);
    }
    let (s, t) = (a + u * (b - a), c + v * (d - c));
    let gap = ci.point_at(s).distance_to(&cj.point_at(t));
    if gap > radius {
        return Some(None);
    }
    // The distance from a point of ci to the segment cj (the overlap
    // rule below).
    let dist = |x: f64| {
        let p = ci.point_at(x);
        let w = ((p - q0).dot(&d2) / l2).clamp(0.0, 1.0);
        p.distance_to(&(q0 + d2 * w))
    };
    let zone = segment_zone((p0, d1), u, (q0, d2), radius);
    // Within the radius a twentieth of the range either side is a
    // stretch the two share, not a meeting (the sampled search's rule).
    let step = super::meetings::OVERLAP_FRACTION * (b - a);
    let still_on = |x: f64| dist(x.clamp(a.min(b), a.max(b))) <= radius;
    if still_on(s - step) && still_on(s + step) {
        return Some(None);
    }
    Some(Some((s, t, gap, zone, sin_angle)))
}

/// The arc length of segment `p0 + u d1` about fraction `u` within
/// `radius` of segment `q0 + v d2`, the larger of its two sides, `u`
/// itself within it (the closed form of the sampled search's zone).
/// The distance from a point moving along the first to the second is
/// convex, and quadratic in each piece the moving foot passes (before
/// the second's start, along it, past its end): along each side the
/// zone ends at the first piece's exit root, or at the first segment's
/// end.
fn segment_zone(
    (p0, d1): (Point3, popkern_core::geometry::Vec3),
    u: f64,
    (q0, d2): (Point3, popkern_core::geometry::Vec3),
    radius: f64,
) -> f64 {
    use popkern_core::geometry::Vec3;
    let (l1, l2) = (d1.dot(&d1), d2.dot(&d2));
    let x = p0 + d1 * u;
    let q1 = q0 + d2;
    let r2 = radius * radius;
    // The larger root of |c + t e|² = r² (t from 0, starting inside),
    // or None when the distance does not grow along e.
    let exit = |c: Vec3, e: Vec3| {
        let (ee, ce, cc) = (e.dot(&e), c.dot(&e), c.dot(&c));
        (ee > 0.0).then(|| (-ce + (ce * ce - ee * (cc - r2)).max(0.0).sqrt()) / ee)
    };
    let mut widest: f64 = 0.0;
    for dir in [1.0_f64, -1.0] {
        let room = if dir > 0.0 { 1.0 - u } else { u };
        let e = d1 * dir;
        // The foot's fraction on the second segment, linear in t.
        let v0 = (x - q0).dot(&d2) / l2;
        let dv = e.dot(&d2) / l2;
        let mut cuts: Vec<f64> = [-v0 / dv, (1.0 - v0) / dv]
            .into_iter()
            .filter(|t| t.is_finite() && *t > 0.0 && *t < room)
            .collect();
        cuts.sort_by(f64::total_cmp);
        cuts.push(room);
        let mut start = 0.0;
        let mut reached = room;
        for end in cuts {
            let mid = 0.5 * (start + end);
            let v = v0 + dv * mid;
            let (c, e_piece) = if v < 0.0 {
                (x - q0, e)
            } else if v > 1.0 {
                (x - q1, e)
            } else {
                let off = x - q0;
                (off - d2 * (off.dot(&d2) / l2), e - d2 * dv)
            };
            if let Some(t) = exit(c, e_piece)
                && t <= end
            {
                reached = t.max(start);
                break;
            }
            start = end;
        }
        widest = widest.max(reached * l1.sqrt());
    }
    widest
}

/// The length of a record between parameters `s0` and `s1`: exact for a
/// straight record (its point moves linearly), else the curve's own.
pub(super) fn length_between(c: &NurbsCurve, s0: f64, s1: f64) -> f64 {
    match straight(c) {
        Some(((a, b), p, q)) if b != a => ((s1 - s0) / (b - a)).abs() * p.distance_to(&q),
        _ => c.length_between(s0, s1),
    }
}

/// The fractions `(u, v)` along segments `p0 + u d1` and `q0 + v d2`,
/// each in [0, 1], of their closest points (the segments not parallel).
fn closest_fractions(
    p0: Point3,
    d1: popkern_core::geometry::Vec3,
    q0: Point3,
    d2: popkern_core::geometry::Vec3,
) -> (f64, f64) {
    let r = p0 - q0;
    let (a, e, f) = (d1.dot(&d1), d2.dot(&d2), d2.dot(&r));
    let (b, c) = (d1.dot(&d2), d1.dot(&r));
    let denom = a * e - b * b;
    let mut u = ((b * f - c * e) / denom).clamp(0.0, 1.0);
    let mut v = (b * u + f) / e;
    if v < 0.0 {
        v = 0.0;
        u = (-c / a).clamp(0.0, 1.0);
    } else if v > 1.0 {
        v = 1.0;
        u = ((b - c) / a).clamp(0.0, 1.0);
    }
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segment(p: Point3, q: Point3) -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 2.0, 2.0],
            control_points: vec![p, q],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn crossing_segments_meet_where_they_pass_closest() {
        let x = segment(Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let y = segment(Point3::new(0.5, -1.0, 4e-7), Point3::new(0.5, 1.0, 4e-7));
        let (s, t, gap, zone, sin) = line_meeting(&x, &y, 1e-6)
            .expect("both straight")
            .expect("they meet");
        assert!(
            (s - 1.5).abs() < 1e-12 && (t - 1.0).abs() < 1e-12,
            "{s} {t}"
        );
        assert!((gap - 4e-7).abs() < 1e-15, "{gap}");
        assert!((sin - 1.0).abs() < 1e-12);
        // Within 1e-6 of the other segment over sqrt(1e-12 - 1.6e-13).
        assert!((zone - (1e-12_f64 - 1.6e-13).sqrt()).abs() < 1e-9, "{zone}");
    }

    /// The closed-form zone against the doubling-and-bisection search it
    /// replaces, over crossings, shallow crossings, T-junctions at and
    /// near either end, and zones cut short by the first segment's end.
    #[test]
    fn the_segment_zone_is_the_searched_zone() {
        let x = segment(Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let r = 1e-6;
        let cases = [
            (Point3::new(0.5, -1.0, 4e-7), Point3::new(0.5, 1.0, 4e-7)),
            (Point3::new(-0.7, -1e-4, 3e-7), Point3::new(0.9, 1e-4, 3e-7)),
            (Point3::new(0.2, 0.0, 5e-7), Point3::new(0.2, 1.0, 5e-7)),
            (Point3::new(0.2, 3e-7, 0.0), Point3::new(0.2, 1.0, 0.0)),
            (Point3::new(0.2, -3e-7, 0.0), Point3::new(0.2, 1.0, 0.0)),
            (
                Point3::new(1.0 - 2e-7, -1.0, 0.0),
                Point3::new(1.0 - 2e-7, 1.0, 0.0),
            ),
            (Point3::new(-1.0, -1e-7, 0.0), Point3::new(0.0, 1e-7, 0.0)),
        ];
        let mut measured = 0;
        for (q0, q1) in cases {
            let y = segment(q0, q1);
            let Some(Some((s, _, _, zone, _))) = line_meeting(&x, &y, r) else {
                continue;
            };
            let d2 = q1 - q0;
            let dist = |t: f64| {
                let p = x.point_at(t);
                let w = ((p - q0).dot(&d2) / d2.dot(&d2)).clamp(0.0, 1.0);
                p.distance_to(&(q0 + d2 * w))
            };
            let searched = super::super::meetings_foot::zone_along(&x, s, (0.0, 2.0), dist, r);
            measured += 1;
            assert!(
                (zone - searched).abs() <= 1e-9 * searched.max(1e-6),
                "{q0:?}-{q1:?}: closed {zone:e}, searched {searched:e}"
            );
        }
        assert!(measured >= 5, "{measured} cases measured");
    }

    #[test]
    fn segments_apart_or_parallel_are_decided_or_left_to_the_search() {
        let x = segment(Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let apart = segment(Point3::new(0.5, -1.0, 1e-3), Point3::new(0.5, 1.0, 1e-3));
        assert_eq!(line_meeting(&x, &apart, 1e-6), Some(None));
        let parallel = segment(Point3::new(-1.0, 1e-7, 0.0), Point3::new(1.0, 1e-7, 0.0));
        assert_eq!(line_meeting(&x, &parallel, 1e-6), None);
        let parallel_apart = segment(Point3::new(3.0, 1e-3, 0.0), Point3::new(5.0, 1e-3, 0.0));
        assert_eq!(line_meeting(&x, &parallel_apart, 1e-6), Some(None));
        let mut curved = segment(Point3::new(0.5, -1.0, 0.0), Point3::new(0.5, 1.0, 0.0));
        curved.weights = vec![1.0, 2.0];
        assert_eq!(line_meeting(&x, &curved, 1e-6), None);
    }
}
