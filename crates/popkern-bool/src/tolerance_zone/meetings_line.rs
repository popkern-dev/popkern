//! Meetings of two straight records in closed form (PK-1161).
//!
//! Most records of planar booleans are segments, degree-1 curves of two
//! control points. Two segments meet at most once, where they pass
//! closest, which the segment-segment closest-point solution gives
//! exactly; the sampled search (129 stations each, a nearest-station
//! pass and a golden-section refinement per local minimum) cost a box
//! union 0.17 ms of its 1.3, and a 50-face prism 0.26 of 2.4. Parallel
//! segments, which meet along a stretch or not at all, are left to the
//! sampled search and its overlap rule.

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
    // The distance from a point of ci to the segment cj.
    let dist = |x: f64| {
        let p = ci.point_at(x);
        let w = ((p - q0).dot(&d2) / l2).clamp(0.0, 1.0);
        p.distance_to(&(q0 + d2 * w))
    };
    let zone = super::meetings_foot::zone_along(ci, s, (a, b), dist, radius);
    // Within the radius a twentieth of the range either side is a
    // stretch the two share, not a meeting (the sampled search's rule).
    let step = super::meetings::OVERLAP_FRACTION * (b - a);
    let still_on = |x: f64| dist(x.clamp(a.min(b), a.max(b))) <= radius;
    if still_on(s - step) && still_on(s + step) {
        return Some(None);
    }
    Some(Some((s, t, gap, zone, sin_angle)))
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
