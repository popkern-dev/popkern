//! Line-line intersection algorithm.

use crate::geometry::intersection::CurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Intersect two lines in 3D space.
///
/// Lines are parameterized as `P1 + t*D1` and `P2 + s*D2`.
/// Handles parallel/coincident and skew cases.
#[allow(dead_code)]
pub(crate) fn intersect_line_line(
    p1: &Point3,
    d1: &Vec3,
    t1_range: (f64, f64),
    p2: &Point3,
    d2: &Vec3,
    t2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let cross = d1.cross(d2);
    let cross_len_sq = cross.length_squared();

    // Epsilon for detecting parallel directions: compare cross product magnitude
    // against the product of the input direction magnitudes.
    let parallel_eps = f64::EPSILON * d1.length_squared() * d2.length_squared();

    if cross_len_sq <= parallel_eps {
        // Lines are parallel (or one/both directions are zero-length — treat as coincident
        // check). Compute distance between the lines.
        let w = *p2 - *p1; // vector from P1 to P2
        let d1_len_sq = d1.length_squared();

        let dist = if d1_len_sq < f64::EPSILON {
            // Degenerate direction — fall back to point distance
            w.length()
        } else {
            // Project w perpendicular to d1
            let proj = w.dot(d1) / d1_len_sq;
            let perp = w - *d1 * proj;
            perp.length()
        };

        if dist >= tolerance.linear {
            // Parallel, non-coincident
            return Vec::new();
        }

        // Coincident lines — compute overlap in curve 1's parameter space.
        // Project each endpoint of line 2's range onto line 1's parameter axis.
        let d1_len_sq = d1.length_squared();
        if d1_len_sq < f64::EPSILON {
            return Vec::new();
        }

        let q_start = *p2 + *d2 * t2_range.0;
        let q_end = *p2 + *d2 * t2_range.1;

        let t_start_on_1 = (q_start - *p1).dot(d1) / d1_len_sq;
        let t_end_on_1 = (q_end - *p1).dot(d1) / d1_len_sq;

        // Interval of line 2's range mapped to line 1's parameter space
        let t2_lo = t_start_on_1.min(t_end_on_1);
        let t2_hi = t_start_on_1.max(t_end_on_1);

        // Intersect with t1_range
        let overlap_lo = t1_range.0.max(t2_lo);
        let overlap_hi = t1_range.1.min(t2_hi);

        if overlap_lo > overlap_hi {
            return Vec::new();
        }

        // Map overlap back to line 2's parameter space.
        // On line 2, the parameter corresponding to a point on line 1 at t1 is:
        //   s = (P1 + t1*D1 - P2) · D2 / (D2·D2)
        let d2_len_sq = d2.length_squared();
        let s_for_t = |t1: f64| -> f64 {
            if d2_len_sq < f64::EPSILON {
                t2_range.0
            } else {
                ((*p1 + *d1 * t1) - *p2).dot(d2) / d2_len_sq
            }
        };

        let s_lo = s_for_t(overlap_lo).min(s_for_t(overlap_hi));
        let s_hi = s_for_t(overlap_lo).max(s_for_t(overlap_hi));

        // Clamp to t2_range
        let s_lo = s_lo.max(t2_range.0);
        let s_hi = s_hi.min(t2_range.1);

        return vec![CurveIntersection::Overlap {
            t1_range: (overlap_lo, overlap_hi),
            t2_range: (s_lo, s_hi),
        }];
    }

    // Non-parallel: solve the 2×2 system for closest approach.
    //   D1·D1 * t - D1·D2 * s = D1·(P2-P1)
    //   D2·D1 * t - D2·D2 * s = D2·(P2-P1)
    let w = *p2 - *p1;
    let a = d1.dot(d1); // D1·D1
    let b = d1.dot(d2); // D1·D2
    let c = d2.dot(d2); // D2·D2
    let d = d1.dot(&w); // D1·(P2-P1)
    let e = d2.dot(&w); // D2·(P2-P1)

    // det = a*c - b*b  (same as |D1×D2|² for unit vectors, but works for
    // non-unit directions too)
    let det = a * c - b * b;

    if det.abs() < f64::EPSILON {
        // Shouldn't reach here given the cross-product check above, but be safe.
        return Vec::new();
    }

    let t = (d * c - e * b) / det;
    let s = (d * b - e * a) / det;

    // Check parameter ranges
    if t < t1_range.0 || t > t1_range.1 {
        return Vec::new();
    }
    if s < t2_range.0 || s > t2_range.1 {
        return Vec::new();
    }

    // Compute closest approach points and check distance
    let point_on_1 = *p1 + *d1 * t;
    let point_on_2 = *p2 + *d2 * s;
    let gap = (point_on_1 - point_on_2).length();

    if gap >= tolerance.linear {
        return Vec::new();
    }

    // Use midpoint as the intersection point
    let point = Point3::new(
        (point_on_1.x + point_on_2.x) * 0.5,
        (point_on_1.y + point_on_2.y) * 0.5,
        (point_on_1.z + point_on_2.z) * 0.5,
    );

    vec![CurveIntersection::Point {
        t1: t,
        t2: s,
        point,
    }]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perpendicular_lines_one_point() {
        let r = intersect_line_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 2.0),
            &Point3::new(0.5, -1.0, 0.0),
            &Vec3::new(0.0, 1.0, 0.0),
            (0.0, 2.0),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 1);
        if let CurveIntersection::Point { t1, t2, point } = &r[0] {
            assert!((*t1 - 0.5).abs() < 1e-10);
            assert!((*t2 - 1.0).abs() < 1e-10);
            assert!((point.x - 0.5).abs() < 1e-10);
        } else {
            panic!("Expected Point");
        }
    }

    #[test]
    fn parallel_non_coincident() {
        let r = intersect_line_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 1.0),
            &Point3::new(0.0, 1.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }

    #[test]
    fn coincident_lines_overlap() {
        let r = intersect_line_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 1.0),
            &Point3::new(0.5, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 1);
        assert!(matches!(&r[0], CurveIntersection::Overlap { .. }));
    }

    #[test]
    fn skew_within_tolerance() {
        let r = intersect_line_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 1.0),
            &Point3::new(0.5, 0.0, 1e-8),
            &Vec3::new(0.0, 1.0, 0.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 1);
    }

    #[test]
    fn skew_outside_tolerance() {
        let r = intersect_line_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 1.0),
            &Point3::new(0.5, 0.0, 1.0),
            &Vec3::new(0.0, 1.0, 0.0),
            (0.0, 1.0),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }

    #[test]
    fn intersection_outside_range() {
        let r = intersect_line_line(
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 0.3),
            &Point3::new(0.5, -1.0, 0.0),
            &Vec3::new(0.0, 1.0, 0.0),
            (0.0, 2.0),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }
}
