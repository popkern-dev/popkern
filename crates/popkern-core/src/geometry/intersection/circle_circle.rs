//! Circle-circle intersection algorithm.

#![allow(dead_code)]

use crate::geometry::analytic_curves::circle_basis;
use crate::geometry::intersection::CurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Compute the angle parameter for `point` on a circle defined by
/// `center`, basis vectors `(u, v)`, and `radius`.
fn angle_on_circle(point: &Point3, center: &Point3, u: &Vec3, v: &Vec3) -> f64 {
    let d = *point - *center;
    let x = d.dot(u);
    let y = d.dot(v);
    let angle = y.atan2(x);
    if angle < 0.0 {
        angle + std::f64::consts::TAU
    } else {
        angle
    }
}

/// Returns true if `angle` is within the range `[lo, hi]` (mod TAU),
/// with a small tolerance to handle boundary cases.
fn angle_in_range(angle: f64, range: (f64, f64), tolerance: f64) -> bool {
    let tau = std::f64::consts::TAU;
    let (lo, hi) = range;

    // Full circle case
    if (hi - lo - tau).abs() < tolerance {
        return true;
    }

    // Normalize angle into [0, TAU)
    let mut a = angle % tau;
    if a < 0.0 {
        a += tau;
    }
    let mut lo_n = lo % tau;
    if lo_n < 0.0 {
        lo_n += tau;
    }
    let mut hi_n = hi % tau;
    if hi_n < 0.0 {
        hi_n += tau;
    }

    if lo_n <= hi_n {
        a >= lo_n - tolerance && a <= hi_n + tolerance
    } else {
        // Range wraps around TAU
        a >= lo_n - tolerance || a <= hi_n + tolerance
    }
}

/// Intersect two circles in 3D space.
///
/// Each circle is defined by center, normal (plane), radius, and an angle range.
/// Returns intersection points and/or overlap regions.
#[allow(clippy::too_many_arguments)]
pub(crate) fn intersect_circle_circle(
    c1_center: &Point3,
    c1_normal: &Vec3,
    c1_radius: f64,
    c1_angle_range: (f64, f64),
    c2_center: &Point3,
    c2_normal: &Vec3,
    c2_radius: f64,
    c2_angle_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let tol = tolerance.linear();

    // Normalize normals
    let n1 = match c1_normal.normalized() {
        Some(n) => n,
        None => return Vec::new(),
    };
    let n2 = match c2_normal.normalized() {
        Some(n) => n,
        None => return Vec::new(),
    };

    // Check if normals are parallel (coplanar case)
    let cross = n1.cross(&n2);
    let cross_len = cross.length();

    if cross_len < tolerance.angular() {
        // Normals are parallel — check if circles are coplanar
        // Distance from c2_center to the plane of circle 1
        let delta = *c2_center - *c1_center;
        let plane_dist = delta.dot(&n1).abs();

        if plane_dist > tol {
            // Circles are in parallel but offset planes — no intersection
            return Vec::new();
        }

        intersect_coplanar(
            c1_center,
            &n1,
            c1_radius,
            c1_angle_range,
            c2_center,
            &n2,
            c2_radius,
            c2_angle_range,
            tolerance,
        )
    } else {
        intersect_non_coplanar(
            c1_center,
            &n1,
            c1_radius,
            c1_angle_range,
            c2_center,
            &n2,
            c2_radius,
            c2_angle_range,
            tolerance,
        )
    }
}

/// Coplanar circle-circle intersection.
#[allow(clippy::too_many_arguments)]
fn intersect_coplanar(
    c1_center: &Point3,
    c1_normal: &Vec3,
    c1_radius: f64,
    c1_angle_range: (f64, f64),
    c2_center: &Point3,
    _c2_normal: &Vec3,
    c2_radius: f64,
    c2_angle_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let tol = tolerance.linear();
    let delta = *c2_center - *c1_center;
    let d = delta.length();

    let r1 = c1_radius;
    let r2 = c2_radius;

    // Case: identical (coincident) circles
    if d < tol && (r1 - r2).abs() < tol {
        // Overlap — the circles are identical
        // Compute the overlap of the angle ranges
        let t1_start = c1_angle_range.0;
        let t1_end = c1_angle_range.1;
        let t2_start = c2_angle_range.0;
        let t2_end = c2_angle_range.1;
        return vec![CurveIntersection::Overlap {
            t1_range: (t1_start, t1_end),
            t2_range: (t2_start, t2_end),
        }];
    }

    // Case: one circle inside the other (concentric with different radii)
    if d + tol < (r1 - r2).abs() {
        return Vec::new();
    }

    // Case: disjoint — too far apart
    if d > r1 + r2 + tol {
        return Vec::new();
    }

    let (u1, v1) = circle_basis(*c1_normal);
    let (u2, v2) = circle_basis(*c1_normal); // Use same normal for c2 since coplanar

    // Tangent or two-point intersection
    let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let h_sq = r1 * r1 - a * a;

    // Direction from c1 to c2
    let dir = match delta.normalized() {
        Some(dir) => dir,
        None => return Vec::new(), // Centers coincide but radii differ — handled above
    };

    let p_mid = *c1_center + dir * a;

    if h_sq.abs() < tol * tol || h_sq < 0.0 {
        // Tangent — single intersection point
        let t1 = angle_on_circle(&p_mid, c1_center, &u1, &v1);
        let t2 = angle_on_circle(&p_mid, c2_center, &u2, &v2);

        if angle_in_range(t1, c1_angle_range, tolerance.parametric_base())
            && angle_in_range(t2, c2_angle_range, tolerance.parametric_base())
        {
            return vec![CurveIntersection::Point {
                point: p_mid,
                t1,
                t2,
            }];
        }
        return Vec::new();
    }

    // Two intersection points
    let h = h_sq.sqrt();
    // Perpendicular direction in the plane
    let perp = c1_normal.cross(&dir);
    let perp = match perp.normalized() {
        Some(p) => p,
        None => return Vec::new(),
    };

    let offset = perp * h;
    let p1 = p_mid + offset;
    let p2 = p_mid - offset;

    let mut results = Vec::new();

    for pt in &[p1, p2] {
        let t1 = angle_on_circle(pt, c1_center, &u1, &v1);
        let t2 = angle_on_circle(pt, c2_center, &u2, &v2);

        if angle_in_range(t1, c1_angle_range, tolerance.parametric_base())
            && angle_in_range(t2, c2_angle_range, tolerance.parametric_base())
        {
            results.push(CurveIntersection::Point { point: *pt, t1, t2 });
        }
    }

    results
}

/// Non-coplanar circle-circle intersection.
///
/// Intersects the two planes to get a line, then finds where that line
/// intersects each circle, and matches coincident points.
#[allow(clippy::too_many_arguments)]
fn intersect_non_coplanar(
    c1_center: &Point3,
    c1_normal: &Vec3,
    c1_radius: f64,
    c1_angle_range: (f64, f64),
    c2_center: &Point3,
    c2_normal: &Vec3,
    c2_radius: f64,
    c2_angle_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let tol = tolerance.linear();

    // Find the line of intersection of the two planes.
    // Plane 1: n1 . (P - C1) = 0
    // Plane 2: n2 . (P - C2) = 0
    let line_dir = c1_normal.cross(c2_normal);
    let line_dir_n = match line_dir.normalized() {
        Some(d) => d,
        None => return Vec::new(), // Parallel planes (shouldn't reach here)
    };

    // Find a point on the line of intersection.
    // We solve: n1 . (P - C1) = 0, n2 . (P - C2) = 0
    // Let P = C1 + s * n1 + t * n2 + u * line_dir
    // We can set u = 0 and solve for s, t.
    let d1 = 0.0; // n1 . (C1 - C1) = 0
    let d2 = c2_normal.dot(&(*c1_center - *c2_center)); // n2 . (C1 - C2)

    let n1n1 = c1_normal.dot(c1_normal);
    let n1n2 = c1_normal.dot(c2_normal);
    let n2n2 = c2_normal.dot(c2_normal);

    let denom = n1n1 * n2n2 - n1n2 * n1n2;
    if denom.abs() < 1e-30 {
        return Vec::new();
    }

    let s = (d1 * n2n2 - d2 * n1n2) / denom;
    let t = (d2 * n1n1 - d1 * n1n2) / denom;
    // Negate t because d2 = n2.(C1-C2), and we want P on both planes
    let line_point = *c1_center + *c1_normal * s + *c2_normal * (-t);

    // Now intersect this line with each circle
    let pts1 = intersect_line_with_circle(
        &line_point,
        &line_dir_n,
        c1_center,
        c1_normal,
        c1_radius,
        tolerance,
    );
    let pts2 = intersect_line_with_circle(
        &line_point,
        &line_dir_n,
        c2_center,
        c2_normal,
        c2_radius,
        tolerance,
    );

    let (u1, v1) = circle_basis(*c1_normal);
    let (u2, v2) = circle_basis(*c2_normal);

    let mut results = Vec::new();

    for p1 in &pts1 {
        for p2 in &pts2 {
            if p1.distance_to(p2) < tol {
                // Use the midpoint for better accuracy
                let mid = Point3::new(
                    (p1.x + p2.x) * 0.5,
                    (p1.y + p2.y) * 0.5,
                    (p1.z + p2.z) * 0.5,
                );

                let t1 = angle_on_circle(&mid, c1_center, &u1, &v1);
                let t2 = angle_on_circle(&mid, c2_center, &u2, &v2);

                if angle_in_range(t1, c1_angle_range, tolerance.parametric_base())
                    && angle_in_range(t2, c2_angle_range, tolerance.parametric_base())
                {
                    results.push(CurveIntersection::Point { point: mid, t1, t2 });
                }
            }
        }
    }

    results
}

/// Intersect a line (point + direction) with a circle in 3D.
///
/// The circle is defined by center, normal, and radius.
/// Returns 0, 1, or 2 intersection points.
fn intersect_line_with_circle(
    line_point: &Point3,
    line_dir: &Vec3,
    circle_center: &Point3,
    circle_normal: &Vec3,
    circle_radius: f64,
    tolerance: &Tolerance,
) -> Vec<Point3> {
    // Project everything into the plane of the circle.
    // The line must lie in (or near) the circle's plane for there to be intersections.

    let (u, v) = circle_basis(*circle_normal);

    // Express the line in 2D coordinates in the circle's plane
    let delta = *line_point - *circle_center;

    // Project line origin and direction into the plane
    let ox = delta.dot(&u);
    let oy = delta.dot(&v);
    let dx = line_dir.dot(&u);
    let dy = line_dir.dot(&v);

    // Solve: (ox + t*dx)^2 + (oy + t*dy)^2 = R^2
    let a = dx * dx + dy * dy;
    let b = 2.0 * (ox * dx + oy * dy);
    let c = ox * ox + oy * oy - circle_radius * circle_radius;

    if a.abs() < 1e-30 {
        return Vec::new();
    }

    let disc = b * b - 4.0 * a * c;

    if disc < -tolerance.solver() {
        return Vec::new();
    }

    let disc_sqrt = if disc < 0.0 { 0.0 } else { disc.sqrt() };

    let mut points = Vec::new();

    let t1 = (-b + disc_sqrt) / (2.0 * a);
    let t2 = (-b - disc_sqrt) / (2.0 * a);

    if (t1 - t2).abs() < tolerance.solver() {
        // Single tangent point
        let p = *line_point + *line_dir * t1;
        points.push(p);
    } else {
        let p1 = *line_point + *line_dir * t1;
        let p2 = *line_point + *line_dir * t2;
        points.push(p1);
        points.push(p2);
    }

    points
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    fn cc(
        c1: Point3,
        n1: Vec3,
        r1: f64,
        a1: (f64, f64),
        c2: Point3,
        n2: Vec3,
        r2: f64,
        a2: (f64, f64),
    ) -> Vec<CurveIntersection> {
        intersect_circle_circle(&c1, &n1, r1, a1, &c2, &n2, r2, a2, &Tolerance::default())
    }
    fn zn() -> Vec3 {
        Vec3::new(0.0, 0.0, 1.0)
    }
    fn o() -> Point3 {
        Point3::new(0.0, 0.0, 0.0)
    }
    fn full() -> (f64, f64) {
        (0.0, TAU)
    }

    #[test]
    fn coplanar_cases() {
        // Two points
        let r = cc(
            o(),
            zn(),
            1.0,
            full(),
            Point3::new(1.0, 0.0, 0.0),
            zn(),
            1.0,
            full(),
        );
        assert_eq!(r.len(), 2);
        // Tangent
        assert_eq!(
            cc(
                o(),
                zn(),
                1.0,
                full(),
                Point3::new(2.0, 0.0, 0.0),
                zn(),
                1.0,
                full()
            )
            .len(),
            1
        );
        // Disjoint
        assert!(
            cc(
                o(),
                zn(),
                1.0,
                full(),
                Point3::new(5.0, 0.0, 0.0),
                zn(),
                1.0,
                full()
            )
            .is_empty()
        );
        // Identical => overlap
        let r = cc(o(), zn(), 1.0, full(), o(), zn(), 1.0, full());
        assert!(matches!(&r[0], CurveIntersection::Overlap { .. }));
        // Concentric different radii
        assert!(cc(o(), zn(), 1.0, full(), o(), zn(), 2.0, full()).is_empty());
    }

    #[test]
    fn non_coplanar_cases() {
        let xn = Vec3::new(1.0, 0.0, 0.0);
        assert_eq!(cc(o(), zn(), 1.0, full(), o(), xn, 1.0, full()).len(), 2);
        // Disjoint non-coplanar
        assert!(
            cc(
                o(),
                zn(),
                1.0,
                full(),
                Point3::new(10.0, 10.0, 10.0),
                xn,
                1.0,
                full()
            )
            .is_empty()
        );
        // Restricted arc range
        assert!(cc(o(), zn(), 1.0, (0.0, 1.0), o(), xn, 1.0, full()).len() <= 2);
    }

    #[test]
    fn degenerate_and_offset_cases() {
        let z = Vec3::new(0.0, 0.0, 0.0);
        let pz = Point3::new(0.0, 0.0, 5.0);
        assert!(cc(o(), z, 1.0, full(), o(), zn(), 1.0, full()).is_empty());
        assert!(cc(o(), zn(), 1.0, full(), o(), z, 1.0, full()).is_empty());
        assert!(cc(o(), zn(), 1.0, full(), pz, zn(), 1.0, full()).is_empty());
    }

    #[test]
    fn angle_in_range_cases() {
        assert!(angle_in_range(-0.1, (0.0, TAU), 0.2));
        assert!(angle_in_range(0.1, (-0.1, 1.0), 0.01));
        // Wrapping range
        assert!(angle_in_range(5.5, (5.0, 1.0), 0.01));
        assert!(angle_in_range(0.5, (5.0, 1.0), 0.01));
        assert!(!angle_in_range(3.0, (5.0, 1.0), 0.01));
    }
}
