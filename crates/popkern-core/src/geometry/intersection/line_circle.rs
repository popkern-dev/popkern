//! Line-circle intersection algorithm.

use crate::geometry::analytic_curves::circle_basis;
use crate::geometry::intersection::CurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Intersect a line with a circle in 3D space.
/// Circle is defined by center, normal (to its plane), and radius.
/// angle_range is the valid angle range for the circle/arc.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_circle(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_line_range: (f64, f64),
    circle_center: &Point3,
    circle_normal: &Vec3,
    circle_radius: f64,
    angle_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let d_dot_n = line_dir.dot(circle_normal);

    if d_dot_n.abs() < tolerance.linear() {
        // Line is parallel to the circle's plane.
        // Check if line lies in the plane.
        let cp = *line_origin - *circle_center;
        let dist_to_plane = cp.dot(circle_normal).abs();
        if dist_to_plane > tolerance.linear() {
            // Parallel but not in plane — no intersection.
            return Vec::new();
        }
        // Line lies in the plane of the circle — solve in 2D.
        in_plane_intersections(
            line_origin,
            line_dir,
            t_line_range,
            circle_center,
            circle_normal,
            circle_radius,
            angle_range,
            tolerance,
        )
    } else {
        // Line pierces the plane at parameter t_pierce.
        let cp = *circle_center - *line_origin;
        let t_pierce = cp.dot(circle_normal) / d_dot_n;

        if t_pierce < t_line_range.0 - tolerance.linear()
            || t_pierce > t_line_range.1 + tolerance.linear()
        {
            return Vec::new();
        }

        let pierce_pt = *line_origin + *line_dir * t_pierce;
        let dist_from_center = (pierce_pt - *circle_center).length();

        if (dist_from_center - circle_radius).abs() > tolerance.linear() {
            return Vec::new();
        }

        // Compute circle angle for the pierce point.
        let (u, v) = circle_basis(*circle_normal);
        let local = pierce_pt - *circle_center;
        let local_u = local.dot(&u);
        let local_v = local.dot(&v);
        let angle = local_v.atan2(local_u);
        let angle = normalise_angle(angle, angle_range);

        if !angle_in_range(angle, angle_range, tolerance) {
            return Vec::new();
        }

        vec![CurveIntersection::Point {
            point: pierce_pt,
            t1: t_pierce,
            t2: angle,
        }]
    }
}

/// Solve line-circle intersection when the line lies in the circle's plane.
#[allow(clippy::too_many_arguments)]
fn in_plane_intersections(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_line_range: (f64, f64),
    circle_center: &Point3,
    circle_normal: &Vec3,
    circle_radius: f64,
    angle_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    // |line_origin + t*line_dir - circle_center|^2 = circle_radius^2
    // Let w = line_origin - circle_center
    // |w + t*D|^2 = R^2
    // t^2 * |D|^2 + 2t*(w·D) + |w|^2 - R^2 = 0
    let w = *line_origin - *circle_center;
    let a = line_dir.dot(line_dir);
    let b = 2.0 * w.dot(line_dir);
    let c = w.dot(&w) - circle_radius * circle_radius;

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < -(tolerance.linear() * tolerance.linear() * a) {
        // No real intersection.
        return Vec::new();
    }

    let (u, v) = circle_basis(*circle_normal);

    let mut results = Vec::new();

    if discriminant <= tolerance.linear() * tolerance.linear() * a {
        // Tangent — one solution.
        let t = -b / (2.0 * a);
        if let Some(isect) = make_intersection(
            line_origin,
            line_dir,
            t,
            t_line_range,
            circle_center,
            &u,
            &v,
            angle_range,
            tolerance,
        ) {
            results.push(isect);
        }
    } else {
        // Two solutions.
        let sqrt_disc = discriminant.sqrt();
        let t1 = (-b - sqrt_disc) / (2.0 * a);
        let t2 = (-b + sqrt_disc) / (2.0 * a);
        if let Some(isect) = make_intersection(
            line_origin,
            line_dir,
            t1,
            t_line_range,
            circle_center,
            &u,
            &v,
            angle_range,
            tolerance,
        ) {
            results.push(isect);
        }
        if let Some(isect) = make_intersection(
            line_origin,
            line_dir,
            t2,
            t_line_range,
            circle_center,
            &u,
            &v,
            angle_range,
            tolerance,
        ) {
            results.push(isect);
        }
    }

    results
}

/// Build a CurveIntersection for a candidate parameter t, filtering by t_line_range and angle_range.
#[allow(clippy::too_many_arguments)]
fn make_intersection(
    line_origin: &Point3,
    line_dir: &Vec3,
    t: f64,
    t_line_range: (f64, f64),
    circle_center: &Point3,
    u: &Vec3,
    v: &Vec3,
    angle_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<CurveIntersection> {
    if t < t_line_range.0 - tolerance.linear() || t > t_line_range.1 + tolerance.linear() {
        return None;
    }

    let pt = *line_origin + *line_dir * t;
    let local = pt - *circle_center;
    let local_u = local.dot(u);
    let local_v = local.dot(v);
    let angle = local_v.atan2(local_u);
    let angle = normalise_angle(angle, angle_range);

    if !angle_in_range(angle, angle_range, tolerance) {
        return None;
    }

    Some(CurveIntersection::Point {
        point: pt,
        t1: t,
        t2: angle,
    })
}

/// Shift angle (from atan2) into the range [angle_range.0, angle_range.0 + TAU) so that
/// a full circle [0, TAU) accepts all angles.
fn normalise_angle(angle: f64, angle_range: (f64, f64)) -> f64 {
    let tau = std::f64::consts::TAU;
    let mut a = angle;
    // Bring a to within [angle_range.0, angle_range.0 + TAU)
    while a < angle_range.0 - f64::EPSILON {
        a += tau;
    }
    while a >= angle_range.0 + tau - f64::EPSILON {
        a -= tau;
    }
    a
}

/// Returns true if `angle` lies within [angle_range.0, angle_range.1] (inclusive within tolerance).
fn angle_in_range(angle: f64, angle_range: (f64, f64), tolerance: &Tolerance) -> bool {
    angle >= angle_range.0 - tolerance.angular() && angle <= angle_range.1 + tolerance.angular()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn secant_two_points() {
        let r = intersect_line_circle(
            &Point3::new(-2.0, 0.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-2.0, 4.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn tangent_one_point() {
        let r = intersect_line_circle(
            &Point3::new(-1.0, 1.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 2.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 1);
    }

    #[test]
    fn miss_no_intersection() {
        let r = intersect_line_circle(
            &Point3::new(-1.0, 3.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (0.0, 2.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }

    #[test]
    fn line_pierces_plane_at_circle() {
        // Line along Z axis, circle in XY plane at radius 0 from origin — pierces at center, not on circle
        let r = intersect_line_circle(
            &Point3::new(0.0, 0.0, -1.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (0.0, 2.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert!(r.is_empty()); // pierces center, not circumference
    }

    #[test]
    fn line_pierces_plane_on_circle() {
        // Line along Z axis at x=1, circle of radius 1 centered at origin in XY plane
        let r = intersect_line_circle(
            &Point3::new(1.0, 0.0, -1.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (0.0, 2.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert_eq!(r.len(), 1);
    }

    #[test]
    fn line_parallel_not_in_plane() {
        // Line parallel to circle plane but not in it
        let r = intersect_line_circle(
            &Point3::new(0.0, 0.0, 5.0), // 5 above plane
            &Vec3::new(1.0, 0.0, 0.0),   // horizontal
            (0.0, 2.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }

    #[test]
    fn line_pierces_outside_t_range() {
        // Line pierces the plane at t_pierce = 1.0, but t range is (0, 0.5)
        let r = intersect_line_circle(
            &Point3::new(1.0, 0.0, -1.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (0.0, 0.5), // pierce at t=1 is outside this range
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }

    #[test]
    fn line_pierces_but_angle_outside_range() {
        // Line pierces at (0,1,0) which is angle=PI/2, but arc is [PI, TAU]
        let r = intersect_line_circle(
            &Point3::new(0.0, 1.0, -1.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (0.0, 2.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (std::f64::consts::PI, std::f64::consts::TAU - 0.01), // angle PI/2 is outside
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }

    #[test]
    fn normalise_angle_wraps() {
        // Test that normalise_angle wraps correctly
        let tau = std::f64::consts::TAU;
        let a = normalise_angle(-0.5, (0.0, tau));
        assert!(a >= 0.0 - f64::EPSILON);
        assert!(a < tau + f64::EPSILON);
    }

    #[test]
    fn in_plane_tangent_outside_t_range() {
        // Line tangent to circle but t parameter is outside the line range
        let r = intersect_line_circle(
            &Point3::new(-1.0, 1.0, 0.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (5.0, 10.0), // very restricted range far from intersection
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            1.0,
            (0.0, std::f64::consts::TAU),
            &Tolerance::default(),
        );
        assert!(r.is_empty());
    }
}
