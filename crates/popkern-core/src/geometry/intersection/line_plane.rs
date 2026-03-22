//! Line-plane intersection algorithm.

use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Intersect a line `P + t*D` with an infinite plane then filter to the finite
/// patch defined by `u_range` and `v_range`.
///
/// The plane is described by an origin `O`, unit normal `N`, and two axis
/// vectors `u_axis` / `v_axis` that span the surface parameter space.
///
/// Returns:
/// - An empty `Vec` when there is no intersection within the given ranges.
/// - A single `Point` result when the line pierces the patch.
/// - A single `Overlap` result when the line lies inside the plane (within
///   `tolerance.linear`) — the full `t_range` is returned because the line is
///   coincident with the surface everywhere along its parameter domain.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn intersect_line_plane(
    line_origin: &Point3,
    line_dir: &Vec3,
    t_range: (f64, f64),
    plane_origin: &Point3,
    plane_normal: &Vec3,
    plane_u_axis: &Vec3,
    plane_v_axis: &Vec3,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    let denom = line_dir.dot(plane_normal);

    // Parallel (or anti-parallel) case: |D·N| < epsilon.
    // Use the same relative epsilon pattern as line_line: compare against the
    // product of the direction magnitudes so that the test is scale-invariant.
    let parallel_eps = f64::EPSILON * line_dir.length_squared() * plane_normal.length_squared();

    if denom * denom <= parallel_eps {
        // Measure the signed distance from the line origin to the plane.
        let offset = *plane_origin - *line_origin; // O - P
        let dist = offset.dot(plane_normal).abs();

        if dist < tolerance.linear {
            // Line lies in the plane — overlap over the full t_range.
            return vec![SurfaceCurveIntersection::Overlap { t_range }];
        }

        // Parallel, no intersection.
        return Vec::new();
    }

    // General case: solve t = (O - P)·N / (D·N).
    let offset = *plane_origin - *line_origin; // O - P
    let t = offset.dot(plane_normal) / denom;

    // Check t against the line's parameter range.
    if t < t_range.0 || t > t_range.1 {
        return Vec::new();
    }

    // Compute the intersection point Q = P + t*D.
    let q = *line_origin + *line_dir * t;

    // Project Q onto the plane's (u, v) axes.
    let q_minus_o = q - *plane_origin;
    let u = q_minus_o.dot(plane_u_axis);
    let v = q_minus_o.dot(plane_v_axis);

    // Filter to the patch's parameter ranges.
    if u < u_range.0 || u > u_range.1 || v < v_range.0 || v > v_range.1 {
        return Vec::new();
    }

    vec![SurfaceCurveIntersection::Point { point: q, t, u, v }]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xy_plane() -> (Point3, Vec3, Vec3, Vec3) {
        (
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        )
    }

    // Helper: intersect a line with the XY plane patch [(-5,5), (-5,5)].
    fn intersect_xy(
        lo: Point3,
        ld: Vec3,
        t_range: (f64, f64),
    ) -> Vec<SurfaceCurveIntersection> {
        let (po, pn, pu, pv) = xy_plane();
        intersect_line_plane(
            &lo, &ld, t_range,
            &po, &pn, &pu, &pv,
            (-5.0, 5.0), (-5.0, 5.0),
            &Tolerance::default(),
        )
    }

    // -----------------------------------------------------------------------
    // Test 1: perpendicular line — one point with correct u, v
    // -----------------------------------------------------------------------
    #[test]
    fn perpendicular_line_one_point() {
        // Line along +Z, starting at (2, 3, -1), going in direction (0, 0, 1).
        // Should hit the XY plane at t = 1.0, point (2, 3, 0), u=2, v=3.
        let results = intersect_xy(
            Point3::new(2.0, 3.0, -1.0),
            Vec3::new(0.0, 0.0, 1.0),
            (0.0, 10.0),
        );
        assert_eq!(results.len(), 1, "expected exactly one intersection");
        if let SurfaceCurveIntersection::Point { point, t, u, v } = &results[0] {
            assert!((*t - 1.0).abs() < 1e-10, "t should be 1.0, got {}", t);
            assert!((point.x - 2.0).abs() < 1e-10);
            assert!((point.y - 3.0).abs() < 1e-10);
            assert!(point.z.abs() < 1e-10);
            assert!((*u - 2.0).abs() < 1e-10, "u should be 2.0, got {}", u);
            assert!((*v - 3.0).abs() < 1e-10, "v should be 3.0, got {}", v);
        } else {
            panic!("expected Point intersection, got {:?}", results[0]);
        }
    }

    // -----------------------------------------------------------------------
    // Test 2: parallel miss — zero intersections
    // -----------------------------------------------------------------------
    #[test]
    fn parallel_miss_no_intersection() {
        // Line parallel to XY plane but 1 unit above it.
        let results = intersect_xy(
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            (0.0, 10.0),
        );
        assert!(results.is_empty(), "parallel miss should yield no intersections");
    }

    // -----------------------------------------------------------------------
    // Test 3: line in the plane — overlap
    // -----------------------------------------------------------------------
    #[test]
    fn line_in_plane_overlap() {
        // Line lying exactly in the XY plane (z=0 throughout).
        let results = intersect_xy(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            (0.0, 5.0),
        );
        assert_eq!(results.len(), 1, "line in plane should give one overlap result");
        if let SurfaceCurveIntersection::Overlap { t_range } = &results[0] {
            assert!((t_range.0 - 0.0).abs() < 1e-10);
            assert!((t_range.1 - 5.0).abs() < 1e-10);
        } else {
            panic!("expected Overlap, got {:?}", results[0]);
        }
    }

    // -----------------------------------------------------------------------
    // Test 4: oblique line — one intersection
    // -----------------------------------------------------------------------
    #[test]
    fn oblique_line_one_intersection() {
        // Line from (0, 0, -2) in direction (1, 1, 2): hits z=0 at t=1,
        // point (1, 1, 0), u=1, v=1.
        let results = intersect_xy(
            Point3::new(0.0, 0.0, -2.0),
            Vec3::new(1.0, 1.0, 2.0),
            (0.0, 10.0),
        );
        assert_eq!(results.len(), 1);
        if let SurfaceCurveIntersection::Point { point, t, u, v } = &results[0] {
            assert!((*t - 1.0).abs() < 1e-10, "t={}", t);
            assert!((point.x - 1.0).abs() < 1e-10);
            assert!((point.y - 1.0).abs() < 1e-10);
            assert!(point.z.abs() < 1e-10);
            assert!((*u - 1.0).abs() < 1e-10);
            assert!((*v - 1.0).abs() < 1e-10);
        } else {
            panic!("expected Point, got {:?}", results[0]);
        }
    }

    // -----------------------------------------------------------------------
    // Test 5: intersection outside u/v range — filtered to zero
    // -----------------------------------------------------------------------
    #[test]
    fn outside_uv_range_filtered() {
        // Line hits z=0 at (10, 0, 0), which is outside u_range (-5, 5).
        let (po, pn, pu, pv) = xy_plane();
        let results = intersect_line_plane(
            &Point3::new(10.0, 0.0, -1.0),
            &Vec3::new(0.0, 0.0, 1.0),
            (0.0, 10.0),
            &po, &pn, &pu, &pv,
            (-5.0, 5.0),
            (-5.0, 5.0),
            &Tolerance::default(),
        );
        assert!(results.is_empty(), "intersection outside uv range should be filtered");
    }

    // -----------------------------------------------------------------------
    // Test 6: intersection outside t_range — filtered to zero
    // -----------------------------------------------------------------------
    #[test]
    fn outside_t_range_filtered() {
        // Line hits z=0 at t=1 but t_range only covers (2, 10).
        let results = intersect_xy(
            Point3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, 0.0, 1.0),
            (2.0, 10.0),
        );
        assert!(results.is_empty(), "intersection outside t_range should be filtered");
    }
}
