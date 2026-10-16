//! Helper functions for center curve computation.
//!
//! Provides surface projection, footprint clipping, NURBS fitting,
//! and parameter range clamping used by the [`super::center_curve`] module.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};

/// Clip footprint points to the face's UV parameter bounds.
///
/// Projects each 3D footprint point into the surface UV domain and retains
/// only those points whose UV coordinates fall within (or very close to)
/// `u_range` and `v_range`. This prevents trim curves from extending
/// beyond face boundaries, which would cause `split_face` to receive
/// mismatched anchor vertices.
///
/// When a consecutive pair of points straddles the boundary (one inside,
/// one outside), a boundary intersection point is inserted so that the
/// trimmed footprint ends precisely at the face edge. After clipping,
/// endpoints are snapped to the nearest face boundary if they lie within
/// a small tolerance.
///
/// If clipping would remove all points, the original footprint is returned
/// unchanged to avoid degenerate results.
pub fn clip_footprint_to_bounds(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    footprint: &[Point3],
) -> Vec<Point3> {
    if footprint.len() < 2 {
        return footprint.to_vec();
    }

    // Boundary-inclusion tolerance (generous for offset-surface noise).
    let eps = 1e-4;

    let in_bounds = |u: f64, v: f64| -> bool {
        u >= u_range.0 - eps && u <= u_range.1 + eps && v >= v_range.0 - eps && v <= v_range.1 + eps
    };

    // Branch-aware projection for EVERY periodic direction: the raw
    // projection lands angles in arbitrary branches and the bounds test
    // then DROPS valid samples. Measured twice: boss/plate seam lost
    // 2/20 on a cylinder (fixed then cylinder-only), cone base lost
    // 18/20 against its (-2π, 0) range (PK-264, "2 vs 20").
    let uvs: Vec<Option<(f64, f64)>> = footprint
        .iter()
        .map(|pt| project_to_face_branch(surface, u_range, v_range, pt))
        .collect();

    let mut clipped: Vec<Point3> = Vec::with_capacity(footprint.len());

    for i in 0..footprint.len() {
        let cur_in = uvs[i].is_none_or(|(u, v)| in_bounds(u, v));

        // When transitioning from outside to inside, insert a boundary
        // interpolation point.
        if i > 0 {
            let prev_in = uvs[i - 1].is_none_or(|(u, v)| in_bounds(u, v));
            if cur_in != prev_in
                && let Some(boundary_pt) = interpolate_boundary_point(
                    &footprint[i - 1],
                    &footprint[i],
                    surface,
                    u_range,
                    v_range,
                )
            {
                clipped.push(boundary_pt);
            }
        }

        if cur_in {
            clipped.push(footprint[i]);
        }
    }

    // Never return an empty footprint — fall back to unclipped.
    if clipped.len() < 2 {
        return footprint.to_vec();
    }

    // Snap endpoints to face boundary edges if they are close.
    // SAFETY: clipped.len() >= 2 guaranteed by the check above.
    let n = clipped.len();
    snap_endpoint_to_boundary(&mut clipped[0], surface, u_range, v_range);
    snap_endpoint_to_boundary(&mut clipped[n - 1], surface, u_range, v_range);

    clipped
}

/// Project a 3D point into a face's UV chart, wrapping every periodic
/// parameter into the branch nearest the face range's midpoint (the
/// raw projection is branch-arbitrary on periodic surfaces).
pub fn project_to_face_branch(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    point: &Point3,
) -> Option<(f64, f64)> {
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, point)?;
    let wrap = |x: f64, range: (f64, f64), period: Option<f64>| match period {
        Some(p) => x - p * ((x - 0.5 * (range.0 + range.1)) / p).round(),
        None => x,
    };
    let (u_period, v_period) = surface.periodicity();
    Some((wrap(u, u_range, u_period), wrap(v, v_range, v_period)))
}

/// Interpolate a point on the face boundary between an inside and outside
/// footprint point using binary search in 3D.
fn interpolate_boundary_point(
    pt_a: &Point3,
    pt_b: &Point3,
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Option<Point3> {
    let eps = 1e-4;
    let in_bounds = |u: f64, v: f64| -> bool {
        u >= u_range.0 - eps && u <= u_range.1 + eps && v >= v_range.0 - eps && v <= v_range.1 + eps
    };

    // Binary search for the boundary crossing.
    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    let a_in = project_to_face_branch(surface, u_range, v_range, pt_a)
        .is_none_or(|(u, v)| in_bounds(u, v));

    for _ in 0..20 {
        let mid = (lo + hi) * 0.5;
        let pt = Point3::new(
            pt_a.x + (pt_b.x - pt_a.x) * mid,
            pt_a.y + (pt_b.y - pt_a.y) * mid,
            pt_a.z + (pt_b.z - pt_a.z) * mid,
        );
        let mid_in = project_to_face_branch(surface, u_range, v_range, &pt)
            .is_none_or(|(u, v)| in_bounds(u, v));
        if mid_in == a_in {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let t = (lo + hi) * 0.5;
    Some(Point3::new(
        pt_a.x + (pt_b.x - pt_a.x) * t,
        pt_a.y + (pt_b.y - pt_a.y) * t,
        pt_a.z + (pt_b.z - pt_a.z) * t,
    ))
}

/// Snap a footprint endpoint to the nearest face boundary if its UV
/// projection is within `snap_tol` of a boundary edge.
fn snap_endpoint_to_boundary(
    pt: &mut Point3,
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) {
    use popkern_core::geometry::surface::SurfaceEval;

    let snap_tol = 1e-3;
    let (u, v) = match project_to_face_branch(surface, u_range, v_range, pt) {
        Some(uv) => uv,
        None => return,
    };

    // Snap u to boundary if close.
    let u_snapped = if (u - u_range.0).abs() < snap_tol {
        u_range.0
    } else if (u - u_range.1).abs() < snap_tol {
        u_range.1
    } else {
        u
    };

    // Snap v to boundary if close.
    let v_snapped = if (v - v_range.0).abs() < snap_tol {
        v_range.0
    } else if (v - v_range.1).abs() < snap_tol {
        v_range.1
    } else {
        v
    };

    if (u_snapped - u).abs() > 1e-15 || (v_snapped - v).abs() > 1e-15 {
        *pt = surface.point_at(u_snapped, v_snapped);
    }
}

/// Clamp a potentially infinite parameter range to a fallback finite range.
pub fn clamp_finite_range(range: (f64, f64), fallback: (f64, f64)) -> (f64, f64) {
    let lo = if range.0.is_finite() {
        range.0
    } else {
        fallback.0
    };
    let hi = if range.1.is_finite() {
        range.1
    } else {
        fallback.1
    };
    (lo, hi)
}

/// Build a degree-1 NURBS (polyline) through the given points.
///
/// Uses uniform chord-length parameterization. This is the simplest
/// possible interpolation; higher-degree fitting can be added later.
pub fn fit_nurbs_through_points(points: &[Point3]) -> Result<NurbsCurve, String> {
    let n = points.len();
    if n < 2 {
        return Err("need at least 2 points to fit a NURBS curve".into());
    }

    // Chord-length parameterization
    let mut params = vec![0.0f64; n];
    for i in 1..n {
        params[i] = params[i - 1] + points[i].distance_to(&points[i - 1]);
    }
    let total = params[n - 1];
    if total < 1e-30 {
        return Err("degenerate: all points coincide".into());
    }
    for p in &mut params {
        *p /= total;
    }

    // Degree-1 NURBS: knots = clamped uniform, control_points = the points themselves.
    let degree = 1;
    let mut knots = Vec::with_capacity(n + degree + 1);
    knots.push(0.0);
    knots.push(0.0);
    for param in params.iter().take(n - 1).skip(1) {
        knots.push(*param);
    }
    knots.push(1.0);
    knots.push(1.0);

    Ok(NurbsCurve {
        degree,
        knots,
        control_points: points.to_vec(),
        weights: vec![1.0; n],
    })
}

pub use crate::center_curve_projection::project_onto_surface;

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::curve::CurveEval;

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn fit_nurbs_through_points_simple() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ];
        let curve = fit_nurbs_through_points(&pts).unwrap();
        assert_eq!(curve.degree, 1);
        assert_eq!(curve.control_points.len(), 3);

        // Endpoints should match
        let (t0, t1) = curve.parameter_range();
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        assert!(p0.distance_to(&pts[0]) < 1e-10);
        assert!(p1.distance_to(&pts[2]) < 1e-10);
    }

    #[test]
    fn fit_nurbs_too_few_points() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0)];
        assert!(fit_nurbs_through_points(&pts).is_err());
    }

    #[test]
    fn project_onto_plane() {
        let plane = xy_plane();
        let pt = Point3::new(1.0, 2.0, 5.0);
        let proj = project_onto_surface(&plane, (-10.0, 10.0), (-10.0, 10.0), &pt).unwrap();
        assert!((proj.z).abs() < 1e-10, "should project to z=0");
        assert!((proj.x - 1.0).abs() < 1e-10);
        assert!((proj.y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn project_onto_sphere() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: popkern_core::geometry::Vec3::Z,
        };
        let pt = Point3::new(3.0, 0.0, 0.0);
        let proj = project_onto_surface(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::PI),
            &pt,
        )
        .unwrap();
        // Should be on the sphere, closest to (3,0,0) -> (1,0,0)
        let r = (proj.x * proj.x + proj.y * proj.y + proj.z * proj.z).sqrt();
        assert!(
            (r - 1.0).abs() < 0.2,
            "projected point should be on sphere, r={}",
            r
        );
    }

    #[test]
    fn project_onto_cylinder() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let pt = Point3::new(5.0, 0.0, 1.0);
        let proj =
            project_onto_surface(&cyl, (0.0, std::f64::consts::TAU), (-5.0, 5.0), &pt).unwrap();
        // Should be on the cylinder surface (r=2 from z-axis).
        let r = (proj.x * proj.x + proj.y * proj.y).sqrt();
        assert!(
            (r - 2.0).abs() < 0.5,
            "projected point should be near cylinder, r={}",
            r
        );
    }

    #[test]
    fn clip_footprint_all_inside() {
        let plane = xy_plane();
        let footprint = vec![
            Point3::new(0.2, 0.3, 0.0),
            Point3::new(0.5, 0.5, 0.0),
            Point3::new(0.8, 0.7, 0.0),
        ];
        let clipped = clip_footprint_to_bounds(&plane, (0.0, 1.0), (0.0, 1.0), &footprint);
        assert_eq!(clipped.len(), 3, "all points inside, should keep all");
    }

    #[test]
    fn clip_footprint_straddle_boundary() {
        let plane = xy_plane();
        // Points go from inside to outside the u=1.0 boundary.
        let footprint = vec![
            Point3::new(0.5, 0.5, 0.0),
            Point3::new(0.9, 0.5, 0.0),
            Point3::new(1.5, 0.5, 0.0), // outside
            Point3::new(2.0, 0.5, 0.0), // outside
        ];
        let clipped = clip_footprint_to_bounds(&plane, (0.0, 1.0), (0.0, 1.0), &footprint);
        // Should retain interior points + interpolated boundary point.
        assert!(
            clipped.len() >= 2,
            "should have at least 2 points after clipping, got {}",
            clipped.len()
        );
        // The last retained point should be near x=1.0.
        let last = clipped.last().unwrap();
        assert!(
            last.x <= 1.01,
            "last point should be at or near boundary, x={}",
            last.x
        );
    }

    #[test]
    fn clip_footprint_single_point_unchanged() {
        let plane = xy_plane();
        let clipped = clip_footprint_to_bounds(
            &plane,
            (0.0, 1.0),
            (0.0, 1.0),
            &[Point3::new(0.5, 0.5, 0.0)],
        );
        assert_eq!(clipped.len(), 1);
    }

    #[test]
    fn snap_endpoint_near_and_far_from_boundary() {
        let plane = xy_plane();
        let mut near = Point3::new(0.0005, 0.5, 0.0);
        snap_endpoint_to_boundary(&mut near, &plane, (0.0, 1.0), (0.0, 1.0));
        assert!(near.x.abs() < 1e-10, "should snap to u=0, got x={}", near.x);

        let mut far = Point3::new(0.5, 0.5, 0.0);
        let orig = far;
        snap_endpoint_to_boundary(&mut far, &plane, (0.0, 1.0), (0.0, 1.0));
        assert!(
            far.distance_to(&orig) < 1e-10,
            "interior point should not snap"
        );
    }

    #[test]
    fn fit_nurbs_various_point_counts() {
        // Two points
        let pts2 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)];
        let c2 = fit_nurbs_through_points(&pts2).unwrap();
        assert_eq!(c2.degree, 1);
        assert_eq!(c2.control_points.len(), 2);

        // Five points
        let pts5: Vec<Point3> = (0..5)
            .map(|i| Point3::new(i as f64, (i as f64).sin(), 0.0))
            .collect();
        let c5 = fit_nurbs_through_points(&pts5).unwrap();
        assert_eq!(c5.control_points.len(), 5);
        let (t0, t1) = c5.parameter_range();
        assert!(c5.point_at(t0).distance_to(&pts5[0]) < 1e-10);
        assert!(c5.point_at(t1).distance_to(&pts5[4]) < 1e-10);

        // Degenerate (coincident)
        let degen = vec![Point3::new(1.0, 2.0, 3.0), Point3::new(1.0, 2.0, 3.0)];
        assert!(fit_nurbs_through_points(&degen).is_err());
    }

    #[test]
    fn clamp_finite_range_cases() {
        assert_eq!(clamp_finite_range((1.0, 5.0), (-10.0, 10.0)), (1.0, 5.0));
        assert_eq!(
            clamp_finite_range((f64::NEG_INFINITY, f64::INFINITY), (-10.0, 10.0)),
            (-10.0, 10.0)
        );
    }
}
