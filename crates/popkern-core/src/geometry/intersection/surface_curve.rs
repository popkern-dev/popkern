//! Surface-curve intersection types and dispatch.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;
use serde::{Deserialize, Serialize};

/// A single intersection result between a surface and a curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceCurveIntersection {
    /// Curve pierces the surface at a point.
    Point {
        /// The intersection point in 3D space.
        point: Point3,
        /// Parameter value on the curve at the intersection.
        t: f64,
        /// Surface parameter u at the intersection.
        u: f64,
        /// Surface parameter v at the intersection.
        v: f64,
    },
    /// Curve lies on the surface over a parameter interval.
    Overlap {
        /// Parameter range [start, end] on the curve.
        t_range: (f64, f64),
    },
}

/// Result of a surface-curve intersection operation.
#[derive(Debug, Clone)]
pub struct SurfaceCurveIntersectionResult {
    /// All intersections found, sorted by t.
    pub intersections: Vec<SurfaceCurveIntersection>,
}

#[path = "surface_curve_dispatch.rs"]
mod surface_curve_dispatch;

/// Intersect a surface with a curve within the given tolerance.
///
/// Dispatches to specialised analytic solvers for known geometry pairs
/// (Line vs Plane/Cylinder/Cone/Sphere/Torus), falling back to the general
/// NURBS-based solver for everything else.
pub fn intersect_surface_curve(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceCurveIntersectionResult, CoreError> {
    intersect_surface_curve_inner(surface, u_range, v_range, curve, t_range, tolerance, false)
}

/// [`intersect_surface_curve`] with COMPLETE crossing enumeration on
/// multi-wrap NURBS charts (PK-494): the general solver's AABB grid
/// spends a total cell budget on the denser knot axis instead of the
/// tuned 48-per-axis clamp. For parity consumers ONLY — point-in-solid
/// classification rays, where a missed crossing on a >40-turn coil
/// flips an in/out verdict. Imprint/edge-hit discovery must keep the
/// tuned entry: the split stage's junction graph is calibrated to its
/// hit set, and enumerating extra near-graze hits there re-wired the
/// M5 seam-pinch serve into a shell-gap decline (measured, PK-461).
pub fn intersect_surface_curve_dense(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceCurveIntersectionResult, CoreError> {
    intersect_surface_curve_inner(surface, u_range, v_range, curve, t_range, tolerance, true)
}

#[allow(clippy::too_many_arguments)]
fn intersect_surface_curve_inner(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
    dense_grid: bool,
) -> Result<SurfaceCurveIntersectionResult, CoreError> {
    if t_range.0 >= t_range.1 {
        return Err(CoreError::PreconditionFailed(format!(
            "t_range invalid: {} >= {}",
            t_range.0, t_range.1
        )));
    }

    let mut intersections = surface_curve_dispatch::dispatch(
        surface, u_range, v_range, curve, t_range, tolerance, dense_grid,
    );

    // Sort by t.
    intersections.sort_by(|a, b| {
        let ta = match a {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        let tb = match b {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(SurfaceCurveIntersectionResult { intersections })
}

#[cfg(test)]
#[path = "surface_curve_tests.rs"]
mod tests;

// Remaining tests in surface_curve_tests.rs
#[cfg(test)]
mod extra_tests {
    use super::*;
    use crate::geometry::curve::{CurveEval, CurveGeometry};
    use crate::geometry::nurbs_curve::NurbsCurve;
    use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;

    // -----------------------------------------------------------------------
    // dispatch_line_plane: verify the analytic path is used (correct result)
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_line_plane() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, -1.0),
            direction: Vec3::new(0.0, 0.0, 1.0),
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &line,
            (0.0, 2.0),
            &tol,
        );
        let res = result.as_ref().ok();
        assert!(res.is_some());
        let pts: Vec<_> = res
            .into_iter()
            .flat_map(|r| r.intersections.iter())
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert_eq!(pts.len(), 1);
        if let SurfaceCurveIntersection::Point { t, point, .. } = &pts[0] {
            assert!((*t - 1.0).abs() < 1e-10, "t should be 1.0, got {}", t);
            assert!(point.z.abs() < 1e-10, "z should be 0.0, got {}", point.z);
        }
    }

    // -----------------------------------------------------------------------
    // dispatch_line_sphere: verify line-sphere analytic path
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_line_sphere() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(2.0, 0.0, 0.0),
            direction: Vec3::new(-1.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::PI),
            &line,
            (0.0, 4.0),
            &tol,
        );
        let res = result.as_ref().ok();
        assert!(res.is_some());
        let pts: Vec<_> = res
            .into_iter()
            .flat_map(|r| r.intersections.iter())
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert!(
            pts.len() >= 2,
            "Expected at least 2 intersections for line through sphere, got {}",
            pts.len()
        );
    }

    // -----------------------------------------------------------------------
    // dispatch_circle_on_plane: general solver path (circle is not a Line)
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_circle_on_plane() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // Circle in xz-plane pierces z=0 at two points.
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            radius: 1.0,
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &circle,
            (0.0, std::f64::consts::TAU),
            &tol,
        );
        let res = result.as_ref().ok();
        assert!(res.is_some());
        let pts: Vec<_> = res
            .into_iter()
            .flat_map(|r| r.intersections.iter())
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert!(
            pts.len() >= 2,
            "Expected at least 2 intersections for circle piercing plane, got {}",
            pts.len()
        );
    }

    // -----------------------------------------------------------------------
    // dispatch_nurbs_on_cylinder: general solver path (NURBS curve)
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_nurbs_on_cylinder() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        // A NURBS line crossing through the cylinder along the x-axis.
        let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(-2.0, 0.0, 0.5), Point3::new(2.0, 0.0, 0.5)],
            weights: vec![1.0, 1.0],
        });
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &cyl,
            (0.0, std::f64::consts::TAU),
            (-5.0, 5.0),
            &nurbs_line,
            (0.0, 1.0),
            &tol,
        );
        let res = result.as_ref().ok();
        assert!(res.is_some());
        let pts: Vec<_> = res
            .into_iter()
            .flat_map(|r| r.intersections.iter())
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert!(
            pts.len() >= 2,
            "Expected at least 2 intersections for line through cylinder, got {}",
            pts.len()
        );
    }

    // -----------------------------------------------------------------------
    // dispatch_line_cylinder: verify analytic path
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_line_cylinder() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(2.0, 0.0, 0.5),
            direction: Vec3::new(-1.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &cyl,
            (0.0, std::f64::consts::TAU),
            (-5.0, 5.0),
            &line,
            (0.0, 4.0),
            &tol,
        )
        .unwrap();
        let pts: Vec<_> = result
            .intersections
            .iter()
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert!(
            pts.len() >= 2,
            "Expected 2 intersections for line through cylinder, got {}",
            pts.len()
        );
    }

    // -----------------------------------------------------------------------
    // dispatch_line_cone: verify analytic path
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_line_cone() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(2.0, 0.0, 1.0),
            direction: Vec3::new(-1.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &cone,
            (0.0, std::f64::consts::TAU),
            (0.0, 10.0),
            &line,
            (0.0, 4.0),
            &tol,
        )
        .unwrap();
        let pts: Vec<_> = result
            .intersections
            .iter()
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert!(
            !pts.is_empty(),
            "Expected at least 1 intersection for line through cone, got {}",
            pts.len()
        );
    }

    // -----------------------------------------------------------------------
    // dispatch_line_torus: verify analytic path
    // -----------------------------------------------------------------------
    #[test]
    fn dispatch_line_torus() {
        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 0.5,
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(5.0, 0.0, 0.0),
            direction: Vec3::new(-1.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &torus,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::TAU),
            &line,
            (0.0, 10.0),
            &tol,
        )
        .unwrap();
        let pts: Vec<_> = result
            .intersections
            .iter()
            .filter(|i| matches!(i, SurfaceCurveIntersection::Point { .. }))
            .collect();
        assert!(
            pts.len() >= 2,
            "Expected at least 2 intersections for line through torus, got {}",
            pts.len()
        );
    }

    // -----------------------------------------------------------------------
    // invariant_all_points_match: for every Point result, verify
    //   curve.point_at(t) ~ surface.point_at(u, v)
    // -----------------------------------------------------------------------
    #[test]
    fn invariant_all_points_match() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let line = CurveGeometry::Line {
            origin: Point3::new(1.0, 2.0, -1.0),
            direction: Vec3::new(0.0, 0.0, 1.0),
        };
        let tol = Tolerance::default();
        let result = intersect_surface_curve(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &line,
            (0.0, 2.0),
            &tol,
        );
        let res = result.ok();
        assert!(res.is_some());
        for i in &res
            .as_ref()
            .into_iter()
            .flat_map(|r| r.intersections.iter())
            .collect::<Vec<_>>()
        {
            if let SurfaceCurveIntersection::Point { point, t, u, v } = i {
                let curve_pt = line.point_at(*t);
                let surface_pt = plane.point_at(*u, *v);
                assert!(
                    point.distance_to(&curve_pt) < 1e-4,
                    "curve.point_at(t) mismatch"
                );
                assert!(
                    point.distance_to(&surface_pt) < 1e-4,
                    "surface.point_at(u,v) mismatch"
                );
            }
        }
    }
}
