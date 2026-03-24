//! Surface offset computation.

pub mod analytic;
pub mod nurbs;
pub mod self_intersect;

pub use analytic::offset_analytic;
pub use nurbs::offset_nurbs;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

/// Compute a surface offset, dispatching to the analytic path when possible
/// and falling back to NURBS-based adaptive sampling otherwise.
pub fn offset_surface(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, String> {
    match offset_analytic(surface, distance, inward) {
        Ok(s) => Ok(s),
        Err(_) => offset_nurbs(surface, distance, inward, u_range, v_range, tolerance),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn offset_surface_dispatches_analytic_for_plane() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let tol = Tolerance::default();
        let result = offset_surface(&plane, 1.0, false, (-1.0, 1.0), (-1.0, 1.0), &tol).unwrap();
        // Should be analytic Plane, not NURBS
        match result {
            SurfaceGeometry::Plane { origin, .. } => {
                assert!((origin.z - 1.0).abs() < 1e-12);
            }
            _ => panic!("expected analytic Plane dispatch"),
        }
    }

    #[test]
    fn offset_surface_dispatches_nurbs_for_cone() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 1.0),
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        };
        let tol = Tolerance::default();
        let result = offset_surface(&cone, 0.1, false, (0.1, 1.0), (0.1, 1.0), &tol).unwrap();
        match result {
            SurfaceGeometry::Nurbs(_) => {} // expected
            _ => panic!("expected NURBS fallback for cone"),
        }
    }
}
