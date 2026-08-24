//! Surface offset computation.

pub mod analytic;
pub mod nurbs;
pub mod self_intersect;

pub use analytic::{AnalyticOffset, offset_analytic, offset_analytic_classified};
pub use nurbs::{offset_nurbs, offset_nurbs_measured};

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

/// An offset surface together with an honest record of how it was made.
#[derive(Debug, Clone)]
pub struct OffsetSurfaceOutcome {
    /// The offset surface.
    pub surface: SurfaceGeometry,
    /// `None` for a closed-form exact offset; `Some(dev)` for a sampled
    /// NURBS approximation, where `dev` is the measured maximum
    /// mid-sample deviation from the exact offset at the final grid.
    /// Callers decide whether `dev` meets their tolerance — the
    /// approximation itself never silently claims exactness.
    pub approximation_deviation: Option<f64>,
}

/// Compute a surface offset with an honest quality record.
///
/// Dispatches to the closed-form analytic path when the surface class
/// supports one, and falls back to NURBS-based adaptive sampling for
/// classes with no closed form (cone, NURBS, Gregory, polyhedral
/// spline).
///
/// # Errors
///
/// A genuine geometric decline (e.g. an inward offset that would drive
/// a radius through zero) is returned as `Err` and is NOT converted
/// into an approximation — the offset surface does not exist at that
/// distance, and approximating it would serve inverted geometry.
pub fn offset_surface_measured(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<OffsetSurfaceOutcome, String> {
    match offset_analytic_classified(surface, distance, inward)? {
        AnalyticOffset::Exact(s) => Ok(OffsetSurfaceOutcome {
            surface: s,
            approximation_deviation: None,
        }),
        AnalyticOffset::Unsupported => {
            let (s, dev) =
                offset_nurbs_measured(surface, distance, inward, u_range, v_range, tolerance)?;
            Ok(OffsetSurfaceOutcome {
                surface: s,
                approximation_deviation: Some(dev),
            })
        }
    }
}

/// Compute a surface offset, dispatching to the analytic path when possible
/// and falling back to NURBS-based adaptive sampling otherwise.
///
/// Convenience wrapper over [`offset_surface_measured`] for callers that
/// do not gate on approximation quality. Genuine geometric declines
/// (radius through zero) still propagate as `Err`.
pub fn offset_surface(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, String> {
    offset_surface_measured(surface, distance, inward, u_range, v_range, tolerance)
        .map(|outcome| outcome.surface)
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

    #[test]
    fn measured_outcome_marks_exact_and_approximate() {
        let tol = Tolerance::default();
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let exact =
            offset_surface_measured(&plane, 1.0, false, (-1.0, 1.0), (-1.0, 1.0), &tol).unwrap();
        assert!(exact.approximation_deviation.is_none(), "plane is exact");

        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 1.0),
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        };
        let approx =
            offset_surface_measured(&cone, 0.1, false, (0.1, 1.0), (0.1, 1.0), &tol).unwrap();
        let dev = approx
            .approximation_deviation
            .expect("cone offset is an approximation and must say so");
        assert!(dev.is_finite() && dev >= 0.0);
    }

    #[test]
    fn genuine_decline_is_not_approximated_away() {
        // An inward cylinder offset past its radius has no offset surface;
        // the old dispatch silently NURBS-approximated inverted geometry.
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let tol = Tolerance::default();
        let result = offset_surface_measured(
            &cyl,
            2.0,
            true,
            (0.0, std::f64::consts::TAU),
            (0.0, 2.0),
            &tol,
        );
        assert!(result.is_err(), "radius-through-zero must decline");
    }
}
