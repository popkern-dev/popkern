//! Surface offset computation.

pub mod analytic;
pub mod error;
pub mod nurbs;
mod nurbs_fit;
pub mod self_intersect;

pub use analytic::{
    AnalyticOffset, corresponding_offset_window, offset_analytic, offset_analytic_classified,
};
pub use error::OffsetError;
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
/// supports one (plane, cylinder, cone, sphere, torus), and falls back
/// to degree-3 least-squares fitting for classes with no closed form
/// (NURBS, Gregory, polyhedral spline).
///
/// # Errors
///
/// A genuine geometric decline (an offset that would drive a radius
/// through zero, or grow a torus into a spindle) is returned as `Err`
/// and is NOT converted into an approximation — the offset surface does
/// not exist at that distance, and approximating it would serve
/// inverted geometry. On the approximation path, a distance exceeding
/// the local curvature radius inside the window folds the offset map,
/// so fitting is declined ([`OffsetError::OffsetMapFolds`]) — exact-arm
/// classes instead handle folds downstream via trimming and the
/// collapse gate.
pub fn offset_surface_measured(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<OffsetSurfaceOutcome, OffsetError> {
    match offset_analytic_classified(surface, distance, inward)? {
        AnalyticOffset::Exact(s) => Ok(OffsetSurfaceOutcome {
            surface: s,
            approximation_deviation: None,
        }),
        AnalyticOffset::Unsupported => {
            if let Some(max_safe) = self_intersect::check_offset_self_intersection(
                surface, distance, inward, u_range, v_range,
            ) {
                return Err(OffsetError::OffsetMapFolds {
                    requested: distance,
                    max_safe,
                });
            }
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
) -> Result<SurfaceGeometry, OffsetError> {
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
    fn offset_surface_dispatches_exact_for_cone() {
        // The cone gained a closed-form arm (PK-322): same half-angle,
        // apex translated d/sin(half_angle) along the axis — never the
        // NURBS fitting path.
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 1.0),
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        };
        let tol = Tolerance::default();
        let result = offset_surface(&cone, 0.1, false, (0.1, 1.0), (0.1, 1.0), &tol).unwrap();
        match result {
            SurfaceGeometry::Cone { apex, .. } => {
                let shift = 0.1 / std::f64::consts::FRAC_PI_4.sin();
                assert!((apex.z - (1.0 + shift)).abs() < 1e-12, "apex.z {}", apex.z);
            }
            other => panic!("expected exact Cone offset, got {other:?}"),
        }
    }

    fn bump_nurbs() -> SurfaceGeometry {
        use popkern_core::geometry::NurbsSurface;
        SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.3),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 0.5, 0.3),
                Point3::new(0.5, 0.5, 0.6),
                Point3::new(1.0, 0.5, 0.3),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.5, 1.0, 0.3),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 3,
            num_v: 3,
            weights: vec![1.0; 9],
        })
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

        let approx =
            offset_surface_measured(&bump_nurbs(), 0.05, false, (0.0, 1.0), (0.0, 1.0), &tol)
                .unwrap();
        let dev = approx
            .approximation_deviation
            .expect("free-form NURBS offset is an approximation and must say so");
        assert!(dev.is_finite() && dev >= 0.0);
    }

    #[test]
    fn folding_approximation_declines_with_map_folds() {
        // A curved free-form surface offset far past its curvature
        // radius: fitting would approximate a folded offset map, so the
        // dispatch declines with the fold bound instead.
        let tol = Tolerance::default();
        let result =
            offset_surface_measured(&bump_nurbs(), 50.0, true, (0.0, 1.0), (0.0, 1.0), &tol);
        match result {
            Err(OffsetError::OffsetMapFolds {
                requested,
                max_safe,
            }) => {
                assert!((requested - 50.0).abs() < 1e-12);
                assert!(max_safe > 0.0 && max_safe < 50.0, "max_safe {max_safe}");
            }
            other => panic!("expected OffsetMapFolds decline, got {other:?}"),
        }
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
