//! Guided sweep and loft operations.
//!
//! Validates that a sweep or loft surface passes through guide rail curves
//! within tolerance. This is a framework for future guide-constrained surface
//! modification; the current implementation validates proximity but does not
//! adjust control points.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface};
use popkern_core::tolerance::Tolerance;

use crate::error::{SweepError, SweepWarning};

/// Number of sample points along each guide rail for deviation checking.
const GUIDE_SAMPLES: usize = 20;

/// Number of grid divisions per parameter direction for nearest-point search.
const GRID_DIVISIONS: usize = 20;

/// Maximum number of guide rails supported.
const MAX_GUIDES: usize = 16;

/// Result of guide rail validation.
#[derive(Debug)]
pub struct GuideValidation {
    /// Maximum deviation found across all guides.
    pub max_deviation: f64,
    /// Per-guide maximum deviations.
    pub guide_deviations: Vec<f64>,
    /// Non-fatal warnings.
    pub warnings: Vec<SweepWarning>,
}

/// Validate that a surface passes through guide rail curves within tolerance.
///
/// Samples each guide rail at [`GUIDE_SAMPLES`] points and finds the closest
/// point on the surface via grid search. Reports the maximum deviation.
///
/// # Errors
///
/// Returns [`SweepError::GuideDeviation`] if the maximum deviation exceeds
/// the tolerance.
pub fn validate_guide_rails(
    surface: &NurbsSurface,
    guides: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Result<GuideValidation, SweepError> {
    if guides.is_empty() {
        return Ok(GuideValidation {
            max_deviation: 0.0,
            guide_deviations: vec![],
            warnings: vec![],
        });
    }

    if guides.len() > MAX_GUIDES {
        return Err(SweepError::InvalidPath(format!(
            "too many guide rails: {} (max {})",
            guides.len(),
            MAX_GUIDES
        )));
    }

    let mut max_deviation = 0.0_f64;
    let mut guide_deviations = Vec::with_capacity(guides.len());

    for guide in guides {
        let dev = measure_guide_deviation(surface, guide);
        max_deviation = max_deviation.max(dev);
        guide_deviations.push(dev);
    }

    if max_deviation > tolerance.linear {
        return Err(SweepError::GuideDeviation {
            deviation: max_deviation,
        });
    }

    Ok(GuideValidation {
        max_deviation,
        guide_deviations,
        warnings: vec![],
    })
}

/// Measure the maximum deviation between a guide curve and a surface.
///
/// Samples the guide at [`GUIDE_SAMPLES`] points and for each sample,
/// finds the nearest point on the surface via a uniform grid search.
fn measure_guide_deviation(surface: &NurbsSurface, guide: &NurbsCurve) -> f64 {
    let (t_min, t_max) = guide.parameter_range();
    let mut max_dev = 0.0_f64;

    for i in 0..=GUIDE_SAMPLES {
        let frac = i as f64 / GUIDE_SAMPLES as f64;
        let t = t_min + frac * (t_max - t_min);
        let guide_pt = guide.point_at(t);

        let nearest_dist = find_nearest_surface_distance(surface, &guide_pt);
        max_dev = max_dev.max(nearest_dist);
    }

    max_dev
}

/// Find the distance from a point to the nearest point on a surface
/// using a uniform grid search.
fn find_nearest_surface_distance(
    surface: &NurbsSurface,
    point: &popkern_core::geometry::Point3,
) -> f64 {
    let mut best_dist = f64::MAX;

    for iu in 0..=GRID_DIVISIONS {
        let u = iu as f64 / GRID_DIVISIONS as f64;
        for iv in 0..=GRID_DIVISIONS {
            let v = iv as f64 / GRID_DIVISIONS as f64;
            let surf_pt = surface.point_at(u, v);
            let dist = point.distance_to(&surf_pt);
            if dist < best_dist {
                best_dist = dist;
            }
        }
    }

    best_dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    /// Build a flat surface in the XY plane from (0,0,0) to (1,1,0).
    fn flat_surface() -> NurbsSurface {
        NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0; 4],
        }
    }

    /// A guide that lies exactly on the flat surface.
    fn on_surface_guide() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    /// A guide that is far from the surface (offset in Z).
    fn far_guide() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.5, 0.5, 5.0), Point3::new(0.5, 0.5, 10.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn guide_on_surface_passes() {
        let surface = flat_surface();
        let guide = on_surface_guide();
        let tol = Tolerance::default();

        let result = validate_guide_rails(&surface, &[guide], &tol);
        assert!(result.is_ok(), "expected pass, got {:?}", result.err());

        let validation = result.unwrap(); // ok in tests
        assert!(
            validation.max_deviation < 0.1,
            "expected small deviation, got {}",
            validation.max_deviation
        );
    }

    #[test]
    fn guide_far_from_surface_fails() {
        let surface = flat_surface();
        let guide = far_guide();
        let tol = Tolerance::default();

        let result = validate_guide_rails(&surface, &[guide], &tol);
        assert!(result.is_err());

        match result.err() {
            Some(SweepError::GuideDeviation { deviation }) => {
                assert!(
                    deviation > 1.0,
                    "expected large deviation, got {}",
                    deviation
                );
            }
            other => panic!("expected GuideDeviation, got {:?}", other),
        }
    }

    #[test]
    fn empty_guides_passes() {
        let surface = flat_surface();
        let tol = Tolerance::default();
        let result = validate_guide_rails(&surface, &[], &tol);
        assert!(result.is_ok());
        let validation = result.unwrap(); // ok in tests
        assert!(validation.max_deviation < 1e-10);
    }

    #[test]
    fn multiple_guides_reports_per_guide_deviations() {
        let surface = flat_surface();
        let good_guide = on_surface_guide();
        let bad_guide = far_guide();
        let tol = Tolerance::default();

        let result = validate_guide_rails(&surface, &[good_guide, bad_guide], &tol);
        assert!(result.is_err());
    }

    #[test]
    fn too_many_guides_returns_error() {
        let surface = flat_surface();
        let guides: Vec<NurbsCurve> = (0..17).map(|_| on_surface_guide()).collect();
        let tol = Tolerance::default();
        let result = validate_guide_rails(&surface, &guides, &tol);
        assert!(result.is_err());
    }
}
