//! Guided sweep and loft operations.
//!
//! Adjusts a sweep or loft surface to pass through guide rail curves by
//! iteratively shifting the nearest control points toward guide sample
//! points. Falls back to deviation validation when guides already lie
//! within tolerance.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
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

/// Maximum iterations for control point adjustment.
const MAX_ADJUSTMENT_ITERATIONS: usize = 5;

/// Adjust a surface's control points so it passes through guide rail
/// curves within tolerance.
///
/// Uses iterative nearest-control-point shifting: for each guide sample
/// point, finds the closest (u, v) on the surface via grid search, then
/// shifts the nearest control point by a distance-weighted fraction of
/// the deviation vector. Repeats for up to [`MAX_ADJUSTMENT_ITERATIONS`]
/// iterations or until convergence.
///
/// # Errors
///
/// Returns [`SweepError::GuideDeviation`] if the surface cannot be
/// brought within tolerance after all iterations.
pub fn apply_guide_rails(
    surface: &mut NurbsSurface,
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

    // Check if already within tolerance.
    let initial_dev = compute_max_deviation(surface, guides);
    if initial_dev <= tolerance.linear {
        return Ok(GuideValidation {
            max_deviation: initial_dev,
            guide_deviations: compute_per_guide_deviations(surface, guides),
            warnings: vec![],
        });
    }

    // Iteratively adjust control points.
    for _iteration in 0..MAX_ADJUSTMENT_ITERATIONS {
        adjust_control_points_once(surface, guides);

        let current_dev = compute_max_deviation(surface, guides);
        if current_dev <= tolerance.linear {
            return Ok(GuideValidation {
                max_deviation: current_dev,
                guide_deviations: compute_per_guide_deviations(surface, guides),
                warnings: vec![],
            });
        }
    }

    // Final check after all iterations.
    let final_dev = compute_max_deviation(surface, guides);
    if final_dev > tolerance.linear {
        return Err(SweepError::GuideDeviation {
            deviation: final_dev,
        });
    }

    Ok(GuideValidation {
        max_deviation: final_dev,
        guide_deviations: compute_per_guide_deviations(surface, guides),
        warnings: vec![],
    })
}

/// Run one pass of control point adjustment for all guides.
///
/// For each guide sample point, finds the nearest (u, v) on the surface,
/// computes the deviation vector, locates the nearest control point, and
/// shifts it by the deviation weighted by inverse distance.
fn adjust_control_points_once(surface: &mut NurbsSurface, guides: &[NurbsCurve]) {
    for guide in guides {
        let (t_min, t_max) = guide.parameter_range();

        for i in 0..=GUIDE_SAMPLES {
            let frac = i as f64 / GUIDE_SAMPLES as f64;
            let t = t_min + frac * (t_max - t_min);
            let guide_pt = guide.point_at(t);

            // Find nearest (u, v) on the surface.
            let (best_u, best_v) = find_nearest_uv(surface, &guide_pt);
            let surf_pt = surface.point_at(best_u, best_v);

            // Deviation vector from surface to guide.
            let delta = guide_pt - surf_pt;
            let delta_len = delta.length();
            if delta_len < 1e-30 {
                continue;
            }

            // Find the nearest control point to (best_u, best_v).
            let (cp_iu, cp_iv) = find_nearest_control_point(surface, best_u, best_v);
            let cp_idx = cp_iu * surface.num_v + cp_iv;

            // Shift the control point by delta, weighted by the control
            // point's proximity to the parametric location.
            let weight = control_point_weight(surface, cp_iu, cp_iv, best_u, best_v);
            let shift = delta * weight;
            surface.control_points[cp_idx] = surface.control_points[cp_idx] + shift;
        }
    }
}

/// Find the (u, v) parameters of the nearest point on the surface to a
/// given 3D point, using a uniform grid search.
fn find_nearest_uv(surface: &NurbsSurface, point: &Point3) -> (f64, f64) {
    let mut best_u = 0.0;
    let mut best_v = 0.0;
    let mut best_dist = f64::MAX;

    for iu in 0..=GRID_DIVISIONS {
        let u = iu as f64 / GRID_DIVISIONS as f64;
        for iv in 0..=GRID_DIVISIONS {
            let v = iv as f64 / GRID_DIVISIONS as f64;
            let surf_pt = surface.point_at(u, v);
            let dist = point.distance_to(&surf_pt);
            if dist < best_dist {
                best_dist = dist;
                best_u = u;
                best_v = v;
            }
        }
    }

    (best_u, best_v)
}

/// Find the control point grid indices nearest to a (u, v) location.
///
/// Maps u and v into the control point grid using Greville abscissae
/// approximation (uniform spacing assumption).
fn find_nearest_control_point(surface: &NurbsSurface, u: f64, v: f64) -> (usize, usize) {
    let iu = ((u * (surface.num_u - 1) as f64).round() as usize).min(surface.num_u - 1);
    let iv = ((v * (surface.num_v - 1) as f64).round() as usize).min(surface.num_v - 1);
    (iu, iv)
}

/// Compute a distance-based weight for shifting a control point.
///
/// Returns a value in (0, 1] based on how close the control point's
/// approximate parametric location is to the target (u, v). Control
/// points near the target get weight close to 1; distant ones get less.
fn control_point_weight(
    surface: &NurbsSurface,
    cp_iu: usize,
    cp_iv: usize,
    target_u: f64,
    target_v: f64,
) -> f64 {
    let cp_u = if surface.num_u > 1 {
        cp_iu as f64 / (surface.num_u - 1) as f64
    } else {
        0.5
    };
    let cp_v = if surface.num_v > 1 {
        cp_iv as f64 / (surface.num_v - 1) as f64
    } else {
        0.5
    };

    let du = cp_u - target_u;
    let dv = cp_v - target_v;
    let param_dist = (du * du + dv * dv).sqrt();

    // Gaussian-like falloff: weight = exp(-k * dist^2)
    // With k = 10, a control point 0.3 away in parameter space gets ~0.4 weight.
    (-10.0 * param_dist * param_dist).exp()
}

/// Compute maximum deviation across all guides.
fn compute_max_deviation(surface: &NurbsSurface, guides: &[NurbsCurve]) -> f64 {
    guides
        .iter()
        .map(|g| measure_guide_deviation(surface, g))
        .fold(0.0_f64, f64::max)
}

/// Compute per-guide maximum deviations.
fn compute_per_guide_deviations(surface: &NurbsSurface, guides: &[NurbsCurve]) -> Vec<f64> {
    guides
        .iter()
        .map(|g| measure_guide_deviation(surface, g))
        .collect()
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
        let result = validate_guide_rails(&surface, &[on_surface_guide()], &Tolerance::default());
        assert!(result.is_ok(), "expected pass, got {:?}", result.err());
        let validation = result.unwrap(); // ok in tests
        assert!(
            validation.max_deviation < 0.1,
            "got {}",
            validation.max_deviation
        );
    }

    #[test]
    fn guide_far_from_surface_fails() {
        let result = validate_guide_rails(&flat_surface(), &[far_guide()], &Tolerance::default());
        assert!(result.is_err());
        match result.err() {
            Some(SweepError::GuideDeviation { deviation }) => {
                assert!(deviation > 1.0, "expected large deviation, got {deviation}");
            }
            other => panic!("expected GuideDeviation, got {:?}", other),
        }
    }

    #[test]
    fn empty_guides_passes() {
        let result = validate_guide_rails(&flat_surface(), &[], &Tolerance::default());
        assert!(result.is_ok());
        assert!(result.unwrap().max_deviation < 1e-10); // ok in tests
    }

    #[test]
    fn multiple_guides_reports_per_guide_deviations() {
        let result = validate_guide_rails(
            &flat_surface(),
            &[on_surface_guide(), far_guide()],
            &Tolerance::default(),
        );
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

    #[test]
    fn apply_guides_noop_for_empty() {
        let mut surface = flat_surface();
        let tol = Tolerance::default();
        let result = apply_guide_rails(&mut surface, &[], &tol);
        assert!(result.is_ok());
    }

    #[test]
    fn apply_guides_already_within_tolerance() {
        let mut surface = flat_surface();
        let guide = on_surface_guide();
        let tol = Tolerance::default();
        let result = apply_guide_rails(&mut surface, &[guide], &tol);
        assert!(
            result.is_ok(),
            "expected pass for on-surface guide, got {:?}",
            result.err()
        );
    }

    #[test]
    fn apply_guides_adjusts_surface_for_small_offset() {
        // Create a flat surface and a guide offset in Z.
        let mut surface = flat_surface();
        let offset = 0.3;
        let guide = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.25, 0.5, offset),
                Point3::new(0.75, 0.5, offset),
            ],
            weights: vec![1.0, 1.0],
        };

        // Before adjustment, deviation should be ~offset.
        let dev_before = measure_guide_deviation(&surface, &guide);
        assert!(
            dev_before > 0.1,
            "expected deviation >{}, got {dev_before}",
            0.1
        );

        // Apply with a tolerance that forces adjustment but is achievable.
        let tol = Tolerance {
            linear: 0.15,
            ..Tolerance::default()
        };
        let _result = apply_guide_rails(&mut surface, &[guide.clone()], &tol);
        // We may or may not converge within tolerance, but deviation
        // should be reduced from the original.
        let dev_after = measure_guide_deviation(&surface, &guide);
        assert!(
            dev_after < dev_before,
            "deviation not reduced: before={dev_before}, after={dev_after}"
        );
    }

    #[test]
    fn apply_guides_too_many_returns_error() {
        let mut surface = flat_surface();
        let guides: Vec<NurbsCurve> = (0..17).map(|_| on_surface_guide()).collect();
        let tol = Tolerance::default();
        let result = apply_guide_rails(&mut surface, &guides, &tol);
        assert!(result.is_err());
    }
}
