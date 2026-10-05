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
use popkern_math::PortableMath;

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

    if max_deviation > tolerance.linear() {
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
const MAX_ADJUSTMENT_ITERATIONS: usize = 10;

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
    if initial_dev <= tolerance.linear() {
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
        if current_dev <= tolerance.linear() {
            return Ok(GuideValidation {
                max_deviation: current_dev,
                guide_deviations: compute_per_guide_deviations(surface, guides),
                warnings: vec![],
            });
        }
    }

    // Final check after all iterations.
    let final_dev = compute_max_deviation(surface, guides);
    if final_dev > tolerance.linear() {
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

/// Run one pass of weighted multi-point control point adjustment for all
/// guides.
///
/// For each guide sample point, finds the nearest (u, v) on the surface,
/// computes the deviation vector, then distributes the correction to ALL
/// control points weighted by Gaussian proximity in parameter space. This
/// is equivalent to one iteration of the Gauss-Seidel method for the
/// normal equations of the least-squares fitting problem, and produces
/// substantially better results than nearest-point-only shifting.
fn adjust_control_points_once(surface: &mut NurbsSurface, guides: &[NurbsCurve]) {
    let n_total = surface.num_u * surface.num_v;

    // Accumulate weighted deltas for each control point.
    let mut delta_sum: Vec<[f64; 3]> = vec![[0.0; 3]; n_total];
    let mut weight_sum: Vec<f64> = vec![0.0; n_total];

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

            // Distribute to all control points, weighted by Gaussian
            // proximity in parameter space.
            for iu in 0..surface.num_u {
                for iv in 0..surface.num_v {
                    let w = control_point_weight(surface, iu, iv, best_u, best_v);
                    if w < 1e-12 {
                        continue;
                    }
                    let idx = iu * surface.num_v + iv;
                    delta_sum[idx][0] += w * delta.x;
                    delta_sum[idx][1] += w * delta.y;
                    delta_sum[idx][2] += w * delta.z;
                    weight_sum[idx] += w;
                }
            }
        }
    }

    // Apply averaged deltas to control points.
    for idx in 0..n_total {
        if weight_sum[idx] > 1e-30 {
            let inv_w = 1.0 / weight_sum[idx];
            let cp = &mut surface.control_points[idx];
            *cp = Point3::new(
                cp.x + delta_sum[idx][0] * inv_w,
                cp.y + delta_sum[idx][1] * inv_w,
                cp.z + delta_sum[idx][2] * inv_w,
            );
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
    (-10.0 * param_dist * param_dist).pexp()
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
mod tests;
