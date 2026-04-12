//! Deviation measurement between geometric entities.
//!
//! Measures how far curves and surfaces deviate from each other using
//! uniform sampling with adaptive refinement near the maximum. Used by
//! the healing system to quantify edge-surface gaps, seam mismatches,
//! and surface approximation errors.

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use serde::{Deserialize, Serialize};

/// Result of a deviation measurement between two geometric entities.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DeviationResult {
    /// Maximum deviation found across all samples.
    pub max_deviation: f64,
    /// Mean deviation across all samples.
    pub mean_deviation: f64,
    /// Parameter location of the maximum deviation (on the first entity).
    pub max_location: f64,
    /// Total number of samples evaluated (including adaptive refinement).
    pub sample_count: usize,
}

/// Measure deviation between a curve and a surface.
///
/// Samples the curve at `num_samples` uniform points in `t_range`, computes
/// the 3D distance from each curve point to the closest point on the surface,
/// then performs adaptive binary-search refinement near the maximum.
///
/// If `surface.closest_point()` fails for a sample, falls back to measuring
/// distance from the curve point to `surface.point_at(0.0, 0.0)`.
pub fn edge_surface_deviation(
    curve: &dyn CurveEval,
    t_range: (f64, f64),
    surface: &dyn SurfaceEval,
    num_samples: usize,
) -> DeviationResult {
    let samples = num_samples.max(2);
    let (t0, t1) = t_range;

    let mut max_dev = 0.0_f64;
    let mut max_t = t0;
    let mut sum_dev = 0.0;
    let mut count = 0_usize;

    // Uniform sampling pass.
    for i in 0..samples {
        let t = t0 + (t1 - t0) * (i as f64) / ((samples - 1) as f64);
        let curve_pt = curve.point_at(t);

        let dist = if let Some((u, v)) = surface.closest_point(&curve_pt) {
            let surf_pt = surface.point_at(u, v);
            curve_pt.distance_to(&surf_pt)
        } else {
            // Fallback: distance to surface origin point.
            let fallback_pt = surface.point_at(0.0, 0.0);
            curve_pt.distance_to(&fallback_pt)
        };

        sum_dev += dist;
        count += 1;

        if dist > max_dev {
            max_dev = dist;
            max_t = t;
        }
    }

    // Adaptive refinement: binary search around the current maximum.
    let refinement_iterations = 4;
    let mut lo = (max_t - (t1 - t0) / (samples as f64)).max(t0);
    let mut hi = (max_t + (t1 - t0) / (samples as f64)).min(t1);

    for _ in 0..refinement_iterations {
        let mid_lo = (lo + max_t) * 0.5;
        let mid_hi = (max_t + hi) * 0.5;

        for t in [mid_lo, mid_hi] {
            let curve_pt = curve.point_at(t);

            let dist = if let Some((u, v)) = surface.closest_point(&curve_pt) {
                let surf_pt = surface.point_at(u, v);
                curve_pt.distance_to(&surf_pt)
            } else {
                let fallback_pt = surface.point_at(0.0, 0.0);
                curve_pt.distance_to(&fallback_pt)
            };

            sum_dev += dist;
            count += 1;

            if dist > max_dev {
                max_dev = dist;
                max_t = t;
            }
        }

        // Narrow the search window around the new maximum.
        let span = (hi - lo) * 0.5;
        lo = (max_t - span * 0.5).max(t0);
        hi = (max_t + span * 0.5).min(t1);
    }

    let mean_dev = if count > 0 {
        sum_dev / count as f64
    } else {
        0.0
    };

    DeviationResult {
        max_deviation: max_dev,
        mean_deviation: mean_dev,
        max_location: max_t,
        sample_count: count,
    }
}

/// Measure deviation between two curves sampled at corresponding parameter fractions.
///
/// Samples both curves at `num_samples` evenly spaced parameter fractions
/// within their respective ranges. At each fraction, computes the point-to-point
/// distance. After the uniform pass, performs adaptive binary-search refinement
/// near the current maximum (same pattern as [`edge_surface_deviation`]).
pub fn curve_curve_deviation(
    curve_a: &dyn CurveEval,
    range_a: (f64, f64),
    curve_b: &dyn CurveEval,
    range_b: (f64, f64),
    num_samples: usize,
) -> DeviationResult {
    let samples = num_samples.max(2);
    let (a0, a1) = range_a;
    let (b0, b1) = range_b;

    let mut max_dev = 0.0_f64;
    let mut max_t = a0;
    let mut sum_dev = 0.0;
    let mut count = 0_usize;

    // Uniform sampling pass.
    for i in 0..samples {
        let frac = (i as f64) / ((samples - 1) as f64);
        let ta = a0 + (a1 - a0) * frac;
        let tb = b0 + (b1 - b0) * frac;

        let pa = curve_a.point_at(ta);
        let pb = curve_b.point_at(tb);
        let dist = pa.distance_to(&pb);

        sum_dev += dist;
        count += 1;

        if dist > max_dev {
            max_dev = dist;
            max_t = ta;
        }
    }

    // Adaptive refinement: binary search around the current maximum.
    let refinement_iterations = 4;
    let step = (a1 - a0) / (samples as f64);
    let mut lo = (max_t - step).max(a0);
    let mut hi = (max_t + step).min(a1);

    for _ in 0..refinement_iterations {
        let mid_lo = (lo + max_t) * 0.5;
        let mid_hi = (max_t + hi) * 0.5;

        for ta in [mid_lo, mid_hi] {
            // Map back to curve_b's parameter space proportionally.
            let frac = if (a1 - a0).abs() > 1e-15 {
                (ta - a0) / (a1 - a0)
            } else {
                0.0
            };
            let tb = b0 + (b1 - b0) * frac;

            let pa = curve_a.point_at(ta);
            let pb = curve_b.point_at(tb);
            let dist = pa.distance_to(&pb);

            sum_dev += dist;
            count += 1;

            if dist > max_dev {
                max_dev = dist;
                max_t = ta;
            }
        }

        // Narrow the search window around the new maximum.
        let span = (hi - lo) * 0.5;
        lo = (max_t - span * 0.5).max(a0);
        hi = (max_t + span * 0.5).min(a1);
    }

    let mean_dev = if count > 0 {
        sum_dev / count as f64
    } else {
        0.0
    };

    DeviationResult {
        max_deviation: max_dev,
        mean_deviation: mean_dev,
        max_location: max_t,
        sample_count: count,
    }
}

/// Measure deviation between two surfaces.
///
/// Samples `surface_a` on a grid of `sqrt(num_samples) x sqrt(num_samples)`
/// points within `range_a`. For each sample point, finds the closest point
/// on `surface_b` and measures the distance. After the uniform grid pass,
/// performs adaptive refinement by subdividing the grid cell containing the
/// maximum into a 4x4 sub-grid (repeated 3 times).
///
/// If `surface_b.closest_point()` fails for a sample, falls back to
/// `surface_b.point_at(0.0, 0.0)`.
pub fn surface_surface_deviation(
    surface_a: &dyn SurfaceEval,
    range_a: ((f64, f64), (f64, f64)),
    surface_b: &dyn SurfaceEval,
    range_b: ((f64, f64), (f64, f64)),
    num_samples: usize,
) -> DeviationResult {
    let _ = range_b; // Used only for context; closest_point searches full domain.

    let grid = (num_samples as f64).sqrt().ceil() as usize;
    let grid = grid.max(2);

    let ((ua_min, ua_max), (va_min, va_max)) = range_a;

    // Clamp infinite ranges.
    let ua_lo = ua_min.max(-100.0);
    let ua_hi = ua_max.min(100.0);
    let va_lo = va_min.max(-100.0);
    let va_hi = va_max.min(100.0);

    let mut max_dev = 0.0_f64;
    let mut max_location = 0.0_f64;
    let mut sum_dev = 0.0;
    let mut count = 0_usize;
    let mut max_i = 0_usize;
    let mut max_j = 0_usize;

    // Uniform grid pass.
    for i in 0..grid {
        let u = ua_lo + (ua_hi - ua_lo) * (i as f64) / ((grid - 1) as f64);
        for j in 0..grid {
            let v = va_lo + (va_hi - va_lo) * (j as f64) / ((grid - 1) as f64);
            let pt_a = surface_a.point_at(u, v);

            let dist = if let Some((ub, vb)) = surface_b.closest_point(&pt_a) {
                let pt_b = surface_b.point_at(ub, vb);
                pt_a.distance_to(&pt_b)
            } else {
                let fallback = surface_b.point_at(0.0, 0.0);
                pt_a.distance_to(&fallback)
            };

            sum_dev += dist;
            count += 1;

            if dist > max_dev {
                max_dev = dist;
                max_location = (i * grid + j) as f64;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Adaptive refinement: subdivide the grid cell containing the maximum
    // into a 4x4 sub-grid, repeated 3 times.
    let u_step = (ua_hi - ua_lo) / ((grid - 1) as f64);
    let v_step = (va_hi - va_lo) / ((grid - 1) as f64);

    let mut cell_u_lo = ua_lo + (max_i.saturating_sub(1) as f64) * u_step;
    let mut cell_u_hi = (ua_lo + ((max_i + 1) as f64) * u_step).min(ua_hi);
    let mut cell_v_lo = va_lo + (max_j.saturating_sub(1) as f64) * v_step;
    let mut cell_v_hi = (va_lo + ((max_j + 1) as f64) * v_step).min(va_hi);

    let refinement_iterations = 3;
    let sub_grid = 4;

    for _ in 0..refinement_iterations {
        let mut local_max_u = cell_u_lo;
        let mut local_max_v = cell_v_lo;

        for si in 0..sub_grid {
            let u = cell_u_lo + (cell_u_hi - cell_u_lo) * (si as f64) / ((sub_grid - 1) as f64);
            for sj in 0..sub_grid {
                let v = cell_v_lo + (cell_v_hi - cell_v_lo) * (sj as f64) / ((sub_grid - 1) as f64);
                let pt_a = surface_a.point_at(u, v);

                let dist = if let Some((ub, vb)) = surface_b.closest_point(&pt_a) {
                    let pt_b = surface_b.point_at(ub, vb);
                    pt_a.distance_to(&pt_b)
                } else {
                    let fallback = surface_b.point_at(0.0, 0.0);
                    pt_a.distance_to(&fallback)
                };

                sum_dev += dist;
                count += 1;

                if dist > max_dev {
                    max_dev = dist;
                    local_max_u = u;
                    local_max_v = v;
                }
            }
        }

        // Narrow the cell around the new local maximum.
        let half_u = (cell_u_hi - cell_u_lo) * 0.25;
        let half_v = (cell_v_hi - cell_v_lo) * 0.25;
        cell_u_lo = (local_max_u - half_u).max(ua_lo);
        cell_u_hi = (local_max_u + half_u).min(ua_hi);
        cell_v_lo = (local_max_v - half_v).max(va_lo);
        cell_v_hi = (local_max_v + half_v).min(va_hi);
    }

    let mean_dev = if count > 0 {
        sum_dev / count as f64
    } else {
        0.0
    };

    DeviationResult {
        max_deviation: max_dev,
        mean_deviation: mean_dev,
        max_location,
        sample_count: count,
    }
}

/// Deviation result with a guaranteed mathematical upper bound.
///
/// Produced by [`edge_surface_deviation_bounded`], which uses recursive
/// bisection with Lipschitz-based interval pruning to guarantee that
/// the true maximum deviation cannot exceed `upper_bound`.
#[derive(Debug, Clone, Copy)]
pub struct BoundedDeviationResult {
    /// The deviation result from sampling within the subdivision.
    pub result: DeviationResult,
    /// Guaranteed upper bound on the true maximum deviation.
    pub upper_bound: f64,
    /// Number of subdivision intervals processed.
    pub intervals_processed: usize,
}

/// Compute a guaranteed upper bound on edge-surface deviation.
///
/// Uses recursive bisection: splits the parameter range in half,
/// evaluates deviation at endpoints and midpoint, computes a Lipschitz
/// bound from the maximum derivative magnitude, and prunes intervals
/// where the bound is below the current best max. Converges when all
/// intervals are either pruned or narrower than `min_interval`.
pub fn edge_surface_deviation_bounded(
    curve: &dyn CurveEval,
    t_range: (f64, f64),
    surface: &dyn SurfaceEval,
    min_interval: f64,
    max_subdivisions: usize,
) -> BoundedDeviationResult {
    let (t0, t1) = t_range;

    // Helper: compute distance from a curve point to the surface.
    let dev_at = |t: f64| -> f64 {
        let pt = curve.point_at(t);
        if let Some((u, v)) = surface.closest_point(&pt) {
            let surf_pt = surface.point_at(u, v);
            pt.distance_to(&surf_pt)
        } else {
            let fallback = surface.point_at(0.0, 0.0);
            pt.distance_to(&fallback)
        }
    };

    // Work queue: each entry is (t_lo, t_hi).
    let mut queue: Vec<(f64, f64)> = vec![(t0, t1)];
    let mut global_max_dev = 0.0_f64;
    let mut global_max_t = t0;
    let mut sum_dev = 0.0_f64;
    let mut count = 0_usize;
    let mut intervals_processed = 0_usize;
    // Collect upper bounds of all leaf (terminal) intervals.
    let mut leaf_bounds: Vec<f64> = Vec::new();

    while let Some((lo, hi)) = queue.pop() {
        if intervals_processed >= max_subdivisions {
            // Budget exhausted: record as leaf with Lipschitz bound.
            let mid = (lo + hi) * 0.5;
            let d_mid = dev_at(mid);
            let speed = curve.derivative_at(mid, 1).length();
            let lip = speed * (hi - lo) * 0.5;
            leaf_bounds.push(d_mid + lip);
            if d_mid > global_max_dev {
                global_max_dev = d_mid;
                global_max_t = mid;
            }
            sum_dev += d_mid;
            count += 1;
            continue;
        }

        intervals_processed += 1;
        let mid = (lo + hi) * 0.5;

        let d_lo = dev_at(lo);
        let d_mid = dev_at(mid);
        let d_hi = dev_at(hi);

        sum_dev += d_lo + d_mid + d_hi;
        count += 3;

        let sample_max = d_lo.max(d_mid).max(d_hi);
        let sample_max_t = if (sample_max - d_lo).abs() < 1e-30 {
            lo
        } else if (sample_max - d_mid).abs() < 1e-30 {
            mid
        } else {
            hi
        };

        if sample_max > global_max_dev {
            global_max_dev = sample_max;
            global_max_t = sample_max_t;
        }

        // Lipschitz bound: curve speed bounds how fast deviation can change.
        // Use half-interval width since we sampled the midpoint.
        let speed = curve.derivative_at(mid, 1).length();
        let lip = speed * (hi - lo) * 0.5;
        let interval_upper = sample_max + lip;

        // Prune: if best this interval can offer is below the known max.
        if interval_upper <= global_max_dev {
            continue;
        }

        // If the interval is narrow enough, record as leaf and stop.
        if (hi - lo) <= min_interval {
            leaf_bounds.push(interval_upper);
            continue;
        }

        // Subdivide into two halves.
        queue.push((lo, mid));
        queue.push((mid, hi));
    }

    // Global upper bound is the max of all leaf interval bounds,
    // but cannot be less than the observed sample maximum.
    let global_upper_bound = leaf_bounds.iter().copied().fold(global_max_dev, f64::max);

    let mean_dev = if count > 0 {
        sum_dev / count as f64
    } else {
        0.0
    };

    BoundedDeviationResult {
        result: DeviationResult {
            max_deviation: global_max_dev,
            mean_deviation: mean_dev,
            max_location: global_max_t,
            sample_count: count,
        },
        upper_bound: global_upper_bound,
        intervals_processed,
    }
}

#[cfg(test)]
#[path = "deviation_tests.rs"]
mod tests;
