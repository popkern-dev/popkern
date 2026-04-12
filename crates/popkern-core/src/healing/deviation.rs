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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    fn x_axis_line() -> CurveGeometry {
        CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        }
    }

    fn y_axis_line() -> CurveGeometry {
        CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn identical_curves_zero_deviation() {
        let line = x_axis_line();
        let range = (0.0, 10.0);

        let result = curve_curve_deviation(&line, range, &line, range, 50);

        assert!(
            result.max_deviation < 1e-12,
            "identical curves should have ~0 deviation, got {}",
            result.max_deviation
        );
        assert!(
            result.mean_deviation < 1e-12,
            "mean deviation should be ~0 for identical curves, got {}",
            result.mean_deviation
        );
        // 50 uniform + 8 adaptive refinement samples.
        assert!(
            result.sample_count >= 50,
            "sample count should include refinement, got {}",
            result.sample_count
        );
    }

    #[test]
    fn curve_on_plane_zero_deviation() {
        let line = x_axis_line();
        let plane = xy_plane();

        let result = edge_surface_deviation(&line, (0.0, 5.0), &plane, 30);

        assert!(
            result.max_deviation < 1e-6,
            "line on XY plane should have ~0 deviation, got {}",
            result.max_deviation
        );
        assert!(
            result.mean_deviation < 1e-6,
            "mean deviation should be ~0 for line on plane, got {}",
            result.mean_deviation
        );
    }

    #[test]
    fn parallel_curves_known_deviation() {
        let line_a = x_axis_line();
        // Line at z=1, parallel to x_axis_line.
        let line_b = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 1.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let range = (0.0, 5.0);

        let result = curve_curve_deviation(&line_a, range, &line_b, range, 20);

        assert!(
            (result.max_deviation - 1.0).abs() < 1e-10,
            "parallel lines 1 unit apart should have max deviation ~1.0, got {}",
            result.max_deviation
        );
        assert!(
            (result.mean_deviation - 1.0).abs() < 1e-10,
            "mean deviation should be ~1.0, got {}",
            result.mean_deviation
        );
    }

    #[test]
    fn surface_surface_same_plane() {
        let plane = xy_plane();
        let range = ((-1.0, 1.0), (-1.0, 1.0));

        let result = surface_surface_deviation(&plane, range, &plane, range, 25);

        assert!(
            result.max_deviation < 1e-6,
            "same plane should have ~0 deviation, got {}",
            result.max_deviation
        );
        assert!(result.sample_count > 0);
    }

    #[test]
    fn deviation_result_fields() {
        let line_a = x_axis_line();
        let line_b = y_axis_line();

        let result = curve_curve_deviation(&line_a, (0.0, 1.0), &line_b, (0.0, 1.0), 10);

        // At t=0 both are at origin, so min deviation is 0.
        // At t=1 they are sqrt(2) apart.
        assert!(result.max_deviation > 0.0);
        // 10 uniform + 8 adaptive refinement samples.
        assert!(
            result.sample_count >= 10,
            "sample count should include refinement, got {}",
            result.sample_count
        );
    }
}
