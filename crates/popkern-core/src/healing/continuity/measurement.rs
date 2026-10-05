//! Cross-boundary continuity measurement between two surfaces.
//!
//! Samples corresponding parameter pairs along a shared boundary and
//! computes positional (G0), normal-angle (G1), and curvature (G2)
//! deviations to determine the achieved continuity level.

use crate::geometry::surface::SurfaceEval;
use crate::healing::curvature::surface_curvature::gaussian_curvature;
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// A pair of corresponding UV parameters on two surfaces along a shared boundary.
#[derive(Debug, Clone, Copy)]
pub struct BoundarySample {
    /// U parameter on surface A.
    pub u_a: f64,
    /// V parameter on surface A.
    pub v_a: f64,
    /// U parameter on surface B.
    pub u_b: f64,
    /// V parameter on surface B.
    pub v_b: f64,
}

/// Result of cross-boundary continuity measurement.
///
/// Contains per-level statistics (G0 position, G1 normal angle, G2 curvature,
/// G3 curvature rate) and the highest achieved continuity level.
#[derive(Debug, Clone)]
pub struct ContinuityMeasurement {
    /// Highest continuity level that passes the tolerance thresholds.
    pub achieved_level: ContinuityLevel,
    /// Maximum positional deviation (model units).
    pub g0_max: f64,
    /// Mean positional deviation.
    pub g0_mean: f64,
    /// 95th percentile positional deviation.
    pub g0_p95: f64,
    /// Median positional deviation.
    pub g0_median: f64,
    /// Maximum normal angle deviation (radians).
    pub g1_max_angle: f64,
    /// Mean normal angle deviation (radians).
    pub g1_mean_angle: f64,
    /// 95th percentile normal angle deviation (radians).
    pub g1_p95_angle: f64,
    /// Median normal angle deviation (radians).
    pub g1_median_angle: f64,
    /// Maximum Gaussian curvature difference.
    pub g2_max_curvature_diff: f64,
    /// Mean Gaussian curvature difference.
    pub g2_mean_curvature_diff: f64,
    /// Maximum curvature-rate difference (derivative of curvature w.r.t. arc length).
    pub g3_max_curvature_rate_diff: f64,
    /// Mean curvature-rate difference.
    pub g3_mean_curvature_rate_diff: f64,
    /// Number of boundary samples evaluated.
    pub sample_count: usize,
}

/// Curvature-difference threshold for G2 continuity.
///
/// If the maximum absolute difference in Gaussian curvature across all
/// samples is below this value, the boundary is considered G2.
const G2_CURVATURE_THRESHOLD: f64 = 1e-3;

/// Curvature-rate-difference threshold for G3 continuity.
///
/// Curvature rate matching is inherently loose, so the default is 1.0.
const G3_CURVATURE_RATE_THRESHOLD: f64 = 1.0;

/// Compute the 95th-percentile value from a sorted slice.
///
/// Returns 0.0 for empty input.
fn percentile_95(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() as f64) * 0.95).ceil() as usize;
    let idx = idx.min(sorted.len() - 1);
    sorted[idx]
}

/// Compute the median value from a sorted slice.
///
/// Returns 0.0 for empty input. For even-length slices, returns the
/// average of the two middle values.
fn median(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) * 0.5
    }
}

/// Measure cross-boundary continuity between two surfaces along a shared boundary.
///
/// For each [`BoundarySample`] the function evaluates:
///
/// - **G0** (positional): Euclidean distance between `surface_a.point_at(u_a, v_a)`
///   and `surface_b.point_at(u_b, v_b)`.
/// - **G1** (tangent-plane): angle between surface normals via `acos(|n_a . n_b|)`.
///   Samples where either normal has near-zero length are skipped.
/// - **G2** (curvature): absolute difference in Gaussian curvature via
///   [`gaussian_curvature`]. Samples where curvature is undefined (degenerate
///   points) are skipped.
///
/// - **G3** (curvature rate): difference in the derivative of Gaussian curvature
///   w.r.t. arc length, estimated via finite differences between adjacent samples.
///
/// The achieved level is determined by comparing statistics against `tolerance`:
/// - G3 if `g3_max < G3_CURVATURE_RATE_THRESHOLD` **and** G2 passes
/// - G2 if `g2_max_curvature_diff < G2_CURVATURE_THRESHOLD` **and** G1 passes
/// - G1 if `g1_max_angle < tolerance.angular` **and** G0 passes
/// - G0 if `g0_max < tolerance.linear`
pub fn cross_boundary_continuity(
    surface_a: &dyn SurfaceEval,
    surface_b: &dyn SurfaceEval,
    samples: &[BoundarySample],
    tolerance: &Tolerance,
) -> ContinuityMeasurement {
    let mut g0_values = Vec::with_capacity(samples.len());
    let mut g1_values = Vec::with_capacity(samples.len());
    let mut g2_values = Vec::with_capacity(samples.len());

    for sample in samples {
        // G0: positional deviation
        let pa = surface_a.point_at(sample.u_a, sample.v_a);
        let pb = surface_b.point_at(sample.u_b, sample.v_b);
        g0_values.push(pa.distance_to(&pb));

        // G1: normal angle deviation
        let na = surface_a.normal_at(sample.u_a, sample.v_a);
        let nb = surface_b.normal_at(sample.u_b, sample.v_b);
        let len_a = na.length();
        let len_b = nb.length();
        if len_a > 1e-12 && len_b > 1e-12 {
            let cos_angle = na.dot(&nb) / (len_a * len_b);
            // Use |cos| because normals may point in opposite directions
            // (opposite face orientations).
            let angle = cos_angle.abs().clamp(0.0, 1.0).pacos();
            g1_values.push(angle);
        }

        // G2: Gaussian curvature difference
        let ka = gaussian_curvature(surface_a, sample.u_a, sample.v_a);
        let kb = gaussian_curvature(surface_b, sample.u_b, sample.v_b);
        if let (Some(ka), Some(kb)) = (ka, kb) {
            g2_values.push((ka - kb).abs());
        }
    }

    // Compute G0 statistics.
    let g0_max = g0_values.iter().copied().fold(0.0_f64, f64::max);
    let g0_mean = if g0_values.is_empty() {
        0.0
    } else {
        g0_values.iter().sum::<f64>() / g0_values.len() as f64
    };
    let mut g0_sorted = g0_values.clone();
    g0_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let g0_p95 = percentile_95(&g0_sorted);
    let g0_median = median(&g0_sorted);

    // Compute G1 statistics.
    let g1_max_angle = g1_values.iter().copied().fold(0.0_f64, f64::max);
    let g1_mean_angle = if g1_values.is_empty() {
        0.0
    } else {
        g1_values.iter().sum::<f64>() / g1_values.len() as f64
    };
    let mut g1_sorted = g1_values.clone();
    g1_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let g1_p95_angle = percentile_95(&g1_sorted);
    let g1_median_angle = median(&g1_sorted);

    // Compute G2 statistics.
    let g2_max_curvature_diff = g2_values.iter().copied().fold(0.0_f64, f64::max);
    let g2_mean_curvature_diff = if g2_values.is_empty() {
        0.0
    } else {
        g2_values.iter().sum::<f64>() / g2_values.len() as f64
    };

    // Compute G3: curvature rate (d(curvature)/d(arc_length)) difference.
    let mut g3_values = Vec::new();
    if g2_values.len() >= 2 {
        let mut ka_vals = Vec::with_capacity(samples.len());
        let mut kb_vals = Vec::with_capacity(samples.len());
        let mut pts_a = Vec::with_capacity(samples.len());
        for s in samples {
            ka_vals.push(gaussian_curvature(surface_a, s.u_a, s.v_a));
            kb_vals.push(gaussian_curvature(surface_b, s.u_b, s.v_b));
            pts_a.push(surface_a.point_at(s.u_a, s.v_a));
        }
        for i in 0..ka_vals.len() - 1 {
            if let (Some(ka_i), Some(ka_n)) = (ka_vals[i], ka_vals[i + 1])
                && let (Some(kb_i), Some(kb_n)) = (kb_vals[i], kb_vals[i + 1])
            {
                let ds = pts_a[i].distance_to(&pts_a[i + 1]);
                if ds > 1e-15 {
                    let rate_a = (ka_n - ka_i) / ds;
                    let rate_b = (kb_n - kb_i) / ds;
                    g3_values.push((rate_a - rate_b).abs());
                }
            }
        }
    }
    let g3_max_curvature_rate_diff = g3_values.iter().copied().fold(0.0_f64, f64::max);
    let g3_mean_curvature_rate_diff = if g3_values.is_empty() {
        0.0
    } else {
        g3_values.iter().sum::<f64>() / g3_values.len() as f64
    };

    // Determine achieved level (each level subsumes the previous).
    let passes_g0 = g0_max < tolerance.linear();
    let passes_g1 = passes_g0 && g1_max_angle < tolerance.angular();
    let passes_g2 = passes_g1 && g2_max_curvature_diff < G2_CURVATURE_THRESHOLD;
    let passes_g3 = passes_g2 && g3_max_curvature_rate_diff < G3_CURVATURE_RATE_THRESHOLD;

    let achieved_level = if passes_g3 {
        ContinuityLevel::G3
    } else if passes_g2 {
        ContinuityLevel::G2
    } else if passes_g1 {
        ContinuityLevel::G1
    } else {
        // G0 is the base level (no "below G0" variant).
        ContinuityLevel::G0
    };

    ContinuityMeasurement {
        achieved_level,
        g0_max,
        g0_mean,
        g0_p95,
        g0_median,
        g1_max_angle,
        g1_mean_angle,
        g1_p95_angle,
        g1_median_angle,
        g2_max_curvature_diff,
        g2_mean_curvature_diff,
        g3_max_curvature_rate_diff,
        g3_mean_curvature_rate_diff,
        sample_count: samples.len(),
    }
}

/// A single detected internal discontinuity on a surface.
#[derive(Debug, Clone)]
pub struct DiscontinuityRecord {
    /// Parameter location of the discontinuity.
    pub u: f64,
    /// Parameter location of the discontinuity.
    pub v: f64,
    /// The highest continuity level achieved at this location.
    pub achieved_level: ContinuityLevel,
    /// Magnitude of the deviation at the target level.
    pub deviation: f64,
}

/// Detect internal discontinuities on a single surface.
///
/// Samples the surface on a grid and evaluates continuity between adjacent
/// sample points. Reports locations where continuity drops below the
/// `target_level`. For NURBS surfaces this catches knot-line discontinuities;
/// for analytical surfaces it should report nothing.
///
/// This implements Section 2.2(a) of the healing requirements.
pub fn internal_discontinuities(
    surface: &dyn SurfaceEval,
    target_level: ContinuityLevel,
    tolerance: &Tolerance,
    sample_density: usize,
) -> Vec<DiscontinuityRecord> {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let u_lo = u0.max(-100.0);
    let u_hi = u1.min(100.0);
    let v_lo = v0.max(-100.0);
    let v_hi = v1.min(100.0);

    let n = sample_density.max(4);
    let mut records = Vec::new();
    let h_u = (u_hi - u_lo) / (n as f64);
    let h_v = (v_hi - v_lo) / (n as f64);

    // Check continuity in the u-direction along each v-row.
    for j in 0..=n {
        let v = v_lo + (j as f64) * h_v;
        for i in 0..n {
            let u_left = u_lo + (i as f64) * h_u;
            let u_right = u_left + h_u;
            let u_mid = (u_left + u_right) * 0.5;
            let pts = ((u_left, v), (u_right, v), (u_mid, v));
            if let Some(rec) = check_adjacent_continuity(surface, pts, target_level, tolerance) {
                records.push(rec);
            }
        }
    }

    // Check continuity in the v-direction along each u-column.
    for i in 0..=n {
        let u = u_lo + (i as f64) * h_u;
        for j in 0..n {
            let v_bot = v_lo + (j as f64) * h_v;
            let v_top = v_bot + h_v;
            let v_mid = (v_bot + v_top) * 0.5;
            let pts = ((u, v_bot), (u, v_top), (u, v_mid));
            if let Some(rec) = check_adjacent_continuity(surface, pts, target_level, tolerance) {
                records.push(rec);
            }
        }
    }

    records
}

/// Check continuity between two adjacent sample points on the same surface.
///
/// `pts` is `((u1,v1), (u2,v2), (u_mid,v_mid))`. Evaluates normals and
/// curvature at the endpoints. If the values are inconsistent, a
/// discontinuity is reported at the midpoint location.
fn check_adjacent_continuity(
    surface: &dyn SurfaceEval,
    pts: ((f64, f64), (f64, f64), (f64, f64)),
    target_level: ContinuityLevel,
    tolerance: &Tolerance,
) -> Option<DiscontinuityRecord> {
    let ((u1, v1), (u2, v2), (u_mid, v_mid)) = pts;
    let n1 = surface.normal_at(u1, v1);
    let n2 = surface.normal_at(u2, v2);
    let n1_len = n1.length();
    let n2_len = n2.length();

    // Skip singular points.
    if n1_len < 1e-15 || n2_len < 1e-15 {
        return None;
    }

    // G1 check: angle between normals at adjacent points.
    let n1_unit = n1 * (1.0 / n1_len);
    let n2_unit = n2 * (1.0 / n2_len);
    let dot = n1_unit.dot(&n2_unit).clamp(-1.0, 1.0);
    let angle = dot.abs().pacos();

    if target_level >= ContinuityLevel::G1 && angle > tolerance.angular() {
        return Some(DiscontinuityRecord {
            u: u_mid,
            v: v_mid,
            achieved_level: ContinuityLevel::G0,
            deviation: angle,
        });
    }

    // G2 check: curvature difference between adjacent points.
    if target_level >= ContinuityLevel::G2 {
        let k1 = gaussian_curvature(surface, u1, v1);
        let k2 = gaussian_curvature(surface, u2, v2);
        if let (Some(kv1), Some(kv2)) = (k1, k2) {
            let curv_diff = (kv1 - kv2).abs();
            let threshold = 0.1; // curvature difference threshold
            if curv_diff > threshold {
                return Some(DiscontinuityRecord {
                    u: u_mid,
                    v: v_mid,
                    achieved_level: ContinuityLevel::G1,
                    deviation: curv_diff,
                });
            }
        }
    }

    None
}

#[cfg(test)]
#[path = "measurement_tests.rs"]
mod tests;
