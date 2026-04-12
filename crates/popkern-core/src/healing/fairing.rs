//! Curve fairing: smooth curves while preserving large-scale shape.
//!
//! Provides curve smoothing with bounded displacement and wrinkle removal.
//! Fairing uses Laplacian smoothing on NURBS control points to reduce
//! high-frequency oscillations while respecting a maximum displacement
//! constraint from the original control polygon.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;

/// Result of a curve fairing operation.
#[derive(Debug, Clone)]
pub struct FairingResult {
    /// Maximum displacement from the original curve.
    pub max_displacement: f64,
    /// Number of control points adjusted.
    pub points_adjusted: usize,
}

/// Fair a NURBS curve by smoothing its control points.
///
/// Applies Laplacian smoothing to interior control points, moving each
/// toward the average of its neighbors. The `strength` parameter (0.0 to 1.0)
/// controls how much smoothing is applied per iteration. The `max_displacement`
/// bounds how far any control point can move from its original position.
///
/// This is a simple but effective approach for removing small-scale
/// oscillations (wrinkles) while preserving large-scale shape.
///
/// # Arguments
///
/// * `curve` - The NURBS curve to smooth (modified in place).
/// * `iterations` - Number of smoothing passes.
/// * `strength` - Blending factor per iteration, clamped to `[0.0, 1.0]`.
/// * `max_displacement` - Maximum distance any control point may move from
///   its original position.
pub fn fair_nurbs_curve(
    curve: &mut NurbsCurve,
    iterations: usize,
    strength: f64,
    max_displacement: f64,
) -> FairingResult {
    let n = curve.control_points.len();

    // Fewer than 3 control points means no interior points to smooth.
    if n < 3 || iterations == 0 {
        return FairingResult {
            max_displacement: 0.0,
            points_adjusted: 0,
        };
    }

    let strength = strength.clamp(0.0, 1.0);
    let original = curve.control_points.clone();

    for _ in 0..iterations {
        // Snapshot current positions for this iteration so we don't
        // use partially-updated values within the same pass.
        let current = curve.control_points.clone();

        for i in 1..(n - 1) {
            let prev = current[i - 1];
            let next = current[i + 1];
            let avg = Point3::new(
                (prev.x + next.x) * 0.5,
                (prev.y + next.y) * 0.5,
                (prev.z + next.z) * 0.5,
            );

            // Move toward the neighbor average by the given strength.
            let target = Point3::new(
                current[i].x + (avg.x - current[i].x) * strength,
                current[i].y + (avg.y - current[i].y) * strength,
                current[i].z + (avg.z - current[i].z) * strength,
            );

            // Clamp: don't exceed max_displacement from the *original* position.
            let disp = original[i].distance_to(&target);
            if disp <= max_displacement {
                curve.control_points[i] = target;
            } else {
                // Scale back toward the original to stay within the bound.
                let ratio = max_displacement / disp;
                curve.control_points[i] = Point3::new(
                    original[i].x + (target.x - original[i].x) * ratio,
                    original[i].y + (target.y - original[i].y) * ratio,
                    original[i].z + (target.z - original[i].z) * ratio,
                );
            }
        }
    }

    // Compute final statistics.
    let mut actual_max_disp = 0.0_f64;
    let mut adjusted = 0_usize;
    let disp_threshold = 1e-15;

    for (orig, cp) in original
        .iter()
        .zip(curve.control_points.iter())
        .skip(1)
        .take(n - 2)
    {
        let d = orig.distance_to(cp);
        actual_max_disp = actual_max_disp.max(d);
        if d > disp_threshold {
            adjusted += 1;
        }
    }

    FairingResult {
        max_displacement: actual_max_disp,
        points_adjusted: adjusted,
    }
}

/// Detect wrinkles in a NURBS curve — high-frequency oscillations.
///
/// Samples the curve at `num_samples` points and computes the curvature
/// at each. Wrinkles appear as curvature values that are significantly
/// larger than the median curvature, indicating localized oscillation.
///
/// Returns a vec of `(parameter, curvature)` for points where curvature
/// exceeds `median_curvature * threshold_multiplier`.
///
/// # Arguments
///
/// * `curve` - Any curve implementing [`CurveEval`].
/// * `t_range` - Parameter range to sample over `(t_start, t_end)`.
/// * `num_samples` - Number of evenly-spaced sample points.
/// * `threshold_multiplier` - How many times the median curvature a value
///   must exceed to be flagged as a wrinkle.
pub fn detect_wrinkles(
    curve: &dyn CurveEval,
    t_range: (f64, f64),
    num_samples: usize,
    threshold_multiplier: f64,
) -> Vec<(f64, f64)> {
    use crate::healing::curvature::curve_curvature::curve_curvature_at;

    if num_samples < 2 {
        return Vec::new();
    }

    let (t0, t1) = t_range;

    // Sample curvature at evenly-spaced points.
    let mut samples: Vec<(f64, f64)> = Vec::with_capacity(num_samples);
    for i in 0..num_samples {
        let t = t0 + (t1 - t0) * (i as f64) / ((num_samples - 1) as f64);
        if let Some(kappa) = curve_curvature_at(curve, t) {
            samples.push((t, kappa));
        }
    }

    if samples.is_empty() {
        return Vec::new();
    }

    // Compute median curvature.
    let mut curvatures: Vec<f64> = samples.iter().map(|(_, k)| *k).collect();
    curvatures.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if curvatures.len().is_multiple_of(2) {
        let mid = curvatures.len() / 2;
        (curvatures[mid - 1] + curvatures[mid]) * 0.5
    } else {
        curvatures[curvatures.len() / 2]
    };

    let threshold = median * threshold_multiplier;

    // Flag points exceeding the threshold.
    samples
        .into_iter()
        .filter(|(_, kappa)| *kappa > threshold)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::{Point3, Vec3};

    /// Build a simple degree-2 NURBS representing a straight line from (0,0,0) to (1,0,0).
    fn line_nurbs() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    /// Build a 5-point NURBS where interior points lie on a line.
    fn straight_5pt_nurbs() -> NurbsCurve {
        NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.25, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.0),
                Point3::new(0.75, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
        }
    }

    #[test]
    fn fair_straight_line_no_change() {
        let mut curve = straight_5pt_nurbs();
        let original = curve.control_points.clone();
        let result = fair_nurbs_curve(&mut curve, 10, 0.5, 1.0);

        // A perfectly straight control polygon should not be changed by
        // Laplacian smoothing (the average of neighbors equals the point).
        assert!(
            result.max_displacement < 1e-12,
            "straight line should have ~zero displacement, got {}",
            result.max_displacement
        );
        assert_eq!(
            result.points_adjusted, 0,
            "no points should be adjusted on a straight line"
        );

        // Verify control points unchanged.
        for (i, (orig, faired)) in original.iter().zip(curve.control_points.iter()).enumerate() {
            assert!(
                orig.distance_to(faired) < 1e-12,
                "control point {i} should be unchanged"
            );
        }
    }

    #[test]
    fn fair_too_few_control_points_is_noop() {
        let mut curve = line_nurbs();
        let result = fair_nurbs_curve(&mut curve, 10, 0.5, 1.0);
        assert!(
            result.max_displacement < f64::EPSILON,
            "2-point curve has no interior points to smooth"
        );
        assert_eq!(result.points_adjusted, 0);
    }

    #[test]
    fn fair_zero_iterations_is_noop() {
        let mut curve = straight_5pt_nurbs();
        let result = fair_nurbs_curve(&mut curve, 0, 0.5, 1.0);
        assert!(result.max_displacement < f64::EPSILON);
        assert_eq!(result.points_adjusted, 0);
    }

    #[test]
    fn fair_respects_max_displacement() {
        // Create a curve with a deliberate wrinkle in the middle.
        let mut curve = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.25, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0), // wrinkle
                Point3::new(0.75, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
        };

        let max_disp = 0.1;
        let result = fair_nurbs_curve(&mut curve, 100, 0.5, max_disp);

        assert!(
            result.max_displacement <= max_disp + 1e-12,
            "displacement {} should not exceed max {}",
            result.max_displacement,
            max_disp
        );
    }

    #[test]
    fn detect_wrinkles_on_smooth_circle_finds_none() {
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let wrinkles = detect_wrinkles(&circle, (0.0, std::f64::consts::TAU), 50, 3.0);

        // A circle has constant curvature, so no point should exceed
        // 3x the median (which equals the constant curvature).
        assert!(
            wrinkles.is_empty(),
            "smooth circle should have no wrinkles, got {} flagged points",
            wrinkles.len()
        );
    }

    #[test]
    fn detect_wrinkles_on_line_finds_none() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let wrinkles = detect_wrinkles(&line, (0.0, 1.0), 20, 2.0);

        // A line has zero curvature everywhere; no wrinkles possible.
        assert!(wrinkles.is_empty(), "straight line should have no wrinkles");
    }

    #[test]
    fn detect_wrinkles_too_few_samples() {
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let wrinkles = detect_wrinkles(&circle, (0.0, 1.0), 1, 2.0);
        assert!(
            wrinkles.is_empty(),
            "fewer than 2 samples should return empty"
        );
    }

    #[test]
    fn fair_strength_clamped() {
        // Strength > 1.0 should be clamped to 1.0 and not panic.
        let mut curve = straight_5pt_nurbs();
        let result = fair_nurbs_curve(&mut curve, 5, 2.0, 1.0);
        assert!(result.max_displacement < 1e-12);

        // Strength < 0.0 should be clamped to 0.0 and produce no change.
        let mut curve2 = straight_5pt_nurbs();
        let result2 = fair_nurbs_curve(&mut curve2, 5, -1.0, 1.0);
        assert!(result2.max_displacement < 1e-12);
        assert_eq!(result2.points_adjusted, 0);
    }
}
