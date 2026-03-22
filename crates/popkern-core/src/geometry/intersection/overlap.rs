//! Curve overlap (coincident interval) detection.

#![allow(dead_code)]

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::CurveIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;

/// Number of samples on each curve for the coarse scan.
const N_SAMPLES: usize = 50;

/// Minimum consecutive in-tolerance samples to declare an overlap.
const MIN_CONSECUTIVE: usize = 3;

/// A sample point pairing between two curves.
struct Sample {
    t1: f64,
    t2: f64,
    dist: f64,
}

/// Find the parameter on `curve2` within `t2_range` that minimises the
/// distance to `target`, starting from a coarse scan then refining with a
/// few Newton steps on f(t) = |curve2(t) - target|².
fn closest_param_on_curve(
    curve2: &NurbsCurve,
    t2_range: (f64, f64),
    target: &crate::geometry::Point3,
) -> (f64, f64) {
    let (t2_lo, t2_hi) = t2_range;

    // Coarse scan — find the best sample.
    let mut best_t = t2_lo;
    let mut best_dist2 = f64::MAX;

    for k in 0..N_SAMPLES {
        let t = t2_lo + (t2_hi - t2_lo) * (k as f64) / (N_SAMPLES as f64 - 1.0);
        let p = curve2.point_at(t);
        let dx = p.x - target.x;
        let dy = p.y - target.y;
        let dz = p.z - target.z;
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < best_dist2 {
            best_dist2 = d2;
            best_t = t;
        }
    }

    // Newton refinement: minimise f(t) = |curve2(t) - target|²
    // f'(t) = 2 (curve2(t) - target) · curve2'(t)
    // f''(t) ≈ 2 |curve2'(t)|²  (dominant term)
    let mut t = best_t;
    for _ in 0..8 {
        let p = curve2.point_at(t);
        let dp = curve2.derivative_at(t, 1);
        let diff_x = p.x - target.x;
        let diff_y = p.y - target.y;
        let diff_z = p.z - target.z;
        let f_prime = 2.0 * (diff_x * dp.x + diff_y * dp.y + diff_z * dp.z);
        let f_double = 2.0 * (dp.x * dp.x + dp.y * dp.y + dp.z * dp.z);
        if f_double.abs() < 1e-30 {
            break;
        }
        let step = f_prime / f_double;
        t = (t - step).clamp(t2_lo, t2_hi);
    }

    let p = curve2.point_at(t);
    let dist = target.distance_to(&p);
    (t, dist)
}

/// Detect overlapping segments between two NURBS curves.
/// Returns a list of [`CurveIntersection::Overlap`] entries.
pub(crate) fn detect_overlaps(
    curve1: &NurbsCurve,
    t1_range: (f64, f64),
    curve2: &NurbsCurve,
    t2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let (t1_lo, t1_hi) = t1_range;

    // Sample curve1 and find closest point on curve2 for each sample.
    let mut samples: Vec<Sample> = Vec::with_capacity(N_SAMPLES);
    for i in 0..N_SAMPLES {
        let t1 = t1_lo + (t1_hi - t1_lo) * (i as f64) / (N_SAMPLES as f64 - 1.0);
        let p1 = curve1.point_at(t1);
        let (t2, dist) = closest_param_on_curve(curve2, t2_range, &p1);
        samples.push(Sample { t1, t2, dist });
    }

    // Collect contiguous runs where dist < tolerance.linear.
    let tol = tolerance.linear;
    let mut overlaps: Vec<CurveIntersection> = Vec::new();

    let mut run_start: Option<usize> = None;

    for i in 0..samples.len() {
        let in_tol = samples[i].dist < tol;
        match (in_tol, run_start) {
            (true, None) => run_start = Some(i),
            (false, Some(start)) => {
                let end = i - 1;
                if end - start + 1 >= MIN_CONSECUTIVE {
                    overlaps.push(build_overlap(&samples[start..=end]));
                }
                run_start = None;
            }
            _ => {}
        }
    }
    // Handle run that reaches the end.
    if let Some(start) = run_start {
        let end = samples.len() - 1;
        if end - start + 1 >= MIN_CONSECUTIVE {
            overlaps.push(build_overlap(&samples[start..=end]));
        }
    }

    overlaps
}

/// Build a [`CurveIntersection::Overlap`] from a contiguous slice of in-tolerance samples.
fn build_overlap(run: &[Sample]) -> CurveIntersection {
    let t1_start = run.first().map(|s| s.t1).unwrap_or(0.0);
    let t1_end = run.last().map(|s| s.t1).unwrap_or(0.0);
    let t2_start = run.first().map(|s| s.t2).unwrap_or(0.0);
    let t2_end = run.last().map(|s| s.t2).unwrap_or(0.0);
    CurveIntersection::Overlap {
        t1_range: (t1_start, t1_end),
        t2_range: (t2_start, t2_end),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;

    #[test]
    fn detect_full_overlap() {
        let c1 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let c2 = c1.clone();
        let overlaps = detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
        assert_eq!(overlaps.len(), 1);
        if let CurveIntersection::Overlap {
            t1_range,
            t2_range: _,
        } = &overlaps[0]
        {
            assert!((t1_range.0 - 0.0).abs() < 0.1);
            assert!((t1_range.1 - 1.0).abs() < 0.1);
        } else {
            panic!("Expected Overlap");
        }
    }

    #[test]
    fn detect_partial_overlap() {
        let c1 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let c2 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let overlaps = detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
        assert_eq!(overlaps.len(), 1);
    }

    #[test]
    fn no_overlap_different_curves() {
        let c1 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let c2 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 5.0, 0.0), Point3::new(1.0, 5.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let overlaps = detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
        assert!(overlaps.is_empty());
    }
}
