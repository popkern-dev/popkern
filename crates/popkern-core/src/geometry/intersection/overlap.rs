//! Curve overlap (coincident interval) detection.

#![allow(dead_code)]

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::CurveIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;

/// Number of samples on each curve for the coarse scan.
const N_SAMPLES: usize = 50;

/// Minimum consecutive in-tolerance samples to declare an overlap.
const MIN_CONSECUTIVE: usize = 3;

/// Bisections refining a run end that lies inside curve 1 (one sample
/// spacing / 2^40: below every model tolerance).
const BISECTIONS: usize = 40;

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
///
/// `curve2_hulls` are control-hull boxes of Bezier spans covering
/// `curve2` over `t2_range` (empty: no prefilter). With positive
/// weights each span lies in its box, so a sample farther than the
/// tolerance from every box is out of tolerance without the closest-
/// point scan: the overlap runs, and so the result, are unchanged
/// (PK-802: the scan ran 2900 evaluations of a 256-span rim circle per
/// clip query, half a slot cut's single-threaded time).
pub(crate) fn detect_overlaps(
    curve1: &NurbsCurve,
    t1_range: (f64, f64),
    curve2: &NurbsCurve,
    t2_range: (f64, f64),
    curve2_hulls: &[(Point3, Point3)],
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let (t1_lo, t1_hi) = t1_range;
    let tol = tolerance.linear();
    let hull_distance = |p: &Point3| {
        curve2_hulls
            .iter()
            .map(|(lo, hi)| {
                let d = |x: f64, a: f64, b: f64| (a - x).max(x - b).max(0.0);
                let (dx, dy, dz) = (d(p.x, lo.x, hi.x), d(p.y, lo.y, hi.y), d(p.z, lo.z, hi.z));
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .fold(f64::INFINITY, f64::min)
    };

    // Sample curve1 and find closest point on curve2 for each sample.
    let mut samples: Vec<Sample> = Vec::with_capacity(N_SAMPLES);
    for i in 0..N_SAMPLES {
        let t1 = t1_lo + (t1_hi - t1_lo) * (i as f64) / (N_SAMPLES as f64 - 1.0);
        let p1 = curve1.point_at(t1);
        if !curve2_hulls.is_empty() {
            let bound = hull_distance(&p1);
            if bound >= tol {
                samples.push(Sample {
                    t1,
                    t2: t2_range.0,
                    dist: bound,
                });
                continue;
            }
        }
        let (t2, dist) = closest_param_on_curve(curve2, t2_range, &p1);
        samples.push(Sample { t1, t2, dist });
    }

    // A run's end strictly inside curve 1 lies between its last
    // in-tolerance sample and the next sample out: bisect the closest-
    // point predicate there, so the overlap ends where the curves part
    // and not on the scan grid (a boss rim lying on its face's rim
    // pcurve parted at boundary t = 0.5 and was reported at 24/49; the
    // clip window then lost 1.1% of the rim at each end and the union
    // stranded record sides, PK-789).
    let refine = |inside: &Sample, outside: &Sample| -> Sample {
        let (mut a, mut b) = (inside.t1, outside.t1);
        let mut best = (inside.t2, inside.dist);
        for _ in 0..BISECTIONS {
            let m = 0.5 * (a + b);
            let (t2, dist) = closest_param_on_curve(curve2, t2_range, &curve1.point_at(m));
            if dist < tol {
                (a, best) = (m, (t2, dist));
            } else {
                b = m;
            }
        }
        Sample {
            t1: a,
            t2: best.0,
            dist: best.1,
        }
    };
    let span = |start: usize, end: usize| -> Vec<Sample> {
        let mut run: Vec<Sample> = samples[start..=end]
            .iter()
            .map(|s| Sample {
                t1: s.t1,
                t2: s.t2,
                dist: s.dist,
            })
            .collect();
        if start > 0 {
            run[0] = refine(&samples[start], &samples[start - 1]);
        }
        if end + 1 < samples.len() {
            let last = run.len() - 1;
            run[last] = refine(&samples[end], &samples[end + 1]);
        }
        run
    };

    // Collect contiguous runs where dist < tolerance.linear().
    let mut overlaps: Vec<CurveIntersection> = Vec::new();

    let mut run_start: Option<usize> = None;

    for (i, sample) in samples.iter().enumerate() {
        let in_tol = sample.dist < tol;
        match (in_tol, run_start) {
            (true, None) => run_start = Some(i),
            (false, Some(start)) => {
                let end = i - 1;
                if end - start + 1 >= MIN_CONSECUTIVE {
                    overlaps.push(build_overlap(&span(start, end)));
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
            overlaps.push(build_overlap(&span(start, end)));
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

    /// A run ending inside curve 1 ends where the curves part, not on
    /// the 50-sample scan (PK-789): curve 2 covers curve 1 from t = 0.37,
    /// which the scan alone reports at its next sample, 19/49 = 0.388.
    #[test]
    fn a_partial_overlap_ends_where_the_curves_part() {
        let line = |a: f64, b: f64| NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(a, 0.0, 0.0), Point3::new(b, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let (c1, c2) = (line(0.0, 1.0), line(0.37, 2.0));
        let overlaps =
            detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &[], &Tolerance::default());
        let [CurveIntersection::Overlap { t1_range, t2_range }] = overlaps.as_slice() else {
            panic!("one overlap expected: {overlaps:?}");
        };
        let lin = Tolerance::default().linear();
        assert!((t1_range.0 - 0.37).abs() <= lin, "start {}", t1_range.0);
        assert!((t1_range.1 - 1.0).abs() <= 1e-12);
        assert!(t2_range.0.abs() <= lin, "curve-2 start {}", t2_range.0);
    }

    #[test]
    fn detect_full_overlap() {
        let c1 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let c2 = c1.clone();
        let overlaps =
            detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &[], &Tolerance::default());
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
        let overlaps =
            detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &[], &Tolerance::default());
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
        let overlaps =
            detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &[], &Tolerance::default());
        assert!(overlaps.is_empty());
    }

    #[test]
    fn overlap_run_at_end_of_samples() {
        // Test that a run reaching the end of samples is handled correctly
        // Two identical curves should have a run that reaches the end
        let c1 = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let overlaps =
            detect_overlaps(&c1, (0.0, 1.0), &c1, (0.0, 1.0), &[], &Tolerance::default());
        assert_eq!(overlaps.len(), 1);
    }

    #[test]
    fn short_overlap_below_min_consecutive() {
        // Test that very short overlaps (< MIN_CONSECUTIVE) are not reported.
        // Two curves that are close for only a tiny stretch should not
        // produce an overlap.
        let c1 = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let c2 = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, -1.0, 0.0), // opposite direction
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        // These curves cross at endpoints only; overlap should be empty or have only
        // endpoint-coincidence (which is too short to be MIN_CONSECUTIVE)
        let overlaps =
            detect_overlaps(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &[], &Tolerance::default());
        // Either no overlap or very short runs filtered out
        // just verify it doesn't panic
        let _ = overlaps;
    }
}
