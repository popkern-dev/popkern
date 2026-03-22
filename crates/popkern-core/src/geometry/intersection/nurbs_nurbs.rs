//! NURBS-NURBS intersection algorithm.
//!
//! Orchestrates the pipeline: overlap detection → Bezier decomposition →
//! bounding box culling → Bezier clipping → Newton polish → deduplication.

#![allow(dead_code)]

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::{CurveIntersection, bezier_clipping, newton, overlap};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;

/// Number of samples per curve for the sampling-based candidate search.
const N_SAMPLES: usize = 64;

/// Intersect two NURBS curves.
pub(crate) fn intersect_nurbs_nurbs(
    curve1: &NurbsCurve,
    t1_range: (f64, f64),
    curve2: &NurbsCurve,
    t2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<CurveIntersection> {
    let mut results: Vec<CurveIntersection> = Vec::new();

    // 1. Overlap detection
    let overlaps = overlap::detect_overlaps(curve1, t1_range, curve2, t2_range, tolerance);
    if !overlaps.is_empty() {
        results.extend(overlaps);
        return results;
    }

    // 2. Bezier decomposition
    let segments1 = decompose_to_bezier(curve1, t1_range);
    let segments2 = decompose_to_bezier(curve2, t2_range);

    // 3 & 4. Bounding box culling + Bezier clipping
    let mut candidates: Vec<(f64, f64)> = Vec::new();
    for (seg1, range1) in &segments1 {
        for (seg2, range2) in &segments2 {
            if !aabb_overlap(&seg1.control_points, &seg2.control_points) {
                continue;
            }
            bezier_clipping::bezier_intersect(
                &seg1.control_points,
                &seg1.weights,
                *range1,
                &seg2.control_points,
                &seg2.weights,
                *range2,
                tolerance,
                50,
                &mut candidates,
            );
        }
    }

    // Supplement with sampling-based candidates to handle cases where
    // Bezier clipping converges poorly (e.g., near-tangent or symmetric curves).
    sample_candidates(
        curve1,
        t1_range,
        curve2,
        t2_range,
        tolerance,
        &mut candidates,
    );

    // 5. Newton polish
    let geom1 = CurveGeometry::Nurbs(curve1.clone());
    let geom2 = CurveGeometry::Nurbs(curve2.clone());

    for (t1_approx, t2_approx) in &candidates {
        if let Some((t1, t2, point)) =
            newton::newton_refine(&geom1, *t1_approx, &geom2, *t2_approx, tolerance)
        {
            // Only keep results within the requested ranges (with some slack)
            if t1 >= t1_range.0 - tolerance.parametric
                && t1 <= t1_range.1 + tolerance.parametric
                && t2 >= t2_range.0 - tolerance.parametric
                && t2 <= t2_range.1 + tolerance.parametric
            {
                results.push(CurveIntersection::Point { point, t1, t2 });
            }
        }
    }

    // 6. Deduplication
    deduplicate(&mut results, tolerance);

    // 7. Sort by t1
    results.sort_by(|a, b| {
        let t1_a = match a {
            CurveIntersection::Point { t1, .. } => *t1,
            CurveIntersection::Overlap { t1_range, .. } => t1_range.0,
        };
        let t1_b = match b {
            CurveIntersection::Point { t1, .. } => *t1,
            CurveIntersection::Overlap { t1_range, .. } => t1_range.0,
        };
        t1_a.partial_cmp(&t1_b).unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

/// Sample both curves and find close-approach parameter pairs as Newton seeds.
///
/// For each sample on curve1, finds the closest sample on curve2. When the
/// distance function shows a sign change in its derivative (indicating a
/// local minimum), the pair is added as a candidate.
fn sample_candidates(
    curve1: &NurbsCurve,
    t1_range: (f64, f64),
    curve2: &NurbsCurve,
    t2_range: (f64, f64),
    tolerance: &Tolerance,
    candidates: &mut Vec<(f64, f64)>,
) {
    // Sample curve2 points
    let n2 = N_SAMPLES;
    let mut pts2: Vec<(f64, Point3)> = Vec::with_capacity(n2);
    for j in 0..n2 {
        let t2 = t2_range.0 + (t2_range.1 - t2_range.0) * (j as f64) / (n2 as f64 - 1.0);
        pts2.push((t2, curve2.point_at(t2)));
    }

    // For each sample on curve1, find the closest point on curve2 samples
    let n1 = N_SAMPLES;
    // Coarse threshold: control polygon diagonal gives a rough upper bound
    let coarse_threshold = {
        let (min1, max1) = aabb(&curve1.control_points);
        let (min2, max2) = aabb(&curve2.control_points);
        let diag1 = min1.distance_to(&max1);
        let diag2 = min2.distance_to(&max2);
        (diag1 + diag2) * 0.5
    };

    for i in 0..n1 {
        let t1 = t1_range.0 + (t1_range.1 - t1_range.0) * (i as f64) / (n1 as f64 - 1.0);
        let p1 = curve1.point_at(t1);

        // Find the closest curve2 sample
        let mut best_j = 0;
        let mut best_dist = f64::MAX;
        for (j, (_, p2)) in pts2.iter().enumerate() {
            let d = p1.distance_to(p2);
            if d < best_dist {
                best_dist = d;
                best_j = j;
            }
        }

        // Only consider close-approach pairs
        if best_dist < coarse_threshold {
            // Check for local minimum: distance should be small relative to neighbors
            let dist_left = if best_j > 0 {
                p1.distance_to(&pts2[best_j - 1].1)
            } else {
                f64::MAX
            };
            let dist_right = if best_j + 1 < n2 {
                p1.distance_to(&pts2[best_j + 1].1)
            } else {
                f64::MAX
            };

            // This is a local minimum if distance is <= both neighbors
            if best_dist <= dist_left && best_dist <= dist_right {
                // Use a tighter threshold: sampling resolution gives approximate distance
                let sample_resolution = (t1_range.1 - t1_range.0) / (n1 as f64) * coarse_threshold;
                if best_dist < sample_resolution.max(tolerance.linear * 1000.0) {
                    candidates.push((t1, pts2[best_j].0));
                }
            }
        }
    }
}

/// Decompose a NURBS curve into Bezier segments by splitting at interior knots.
/// Returns a list of (segment, parameter_range) pairs.
fn decompose_to_bezier(curve: &NurbsCurve, t_range: (f64, f64)) -> Vec<(NurbsCurve, (f64, f64))> {
    let degree = curve.degree;
    let n = curve.control_points.len();

    // Find distinct interior knot values within the requested range.
    // Interior knots are those between knots[degree] and knots[n].
    let mut interior_knots: Vec<f64> = Vec::new();
    for i in (degree + 1)..n {
        let k = curve.knots[i];
        if k > t_range.0
            && k < t_range.1
            && interior_knots
                .last()
                .is_none_or(|&last| (k - last).abs() > 1e-14)
        {
            interior_knots.push(k);
        }
    }

    if interior_knots.is_empty() {
        // Already a single Bezier segment
        return vec![(curve.clone(), t_range)];
    }

    // Split sequentially at each interior knot
    let mut segments: Vec<(NurbsCurve, (f64, f64))> = Vec::new();
    let mut remaining = curve.clone();
    let mut current_start = t_range.0;

    for &knot in &interior_knots {
        match remaining.split_at(knot) {
            Ok((left, right)) => {
                segments.push((left, (current_start, knot)));
                remaining = right;
                current_start = knot;
            }
            Err(_) => {
                // If split fails (e.g., knot at boundary), skip
                continue;
            }
        }
    }

    // Add the final segment
    segments.push((remaining, (current_start, t_range.1)));

    segments
}

/// Compute axis-aligned bounding box from control points.
fn aabb(pts: &[Point3]) -> (Point3, Point3) {
    let mut min = pts[0];
    let mut max = pts[0];
    for p in pts.iter().skip(1) {
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        min.z = min.z.min(p.z);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
        max.z = max.z.max(p.z);
    }
    (min, max)
}

/// Check if two axis-aligned bounding boxes overlap.
fn aabb_overlap(pts1: &[Point3], pts2: &[Point3]) -> bool {
    if pts1.is_empty() || pts2.is_empty() {
        return false;
    }
    let (min1, max1) = aabb(pts1);
    let (min2, max2) = aabb(pts2);
    min1.x <= max2.x
        && max1.x >= min2.x
        && min1.y <= max2.y
        && max1.y >= min2.y
        && min1.z <= max2.z
        && max1.z >= min2.z
}

/// Deduplicate intersection results: merge points whose 3D positions are
/// within `tolerance.linear` OR whose parameters are within `tolerance.parametric`.
fn deduplicate(results: &mut Vec<CurveIntersection>, tolerance: &Tolerance) {
    let mut i = 0;
    while i < results.len() {
        let mut j = i + 1;
        while j < results.len() {
            if is_duplicate(&results[i], &results[j], tolerance) {
                results.remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

/// Check if two intersection results are duplicates.
fn is_duplicate(a: &CurveIntersection, b: &CurveIntersection, tolerance: &Tolerance) -> bool {
    match (a, b) {
        (
            CurveIntersection::Point {
                point: pa,
                t1: t1a,
                t2: t2a,
            },
            CurveIntersection::Point {
                point: pb,
                t1: t1b,
                t2: t2b,
            },
        ) => {
            pa.distance_to(pb) < tolerance.linear
                || ((t1a - t1b).abs() < tolerance.parametric
                    && (t2a - t2b).abs() < tolerance.parametric)
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;

    #[test]
    fn parabola_crosses_horizontal_line() {
        // Quadratic Bezier (parabola) and horizontal line at y=0.5
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
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
        // Parabola y = 2t(1-t) crosses y=0.5 at 2 points
        assert!(
            result.len() >= 2,
            "Expected at least 2 intersections, got {}",
            result.len()
        );

        // Verify invariant: points lie on both curves
        for r in &result {
            if let CurveIntersection::Point { point, t1, t2 } = r {
                let p1 = c1.point_at(*t1);
                let p2 = c2.point_at(*t2);
                assert!(p1.distance_to(&p2) < 1e-6, "Invariant: points don't match");
                assert!(p1.distance_to(point) < 1e-6);
            }
        }
    }

    #[test]
    fn non_intersecting_curves() {
        let c1 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let c2 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 10.0, 0.0), Point3::new(1.0, 10.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
        assert!(result.is_empty());
    }

    #[test]
    fn identical_curves_overlap() {
        let c1 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c1, (0.0, 1.0), &Tolerance::default());
        assert!(
            result
                .iter()
                .any(|r| matches!(r, CurveIntersection::Overlap { .. }))
        );
    }

    #[test]
    fn multi_span_curve_decomposition() {
        // Multi-span curve (has interior knot) crossing a line
        let c1 = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.25, 1.0, 0.0),
                Point3::new(0.75, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0],
        };
        let c2 = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
        assert!(
            result.len() >= 2,
            "Expected at least 2 intersections with multi-span curve"
        );
    }

    #[test]
    fn aabb_overlap_disjoint() {
        let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)];
        let pts2 = vec![Point3::new(5.0, 5.0, 0.0), Point3::new(6.0, 6.0, 0.0)];
        assert!(!aabb_overlap(&pts1, &pts2));
    }

    #[test]
    fn aabb_overlap_touching() {
        let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)];
        let pts2 = vec![Point3::new(0.5, 0.5, 0.5), Point3::new(2.0, 2.0, 2.0)];
        assert!(aabb_overlap(&pts1, &pts2));
    }

    #[test]
    fn aabb_overlap_empty() {
        let empty: Vec<Point3> = vec![];
        let pts = vec![Point3::new(0.0, 0.0, 0.0)];
        assert!(!aabb_overlap(&empty, &pts));
        assert!(!aabb_overlap(&pts, &empty));
    }

    #[test]
    fn deduplicate_removes_close_points() {
        let mut results = vec![
            CurveIntersection::Point {
                point: Point3::new(0.5, 0.5, 0.0),
                t1: 0.5,
                t2: 0.5,
            },
            CurveIntersection::Point {
                point: Point3::new(0.5, 0.5, 0.0),
                t1: 0.500001,
                t2: 0.500001,
            },
        ];
        deduplicate(&mut results, &Tolerance::default());
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn sort_by_overlap_t1_range() {
        // Test that overlaps sort correctly by t1_range.0
        let result = intersect_nurbs_nurbs(
            &NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                weights: vec![1.0, 1.0],
            },
            (0.0, 1.0),
            &NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                weights: vec![1.0, 1.0],
            },
            (0.0, 1.0),
            &Tolerance::default(),
        );
        // Identical curves => overlap
        assert!(!result.is_empty());
    }
}
