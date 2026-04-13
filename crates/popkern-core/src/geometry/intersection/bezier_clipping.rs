//! Bezier clipping algorithm for curve intersection.

#![allow(dead_code)]

use crate::geometry::Point3;
use crate::tolerance::Tolerance;

/// Clip the parameter range of curve_b against the fat line of curve_a.
/// Returns the clipped parameter range of curve_b, or None if no overlap.
pub(crate) fn clip_against_fat_line(
    a_points: &[Point3],
    b_points: &[Point3],
    b_range: (f64, f64),
) -> Option<(f64, f64)> {
    if a_points.len() < 2 || b_points.is_empty() {
        return None;
    }

    let a_first = a_points[0];
    let a_last = a_points[a_points.len() - 1];

    // Fat line direction: line through A[0] and A[last]
    let line_dir = a_last - a_first;
    let line_len_sq = line_dir.length_squared();

    // Degenerate A (all control points coincide)
    if line_len_sq < 1e-30 {
        return None;
    }

    // Compute a normal perpendicular to the fat line
    let normal = compute_fat_line_normal(&line_dir, a_points, b_points)?;

    // Signed distance of a point from the fat line
    let signed_distance = |p: Point3| -> f64 { (p - a_first).dot(&normal) };

    // Compute d_max: max absolute distance of A's control points from the fat line
    let mut d_max: f64 = 0.0;
    for &ap in a_points {
        let d = signed_distance(ap).abs();
        if d > d_max {
            d_max = d;
        }
    }

    // For linear A, d_max is 0. Use a small epsilon to allow crossing detection.
    if d_max < 1e-14 {
        d_max = 1e-14;
    }

    // Compute signed distances of B's control points from the fat line.
    // These form Bernstein coefficients of a distance polynomial D(t).
    let n = b_points.len();
    let distances: Vec<f64> = b_points.iter().map(|&p| signed_distance(p)).collect();

    // Build the Bernstein polygon: points (i/(n-1), distances[i])
    // Find the t-range where the polygon stays within [-d_max, +d_max].
    let mut t_min = 1.0_f64;
    let mut t_max = 0.0_f64;

    let n_f = if n > 1 { (n - 1) as f64 } else { 1.0 };

    // Check each vertex of the Bernstein polygon
    for (i, &d) in distances.iter().enumerate() {
        let t = if n == 1 { 0.5 } else { i as f64 / n_f };
        if d.abs() <= d_max {
            t_min = t_min.min(t);
            t_max = t_max.max(t);
        }
    }

    // Check each edge of the Bernstein polygon for crossings with +/- d_max
    for i in 0..n.saturating_sub(1) {
        let t0 = i as f64 / n_f;
        let t1 = (i + 1) as f64 / n_f;
        let d0 = distances[i];
        let d1 = distances[i + 1];

        // Check crossing with +d_max
        if (d0 - d_max) * (d1 - d_max) < 0.0 {
            let alpha = (d_max - d0) / (d1 - d0);
            let t_cross = t0 + alpha * (t1 - t0);
            t_min = t_min.min(t_cross);
            t_max = t_max.max(t_cross);
        }

        // Check crossing with -d_max
        if (d0 + d_max) * (d1 + d_max) < 0.0 {
            let alpha = (-d_max - d0) / (d1 - d0);
            let t_cross = t0 + alpha * (t1 - t0);
            t_min = t_min.min(t_cross);
            t_max = t_max.max(t_cross);
        }
    }

    if t_min > t_max {
        return None;
    }

    // Clamp to [0, 1]
    t_min = t_min.max(0.0);
    t_max = t_max.min(1.0);

    // Map from [0,1] back to b_range
    let b_len = b_range.1 - b_range.0;
    let result_min = b_range.0 + t_min * b_len;
    let result_max = b_range.0 + t_max * b_len;

    Some((result_min, result_max))
}

/// Compute the normal direction for the fat line.
fn compute_fat_line_normal(
    line_dir: &crate::geometry::Vec3,
    a_points: &[Point3],
    b_points: &[Point3],
) -> Option<crate::geometry::Vec3> {
    let a_first = a_points[0];

    let mut best_cross = crate::geometry::Vec3::zero();
    let mut best_cross_len_sq = 0.0_f64;

    for points in [a_points, b_points] {
        for &p in points {
            let v = p - a_first;
            let c = line_dir.cross(&v);
            let len_sq = c.length_squared();
            if len_sq > best_cross_len_sq {
                best_cross = c;
                best_cross_len_sq = len_sq;
            }
        }
    }

    if best_cross_len_sq < 1e-30 {
        let candidates = [
            crate::geometry::Vec3::new(1.0, 0.0, 0.0),
            crate::geometry::Vec3::new(0.0, 1.0, 0.0),
            crate::geometry::Vec3::new(0.0, 0.0, 1.0),
        ];
        for c in &candidates {
            let cross = line_dir.cross(c);
            if cross.length_squared() > 1e-20 {
                return cross.normalized();
            }
        }
        return None;
    }

    let normal = line_dir.cross(&best_cross);
    if normal.length_squared() < 1e-30 {
        return best_cross.normalized();
    }
    normal.normalized()
}

/// De Casteljau subdivision: extract control points for the sub-curve
/// on [t_lo, t_hi] from the original Bezier defined on [0, 1].
fn subdivide_bezier(points: &[Point3], t_lo: f64, t_hi: f64) -> Vec<Point3> {
    if points.len() <= 1 {
        return points.to_vec();
    }
    if (t_hi - t_lo) > 1.0 - 1e-14 && t_lo < 1e-14 {
        return points.to_vec();
    }

    // Split at t_hi to get [0, t_hi] portion (left half)
    let left = de_casteljau_left(points, t_hi);
    // Now reparameterise: [t_lo, t_hi] maps to [t_lo/t_hi, 1] on the left curve
    if t_hi < 1e-30 {
        return points.to_vec();
    }
    let s = t_lo / t_hi;
    if s < 1e-14 {
        return left;
    }
    // Take [s, 1] of the left curve = right portion at s
    de_casteljau_right(&left, s)
}

/// De Casteljau: return control points for the LEFT portion [0, t].
/// The left portion collects the first point at each reduction level.
fn de_casteljau_left(points: &[Point3], t: f64) -> Vec<Point3> {
    let n = points.len();
    let mut work: Vec<Point3> = points.to_vec();
    let mut left = Vec::with_capacity(n);
    left.push(work[0]);

    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp_point(work[i], work[i + 1], t);
        }
        left.push(work[0]);
    }
    left
}

/// De Casteljau: return control points for the RIGHT portion [t, 1].
/// The right portion collects the last point at each reduction level.
fn de_casteljau_right(points: &[Point3], t: f64) -> Vec<Point3> {
    let n = points.len();
    let mut work: Vec<Point3> = points.to_vec();
    let mut right = Vec::with_capacity(n);
    right.push(work[n - 1]);

    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp_point(work[i], work[i + 1], t);
        }
        right.push(work[n - level - 1]);
    }
    right.reverse();
    right
}

fn lerp_point(a: Point3, b: Point3, t: f64) -> Point3 {
    Point3::new(
        a.x + (b.x - a.x) * t,
        a.y + (b.y - a.y) * t,
        a.z + (b.z - a.z) * t,
    )
}

/// Find all intersection parameter pairs between two Bezier segments.
/// Results are approximate (t1, t2) pairs that should be refined by Newton.
#[allow(clippy::too_many_arguments)]
pub(crate) fn bezier_intersect(
    a_points: &[Point3],
    a_weights: &[f64],
    a_range: (f64, f64),
    b_points: &[Point3],
    b_weights: &[f64],
    b_range: (f64, f64),
    tolerance: &Tolerance,
    max_depth: usize,
    results: &mut Vec<(f64, f64)>,
) {
    // Convert rational control points to Cartesian by dividing by weights
    let a_cart = to_cartesian(a_points, a_weights);
    let b_cart = to_cartesian(b_points, b_weights);

    bezier_intersect_inner(
        &a_cart, a_range, &b_cart, b_range, tolerance, max_depth, results,
    );
}

fn to_cartesian(points: &[Point3], weights: &[f64]) -> Vec<Point3> {
    points
        .iter()
        .zip(weights.iter())
        .map(|(&p, &w)| {
            if (w - 1.0).abs() < f64::EPSILON {
                p
            } else {
                Point3::new(p.x / w, p.y / w, p.z / w)
            }
        })
        .collect()
}

/// Maximum number of intersection candidates to prevent runaway recursion.
const MAX_CANDIDATES: usize = 100;

/// Evaluate the midpoint of a Bezier curve using De Casteljau.
fn evaluate_bezier_midpoint(points: &[Point3]) -> Point3 {
    if points.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let mut work = points.to_vec();
    let n = work.len();
    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp_point(work[i], work[i + 1], 0.5);
        }
    }
    work[0]
}

#[allow(clippy::too_many_arguments)]
fn bezier_intersect_inner(
    a_points: &[Point3],
    a_range: (f64, f64),
    b_points: &[Point3],
    b_range: (f64, f64),
    tolerance: &Tolerance,
    max_depth: usize,
    results: &mut Vec<(f64, f64)>,
) {
    bezier_intersect_inner_with_stalls(
        a_points, a_range, b_points, b_range, tolerance, max_depth, results, 0,
    );
}

#[allow(clippy::too_many_arguments)]
fn bezier_intersect_inner_with_stalls(
    a_points: &[Point3],
    a_range: (f64, f64),
    b_points: &[Point3],
    b_range: (f64, f64),
    tolerance: &Tolerance,
    max_depth: usize,
    results: &mut Vec<(f64, f64)>,
    consecutive_stalls: usize,
) {
    // Safety: cap results to prevent runaway recursion
    if results.len() >= MAX_CANDIDATES {
        return;
    }

    let a_span = a_range.1 - a_range.0;
    let b_span = b_range.1 - b_range.0;

    // Base case: both ranges small enough
    if a_span < tolerance.parametric_base() && b_span < tolerance.parametric_base() {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        results.push((t1, t2));
        return;
    }

    // Asymmetric convergence: one parameter has converged but the other has not.
    // This happens in near-tangent cases where the clipping isolates A but B remains
    // wide. Report the intersection at the A midpoint; Newton polishing downstream
    // will refine the B parameter.
    if a_span < tolerance.parametric_base() {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        results.push((t1, t2));
        return;
    }
    if b_span < tolerance.parametric_base() {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        results.push((t1, t2));
        return;
    }

    // Max depth fail-safe
    if max_depth == 0 {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        results.push((t1, t2));
        return;
    }

    // Tangent-specific termination: when intervals are small but not yet
    // at parametric tolerance, check for closest approach and report if
    // the curves are within linear tolerance (handles near-tangent cases).
    if a_span < tolerance.parametric_base() * 10.0 && b_span < tolerance.parametric_base() * 10.0 {
        let a_sub_mid = subdivide_bezier(a_points, a_range.0, a_range.1);
        let b_sub_mid = subdivide_bezier(b_points, b_range.0, b_range.1);
        let pa = evaluate_bezier_midpoint(&a_sub_mid);
        let pb = evaluate_bezier_midpoint(&b_sub_mid);
        let lin_sq = tolerance.linear() * tolerance.linear();
        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (pa.x, pa.y, pa.z),
            (pb.x, pb.y, pb.z),
            lin_sq,
        ) {
            let t1 = (a_range.0 + a_range.1) * 0.5;
            let t2 = (b_range.0 + b_range.1) * 0.5;
            results.push((t1, t2));
            return;
        }
    }

    // Subdivide both curves to get tight control polygons for current ranges.
    let a_sub = subdivide_bezier(a_points, a_range.0, a_range.1);
    let b_sub = subdivide_bezier(b_points, b_range.0, b_range.1);

    // Clip B against A's fat line
    let clipped = clip_against_fat_line(&a_sub, &b_sub, b_range);
    let clipped_b_range = match clipped {
        None => return, // no intersection
        Some(r) => r,
    };

    let clipped_b_span = clipped_b_range.1 - clipped_b_range.0;
    let b_reduction = if b_span > 1e-30 {
        clipped_b_span / b_span
    } else {
        1.0
    };

    // Also clip A against B's fat line to get the best reduction from either direction.
    // This is critical for near-tangent cases where one direction gives poor reduction
    // but the other direction (e.g., clipping a wide curve against a narrow line) gives
    // excellent reduction.
    let b_sub2 = subdivide_bezier(b_points, clipped_b_range.0, clipped_b_range.1);
    let clipped_a = clip_against_fat_line(&b_sub2, &a_sub, a_range);
    let clipped_a_range = match clipped_a {
        None => return,
        Some(r) => r,
    };

    let clipped_a_span = clipped_a_range.1 - clipped_a_range.0;
    let a_reduction = if a_span > 1e-30 {
        clipped_a_span / a_span
    } else {
        1.0
    };

    // Use the best reduction from either clipping direction.
    let best_reduction = b_reduction.min(a_reduction);

    // Track consecutive stalls at this recursion level.
    let new_stalls = if best_reduction > 0.8 {
        consecutive_stalls + 1
    } else {
        0
    };

    // After 3 consecutive stalls (near-tangent convergence failure), bisect
    // BOTH curves and recurse on all 4 sub-pairs for aggressive resolution.
    if new_stalls >= 3 {
        let a_mid = (clipped_a_range.0 + clipped_a_range.1) * 0.5;
        let b_mid = (clipped_b_range.0 + clipped_b_range.1) * 0.5;
        // All 4 sub-pairs
        for &a_sub_range in &[(clipped_a_range.0, a_mid), (a_mid, clipped_a_range.1)] {
            for &b_sub_range in &[(clipped_b_range.0, b_mid), (b_mid, clipped_b_range.1)] {
                bezier_intersect_inner_with_stalls(
                    a_points,
                    a_sub_range,
                    b_points,
                    b_sub_range,
                    tolerance,
                    max_depth - 1,
                    results,
                    0, // reset stall counter after bisection
                );
            }
        }
        return;
    }

    if best_reduction < 0.8 {
        // Good convergence from at least one direction: recurse with both clipped ranges.
        bezier_intersect_inner_with_stalls(
            a_points,
            clipped_a_range,
            b_points,
            clipped_b_range,
            tolerance,
            max_depth - 1,
            results,
            0, // good convergence resets stall counter
        );
    } else {
        // Slow convergence (near-tangent case): bisect the curve with the longer range
        if clipped_a_span >= clipped_b_span {
            let a_mid = (clipped_a_range.0 + clipped_a_range.1) * 0.5;
            bezier_intersect_inner_with_stalls(
                a_points,
                (clipped_a_range.0, a_mid),
                b_points,
                clipped_b_range,
                tolerance,
                max_depth - 1,
                results,
                new_stalls,
            );
            bezier_intersect_inner_with_stalls(
                a_points,
                (a_mid, clipped_a_range.1),
                b_points,
                clipped_b_range,
                tolerance,
                max_depth - 1,
                results,
                new_stalls,
            );
        } else {
            let b_mid = (clipped_b_range.0 + clipped_b_range.1) * 0.5;
            bezier_intersect_inner_with_stalls(
                a_points,
                clipped_a_range,
                b_points,
                (clipped_b_range.0, b_mid),
                tolerance,
                max_depth - 1,
                results,
                new_stalls,
            );
            bezier_intersect_inner_with_stalls(
                a_points,
                clipped_a_range,
                b_points,
                (b_mid, clipped_b_range.1),
                tolerance,
                max_depth - 1,
                results,
                new_stalls,
            );
        }
    }
}

#[cfg(test)]
#[path = "bezier_clipping_tests.rs"]
mod tests;
