//! Bezier clipping algorithm for curve intersection.
#![allow(dead_code)]

use crate::geometry::Point3;
use crate::tolerance::Tolerance;

#[cfg(test)]
pub(crate) use super::rational_bezier::compute_fat_line_normal;
pub(crate) use super::rational_bezier::{
    HPoint, clip_against_fat_line, evaluate_bezier_midpoint_h, project, subdivide_bezier_h,
    to_homogeneous,
};

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

pub(crate) fn lerp_point(a: Point3, b: Point3, t: f64) -> Point3 {
    Point3::new(
        a.x + (b.x - a.x) * t,
        a.y + (b.y - a.y) * t,
        a.z + (b.z - a.z) * t,
    )
}

/// Find all intersection parameter pairs between two Bezier segments.
/// Results are approximate (t1, t2) pairs in the caller's a_range /
/// b_range coordinates, to be refined by Newton.
///
/// The recursion runs in each segment's LOCAL [0, 1] domain — the
/// subdivision helpers interpret ranges as fractions of the control
/// polygon, so feeding them absolute knot values (a quarter arc's
/// [0.25, 0.5], a trimmed line's [17, 33]) extrapolates garbage
/// polygons. Results are mapped back to the absolute ranges here.
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
    let a_h = to_homogeneous(a_points, a_weights);
    let b_h = to_homogeneous(b_points, b_weights);

    // Local-domain candidates of THIS pair only: the dedupe and the
    // candidate cap compare like with like (earlier pairs' results are
    // already mapped to absolute ranges).
    let mut local: Vec<(f64, f64)> = Vec::new();
    let mut budget = MAX_NODES;
    bezier_intersect_inner(
        &a_h,
        (0.0, 1.0),
        &b_h,
        (0.0, 1.0),
        tolerance,
        max_depth,
        &mut local,
        &mut budget,
    );
    results.extend(local.into_iter().map(|r| {
        (
            a_range.0 + r.0 * (a_range.1 - a_range.0),
            b_range.0 + r.1 * (b_range.1 - b_range.0),
        )
    }));
}

/// Maximum number of intersection candidates to prevent runaway recursion.
const MAX_CANDIDATES: usize = 100;

/// Global subdivision-node budget per `bezier_intersect` call. The
/// stall path 4-way bisects under only a DEPTH cap, so a persistent
/// tangency — curves hugging inside the fat line for a stretch while
/// staying outside linear tolerance — visits O(4^depth) nodes and
/// emits zero candidates for MAX_CANDIDATES to cap (PK-295, measured:
/// a spiric pcurve tangential to a chart boundary line spun one
/// boolean for 25 minutes). Honest configurations use well under a
/// thousand nodes; on exhaustion the node reports its midpoint like
/// the depth fail-safe and lets Newton polishing downstream decide.
const MAX_NODES: usize = 20_000;

/// Record a candidate unless one already sits within `100 ×
/// parametric_base` of it in both parameters. The recursion revisits a
/// root many times at its noise floor (stall bisection splits four ways
/// around it); counting those revisits against `MAX_CANDIDATES` starved
/// every LATER root of the segment pair (PK-739, measured: a square
/// side crossing a circle twice — ~100 copies of the first crossing,
/// the second never searched). Newton polishes the survivor.
fn push_candidate(results: &mut Vec<(f64, f64)>, c: (f64, f64), tolerance: &Tolerance) {
    let eps = tolerance.parametric_base() * 100.0;
    if results
        .iter()
        .any(|r| (r.0 - c.0).abs() < eps && (r.1 - c.1).abs() < eps)
    {
        return;
    }
    results.push(c);
}

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
    a_points: &[HPoint],
    a_range: (f64, f64),
    b_points: &[HPoint],
    b_range: (f64, f64),
    tolerance: &Tolerance,
    max_depth: usize,
    results: &mut Vec<(f64, f64)>,
    budget: &mut usize,
) {
    bezier_intersect_inner_with_stalls(
        a_points, a_range, b_points, b_range, tolerance, max_depth, results, 0, budget,
    );
}

#[allow(clippy::too_many_arguments)]
fn bezier_intersect_inner_with_stalls(
    a_points: &[HPoint],
    a_range: (f64, f64),
    b_points: &[HPoint],
    b_range: (f64, f64),
    tolerance: &Tolerance,
    max_depth: usize,
    results: &mut Vec<(f64, f64)>,
    consecutive_stalls: usize,
    budget: &mut usize,
) {
    // Safety: cap results to prevent runaway recursion
    if results.len() >= MAX_CANDIDATES {
        return;
    }
    // Node-budget fail-safe (see MAX_NODES): report the midpoint —
    // exactly the depth-exhaustion contract — instead of recursing.
    if *budget == 0 {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        push_candidate(results, (t1, t2), tolerance);
        return;
    }
    *budget -= 1;

    let a_span = a_range.1 - a_range.0;
    let b_span = b_range.1 - b_range.0;

    // Base case: both ranges small enough
    if a_span < tolerance.parametric_base() && b_span < tolerance.parametric_base() {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        push_candidate(results, (t1, t2), tolerance);
        return;
    }

    // Asymmetric convergence: one parameter has converged but the other has not.
    // This happens in near-tangent cases where the clipping isolates A but B remains
    // wide. Report the intersection at the A midpoint; Newton polishing downstream
    // will refine the B parameter.
    if a_span < tolerance.parametric_base() {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        push_candidate(results, (t1, t2), tolerance);
        return;
    }
    if b_span < tolerance.parametric_base() {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        push_candidate(results, (t1, t2), tolerance);
        return;
    }

    // Max depth fail-safe
    if max_depth == 0 {
        let t1 = (a_range.0 + a_range.1) * 0.5;
        let t2 = (b_range.0 + b_range.1) * 0.5;
        push_candidate(results, (t1, t2), tolerance);
        return;
    }

    // Tangent-specific termination: when intervals are small but not yet
    // at parametric tolerance, check for closest approach and report if
    // the curves are within linear tolerance (handles near-tangent cases).
    if a_span < tolerance.parametric_base() * 10.0 && b_span < tolerance.parametric_base() * 10.0 {
        let a_sub_mid = subdivide_bezier_h(a_points, a_range.0, a_range.1);
        let b_sub_mid = subdivide_bezier_h(b_points, b_range.0, b_range.1);
        let pa = evaluate_bezier_midpoint_h(&a_sub_mid);
        let pb = evaluate_bezier_midpoint_h(&b_sub_mid);
        let lin_sq = tolerance.linear() * tolerance.linear();
        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (pa.x, pa.y, pa.z),
            (pb.x, pb.y, pb.z),
            lin_sq,
        ) {
            let t1 = (a_range.0 + a_range.1) * 0.5;
            let t2 = (b_range.0 + b_range.1) * 0.5;
            push_candidate(results, (t1, t2), tolerance);
            return;
        }
    }

    // Subdivide both curves to get tight control polygons for current ranges.
    let a_sub_h = subdivide_bezier_h(a_points, a_range.0, a_range.1);
    let a_sub = project(&a_sub_h);
    let b_sub_h = subdivide_bezier_h(b_points, b_range.0, b_range.1);

    // Clip B against A's fat line
    let clipped = clip_against_fat_line(&a_sub, &b_sub_h, b_range);
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
    let b_sub2 = project(&subdivide_bezier_h(
        b_points,
        clipped_b_range.0,
        clipped_b_range.1,
    ));
    let clipped_a = clip_against_fat_line(&b_sub2, &a_sub_h, a_range);
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
                    budget,
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
            budget,
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
                budget,
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
                budget,
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
                budget,
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
                budget,
            );
        }
    }
}

#[cfg(test)]
#[path = "bezier_clipping_tests.rs"]
mod tests;
