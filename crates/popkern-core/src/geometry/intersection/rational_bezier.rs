//! Rational (homogeneous) Bezier machinery for the clipping
//! intersection solver: weighted fat-line clipping and exact
//! rational de Casteljau subdivision.
//!
//! Split from `bezier_clipping.rs` to keep it under the structural
//! 500-line cap.
#![allow(dead_code)]

use crate::geometry::Point3;

use super::bezier_clipping::lerp_point;

/// Homogeneous control point (w·x, w·y, w·z, w). Rational Bezier
/// subdivision is exact only in this space: de Casteljau on the 4D
/// coordinates, then projection, yields the true control polygon of
/// the sub-curve. Subdividing Cartesian points directly is wrong for
/// any weight ≠ 1 (the projected polygon no longer bounds the curve,
/// so fat-line clipping discards true crossings on circle arcs).
#[derive(Clone, Copy)]
pub(crate) struct HPoint {
    pub(crate) wp: Point3,
    pub(crate) w: f64,
}

pub(crate) fn to_homogeneous(points: &[Point3], weights: &[f64]) -> Vec<HPoint> {
    points
        .iter()
        .zip(weights.iter())
        .map(|(&p, &w)| HPoint {
            wp: Point3::new(p.x * w, p.y * w, p.z * w),
            w,
        })
        .collect()
}

/// Project a homogeneous polygon to Cartesian. For positive weights
/// the rational curve lies in the convex hull of these points, which
/// is what the fat-line bound requires.
pub(crate) fn project(points: &[HPoint]) -> Vec<Point3> {
    points
        .iter()
        .map(|h| {
            if h.w.abs() < 1e-300 {
                h.wp
            } else {
                Point3::new(h.wp.x / h.w, h.wp.y / h.w, h.wp.z / h.w)
            }
        })
        .collect()
}

pub(crate) fn lerp_h(a: HPoint, b: HPoint, t: f64) -> HPoint {
    HPoint {
        wp: lerp_point(a.wp, b.wp, t),
        w: a.w + (b.w - a.w) * t,
    }
}

/// De Casteljau subdivision in homogeneous space: control points of
/// the sub-curve on [t_lo, t_hi], where t_lo/t_hi are fractions of the
/// polygon's own [0, 1] domain.
pub(crate) fn subdivide_bezier_h(points: &[HPoint], t_lo: f64, t_hi: f64) -> Vec<HPoint> {
    if points.len() <= 1 {
        return points.to_vec();
    }
    if (t_hi - t_lo) > 1.0 - 1e-14 && t_lo < 1e-14 {
        return points.to_vec();
    }
    if t_hi < 1e-30 {
        return points.to_vec();
    }
    // Left portion [0, t_hi]: first point of each reduction level.
    let n = points.len();
    let mut work: Vec<HPoint> = points.to_vec();
    let mut left = Vec::with_capacity(n);
    left.push(work[0]);
    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp_h(work[i], work[i + 1], t_hi);
        }
        left.push(work[0]);
    }
    let s = t_lo / t_hi;
    if s < 1e-14 {
        return left;
    }
    // Right portion [s, 1] of the left curve: last point of each level.
    let mut work = left;
    let mut right = Vec::with_capacity(n);
    right.push(work[n - 1]);
    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp_h(work[i], work[i + 1], s);
        }
        right.push(work[n - level - 1]);
    }
    right.reverse();
    right
}

/// Evaluate the midpoint of a homogeneous Bezier polygon.
pub(crate) fn evaluate_bezier_midpoint_h(points: &[HPoint]) -> Point3 {
    if points.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let mut work = points.to_vec();
    let n = work.len();
    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp_h(work[i], work[i + 1], 0.5);
        }
    }
    let h = work[0];
    if h.w.abs() < 1e-300 {
        h.wp
    } else {
        Point3::new(h.wp.x / h.w, h.wp.y / h.w, h.wp.z / h.w)
    }
}

/// Clip the parameter range of curve_b against the fat line of curve_a.
/// Returns the clipped parameter range of curve_b, or None if no overlap.
///
/// B is given HOMOGENEOUSLY: its signed-distance function is the
/// rational d(t) = Σ Bᵢ(t)·wᵢ·dᵢ / Σ Bᵢ(t)·wᵢ, so band membership
/// |d(t)| ≤ d_max is equivalent (for positive weights) to
/// f(t) = Σ Bᵢ·wᵢ·(dᵢ − d_max) ≤ 0  ∧  g(t) = Σ Bᵢ·wᵢ·(dᵢ + d_max) ≥ 0,
/// and the polygon crossing tests must use the WEIGHTED coefficients
/// wᵢ·dᵢ. Using raw dᵢ (the polynomial formula) clips off true roots
/// on rational arcs.
/// Is a collapsed (sub-epsilon) A segment within convergence reach of
/// B's control polygon? The band is absolute and generous relative
/// to convergence scales (the measured converged collapse sat 5e-7
/// from B) while far below the flush-graze separations whose
/// exploration strands tuned welds.
fn collapsed_segment_near(a: Point3, b_cart: &[Point3]) -> bool {
    const NEAR: f64 = 1.0e-5;
    if b_cart.len() == 1 {
        return b_cart[0].distance_to(&a) < NEAR;
    }
    b_cart.windows(2).any(|seg| {
        let ab = seg[1] - seg[0];
        let len_sq = ab.length_squared();
        let t = if len_sq > 1e-30 {
            ((a - seg[0]).dot(&ab) / len_sq).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (seg[0] + ab * t).distance_to(&a) < NEAR
    })
}

pub(crate) fn clip_against_fat_line(
    a_points: &[Point3],
    b_points: &[HPoint],
    b_range: (f64, f64),
) -> Option<(f64, f64)> {
    if a_points.len() < 2 || b_points.is_empty() {
        return None;
    }

    let b_cart = project(b_points);

    let a_first = a_points[0];
    let a_last = a_points[a_points.len() - 1];

    // Fat line direction: line through A[0] and A[last]
    let line_dir = a_last - a_first;
    let line_len_sq = line_dir.length_squared();

    // Degenerate A (all control points coincide): no fat line can be
    // established. That is NOT automatically a proof of disjointness
    // — near a converged transversal crossing the A segment
    // legitimately collapses below any absolute epsilon while B
    // still spans the root (measured: two rational circle pcurves
    // crossing at (−0.187, 0.066); the discarded branch was the
    // crossing itself and the cap face's clip saw one crossing
    // instead of two). But a blanket "no reduction" answer explores
    // every far-away collapse too, and the extra depth-exhaustion
    // candidates it mints strand tuned flush welds (measured: the
    // fanned gusset union). Decide by MEASURED proximity: a collapse
    // NEAR B's control polygon keeps converging; a far one is a
    // genuine exclusion.
    if line_len_sq < 1e-30 {
        if collapsed_segment_near(a_points[0], &b_cart) {
            return Some(b_range);
        }
        return None;
    }

    // Compute a normal perpendicular to the fat line
    // A missing normal is the same converged-collapse class as the
    // degenerate chord above — its cross products fell below the
    // helper's absolute thresholds — judged by the same proximity
    // rule.
    let Some(normal) = compute_fat_line_normal(&line_dir, a_points, &b_cart) else {
        if collapsed_segment_near(a_points[0], &b_cart) {
            return Some(b_range);
        }
        return None;
    };

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

    // Signed distances of B's projected control points, and the
    // rational band coefficients fᵢ = wᵢ(dᵢ − d_max), gᵢ = wᵢ(dᵢ + d_max).
    let n = b_points.len();
    let distances: Vec<f64> = b_cart.iter().map(|&p| signed_distance(p)).collect();
    let band: Vec<(f64, f64)> = distances
        .iter()
        .zip(b_points.iter())
        .map(|(&d, h)| {
            let w = if h.w.abs() < 1e-300 { 1.0 } else { h.w.abs() };
            (w * (d - d_max), w * (d + d_max))
        })
        .collect();

    // Build the Bernstein polygon over (i/(n-1), coefficient).
    // Find the t-range where f stays ≤ 0 and g stays ≥ 0.
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

    // Zero crossings of the coefficient polygons' CONVEX HULLS
    // (f = 0 is the +d_max band edge, g = 0 the −d_max band edge). The
    // distance function lies in the hull of its Bernstein points, not
    // on their polygon: consecutive edges alone miss the hull's chords
    // (first→last, …) and can return an interval that EXCLUDES the root
    // (PK-739, measured: a fin side crossing a rational quarter circle
    // was clipped away and the top cap of a 16-fin ring lost a crossing).
    // A hull ∩ line is the span of all pairwise chord crossings.
    for i in 0..n {
        for j in i + 1..n {
            let t0 = i as f64 / n_f;
            let t1 = j as f64 / n_f;
            let (f0, g0) = band[i];
            let (f1, g1) = band[j];

            if f0 * f1 < 0.0 {
                let alpha = f0 / (f0 - f1);
                let t_cross = t0 + alpha * (t1 - t0);
                t_min = t_min.min(t_cross);
                t_max = t_max.max(t_cross);
            }

            if g0 * g1 < 0.0 {
                let alpha = g0 / (g0 - g1);
                let t_cross = t0 + alpha * (t1 - t0);
                t_min = t_min.min(t_cross);
                t_max = t_max.max(t_cross);
            }
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
pub(crate) fn compute_fat_line_normal(
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
