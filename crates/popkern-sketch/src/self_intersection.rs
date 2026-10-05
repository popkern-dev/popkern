// Items are pub(crate) for use by the pipeline module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! Self-intersection detection for traced sketch loops.
//!
//! Checks each pair of non-adjacent edges within a loop for geometric
//! crossing. Lines use parametric intersection; arcs and curves are
//! approximated as polylines with 20 sample points per edge.

use crate::error::SketchError;
use crate::loop_tracer::{DirectedEdge, TracedLoop};
use popkern_constraint::{ElementHandle, Sketch};
use popkern_core::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Number of sample points for polyline approximation of curves.
const CURVE_SAMPLES: usize = 20;

/// Check all traced loops for self-intersecting edges.
///
/// For each loop, checks every pair of non-adjacent edges for crossing.
/// Non-adjacent means the two edges share no node in the loop.
///
/// Returns `Ok(())` if no crossings are found, or
/// `Err(SketchError::SelfIntersectingProfile(x, y))` at the first
/// crossing detected.
pub(crate) fn check_self_intersection(
    loops: &[TracedLoop],
    sketch: &Sketch,
    tolerance: &Tolerance,
) -> Result<(), SketchError> {
    for traced in loops {
        check_loop_self_intersection(traced, sketch, tolerance)?;
    }
    Ok(())
}

/// Check a single loop for self-intersecting edge pairs.
fn check_loop_self_intersection(
    traced: &TracedLoop,
    sketch: &Sketch,
    tolerance: &Tolerance,
) -> Result<(), SketchError> {
    let edges = &traced.edges;
    let n = edges.len();
    if n < 3 {
        return Ok(());
    }

    // Polyline-approximate each edge for intersection testing.
    let polylines: Vec<Vec<(f64, f64)>> = edges.iter().map(|e| sample_edge(e, sketch)).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            // Skip adjacent edges (they share a vertex naturally).
            if are_adjacent(i, j, n) {
                continue;
            }

            if let Some((x, y)) = check_polyline_crossing(&polylines[i], &polylines[j], tolerance) {
                return Err(SketchError::SelfIntersectingProfile(x, y));
            }
        }
    }

    Ok(())
}

/// Return `true` if edges at indices `i` and `j` are adjacent in a
/// loop of length `n` (i.e., they share a node).
fn are_adjacent(i: usize, j: usize, n: usize) -> bool {
    let diff = i.abs_diff(j);
    diff == 1 || diff == n - 1
}

/// Sample an edge as a polyline of 2D points.
///
/// For lines, returns just the two endpoints.
/// For arcs and other curves, returns `CURVE_SAMPLES + 1` points
/// along the curve.
fn sample_edge(edge: &DirectedEdge, sketch: &Sketch) -> Vec<(f64, f64)> {
    match &edge.element {
        ElementHandle::Line(_) => {
            // Lines need only their two endpoints.
            vec![edge.start, edge.end]
        }
        ElementHandle::Arc(h) => {
            if let Some((cx, cy, r, start_angle, end_angle)) = sketch.arc_params(*h) {
                let (a0, a1) = if edge.reversed {
                    (end_angle, start_angle)
                } else {
                    (start_angle, end_angle)
                };
                sample_arc(cx, cy, r, a0, a1)
            } else {
                vec![edge.start, edge.end]
            }
        }
        ElementHandle::Spline(h) => {
            if let Some(knots) = sketch.spline_knots(*h) {
                let degree = sketch.spline_degree(*h).unwrap_or(3);
                let n_cp = sketch
                    .spline_control_points(*h)
                    .map(|c| c.len())
                    .unwrap_or(0);
                if n_cp == 0 || knots.is_empty() {
                    return vec![edge.start, edge.end];
                }
                let t_start = knots[degree.min(knots.len() - 1)];
                let t_end = knots[n_cp.min(knots.len() - 1)];
                let (t0, t1) = if edge.reversed {
                    (t_end, t_start)
                } else {
                    (t_start, t_end)
                };
                sample_spline(sketch, *h, t0, t1)
            } else {
                vec![edge.start, edge.end]
            }
        }
        // Ellipse, Parabola, Hyperbola: fall back to straight-line approximation.
        _ => vec![edge.start, edge.end],
    }
}

/// Sample an arc as a polyline.
fn sample_arc(cx: f64, cy: f64, r: f64, a0: f64, a1: f64) -> Vec<(f64, f64)> {
    let mut pts = Vec::with_capacity(CURVE_SAMPLES + 1);
    for i in 0..=CURVE_SAMPLES {
        let t = i as f64 / CURVE_SAMPLES as f64;
        let angle = a0 + t * (a1 - a0);
        pts.push((cx + r * angle.pcos(), cy + r * angle.psin()));
    }
    pts
}

/// Sample a spline as a polyline.
fn sample_spline(
    sketch: &Sketch,
    handle: popkern_constraint::SplineHandle,
    t0: f64,
    t1: f64,
) -> Vec<(f64, f64)> {
    let mut pts = Vec::with_capacity(CURVE_SAMPLES + 1);
    for i in 0..=CURVE_SAMPLES {
        let t = t0 + (i as f64 / CURVE_SAMPLES as f64) * (t1 - t0);
        if let Some(pt) = sketch.spline_point_at(handle, t) {
            pts.push(pt);
        }
    }
    if pts.len() < 2 {
        // Fallback: shouldn't happen, but return at least 2 points.
        return vec![(0.0, 0.0), (0.0, 0.0)];
    }
    pts
}

/// Check two polylines for crossing. Returns the approximate crossing
/// point `(x, y)` if found.
fn check_polyline_crossing(
    poly_a: &[(f64, f64)],
    poly_b: &[(f64, f64)],
    tolerance: &Tolerance,
) -> Option<(f64, f64)> {
    for i in 0..poly_a.len().saturating_sub(1) {
        let a1 = poly_a[i];
        let a2 = poly_a[i + 1];
        for j in 0..poly_b.len().saturating_sub(1) {
            let b1 = poly_b[j];
            let b2 = poly_b[j + 1];
            if let Some(pt) = segment_intersection(a1, a2, b1, b2, tolerance) {
                return Some(pt);
            }
        }
    }
    None
}

/// Parametric segment-segment intersection.
///
/// Returns the intersection point if the two segments cross at an
/// interior point (exclusive of endpoints within tolerance).
fn segment_intersection(
    p1: (f64, f64),
    p2: (f64, f64),
    p3: (f64, f64),
    p4: (f64, f64),
    tolerance: &Tolerance,
) -> Option<(f64, f64)> {
    let d_a = (p2.0 - p1.0, p2.1 - p1.1);
    let d_b = (p4.0 - p3.0, p4.1 - p3.1);

    // 2D cross product of the two direction vectors.
    let cross = d_a.0 * d_b.1 - d_a.1 * d_b.0;
    if cross.abs() < 1e-15 {
        // Parallel or collinear -- skip.
        return None;
    }

    let d13 = (p3.0 - p1.0, p3.1 - p1.1);
    let t = (d13.0 * d_b.1 - d13.1 * d_b.0) / cross;
    let u = (d13.0 * d_a.1 - d13.1 * d_a.0) / cross;

    // Exclusive of endpoints: use tolerance-based margin.
    let margin = tolerance.linear() * 0.01; // small margin to avoid endpoint noise
    let lo = margin;
    let hi = 1.0 - margin;

    if t > lo && t < hi && u > lo && u < hi {
        let x = p1.0 + t * d_a.0;
        let y = p1.1 + t * d_a.1;
        Some((x, y))
    } else {
        None
    }
}

#[cfg(test)]
#[path = "self_intersection_tests.rs"]
mod tests;
