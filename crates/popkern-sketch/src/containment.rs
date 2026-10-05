// Items are pub(crate) for use by the nesting module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! Point-in-loop containment testing using ray casting.
//!
//! Casts a ray from a test point in the +X direction and counts
//! intersections with each edge of a closed loop. Lines use exact
//! ray-segment intersection; arcs and curves are approximated as
//! polylines with 20 sample points per edge.

use crate::loop_tracer::{DirectedEdge, TracedLoop};
use popkern_constraint::{ElementHandle, Sketch};
use popkern_math::PortableMath;

/// Number of sample points for polyline approximation of curves.
const CURVE_SAMPLES: usize = 20;

/// Test if a point is inside a closed loop using ray casting.
///
/// Casts a ray from `point` in the +X direction and counts edge
/// intersections. An odd count means the point is inside.
pub(crate) fn point_in_loop(point: (f64, f64), loop_: &TracedLoop, sketch: &Sketch) -> bool {
    let mut crossings = 0;
    for edge in &loop_.edges {
        crossings += count_ray_crossings(point, edge, sketch);
    }
    crossings % 2 == 1
}

/// Count the number of times a +X ray from `point` crosses a single edge.
fn count_ray_crossings(point: (f64, f64), edge: &DirectedEdge, sketch: &Sketch) -> usize {
    let segments = edge_segments(edge, sketch);
    let mut count = 0;
    for seg in segments.windows(2) {
        if ray_crosses_segment(point, seg[0], seg[1]) {
            count += 1;
        }
    }
    count
}

/// Break an edge into polyline segments (list of points).
fn edge_segments(edge: &DirectedEdge, sketch: &Sketch) -> Vec<(f64, f64)> {
    match &edge.element {
        ElementHandle::Line(_) => {
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
        // Ellipse, Parabola, Hyperbola: straight-line fallback.
        _ => vec![edge.start, edge.end],
    }
}

/// Sample an arc as a polyline of `CURVE_SAMPLES + 1` points.
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
        return vec![(0.0, 0.0), (0.0, 0.0)];
    }
    pts
}

/// Test if a +X ray from `point` crosses the segment from `a` to `b`.
///
/// Uses the standard ray-casting rule: the ray crosses if exactly one
/// endpoint is above the ray's y-coordinate and the x-intersection is
/// to the right of the test point.
fn ray_crosses_segment(point: (f64, f64), a: (f64, f64), b: (f64, f64)) -> bool {
    let (px, py) = point;

    // One endpoint must be above and the other at or below (exclusive test
    // to avoid double-counting shared vertices).
    let a_above = a.1 > py;
    let b_above = b.1 > py;
    if a_above == b_above {
        return false;
    }

    // Compute the x-coordinate where the segment crosses y = py.
    let dy = b.1 - a.1;
    if dy.abs() < 1e-15 {
        return false;
    }
    let t = (py - a.1) / dy;
    let x_cross = a.0 + t * (b.0 - a.0);

    x_cross > px
}

#[cfg(test)]
#[path = "containment_tests.rs"]
mod tests;
