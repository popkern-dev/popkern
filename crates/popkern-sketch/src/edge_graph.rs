// Items are pub(crate) for use by future modules (loop tracer, etc.).
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! Edge graph construction from solved sketch geometry.
//!
//! Converts non-construction sketch elements into a node/edge graph
//! suitable for loop tracing. Degenerate edges (zero-length) are
//! filtered and reported as warnings.

use crate::error::SketchWarning;
use popkern_constraint::ElementHandle;
use popkern_constraint::Sketch;
use popkern_core::tolerance::Tolerance;
use popkern_math::PortableMath;

/// A graph of 2D nodes and edges extracted from a solved sketch.
///
/// Nodes are merged endpoint positions (within tolerance), and edges
/// are the non-construction sketch elements connecting them.
#[derive(Debug)]
pub(crate) struct EdgeGraph {
    /// Merged node positions in 2D.
    pub nodes: Vec<(f64, f64)>,
    /// Edges connecting nodes.
    pub edges: Vec<GraphEdge>,
}

/// A single edge in the graph, referencing back to the sketch element
/// that produced it.
#[derive(Debug)]
pub(crate) struct GraphEdge {
    /// Handle to the source sketch element.
    pub element: ElementHandle,
    /// Index into [`EdgeGraph::nodes`] for the start endpoint.
    pub start_node: usize,
    /// Index into [`EdgeGraph::nodes`] for the end endpoint.
    pub end_node: usize,
    /// Exact 2D position of the start endpoint (before merging).
    pub start_pos: (f64, f64),
    /// Exact 2D position of the end endpoint (before merging).
    pub end_pos: (f64, f64),
}

/// Build an edge graph from a solved sketch.
///
/// Collects all non-construction elements, computes their 2D endpoints,
/// filters degenerate (zero-length) edges, merges coincident endpoints
/// into shared nodes, and returns the graph along with any warnings.
pub(crate) fn build_edge_graph(
    sketch: &Sketch,
    tolerance: &Tolerance,
) -> (EdgeGraph, Vec<SketchWarning>) {
    let mut warnings = Vec::new();
    let mut nodes: Vec<(f64, f64)> = Vec::new();
    let mut edges: Vec<GraphEdge> = Vec::new();

    let elements = sketch.non_construction_curves();

    for element in elements {
        let endpoints = element_endpoints(sketch, &element);
        let Some((start, end)) = endpoints else {
            continue;
        };

        // Filter degenerate edges, but allow full-circle arcs, full ellipses,
        // and periodic (closed) splines — these have coincident start/end.
        let dx = end.0 - start.0;
        let dy = end.1 - start.1;
        let dist = (dx * dx + dy * dy).sqrt();
        if dist < tolerance.linear() {
            // Check if this is a full-circle arc (angular span > PI).
            let is_full_circle = match &element {
                ElementHandle::Arc(ah) => {
                    if let Some((_, _, _, start_angle, end_angle)) = sketch.arc_params(*ah) {
                        (end_angle - start_angle).abs() > std::f64::consts::PI
                    } else {
                        false
                    }
                }
                _ => false,
            };

            // Check if this is a full ellipse (always a closed curve).
            let is_full_ellipse = matches!(&element, ElementHandle::Ellipse(_));

            // Check if this is a periodic (closed) spline.
            let is_periodic_spline = match &element {
                ElementHandle::Spline(sh) => sketch.spline_periodic(*sh).unwrap_or(false),
                _ => false,
            };

            if is_full_circle {
                // Split the full circle into two semicircular arcs sharing
                // a midpoint at the opposite side of the circle. This avoids
                // self-loop edges which complicate Euler operators.
                if let ElementHandle::Arc(ah) = &element
                    && let Some((cx, cy, r, start_angle, end_angle)) = sketch.arc_params(*ah)
                {
                    let mid_angle = start_angle + (end_angle - start_angle) * 0.5;
                    let mid_pos = (cx + r * mid_angle.pcos(), cy + r * mid_angle.psin());

                    let start_node = find_or_insert_node(&mut nodes, start, tolerance);
                    let mid_node = find_or_insert_node(&mut nodes, mid_pos, tolerance);

                    // First half-arc: start -> midpoint
                    edges.push(GraphEdge {
                        element,
                        start_node,
                        end_node: mid_node,
                        start_pos: start,
                        end_pos: mid_pos,
                    });
                    // Second half-arc: midpoint -> start (closes the circle)
                    edges.push(GraphEdge {
                        element,
                        start_node: mid_node,
                        end_node: start_node,
                        start_pos: mid_pos,
                        end_pos: start,
                    });
                }
            } else if is_full_ellipse {
                // Split the full ellipse into four quarter-ellipses at the
                // major/minor axis endpoints (t=0, PI/2, PI, 3*PI/2). This
                // gives four graph edges forming a closed loop, avoiding the
                // self-loop problem and giving the loop tracer enough nodes.
                if let ElementHandle::Ellipse(eh) = &element
                    && let Some((cx, cy, semi_a, semi_b, rotation)) = sketch.ellipse_params(*eh)
                {
                    let cos_r = rotation.pcos();
                    let sin_r = rotation.psin();
                    // Evaluate ellipse at t=0, PI/2, PI, 3*PI/2.
                    let eval_ellipse = |t: f64| -> (f64, f64) {
                        let lx = semi_a * t.pcos();
                        let ly = semi_b * t.psin();
                        (cx + lx * cos_r - ly * sin_r, cy + lx * sin_r + ly * cos_r)
                    };
                    use std::f64::consts::PI;
                    let p0 = eval_ellipse(0.0);
                    let p1 = eval_ellipse(PI * 0.5);
                    let p2 = eval_ellipse(PI);
                    let p3 = eval_ellipse(PI * 1.5);

                    let n0 = find_or_insert_node(&mut nodes, p0, tolerance);
                    let n1 = find_or_insert_node(&mut nodes, p1, tolerance);
                    let n2 = find_or_insert_node(&mut nodes, p2, tolerance);
                    let n3 = find_or_insert_node(&mut nodes, p3, tolerance);

                    // Four quarter-ellipses: p0->p1, p1->p2, p2->p3, p3->p0
                    for (sn, en, sp, ep) in [
                        (n0, n1, p0, p1),
                        (n1, n2, p1, p2),
                        (n2, n3, p2, p3),
                        (n3, n0, p3, p0),
                    ] {
                        edges.push(GraphEdge {
                            element,
                            start_node: sn,
                            end_node: en,
                            start_pos: sp,
                            end_pos: ep,
                        });
                    }
                }
            } else if is_periodic_spline {
                // Split the periodic spline into four quarter-splines at
                // t=0, t=period/4, t=period/2, t=3*period/4. Using four
                // segments (rather than two halves) ensures the chord polygon
                // has non-zero signed area for reliable outer-boundary detection.
                if let ElementHandle::Spline(sh) = &element {
                    let split_ok = (|| -> Option<()> {
                        let knots = sketch.spline_knots(*sh)?;
                        let t0 = knots[0];
                        let t_end = knots[knots.len() - 1];
                        let period = t_end - t0;
                        let t1 = t0 + period * 0.25;
                        let t2 = t0 + period * 0.5;
                        let t3 = t0 + period * 0.75;
                        let p0 = start; // = spline_point_at(t0)
                        let p1 = sketch.spline_point_at(*sh, t1)?;
                        let p2 = sketch.spline_point_at(*sh, t2)?;
                        let p3 = sketch.spline_point_at(*sh, t3)?;
                        let n0 = find_or_insert_node(&mut nodes, p0, tolerance);
                        let n1 = find_or_insert_node(&mut nodes, p1, tolerance);
                        let n2 = find_or_insert_node(&mut nodes, p2, tolerance);
                        let n3 = find_or_insert_node(&mut nodes, p3, tolerance);
                        for (sn, en, sp, ep) in [
                            (n0, n1, p0, p1),
                            (n1, n2, p1, p2),
                            (n2, n3, p2, p3),
                            (n3, n0, p3, p0),
                        ] {
                            edges.push(GraphEdge {
                                element,
                                start_node: sn,
                                end_node: en,
                                start_pos: sp,
                                end_pos: ep,
                            });
                        }
                        Some(())
                    })();
                    if split_ok.is_none() {
                        warnings.push(SketchWarning::DegenerateEdge);
                    }
                }
            } else {
                warnings.push(SketchWarning::DegenerateEdge);
            }
            continue;
        }

        let start_node = find_or_insert_node(&mut nodes, start, tolerance);
        let end_node = find_or_insert_node(&mut nodes, end, tolerance);

        edges.push(GraphEdge {
            element,
            start_node,
            end_node,
            start_pos: start,
            end_pos: end,
        });
    }

    (EdgeGraph { nodes, edges }, warnings)
}

/// Find an existing node within tolerance, or insert a new one.
fn find_or_insert_node(
    nodes: &mut Vec<(f64, f64)>,
    pos: (f64, f64),
    tolerance: &Tolerance,
) -> usize {
    for (i, node) in nodes.iter().enumerate() {
        let dx = node.0 - pos.0;
        let dy = node.1 - pos.1;
        if (dx * dx + dy * dy).sqrt() < tolerance.linear() {
            return i;
        }
    }
    let idx = nodes.len();
    nodes.push(pos);
    idx
}

/// Compute the 2D endpoints of a sketch element.
///
/// Returns `None` for element types that cannot be converted (e.g.,
/// standalone points) or when required data is missing.
fn element_endpoints(sketch: &Sketch, element: &ElementHandle) -> Option<((f64, f64), (f64, f64))> {
    match element {
        ElementHandle::Line(h) => {
            let (start_pt, end_pt) = sketch.line_endpoints(*h)?;
            let start = sketch.point_position(start_pt)?;
            let end = sketch.point_position(end_pt)?;
            Some((start, end))
        }
        ElementHandle::Arc(h) => {
            let (cx, cy, r, start_angle, end_angle) = sketch.arc_params(*h)?;
            let start = (cx + r * start_angle.pcos(), cy + r * start_angle.psin());
            let end = (cx + r * end_angle.pcos(), cy + r * end_angle.psin());
            Some((start, end))
        }
        ElementHandle::Ellipse(h) => {
            let (cx, cy, semi_a, _semi_b, rotation) = sketch.ellipse_params(*h)?;
            // Ellipse is a full closed curve. Return coincident start/end at
            // the t=0 point on the major axis so the build_edge_graph loop
            // detects dist ≈ 0 and routes into the full-ellipse splitting path.
            let cos_r = rotation.pcos();
            let sin_r = rotation.psin();
            let x0 = cx + semi_a * cos_r;
            let y0 = cy + semi_a * sin_r;
            Some(((x0, y0), (x0, y0)))
        }
        ElementHandle::Parabola(h) => {
            let (vx, vy, focal_dist, rotation) = sketch.parabola_params(*h)?;
            // Parabola in local frame: y^2 = 4*f*x, parameterized as
            // (f*t^2, 2*f*t). Use t in [-1, 1] as a finite segment.
            let cos_r = rotation.pcos();
            let sin_r = rotation.psin();
            let t_start = -1.0;
            let t_end = 1.0;
            let local_start = (focal_dist * t_start * t_start, 2.0 * focal_dist * t_start);
            let local_end = (focal_dist * t_end * t_end, 2.0 * focal_dist * t_end);
            let start = rotate_and_translate(local_start, cos_r, sin_r, vx, vy);
            let end = rotate_and_translate(local_end, cos_r, sin_r, vx, vy);
            Some((start, end))
        }
        ElementHandle::Hyperbola(h) => {
            let (cx, cy, semi_a, semi_b, rotation) = sketch.hyperbola_params(*h)?;
            // Hyperbola parameterized as (semi_a * cosh(t), semi_b * sinh(t)).
            // Use t in [-1, 1] as a finite segment.
            let cos_r = rotation.pcos();
            let sin_r = rotation.psin();
            let t_start: f64 = -1.0;
            let t_end: f64 = 1.0;
            let local_start = (semi_a * t_start.pcosh(), semi_b * t_start.psinh());
            let local_end = (semi_a * t_end.pcosh(), semi_b * t_end.psinh());
            let start = rotate_and_translate(local_start, cos_r, sin_r, cx, cy);
            let end = rotate_and_translate(local_end, cos_r, sin_r, cx, cy);
            Some((start, end))
        }
        ElementHandle::Spline(h) => {
            let knots = sketch.spline_knots(*h)?;
            let degree = sketch.spline_degree(*h)?;
            let n_cp = sketch.spline_control_points(*h)?.len();
            let periodic = sketch.spline_periodic(*h)?;

            if periodic {
                // Periodic spline: endpoints are the same (closed curve).
                let t_start = knots[0];
                let start = sketch.spline_point_at(*h, t_start)?;
                // For a periodic spline, the curve wraps, so use the same
                // point for both endpoints (closed loop).
                Some((start, start))
            } else {
                // Non-periodic: evaluate at first and last knot of the
                // active range.
                let t_start = knots[degree];
                let t_end = knots[n_cp];
                let start = sketch.spline_point_at(*h, t_start)?;
                let end = sketch.spline_point_at(*h, t_end)?;
                Some((start, end))
            }
        }
        ElementHandle::Point(_) => None,
    }
}

/// Rotate a 2D point by angle (given as cos/sin) and translate.
fn rotate_and_translate(local: (f64, f64), cos_r: f64, sin_r: f64, tx: f64, ty: f64) -> (f64, f64) {
    (
        tx + local.0 * cos_r - local.1 * sin_r,
        ty + local.0 * sin_r + local.1 * cos_r,
    )
}

#[cfg(test)]
#[path = "edge_graph_tests.rs"]
mod tests;
