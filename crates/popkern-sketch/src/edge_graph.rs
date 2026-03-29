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

        // Filter degenerate edges.
        let dx = end.0 - start.0;
        let dy = end.1 - start.1;
        let dist = (dx * dx + dy * dy).sqrt();
        if dist < tolerance.linear {
            warnings.push(SketchWarning::DegenerateEdge);
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
        if (dx * dx + dy * dy).sqrt() < tolerance.linear {
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
            let start = (cx + r * start_angle.cos(), cy + r * start_angle.sin());
            let end = (cx + r * end_angle.cos(), cy + r * end_angle.sin());
            Some((start, end))
        }
        ElementHandle::Ellipse(h) => {
            let (cx, cy, semi_a, _semi_b, rotation) = sketch.ellipse_params(*h)?;
            // Ellipse is a full curve; endpoints at t=0 and t=2*PI
            let cos_r = rotation.cos();
            let sin_r = rotation.sin();
            // At t=0: (cx + semi_a * cos_r, cy + semi_a * sin_r)
            let start_x = cx + semi_a * cos_r;
            let start_y = cy + semi_a * sin_r;
            // At t=PI: (cx - semi_a * cos_r, cy - semi_a * sin_r)
            let end_x = cx - semi_a * cos_r;
            let end_y = cy - semi_a * sin_r;
            // Full ellipse is closed, so start == end. Use t=0 for both.
            // For the graph, a closed curve has coincident endpoints.
            Some(((start_x, start_y), (end_x, end_y)))
        }
        ElementHandle::Parabola(h) => {
            let (vx, vy, focal_dist, rotation) = sketch.parabola_params(*h)?;
            // Parabola in local frame: y^2 = 4*f*x, parameterized as
            // (f*t^2, 2*f*t). Use t in [-1, 1] as a finite segment.
            let cos_r = rotation.cos();
            let sin_r = rotation.sin();
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
            let cos_r = rotation.cos();
            let sin_r = rotation.sin();
            let t_start: f64 = -1.0;
            let t_end: f64 = 1.0;
            let local_start = (semi_a * t_start.cosh(), semi_b * t_start.sinh());
            let local_end = (semi_a * t_end.cosh(), semi_b * t_end.sinh());
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
mod tests {
    use super::*;

    #[test]
    fn rectangle_four_lines() {
        let mut sk = Sketch::new();
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(1.0, 0.0);
        let p2 = sk.add_point(1.0, 1.0);
        let p3 = sk.add_point(0.0, 1.0);
        let _l0 = sk.add_line(p0, p1);
        let _l1 = sk.add_line(p1, p2);
        let _l2 = sk.add_line(p2, p3);
        let _l3 = sk.add_line(p3, p0);

        let tol = Tolerance::default();
        let (graph, warnings) = build_edge_graph(&sk, &tol);

        assert!(warnings.is_empty());
        assert_eq!(graph.nodes.len(), 4);
        assert_eq!(graph.edges.len(), 4);
    }

    #[test]
    fn zero_length_line_filtered() {
        let mut sk = Sketch::new();
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(0.0, 0.0); // same position
        let _l = sk.add_line(p0, p1);

        let tol = Tolerance::default();
        let (graph, warnings) = build_edge_graph(&sk, &tol);

        assert_eq!(graph.edges.len(), 0);
        assert_eq!(warnings.len(), 1);
        assert!(matches!(warnings[0], SketchWarning::DegenerateEdge));
    }

    #[test]
    fn two_separate_shapes() {
        let mut sk = Sketch::new();
        // Shape 1: triangle
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(1.0, 0.0);
        let p2 = sk.add_point(0.5, 1.0);
        let _l0 = sk.add_line(p0, p1);
        let _l1 = sk.add_line(p1, p2);
        let _l2 = sk.add_line(p2, p0);

        // Shape 2: separate line segment
        let p3 = sk.add_point(10.0, 10.0);
        let p4 = sk.add_point(11.0, 10.0);
        let _l3 = sk.add_line(p3, p4);

        let tol = Tolerance::default();
        let (graph, warnings) = build_edge_graph(&sk, &tol);

        assert!(warnings.is_empty());
        assert_eq!(graph.edges.len(), 4); // 3 triangle + 1 separate
        assert_eq!(graph.nodes.len(), 5); // 3 triangle + 2 separate
    }

    #[test]
    fn construction_lines_excluded() {
        let mut sk = Sketch::new();
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(1.0, 0.0);
        let p2 = sk.add_point(1.0, 1.0);
        let l0 = sk.add_line(p0, p1);
        let _l1 = sk.add_line(p1, p2);

        // Mark first line as construction
        sk.set_construction(ElementHandle::Line(l0), true);

        let tol = Tolerance::default();
        let (graph, warnings) = build_edge_graph(&sk, &tol);

        assert!(warnings.is_empty());
        assert_eq!(graph.edges.len(), 1);
    }

    #[test]
    fn arc_endpoints_computed_correctly() {
        use std::f64::consts::{FRAC_PI_2, PI};

        let mut sk = Sketch::new();
        // Arc centered at origin, radius 1, from 0 to PI/2
        let _arc = sk.add_arc(0.0, 0.0, 1.0, 0.0, FRAC_PI_2);

        let tol = Tolerance::default();
        let (graph, warnings) = build_edge_graph(&sk, &tol);

        assert!(warnings.is_empty());
        assert_eq!(graph.edges.len(), 1);

        let edge = &graph.edges[0];
        // Start at angle 0: (1, 0)
        assert!((edge.start_pos.0 - 1.0).abs() < 1e-10);
        assert!((edge.start_pos.1 - 0.0).abs() < 1e-10);
        // End at angle PI/2: (0, 1)
        assert!((edge.end_pos.0 - 0.0).abs() < 1e-10);
        assert!((edge.end_pos.1 - 1.0).abs() < 1e-10);

        // Full circle arc (start == end, should be degenerate check)
        let mut sk2 = Sketch::new();
        let _arc2 = sk2.add_arc(0.0, 0.0, 1.0, 0.0, 2.0 * PI);
        let (graph2, _) = build_edge_graph(&sk2, &tol);
        // 2*PI arc: start and end are extremely close but cos(2pi) ~= 1,
        // sin(2pi) ~= 0, same as start. Will be degenerate.
        assert_eq!(graph2.edges.len(), 0);
    }

    #[test]
    fn node_merging_within_tolerance() {
        let mut sk = Sketch::new();
        // Two lines sharing an endpoint, but with slight numerical noise
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(1.0, 0.0);
        let p2 = sk.add_point(1.0 + 1e-9, 0.0); // very close to p1
        let p3 = sk.add_point(2.0, 0.0);
        let _l0 = sk.add_line(p0, p1);
        let _l1 = sk.add_line(p2, p3);

        let tol = Tolerance::default(); // linear = 1e-6
        let (graph, warnings) = build_edge_graph(&sk, &tol);

        assert!(warnings.is_empty());
        assert_eq!(graph.edges.len(), 2);
        // p1 and p2 should merge into the same node
        assert_eq!(graph.nodes.len(), 3);
        // The two edges should share a node
        assert_eq!(graph.edges[0].end_node, graph.edges[1].start_node);
    }
}
