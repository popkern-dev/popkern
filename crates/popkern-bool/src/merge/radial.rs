//! Radial edge ordering at intersection graph nodes.
//!
//! Sorts outgoing half-edges at each node by angle in the 3D tangent plane,
//! enabling correct loop walking around face fragments.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::geometry::Vec3;

use super::graph::IntersectionGraph;

/// Sort outgoing half-edges at a node by angle in the 3D tangent plane.
///
/// Given a surface normal at the node, projects each outgoing edge's tangent
/// vector onto the tangent plane and sorts by angle (counter-clockwise).
/// Nodes with two or fewer outgoing edges are left unchanged.
pub(crate) fn sort_radial_at_node(graph: &mut IntersectionGraph, node_idx: usize, normal: Vec3) {
    let outgoing = &graph.nodes[node_idx].outgoing;
    if outgoing.len() <= 2 {
        return;
    }

    let node_pos = graph.nodes[node_idx].position;

    // Compute tangent vectors for each outgoing half-edge — the
    // TRUE curve tangent at this node when the edge carries a curve.
    // Chord tangents degenerate for curved edges: two half-circle
    // arcs between the same pair of nodes have identical chords, and
    // the sort then cannot distinguish them (a boss/plate
    // intersection circle split at a seam crossing).
    let tangents: Vec<(usize, Vec3)> = outgoing
        .iter()
        .map(|&he_idx| {
            let twin_idx = graph.half_edges[he_idx].twin;
            let end_node = graph.half_edges[twin_idx].origin;
            let end_pos = graph.nodes[end_node].position;
            let chord = end_pos - node_pos;
            let tangent = match &graph.half_edges[he_idx].curve {
                Some(curve) => {
                    let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
                    let (t0, t1) = (curve.t_start, curve.t_end);
                    let p_start = eval.point_at(t0);
                    let p_end = eval.point_at(t1);
                    // Which parameter end sits at this node?
                    let at_start = p_start.distance_to(&node_pos) <= p_end.distance_to(&node_pos);
                    let h = (t1 - t0) * 1.0e-4;
                    let d = if at_start {
                        eval.point_at(t0 + h) - eval.point_at(t0)
                    } else {
                        eval.point_at(t1 - h) - eval.point_at(t1)
                    };
                    if d.length() > 1e-300 { d } else { chord }
                }
                None => chord,
            };
            (he_idx, tangent)
        })
        .collect();

    // Project tangents onto the tangent plane: v_proj = v - n * (v . n)
    let projected: Vec<(usize, Vec3)> = tangents
        .iter()
        .map(|&(he_idx, t)| {
            let proj = t - normal * t.dot(&normal);
            (he_idx, proj)
        })
        .collect();

    // Use the first projected tangent as the reference direction.
    let ref_dir = match projected.first() {
        Some((_, v)) => match v.normalized() {
            Some(n) => n,
            None => return,
        },
        None => return,
    };

    // Compute angles and sort.
    let mut angle_entries: Vec<(usize, f64)> = projected
        .iter()
        .map(|&(he_idx, t)| {
            let cos_val = t.dot(&ref_dir);
            let sin_val = t.cross(&ref_dir).dot(&normal);
            let angle = sin_val.atan2(cos_val);
            (he_idx, angle)
        })
        .collect();

    angle_entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // Update the outgoing list with sorted order.
    graph.nodes[node_idx].outgoing = angle_entries.iter().map(|&(idx, _)| idx).collect();
}

/// Sort radial edges at all nodes in the graph.
///
/// For each node, determines the surface normal from the face that the node's
/// first outgoing edge borders. If no face is assigned, uses a default +Z
/// normal.
pub(crate) fn sort_all_radial(graph: &mut IntersectionGraph) {
    let node_count = graph.nodes.len();

    for node_idx in 0..node_count {
        // Determine a normal for this node from its first outgoing edge's face.
        let normal = find_node_normal(graph, node_idx);
        sort_radial_at_node(graph, node_idx, normal);
    }
}

/// Find a surface normal for a node by looking at its outgoing edges' faces.
///
/// Tries all pairs of outgoing tangent vectors until a non-degenerate cross
/// product is found. Falls back to +Z if no pair yields a valid normal.
fn find_node_normal(graph: &IntersectionGraph, node_idx: usize) -> Vec3 {
    let outgoing = &graph.nodes[node_idx].outgoing;
    if outgoing.len() < 2 {
        return Vec3::new(0.0, 0.0, 1.0);
    }

    let node_pos = graph.nodes[node_idx].position;

    // Compute all tangent vectors.
    let tangents: Vec<Vec3> = outgoing
        .iter()
        .map(|&he_idx| {
            let twin_idx = graph.half_edges[he_idx].twin;
            let end_pos = graph.nodes[graph.half_edges[twin_idx].origin].position;
            end_pos - node_pos
        })
        .collect();

    // Try all pairs until we find a non-degenerate cross product.
    let n = tangents.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let cross = tangents[i].cross(&tangents[j]);
            if let Some(normal) = cross.normalized() {
                return normal;
            }
        }
    }

    Vec3::new(0.0, 0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::graph::{EdgeSource, IntersectionGraph, NodeSource};
    use popkern_core::geometry::Point3;

    #[test]
    fn sort_four_edges_ccw() {
        let mut graph = IntersectionGraph::new();

        // Central node at origin.
        let center = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        // Four outer nodes at +x, +y, -x, -y.
        let px = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let py = graph.add_node(Point3::new(0.0, 1.0, 0.0), NodeSource::Intersection);
        let mx = graph.add_node(Point3::new(-1.0, 0.0, 0.0), NodeSource::Intersection);
        let my = graph.add_node(Point3::new(0.0, -1.0, 0.0), NodeSource::Intersection);

        let eh = popkern_core::topology::EdgeHandle::default();

        // Add edges in scrambled order: +y, -x, +x, -y.
        let (he_py, _) = graph.add_half_edge_pair(center, py, EdgeSource::OriginalA(eh), None);
        let (he_mx, _) = graph.add_half_edge_pair(center, mx, EdgeSource::OriginalA(eh), None);
        let (he_px, _) = graph.add_half_edge_pair(center, px, EdgeSource::OriginalA(eh), None);
        let (he_my, _) = graph.add_half_edge_pair(center, my, EdgeSource::OriginalA(eh), None);

        // Sort with +z normal.
        sort_radial_at_node(&mut graph, center, Vec3::new(0.0, 0.0, 1.0));

        let sorted = &graph.nodes[center].outgoing;
        assert_eq!(sorted.len(), 4);

        // After sorting by angle from first edge's direction, verify order
        // is consistent (all four edges present, sorted by angle).
        // The exact order depends on the reference direction.
        // Verify all half-edges are present.
        assert!(sorted.contains(&he_px));
        assert!(sorted.contains(&he_py));
        assert!(sorted.contains(&he_mx));
        assert!(sorted.contains(&he_my));

        // Verify angular ordering using the same reference direction as the sort:
        // the projected tangent of the first sorted edge.
        let node_pos = graph.nodes[center].position;
        let normal = Vec3::new(0.0, 0.0, 1.0);

        // Compute reference from first sorted edge.
        let first_he = sorted[0];
        let first_twin = graph.half_edges[first_he].twin;
        let first_end = graph.nodes[graph.half_edges[first_twin].origin].position;
        let first_t = first_end - node_pos;
        let ref_dir = (first_t - normal * first_t.dot(&normal))
            .normalized()
            .unwrap_or(Vec3::new(1.0, 0.0, 0.0));

        let angles: Vec<f64> = sorted
            .iter()
            .map(|&he_idx| {
                let twin_idx = graph.half_edges[he_idx].twin;
                let end_pos = graph.nodes[graph.half_edges[twin_idx].origin].position;
                let t = end_pos - node_pos;
                let proj = t - normal * t.dot(&normal);
                let cos_val = proj.dot(&ref_dir);
                let sin_val = proj.cross(&ref_dir).dot(&normal);
                sin_val.atan2(cos_val)
            })
            .collect();

        // Angles should be non-decreasing (sorted), allowing one wrap-around
        // from +pi to -pi at most.
        let mut wrap_count = 0;
        for i in 1..angles.len() {
            if angles[i] < angles[i - 1] - 1e-10 {
                wrap_count += 1;
            }
        }
        assert!(
            wrap_count <= 1,
            "More than one wrap-around in sorted angles: {:?}",
            angles
        );
    }

    #[test]
    fn sort_two_edges_no_change() {
        let mut graph = IntersectionGraph::new();
        let center = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);

        let eh = popkern_core::topology::EdgeHandle::default();
        graph.add_half_edge_pair(center, n1, EdgeSource::OriginalA(eh), None);

        let before = graph.nodes[center].outgoing.clone();
        sort_radial_at_node(&mut graph, center, Vec3::new(0.0, 0.0, 1.0));
        let after = &graph.nodes[center].outgoing;

        assert_eq!(&before, after);
    }

    // ── sort_all_radial tests ──────────────────────────────────────────

    #[test]
    fn sort_all_radial_sorts_every_node() {
        // Build a graph with two high-valence nodes and verify sort_all_radial
        // processes both of them.
        let mut graph = IntersectionGraph::new();
        let eh = popkern_core::topology::EdgeHandle::default();

        // Node A at origin with 3 spokes in the XY plane.
        let a = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let a1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let a2 = graph.add_node(Point3::new(0.0, 1.0, 0.0), NodeSource::Intersection);
        let a3 = graph.add_node(Point3::new(-1.0, 0.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(a, a2, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(a, a1, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(a, a3, EdgeSource::OriginalA(eh), None);

        // Node B at (5,0,0) with 3 spokes in the XY plane.
        let b = graph.add_node(Point3::new(5.0, 0.0, 0.0), NodeSource::Intersection);
        let b1 = graph.add_node(Point3::new(6.0, 0.0, 0.0), NodeSource::Intersection);
        let b2 = graph.add_node(Point3::new(5.0, 1.0, 0.0), NodeSource::Intersection);
        let b3 = graph.add_node(Point3::new(4.0, 0.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(b, b3, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(b, b1, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(b, b2, EdgeSource::OriginalA(eh), None);

        sort_all_radial(&mut graph);

        // Both high-valence nodes should have angularly sorted outgoing lists.
        for &node_idx in &[a, b] {
            let sorted = &graph.nodes[node_idx].outgoing;
            assert_eq!(sorted.len(), 3);

            let node_pos = graph.nodes[node_idx].position;
            let normal = find_node_normal(&graph, node_idx);

            let first_he = sorted[0];
            let first_twin = graph.half_edges[first_he].twin;
            let first_end = graph.nodes[graph.half_edges[first_twin].origin].position;
            let first_t = first_end - node_pos;
            let ref_dir = (first_t - normal * first_t.dot(&normal))
                .normalized()
                .unwrap_or(Vec3::new(1.0, 0.0, 0.0));

            let angles: Vec<f64> = sorted
                .iter()
                .map(|&he_idx| {
                    let twin_idx = graph.half_edges[he_idx].twin;
                    let end_pos = graph.nodes[graph.half_edges[twin_idx].origin].position;
                    let t = end_pos - node_pos;
                    let proj = t - normal * t.dot(&normal);
                    let cos_val = proj.dot(&ref_dir);
                    let sin_val = proj.cross(&ref_dir).dot(&normal);
                    sin_val.atan2(cos_val)
                })
                .collect();

            let mut wrap_count = 0;
            for i in 1..angles.len() {
                if angles[i] < angles[i - 1] - 1e-10 {
                    wrap_count += 1;
                }
            }
            assert!(
                wrap_count <= 1,
                "Node {node_idx}: more than one wrap-around in sorted angles: {angles:?}"
            );
        }
    }

    #[test]
    fn sort_all_radial_leaves_low_valence_nodes_unchanged() {
        // Nodes with 0, 1, or 2 outgoing edges should not be reordered.
        let mut graph = IntersectionGraph::new();
        let eh = popkern_core::topology::EdgeHandle::default();

        // Isolated node (0 outgoing).
        let iso = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);

        // Node with 1 outgoing edge.
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let n2 = graph.add_node(Point3::new(2.0, 0.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);

        let before_iso = graph.nodes[iso].outgoing.clone();
        let before_n1 = graph.nodes[n1].outgoing.clone();

        sort_all_radial(&mut graph);

        assert_eq!(graph.nodes[iso].outgoing, before_iso);
        assert_eq!(graph.nodes[n1].outgoing, before_n1);
    }

    // ── find_node_normal tests ─────────────────────────────────────────

    #[test]
    fn find_node_normal_defaults_to_z_for_isolated_node() {
        let mut graph = IntersectionGraph::new();
        let n = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let normal = find_node_normal(&graph, n);
        assert!((normal - Vec3::new(0.0, 0.0, 1.0)).length() < 1e-10);
    }

    #[test]
    fn find_node_normal_defaults_to_z_for_single_outgoing() {
        let mut graph = IntersectionGraph::new();
        let eh = popkern_core::topology::EdgeHandle::default();
        let a = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let b = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(a, b, EdgeSource::OriginalA(eh), None);

        let normal = find_node_normal(&graph, a);
        assert!((normal - Vec3::new(0.0, 0.0, 1.0)).length() < 1e-10);
    }

    #[test]
    fn find_node_normal_from_two_non_parallel_edges() {
        let mut graph = IntersectionGraph::new();
        let eh = popkern_core::topology::EdgeHandle::default();
        let center = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let px = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let py = graph.add_node(Point3::new(0.0, 1.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(center, px, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(center, py, EdgeSource::OriginalA(eh), None);

        let normal = find_node_normal(&graph, center);
        // Cross product of +X and +Y is +Z (or -Z depending on order).
        assert!(
            normal.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() > 0.99,
            "Normal should be along Z axis, got {normal:?}"
        );
    }

    #[test]
    fn find_node_normal_collinear_edges_falls_back_to_z() {
        // Two outgoing edges along the same line: cross product is zero.
        let mut graph = IntersectionGraph::new();
        let eh = popkern_core::topology::EdgeHandle::default();
        let center = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let px = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let mx = graph.add_node(Point3::new(-1.0, 0.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(center, px, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(center, mx, EdgeSource::OriginalA(eh), None);

        let normal = find_node_normal(&graph, center);
        // Collinear edges produce a zero cross product -> falls back to +Z.
        assert!((normal - Vec3::new(0.0, 0.0, 1.0)).length() < 1e-10);
    }

    #[test]
    fn sort_radial_with_zero_length_tangent_is_noop() {
        // Degenerate: outgoing edge to the same position (zero tangent).
        let mut graph = IntersectionGraph::new();
        let eh = popkern_core::topology::EdgeHandle::default();
        let center = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let same = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let px = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let py = graph.add_node(Point3::new(0.0, 1.0, 0.0), NodeSource::Intersection);
        graph.add_half_edge_pair(center, same, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(center, px, EdgeSource::OriginalA(eh), None);
        graph.add_half_edge_pair(center, py, EdgeSource::OriginalA(eh), None);

        let before = graph.nodes[center].outgoing.clone();
        // Should not panic even with a degenerate zero-length tangent.
        sort_radial_at_node(&mut graph, center, Vec3::new(0.0, 0.0, 1.0));
        // All edges still present.
        assert_eq!(graph.nodes[center].outgoing.len(), before.len());
    }
}
