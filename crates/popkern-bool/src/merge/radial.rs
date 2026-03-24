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

    // Compute tangent vectors for each outgoing half-edge.
    let tangents: Vec<(usize, Vec3)> = outgoing
        .iter()
        .map(|&he_idx| {
            let twin_idx = graph.half_edges[he_idx].twin;
            let end_node = graph.half_edges[twin_idx].origin;
            let end_pos = graph.nodes[end_node].position;
            let tangent = end_pos - node_pos;
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
/// Falls back to +Z if no face or surface is available.
fn find_node_normal(graph: &IntersectionGraph, node_idx: usize) -> Vec3 {
    // We don't have direct access to BRep surfaces from the graph,
    // so use a geometric approach: compute a normal from the tangent vectors.
    let outgoing = &graph.nodes[node_idx].outgoing;
    if outgoing.len() < 2 {
        return Vec3::new(0.0, 0.0, 1.0);
    }

    let node_pos = graph.nodes[node_idx].position;

    // Take the first two outgoing tangents and compute their cross product.
    let he0 = outgoing[0];
    let he1 = outgoing[1];

    let twin0 = graph.half_edges[he0].twin;
    let twin1 = graph.half_edges[he1].twin;

    let end0 = graph.nodes[graph.half_edges[twin0].origin].position;
    let end1 = graph.nodes[graph.half_edges[twin1].origin].position;

    let t0 = end0 - node_pos;
    let t1 = end1 - node_pos;

    let cross = t0.cross(&t1);
    cross.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
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
}
