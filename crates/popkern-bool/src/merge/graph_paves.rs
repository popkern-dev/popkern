//! Splitting inserted records at the T-junctions the vertex table
//! decided (plan §4.3: an interior incidence in a group is a pave on
//! that record).
//!
//! The graph's own pass (`split_shared_edges_at_nodes`) finds the same
//! junctions after the fact, testing every node against every record
//! edge; with the table deciding vertices, each record arrives knowing
//! its paves, and the pass is left only to scenes with common blocks,
//! whose same-domain junctions the table does not yet decide (plan
//! §4.4; measured, suite 22: every junction the table missed was in
//! such a scene).

use popkern_core::geometry::Point3;

use super::graph::{IntersectionGraph, NodeSource};

/// Split each inserted record at its paves: `records` holds a record's
/// forward half-edge and its paves by ascending curve parameter. A
/// pave's node is the graph node within `tolerance` of the pave's
/// anchor, minted there when none is. Returns the splits made.
pub(crate) fn split_records_at_paves(
    graph: &mut IntersectionGraph,
    records: &[(usize, &[(f64, Point3)])],
    tolerance: f64,
) -> usize {
    let mut splits = 0;
    for &(fwd, paves) in records {
        // Each split keeps the head on the parent half-edge and puts
        // the tail on a new pair, so the next pave (further along)
        // lies on the newest tail.
        let mut tail = fwd;
        for &(_, anchor) in paves {
            let node = graph
                .find_node(&anchor, tolerance)
                .unwrap_or_else(|| graph.add_node(anchor, NodeSource::Intersection));
            let twin = graph.half_edges[tail].twin;
            if graph.half_edges[tail].origin == node || graph.half_edges[twin].origin == node {
                continue;
            }
            let next = graph.half_edges.len();
            if super::graph_split_node::split_at_existing_node(graph, tail, node) {
                tail = next;
                splits += 1;
            }
        }
    }
    splits
}
