//! Splitting inserted records at the T-junctions the vertex table
//! decided (plan §4.3: an interior incidence in a group is a pave on
//! that record).
//!
//! The graph's own pass (`split_shared_edges_at_nodes`) finds the same
//! junctions after the fact, testing every node against every record
//! edge; with the table deciding vertices, each record arrives knowing
//! its paves, and the pass is left to the graph's own identity (a
//! declined table, or `POPKERN_VERTEX_TABLE=0`). Common-block scenes
//! once kept it (suite 22: every junction the table missed was in
//! one); the misses were operand edges split at edge-face hits off
//! their face, which the graph now leaves whole
//! (`candidates::interfering_hits`).

use popkern_core::geometry::Point3;

use popkern_core::topology::{BRep, EdgeHandle};

use super::graph::{EdgeSource, IntersectionGraph, NodeSource};

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

/// Split operand edges where records end on them mid-span: each of
/// `paves` names the operand (0 = A), the edge and the vertex. The
/// vertex's node is the graph node within `tolerance` of it (a record
/// end, so it exists); the piece of the edge it lies on, by the graph's
/// own on-edge test, is split there with its pcurves cut (`breps`).
/// Returns the splits made.
pub(crate) fn split_operand_edges_at_paves(
    graph: &mut IntersectionGraph,
    paves: &[(usize, EdgeHandle, Point3)],
    tolerance: f64,
    breps: Option<(&BRep, &BRep, f64)>,
) -> usize {
    let mut splits = 0;
    for &(operand, edge, at) in paves {
        let Some(node) = graph.find_node(&at, tolerance) else {
            continue;
        };
        let pos = graph.nodes[node].position;
        let piece = (0..graph.half_edges.len()).find(|&i| {
            let he = &graph.half_edges[i];
            let of_edge = match he.source {
                EdgeSource::OriginalA(e) | EdgeSource::SplitA { original: e } => {
                    operand == 0 && e == edge
                }
                EdgeSource::OriginalB(e) | EdgeSource::SplitB { original: e } => {
                    operand == 1 && e == edge
                }
                EdgeSource::Shared { .. } => false,
            };
            of_edge
                && i < he.twin
                && he.origin != node
                && graph.half_edges[he.twin].origin != node
                && super::graph_split_shared::node_on_edge(graph, i, pos, tolerance)
        });
        if let Some(i) = piece
            && super::graph_split_node::split_at_existing_node_with_pcurves(graph, i, node, breps)
        {
            splits += 1;
        }
    }
    splits
}
