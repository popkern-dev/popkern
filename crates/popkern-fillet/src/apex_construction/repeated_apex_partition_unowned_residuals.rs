//! Cleanup for unowned construction half-edges left by repeated-apex rewrites.

use super::*;

pub(super) fn mark_unowned_half_edges_construction_residual(brep: &mut BRep) {
    let edges = brep
        .half_edges
        .iter()
        .filter_map(|(_, half_edge)| {
            (half_edge.loop_ == LoopHandle::default()).then_some(half_edge.edge)
        })
        .collect::<std::collections::BTreeSet<_>>();
    for edge in edges {
        if let Some(edge_data) = brep.edges.get_mut(edge) {
            edge_data.role = TopologyRole::ConstructionResidual;
        }
    }
}
