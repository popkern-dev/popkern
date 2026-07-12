//! Cleanup for unowned construction half-edges left by repeated-apex rewrites.

use super::*;

pub(super) fn materialize_unowned_half_edge_residual_faces(
    brep: &mut BRep,
    template_face: FaceHandle,
) -> Result<(), ApexBuildError> {
    let half_edges = brep
        .half_edges
        .iter()
        .filter_map(|(half_edge, data)| {
            // Already-retired records (residual edges and residual
            // faces' loops) are excused by role — collecting them
            // again would mint bigons for prior sweeps' connectors.
            let retired = brep
                .edges
                .get(data.edge)
                .is_some_and(|edge| edge.role == TopologyRole::ConstructionResidual)
                || brep
                    .loops
                    .get(data.loop_)
                    .and_then(|loop_data| brep.faces.get(loop_data.face))
                    .is_some_and(|face| face.role == TopologyRole::ConstructionResidual);
            // Null-loop records are collected even when retired —
            // referential integrity (`validate_topology_references`)
            // has no role exemption; retired ones get homed onto the
            // twin's loop instead of bigoned.
            (data.loop_ == LoopHandle::default()
                || (!retired && !half_edge_is_cycle_member(brep, half_edge)))
            .then_some(half_edge)
        })
        .collect::<std::collections::BTreeSet<_>>();
    let displaced = half_edges.len();
    for half_edge in half_edges {
        let null_loop = brep
            .half_edges
            .get(half_edge)
            .is_some_and(|data| data.loop_ == LoopHandle::default());
        if !null_loop && half_edge_is_cycle_member(brep, half_edge) {
            continue;
        }
        if null_loop {
            // An unwired side of an already-retired record: home it
            // on its twin's loop for referential integrity; the
            // residual role excuses it from cycle checks there.
            let retired_twin_loop = brep.half_edges.get(half_edge).and_then(|data| {
                let edge_retired = brep
                    .edges
                    .get(data.edge)
                    .is_some_and(|edge| edge.role == TopologyRole::ConstructionResidual);
                let twin_loop = brep.half_edges.get(data.twin).map(|twin| twin.loop_)?;
                (edge_retired && brep.loops.contains_key(twin_loop)).then_some(twin_loop)
            });
            if let Some(twin_loop) = retired_twin_loop {
                if let Some(data) = brep.half_edges.get_mut(half_edge) {
                    data.loop_ = twin_loop;
                }
                continue;
            }
        }
        if !null_loop && half_edge_end_vertex(brep, half_edge).is_none() {
            // A displaced record whose twin is already gone cannot
            // anchor a residual bigon; retire it by role so cycle
            // validation and seam scans excuse it like the collected
            // ones.
            if let Some(edge) = brep.half_edges.get(half_edge).map(|data| data.edge)
                && let Some(edge_data) = brep.edges.get_mut(edge)
            {
                edge_data.role = TopologyRole::ConstructionResidual;
            }
            continue;
        }
        materialize_unowned_half_edge_residual_face(brep, template_face, half_edge)?;
    }
    if displaced > 0 && crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/partition] collected {displaced} unowned/displaced \
             half-edge(s) into residual faces"
        );
    }
    Ok(())
}

/// True when the half-edge is reachable by walking `next` from its
/// loop's head. Boundary rewrites re-home only the NEW cycle's
/// half-edges, so a displaced side keeps a stale `loop_` that the
/// loop's cycle no longer traverses (measured on cube_extreme: 12
/// phantom cap-side chord halves counted as cap/support seams). Such
/// records are unowned in fact and must be collected exactly like
/// null-loop ones.
fn half_edge_is_cycle_member(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some(data) = brep.half_edges.get(half_edge) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(data.loop_) else {
        return false;
    };
    let head = loop_data.half_edge;
    let mut cursor = head;
    for _ in 0..brep.half_edges.len() {
        if cursor == half_edge {
            return true;
        }
        let Some(c) = brep.half_edges.get(cursor) else {
            return false;
        };
        cursor = c.next;
        if cursor == head {
            return false;
        }
    }
    false
}

fn materialize_unowned_half_edge_residual_face(
    brep: &mut BRep,
    template_face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> Result<(), ApexBuildError> {
    let start = brep
        .half_edges
        .get(half_edge)
        .map(|data| data.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot materialize stale unowned half-edge \
                 {half_edge:?}"
            ),
        })?;
    let end = half_edge_end_vertex(brep, half_edge).ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition cannot materialize unowned half-edge {half_edge:?}: \
             missing endpoint"
        ),
    })?;
    let (connector_edge, connector, connector_twin) = insert_unwired_edge_pair(brep, end, start);
    brep.edges[connector_edge].role = TopologyRole::ConstructionResidual;
    let residual_face = insert_face_with_outer_loop(brep, template_face, &[half_edge, connector])?;
    brep.faces[residual_face].role = TopologyRole::ConstructionResidual;
    brep.faces[residual_face].surface = None;
    brep.edges[connector_edge].role = TopologyRole::ConstructionResidual;
    // The connector twin closes no cycle, but a final model must be
    // referentially intact (`validate_topology_references`): home it
    // on the residual face's loop, where the ConstructionResidual
    // role excuses it from cycle checks like the connector itself.
    let residual_loop = brep.faces[residual_face].outer_loop;
    if let Some(twin_data) = brep.half_edges.get_mut(connector_twin) {
        twin_data.loop_ = residual_loop;
    }
    Ok(())
}
