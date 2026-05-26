//! Cleanup for unowned construction half-edges left by repeated-apex rewrites.

use super::*;

pub(super) fn materialize_unowned_half_edge_residual_faces(
    brep: &mut BRep,
    template_face: FaceHandle,
) -> Result<(), ApexBuildError> {
    let half_edges = brep
        .half_edges
        .iter()
        .filter_map(|(half_edge, data)| (data.loop_ == LoopHandle::default()).then_some(half_edge))
        .collect::<std::collections::BTreeSet<_>>();
    for half_edge in half_edges {
        if brep
            .half_edges
            .get(half_edge)
            .is_none_or(|data| data.loop_ != LoopHandle::default())
        {
            continue;
        }
        materialize_unowned_half_edge_residual_face(brep, template_face, half_edge)?;
    }
    Ok(())
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
    let (connector_edge, connector, _connector_twin) = insert_unwired_edge_pair(brep, end, start);
    brep.edges[connector_edge].role = TopologyRole::ConstructionResidual;
    let residual_face = insert_face_with_outer_loop(brep, template_face, &[half_edge, connector])?;
    brep.faces[residual_face].role = TopologyRole::ConstructionResidual;
    brep.faces[residual_face].surface = None;
    brep.edges[connector_edge].role = TopologyRole::ConstructionResidual;
    Ok(())
}
