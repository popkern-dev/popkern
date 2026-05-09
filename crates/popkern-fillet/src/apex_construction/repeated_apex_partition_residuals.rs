//! Retired repeated-apex residual topology handling.

use super::*;

pub(super) fn residual_strip_loop_with_duplicated_seam(
    brep: &mut BRep,
    ownership: &FilletEndOwnership,
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    let shared = duplicate_boundary_half_edge_for_residual_strip(
        brep,
        ownership.cylinder_candidate_half_edge,
        ownership.fillet_index,
    )?;
    Ok(vec![
        shared,
        ownership.retired_apex_cell.t_a_to_apex,
        ownership.retired_apex_cell.apex_to_t_b,
    ])
}

pub(super) fn demote_retired_apex_residual_faces(
    brep: &mut BRep,
    cap_face: FaceHandle,
    apex_vertex: VertexHandle,
    retired_spans: &[RetiredApexSpan],
    fillet_ends: &[FilletEndOwnership],
) {
    let mut residual_faces = std::collections::BTreeSet::new();
    for span in retired_spans {
        for half_edge in &span.half_edges {
            if let Some(face) = face_for_half_edge(brep, *half_edge) {
                residual_faces.insert(face);
            }
        }
    }
    for ownership in fillet_ends {
        for half_edge in [
            ownership.retired_apex_cell.t_a_to_apex,
            ownership.retired_apex_cell.apex_to_t_b,
        ] {
            if let Some(face) = face_for_half_edge(brep, half_edge) {
                residual_faces.insert(face);
            }
        }
    }
    let target_cylinder_faces = fillet_ends
        .iter()
        .map(|ownership| ownership.target_cylinder_face)
        .collect::<std::collections::BTreeSet<_>>();
    residual_faces.extend(brep.faces.iter().filter_map(|(face, face_data)| {
        if face_data.surface.is_none() || face == cap_face || target_cylinder_faces.contains(&face)
        {
            return None;
        }
        face_contains_vertex(brep, face, apex_vertex).then_some(face)
    }));
    for face in residual_faces {
        if face == cap_face {
            continue;
        }
        if face_contains_vertex(brep, face, apex_vertex) {
            brep.faces[face].surface = None;
        }
    }
}

pub(super) fn validate_existing_target_cylinder_strip(
    brep: &BRep,
    ownership: &FilletEndOwnership,
    apex_vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(), ApexBuildError> {
    let face = brep
        .faces
        .get(ownership.target_cylinder_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] selected stale target cylinder \
                 face {:?} as the authoritative strip",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        })?;
    let cycle = collect_loop_half_edges(brep, face.outer_loop)?;
    validate_half_edge_cycle(
        brep,
        &cycle,
        "authoritative target-cylinder strip loop",
        ownership.target_cylinder_face,
    )?;
    if cycle.iter().any(|half_edge| {
        brep.half_edges
            .get(*half_edge)
            .is_some_and(|data| data.origin == apex_vertex)
    }) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} still \
                 retains apex vertex {apex_vertex:?} after the cylinder trim rewrite",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        });
    }
    let expected = fillet_surfaces.get(ownership.fillet_index).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] has no fillet surface for \
                 target-cylinder strip validation",
                ownership.fillet_index
            ),
        }
    })?;
    let actual = face
        .surface
        .as_ref()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} has no \
                 production surface after the cylinder trim rewrite",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        })?;
    if std::mem::discriminant(&actual.geometry) != std::mem::discriminant(&expected.surface) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target face {:?} is not the \
                 expected cylinder strip surface",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        });
    }
    Ok(())
}

fn duplicate_boundary_half_edge_for_residual_strip(
    brep: &mut BRep,
    source: HalfEdgeHandle,
    fillet_index: usize,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let source_data =
        brep.half_edges
            .get(source)
            .cloned()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot duplicate \
                     stale residual strip seam {source:?}"
                ),
            })?;
    let curve = brep
        .edges
        .get(source_data.edge)
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot duplicate \
                 residual strip seam {source:?}: source edge has no curve"
            ),
        })?;
    let (edge, half_edge) = insert_boundary_half_edge(brep, source_data.origin);
    brep.edges[edge].curve = Some(curve);
    Ok(half_edge)
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    let loop_handle = brep.half_edges.get(half_edge)?.loop_;
    brep.loops.get(loop_handle).map(|loop_data| loop_data.face)
}

fn face_contains_vertex(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    brep.faces
        .get(face)
        .and_then(|face_data| collect_loop_half_edges(brep, face_data.outer_loop).ok())
        .is_some_and(|cycle| {
            cycle.iter().any(|half_edge| {
                brep.half_edges
                    .get(*half_edge)
                    .is_some_and(|data| data.origin == vertex)
            })
        })
}
