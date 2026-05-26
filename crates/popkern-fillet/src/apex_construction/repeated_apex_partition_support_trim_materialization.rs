use super::partition_support_boundary_claim::claim_existing_support_boundary_for_cylinder_trim;
use super::partition_support_interior_trim_patch::{
    InteriorSupportTrimPatchRequest, author_interior_support_trim_patch_for_cylinder_trim,
};
use super::*;

pub(super) struct GeneratorTrimMaterializationRequest<'a> {
    pub(super) cylinder_face: FaceHandle,
    pub(super) support_face: FaceHandle,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) label: &'static str,
    pub(super) fillet_index: usize,
    pub(super) side_support_bridge_loops: &'a [LoopHandle],
}

pub(super) fn materialized_generator_trim_side(
    brep: &mut BRep,
    request: GeneratorTrimMaterializationRequest<'_>,
) -> Result<MaterializedGeneratorTrimSide, ApexBuildError> {
    let GeneratorTrimMaterializationRequest {
        cylinder_face,
        support_face,
        start,
        end,
        label,
        fillet_index,
        side_support_bridge_loops,
    } = request;
    let curve =
        exact_cylinder_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index)?;
    if let Some((half_edge, actual_support_face)) =
        claim_existing_support_boundary_for_cylinder_trim(
            brep,
            support_face,
            start,
            end,
            curve.clone(),
            label,
            fillet_index,
        )?
    {
        return Ok(MaterializedGeneratorTrimSide {
            half_edge,
            support_face: actual_support_face,
            support_side_owned: true,
            support_source: None,
            authored_boundary_spans: Vec::new(),
        });
    }
    if let Some(claim) = claim_side_support_bridge_boundary_for_cylinder_trim(
        brep,
        SideSupportBridgeTrimClaimRequest {
            support_face,
            bridge_loops: side_support_bridge_loops,
            start,
            end,
            curve: curve.clone(),
            label,
            fillet_index,
        },
    )? {
        return Ok(MaterializedGeneratorTrimSide {
            half_edge: claim.half_edge,
            support_face: claim.support_face,
            support_side_owned: true,
            support_source: None,
            authored_boundary_spans: claim.authored_boundary_spans,
        });
    }
    match materialize_support_face_mef_trim(
        brep,
        support_face,
        start,
        end,
        curve.clone(),
        label,
        fillet_index,
    ) {
        Ok((half_edge, actual_support_face, support_source)) => {
            return Ok(MaterializedGeneratorTrimSide {
                half_edge,
                support_face: actual_support_face,
                support_side_owned: true,
                support_source: Some(support_source),
                authored_boundary_spans: Vec::new(),
            });
        }
        Err(err) => {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} through support-face mef on {support_face:?}: \
                     {err}"
                );
            }
        }
    }
    if let Some(patch) = author_interior_support_trim_patch_for_cylinder_trim(
        brep,
        InteriorSupportTrimPatchRequest {
            support_face,
            start,
            end,
            curve: curve.clone(),
            label,
            fillet_index,
        },
    )? {
        return Ok(patch);
    }
    let (edge, half_edge, _twin) = insert_unwired_edge_pair(brep, start, end);
    brep.edges[edge].curve = Some(curve);
    Ok(MaterializedGeneratorTrimSide {
        half_edge,
        support_face,
        support_side_owned: false,
        support_source: None,
        authored_boundary_spans: Vec::new(),
    })
}

fn materialize_support_face_mef_trim(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<(HalfEdgeHandle, FaceHandle, SupportOwnedTrimSource), ApexBuildError> {
    let rollback = brep.clone();
    let result = (|| {
        let (actual_support_face, support_start, support_end) =
            support_residual_face_split_half_edges(brep, support_face, start, end).ok_or_else(
                || ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                         {label} {start:?}->{end:?} by support-face mef: no residual support \
                         face has both endpoints on one outer loop. {}",
                        support_family_endpoint_diagnostic(brep, support_face, start, end)
                    ),
                },
            )?;
        let result =
            brep.mef(support_start, support_end)
                .map_err(|err| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} by support-face mef on {actual_support_face:?}: \
                     {err:?}"
                    ),
                })?;
        brep.edges[result.edge].curve = Some(curve);
        let support_source = SupportOwnedTrimSource {
            loop_: brep.half_edges[result.half_edge].loop_,
            prev: brep.half_edges[result.half_edge].prev,
            next: brep.half_edges[result.half_edge].next,
        };
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
                 {start:?}->{end:?} by support-face mef on {actual_support_face:?}; cylinder \
                 side {:?}, support side {:?}",
                result.half_edge, result.half_edge_twin
            );
        }
        Ok((result.half_edge, actual_support_face, support_source))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result
}

fn support_family_endpoint_diagnostic(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> String {
    let mut candidate_faces = Vec::new();
    candidate_faces.push(support_face);
    candidate_faces.extend(coplanar_residual_support_faces(brep, support_face));
    let reports = candidate_faces
        .into_iter()
        .map(|face| {
            let start_on_loop = vertex_is_on_face_outer_loop(brep, face, start);
            let end_on_loop = vertex_is_on_face_outer_loop(brep, face, end);
            let start_edge = brep.vertices.get(start).and_then(|vertex| {
                containing_outer_loop_edge_by_curve(brep, face, vertex.position)
            });
            let end_edge = brep.vertices.get(end).and_then(|vertex| {
                containing_outer_loop_edge_by_curve(brep, face, vertex.position)
            });
            format!(
                "{face:?}: start_on_loop={start_on_loop} end_on_loop={end_on_loop} \
                 start_edge={start_edge:?} end_edge={end_edge:?}"
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!("support family endpoint status [{reports}]")
}

fn support_residual_face_split_half_edges(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let mut candidate_faces = Vec::new();
    candidate_faces.push(support_face);
    candidate_faces.extend(coplanar_residual_support_faces(brep, support_face));
    candidate_faces.into_iter().find_map(|candidate_face| {
        support_face_split_half_edges(brep, candidate_face, start, end)
            .map(|(support_start, support_end)| (candidate_face, support_start, support_end))
    })
}

fn support_face_split_half_edges(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let face = brep.faces.get(support_face)?;
    let loop_handle = face.outer_loop;
    let loop_data = brep.loops.get(loop_handle)?;
    let mut start_edges = Vec::new();
    let mut end_edges = Vec::new();
    let mut cursor = loop_data.half_edge;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        if half_edge.origin == start {
            start_edges.push(cursor);
        }
        if half_edge.origin == end {
            end_edges.push(cursor);
        }
        cursor = half_edge.next;
        if cursor == loop_data.half_edge {
            break;
        }
    }
    nearest_directed_loop_pair(brep, loop_handle, &start_edges, &end_edges)
}

fn nearest_directed_loop_pair(
    brep: &BRep,
    loop_handle: LoopHandle,
    starts: &[HalfEdgeHandle],
    ends: &[HalfEdgeHandle],
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    starts
        .iter()
        .flat_map(|start| {
            ends.iter()
                .filter(move |end| *end != start)
                .filter_map(move |end| {
                    directed_loop_distance(brep, loop_handle, *start, *end)
                        .map(|distance| (distance, *start, *end))
                })
        })
        .min_by_key(|(distance, start, end)| (*distance, *start, *end))
        .map(|(_, start, end)| (start, end))
}

fn directed_loop_distance(
    brep: &BRep,
    loop_handle: LoopHandle,
    start: HalfEdgeHandle,
    end: HalfEdgeHandle,
) -> Option<usize> {
    let mut cursor = start;
    for distance in 0..=brep.half_edges.len() {
        if cursor == end {
            return (distance > 0).then_some(distance);
        }
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        cursor = half_edge.next;
    }
    None
}
