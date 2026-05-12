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
) -> Result<
    (
        HalfEdgeHandle,
        FaceHandle,
        bool,
        Option<SupportOwnedTrimSource>,
    ),
    ApexBuildError,
> {
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
        return Ok((half_edge, actual_support_face, true, None));
    }
    if let Some((half_edge, actual_support_face)) =
        claim_side_support_bridge_boundary_for_cylinder_trim(
            brep,
            side_support_bridge_loops,
            start,
            end,
            curve.clone(),
            label,
            fillet_index,
        )?
    {
        return Ok((half_edge, actual_support_face, true, None));
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
            return Ok((half_edge, actual_support_face, true, Some(support_source)));
        }
        Err(err) => {
            if crate::a2_trace::enabled() {
                eprintln!(
                    "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} through support-face mef on {support_face:?}: \
                     {err}"
                );
            }
        }
    }
    let (edge, half_edge, _twin) = insert_unwired_edge_pair(brep, start, end);
    brep.edges[edge].curve = Some(curve);
    Ok((half_edge, support_face, false, None))
}

fn claim_existing_support_boundary_for_cylinder_trim(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<Option<(HalfEdgeHandle, FaceHandle)>, ApexBuildError> {
    let candidate_faces = support_family_faces(brep, support_face);
    let Some((actual_support_face, support_half_edge)) =
        candidate_faces.into_iter().find_map(|candidate_face| {
            support_boundary_half_edge_for_cylinder_side(brep, candidate_face, start, end)
                .map(|half_edge| (candidate_face, half_edge))
        })
    else {
        trace_support_boundary_claim_miss(brep, support_face, start, end, label, fillet_index);
        return Ok(None);
    };
    let rollback = brep.clone();
    let result = (|| {
        let edge = brep
            .half_edges
            .get(support_half_edge)
            .map(|half_edge| half_edge.edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot claim \
                     stale support boundary side {support_half_edge:?} for {label}"
                ),
            })?;
        if !brep.boundary_edges.contains(&edge) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot claim \
                     non-boundary support side {support_half_edge:?} for {label}"
                ),
            });
        }
        if brep
            .half_edges
            .get(support_half_edge)
            .is_some_and(|half_edge| brep.half_edges.contains_key(half_edge.twin))
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot claim \
                     already-owned support side {support_half_edge:?} for {label}"
                ),
            });
        }
        if brep.edges[edge].curve.is_none() {
            brep.edges[edge].curve = Some(curve);
        }
        let half_edge = brep.half_edges.insert(HalfEdge {
            origin: start,
            twin: support_half_edge,
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge,
            loop_: LoopHandle::default(),
            pcurve: None,
        });
        brep.half_edges[support_half_edge].twin = half_edge;
        brep.boundary_edges.remove(&edge);
        if support_split_half_edge_end_vertex(brep, half_edge) != Some(end) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] claimed support \
                     boundary side {support_half_edge:?} for {label} with wrong cylinder \
                     orientation: expected {start:?}->{end:?}, got {half_edge:?}"
                ),
            });
        }
        validate_claimed_support_boundary_owner(
            brep,
            actual_support_face,
            support_half_edge,
            half_edge,
            label,
            fillet_index,
        )?;
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] claimed {label} \
                 {start:?}->{end:?} from support boundary on {actual_support_face:?}; cylinder \
                 side {half_edge:?}, support side {support_half_edge:?}"
            );
        }
        Ok((half_edge, actual_support_face))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result.map(Some)
}

fn trace_support_boundary_claim_miss(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) {
    if !crate::a2_trace::enabled() {
        return;
    }
    let reports = support_family_faces(brep, support_face)
        .into_iter()
        .map(|face| {
            let spans = brep
                .faces
                .get(face)
                .and_then(|face_data| collect_loop_half_edges(brep, face_data.outer_loop).ok())
                .unwrap_or_default()
                .into_iter()
                .filter_map(|half_edge| {
                    let data = brep.half_edges.get(half_edge)?;
                    let span_end = support_split_half_edge_end_vertex(brep, half_edge)?;
                    let touches_request = data.origin == start
                        || data.origin == end
                        || span_end == start
                        || span_end == end;
                    if !touches_request {
                        return None;
                    }
                    let boundary = brep.boundary_edges.contains(&data.edge);
                    let twin_live = brep.half_edges.contains_key(data.twin);
                    Some(format!(
                        "{half_edge:?}:{:?}->{span_end:?} boundary={boundary} twin_live={twin_live}",
                        data.origin
                    ))
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!("{face:?} [{spans}]")
        })
        .collect::<Vec<_>>()
        .join("; ");
    eprintln!(
        "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] found no support boundary \
         side to claim for {label} {start:?}->{end:?} on support family {support_face:?}; \
         candidate spans {reports}"
    );
}

fn support_boundary_half_edge_for_cylinder_side(
    brep: &BRep,
    support_face: FaceHandle,
    cylinder_start: VertexHandle,
    cylinder_end: VertexHandle,
) -> Option<HalfEdgeHandle> {
    let face = brep.faces.get(support_face)?;
    collect_loop_half_edges(brep, face.outer_loop)
        .ok()?
        .into_iter()
        .find(|half_edge| {
            let Some(data) = brep.half_edges.get(*half_edge) else {
                return false;
            };
            data.origin == cylinder_end
                && support_split_half_edge_end_vertex(brep, *half_edge) == Some(cylinder_start)
                && brep.boundary_edges.contains(&data.edge)
                && !brep.half_edges.contains_key(data.twin)
        })
}

fn validate_claimed_support_boundary_owner(
    brep: &BRep,
    support_face: FaceHandle,
    support_half_edge: HalfEdgeHandle,
    cylinder_half_edge: HalfEdgeHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<(), ApexBuildError> {
    let support_loop = brep
        .half_edges
        .get(support_half_edge)
        .map(|half_edge| half_edge.loop_)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed stale support \
                 half-edge {support_half_edge:?} for {label}"
            ),
        })?;
    let support_owner = brep.loops.get(support_loop).map(|loop_data| loop_data.face);
    if support_owner != Some(support_face) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} from \
                 support half-edge {support_half_edge:?}, but its loop owner is {support_owner:?} \
                 instead of {support_face:?}"
            ),
        });
    }
    let Some(cylinder_data) = brep.half_edges.get(cylinder_half_edge) else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed stale cylinder \
                 half-edge {cylinder_half_edge:?} for {label}"
            ),
        });
    };
    if cylinder_data.twin != support_half_edge
        || brep.half_edges[support_half_edge].twin != cylinder_half_edge
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} with \
                 inconsistent twin ownership between {cylinder_half_edge:?} and \
                 {support_half_edge:?}"
            ),
        });
    }
    if brep.boundary_edges.contains(&cylinder_data.edge) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} but \
                 edge {:?} remains marked as boundary",
                cylinder_data.edge
            ),
        });
    }
    Ok(())
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
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
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
