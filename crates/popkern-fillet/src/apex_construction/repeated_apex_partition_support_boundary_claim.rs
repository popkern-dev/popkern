//! Existing support-boundary ownership claims for repeated-apex trim materialization.

use super::*;

pub(super) fn claim_existing_support_boundary_for_cylinder_trim(
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
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] claimed {label} \
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
    if !crate::direct_apex_trace::enabled() {
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
        "[direct-apex-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] found no support boundary \
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
