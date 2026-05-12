//! Claim support-continuation bridge boundary spans for cylinder trims.

use super::*;

pub(super) fn claim_side_support_bridge_boundary_for_cylinder_trim(
    brep: &mut BRep,
    bridge_loops: &[LoopHandle],
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<Option<(HalfEdgeHandle, FaceHandle)>, ApexBuildError> {
    let Some((bridge_face, support_half_edge)) = bridge_loops.iter().find_map(|loop_handle| {
        let loop_data = brep.loops.get(*loop_handle)?;
        collect_loop_half_edges(brep, *loop_handle)
            .ok()?
            .into_iter()
            .find(|half_edge| {
                let Some(data) = brep.half_edges.get(*half_edge) else {
                    return false;
                };
                data.origin == end
                    && support_split_half_edge_end_vertex(brep, *half_edge) == Some(start)
                    && brep.boundary_edges.contains(&data.edge)
                    && !brep.half_edges.contains_key(data.twin)
            })
            .map(|half_edge| (loop_data.face, half_edge))
    }) else {
        trace_side_support_bridge_claim_miss(brep, bridge_loops, start, end, label, fillet_index);
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
                     stale side-support bridge boundary side {support_half_edge:?} for {label}"
                ),
            })?;
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
        validate_side_support_bridge_claim(
            brep,
            bridge_face,
            support_half_edge,
            half_edge,
            label,
            fillet_index,
        )?;
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] claimed {label} \
                 {start:?}->{end:?} from side-support bridge loop on {bridge_face:?}; cylinder \
                 side {half_edge:?}, support side {support_half_edge:?}"
            );
        }
        Ok((half_edge, bridge_face))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result.map(Some)
}

fn trace_side_support_bridge_claim_miss(
    brep: &BRep,
    bridge_loops: &[LoopHandle],
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) {
    if !crate::a2_trace::enabled() {
        return;
    }
    let reports = bridge_loops
        .iter()
        .map(|loop_handle| {
            let face = brep.loops.get(*loop_handle).map(|loop_data| loop_data.face);
            let spans = collect_loop_half_edges(brep, *loop_handle)
                .unwrap_or_default()
                .into_iter()
                .filter_map(|half_edge| {
                    let data = brep.half_edges.get(half_edge)?;
                    let span_end = support_split_half_edge_end_vertex(brep, half_edge)?;
                    let boundary = brep.boundary_edges.contains(&data.edge);
                    let twin_live = brep.half_edges.contains_key(data.twin);
                    Some(format!(
                        "{half_edge:?}:{:?}->{span_end:?} boundary={boundary} twin_live={twin_live}",
                        data.origin
                    ))
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!("{loop_handle:?}/{face:?} [{spans}]")
        })
        .collect::<Vec<_>>()
        .join("; ");
    eprintln!(
        "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] found no side-support bridge \
         boundary side to claim for {label} {start:?}->{end:?}; bridge candidates {reports}"
    );
}

fn validate_side_support_bridge_claim(
    brep: &BRep,
    bridge_face: FaceHandle,
    support_half_edge: HalfEdgeHandle,
    cylinder_half_edge: HalfEdgeHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<(), ApexBuildError> {
    let support_owner = brep
        .half_edges
        .get(support_half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face);
    if support_owner != Some(bridge_face) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} from \
                 side-support bridge half-edge {support_half_edge:?}, but its owner is \
                 {support_owner:?} instead of {bridge_face:?}"
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
        || brep.boundary_edges.contains(&cylinder_data.edge)
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] claimed {label} with \
                 invalid bridge ownership on {cylinder_half_edge:?}/{support_half_edge:?}"
            ),
        });
    }
    Ok(())
}
