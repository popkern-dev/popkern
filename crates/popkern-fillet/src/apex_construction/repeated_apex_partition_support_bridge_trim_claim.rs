//! Claim support-continuation bridge boundary spans for cylinder trims.

#[path = "repeated_apex_partition_support_bridge_path_patch.rs"]
mod side_support_bridge_path_patch;

use super::*;
use side_support_bridge_path_patch::{
    SideSupportBridgePathPatch, SideSupportBridgePathPatchRequest,
    author_side_support_bridge_path_patch_for_cylinder_trim, trace_side_support_bridge_claim_miss,
};

pub(super) fn claim_side_support_bridge_boundary_for_cylinder_trim(
    brep: &mut BRep,
    request: SideSupportBridgeTrimClaimRequest<'_>,
) -> Result<Option<SideSupportBridgeClaim>, ApexBuildError> {
    let SideSupportBridgeTrimClaimRequest {
        support_face,
        bridge_loops,
        start,
        end,
        curve,
        label,
        fillet_index,
    } = request;
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
        if let Some(result) = author_side_support_bridge_path_patch_for_cylinder_trim(
            brep,
            SideSupportBridgePathPatchRequest {
                support_face,
                bridge_loops,
                start,
                end,
                curve,
                label,
                fillet_index,
            },
        )? {
            return Ok(Some(SideSupportBridgeClaim::from_path_patch(result)));
        }
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
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] claimed {label} \
                 {start:?}->{end:?} from side-support bridge loop on {bridge_face:?}; cylinder \
                 side {half_edge:?}, support side {support_half_edge:?}"
            );
        }
        Ok(SideSupportBridgeClaim {
            half_edge,
            support_face: bridge_face,
            authored_boundary_spans: Vec::new(),
        })
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result.map(Some)
}

pub(super) struct SideSupportBridgeClaim {
    pub(super) half_edge: HalfEdgeHandle,
    pub(super) support_face: FaceHandle,
    pub(super) authored_boundary_spans: Vec<HalfEdgeHandle>,
}

impl SideSupportBridgeClaim {
    fn from_path_patch(patch: SideSupportBridgePathPatch) -> Self {
        Self {
            half_edge: patch.cylinder_half_edge,
            support_face: patch.support_patch_face,
            authored_boundary_spans: patch.boundary_spans,
        }
    }
}

pub(super) struct SideSupportBridgeTrimClaimRequest<'a> {
    pub(super) support_face: FaceHandle,
    pub(super) bridge_loops: &'a [LoopHandle],
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) curve: CurveOnEdge,
    pub(super) label: &'static str,
    pub(super) fillet_index: usize,
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
