//! Owner-side span derivation for cylinder-mediated support continuation.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::{
    CylinderMediatedOwnerRepairSpan, CylinderMediatedSupportContinuationSegmentAuthoringPlan,
};

pub(super) fn owner_repair_spans(
    brep: &BRep,
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
) -> Result<Vec<CylinderMediatedOwnerRepairSpan>, CoreError> {
    let mut spans = Vec::new();
    let mut current: Option<OwnerRepairSpanBuilder> = None;
    for &outside_half_edge in &segment.half_edges {
        if !segment_requires_owner_repair(segment, outside_half_edge) {
            push_owner_repair_span(brep, &mut spans, current.take())?;
            continue;
        }
        let twin = segment_twin_half_edge(brep, outside_half_edge)?;
        let owner_loop = segment_owner_loop(brep, twin)?;
        let owner_face = brep.loops[owner_loop].face;
        match current.as_mut() {
            Some(builder) if builder.owner_loop == owner_loop => {
                builder.outside_half_edges.push(outside_half_edge);
                builder.owner_half_edges.push(twin);
            }
            Some(_) => {
                push_owner_repair_span(brep, &mut spans, current.take())?;
                current = Some(OwnerRepairSpanBuilder::new(
                    owner_face,
                    owner_loop,
                    outside_half_edge,
                    twin,
                ));
            }
            None => {
                current = Some(OwnerRepairSpanBuilder::new(
                    owner_face,
                    owner_loop,
                    outside_half_edge,
                    twin,
                ));
            }
        }
    }
    push_owner_repair_span(brep, &mut spans, current)?;
    Ok(spans)
}

fn segment_requires_owner_repair(
    segment: &CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    half_edge: HalfEdgeHandle,
) -> bool {
    segment.cylinder_owned_half_edges.contains(&half_edge)
        || segment.outside_owned_half_edges.contains(&half_edge)
        || segment.support_owned_half_edges.contains(&half_edge)
        || segment.context_owned_half_edges.contains(&half_edge)
}

fn segment_twin_half_edge(
    brep: &BRep,
    outside_half_edge: HalfEdgeHandle,
) -> Result<HalfEdgeHandle, CoreError> {
    let twin = brep
        .half_edges
        .get(outside_half_edge)
        .ok_or_else(|| {
            CoreError::InvalidHandle(format!(
                "owner repair span references stale outside half-edge {outside_half_edge:?}"
            ))
        })?
        .twin;
    if !brep.half_edges.contains_key(twin) {
        return Err(CoreError::PreconditionFailed(format!(
            "owner repair span outside half-edge {outside_half_edge:?} has no live twin"
        )));
    }
    Ok(twin)
}

fn segment_owner_loop(brep: &BRep, twin: HalfEdgeHandle) -> Result<LoopHandle, CoreError> {
    let owner_loop = brep.half_edges[twin].loop_;
    if !brep.loops.contains_key(owner_loop) {
        return Err(CoreError::PreconditionFailed(format!(
            "owner repair span twin half-edge {twin:?} has no live owner loop"
        )));
    }
    Ok(owner_loop)
}

fn push_owner_repair_span(
    brep: &BRep,
    spans: &mut Vec<CylinderMediatedOwnerRepairSpan>,
    builder: Option<OwnerRepairSpanBuilder>,
) -> Result<(), CoreError> {
    let Some(mut builder) = builder else {
        return Ok(());
    };
    builder.owner_half_edges.reverse();
    for pair in builder.owner_half_edges.windows(2) {
        if brep.half_edges[pair[0]].next != pair[1] {
            return Err(CoreError::PreconditionFailed(format!(
                "owner repair span on {:?}/{:?} is not contiguous at {:?}->{:?}",
                builder.owner_face, builder.owner_loop, pair[0], pair[1]
            )));
        }
    }
    spans.push(CylinderMediatedOwnerRepairSpan {
        owner_face: builder.owner_face,
        owner_loop: builder.owner_loop,
        outside_half_edges: builder.outside_half_edges,
        owner_half_edges: builder.owner_half_edges,
    });
    Ok(())
}

struct OwnerRepairSpanBuilder {
    owner_face: FaceHandle,
    owner_loop: LoopHandle,
    outside_half_edges: Vec<HalfEdgeHandle>,
    owner_half_edges: Vec<HalfEdgeHandle>,
}

impl OwnerRepairSpanBuilder {
    fn new(
        owner_face: FaceHandle,
        owner_loop: LoopHandle,
        outside_half_edge: HalfEdgeHandle,
        owner_half_edge: HalfEdgeHandle,
    ) -> Self {
        Self {
            owner_face,
            owner_loop,
            outside_half_edges: vec![outside_half_edge],
            owner_half_edges: vec![owner_half_edge],
        }
    }
}
