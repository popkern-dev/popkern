//! Structured diagnostics for boundary-owned pcurve authoring declines.

use popkern_core::topology::{BRep, HalfEdgeHandle};

use crate::blend::boundary::BoundarySegment;

pub(super) fn boundary_owned_pcurve_decline_diagnostic(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    segment: &BoundarySegment,
    segment_index: usize,
    max_disagreement: f64,
    endpoint_gap: f64,
    tolerance: f64,
) -> crate::types::ApexConstructionDiagnostic {
    let edge = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.edge);
    let support_half_edge = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.twin);
    let support_face = support_half_edge
        .and_then(|support_half_edge| brep.half_edges.get(support_half_edge))
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face)
        .or_else(|| match segment {
            BoundarySegment::FaceEdge { face, .. } => Some(*face),
            BoundarySegment::FilletEnd { .. } => None,
        });
    let code = if endpoint_gap > tolerance {
        "boundary_owned_pcurve_endpoint_mismatch"
    } else {
        "boundary_owned_pcurve_not_on_canonical_boundary"
    };
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "native_gregory_boundary_owned_pcurve".to_string(),
            code: code.to_string(),
            sector_name: format!("cap_boundary.segment[{segment_index}]"),
            edge,
            support_half_edge,
            support_face,
            residual: Some(max_disagreement.max(endpoint_gap)),
            shared_fact_ids: Vec::new(),
        },
    )
}
