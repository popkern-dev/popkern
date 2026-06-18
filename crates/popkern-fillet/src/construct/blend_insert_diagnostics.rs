//! Structured diagnostics for vertex-blend dispatch declines.

use popkern_core::topology::{BRep, FaceHandle};

use crate::blend::BlendBoundary;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

pub(super) fn rational_gregory_dispatch_declined_diagnostic(
    n_sides: usize,
    reason_code: &str,
    detail_code: Option<String>,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: format!("rational_gregory.dispatch.{reason_code}"),
        code: "rational_gregory.dispatch.declined".to_string(),
        detail_code,
        sector_name: format!("n{n_sides}_corner_cap"),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: None,
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
        shared_fact_ids: Vec::new(),
    })
}

pub(super) fn pns_dispatch_declined_diagnostic(
    n_sides: usize,
    reason_code: &str,
    detail_code: Option<String>,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: format!("pns.dispatch.{reason_code}"),
        code: "pns.dispatch.declined".to_string(),
        detail_code,
        sector_name: format!("n{n_sides}_corner_cap"),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: None,
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
        shared_fact_ids: Vec::new(),
    })
}

pub(super) fn rational_gregory_dispatch_precondition_reason(
    brep: &BRep,
    blend_face: FaceHandle,
    n_sides: usize,
    boundary: &BlendBoundary,
) -> String {
    if boundary.segments.len() != 2 * n_sides {
        return format!(
            "blend boundary has {} segments, expected {}",
            boundary.segments.len(),
            2 * n_sides
        );
    }
    if !crate::construct::poke_hex::face_has_2n_outer_half_edges_with(
        brep,
        blend_face,
        n_sides,
        boundary.sphere_fit_substituted,
    ) {
        return format!(
            "blend face outer loop is not dispatch-ready for N-sector Gregory/polar construction \
             (outer_half_edges={:?}, expected {} or {})",
            outer_loop_half_edge_count(brep, blend_face),
            n_sides,
            2 * n_sides
        );
    }
    "unknown rational-Gregory dispatch precondition failure".to_string()
}

fn outer_loop_half_edge_count(brep: &BRep, face: FaceHandle) -> Option<usize> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut cursor = start;
    for count in 1..=64 {
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Some(count);
        }
    }
    None
}
