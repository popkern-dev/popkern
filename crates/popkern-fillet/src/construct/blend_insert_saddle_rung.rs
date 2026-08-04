//! The saddle boundary-ring rung of `insert_vertex_blend` (PK-234).
//! Extracted from `blend_insert.rs` for the structural file budget;
//! see the call site for the rung ordering contract.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::boundary::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};
use crate::vertex_classifier::VertexBlendClassification;

/// Attempt the boundary-driven saddle ring carve. Returns the blend
/// face on success; on decline the BRep is untouched and the reason
/// is recorded for [`saddle_ring_diagnostic`].
pub(super) fn try_saddle_ring_rung(
    brep: &mut BRep,
    vertex: VertexHandle,
    classification: VertexBlendClassification,
    use_rational_gregory: bool,
    n_gon_boundary: Option<&BlendBoundary>,
    fillet_surfaces_pre: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> (Option<FaceHandle>, Option<String>) {
    if !matches!(classification, VertexBlendClassification::OppositeSide) || !use_rational_gregory {
        return (None, None);
    }
    // Temporary containment (PK-234 stage 3): the ring itself carves
    // machine-exact (worst on-surface deviation 4.7e-12), but the
    // downstream post-apply repair passes assume the mitred-era
    // topology and misfire on the new ring (measured: rail-edge
    // curves rewritten to span to the cap's central vertex, support
    // planes refit to NURBS — worst check-9 deviation 8.8e-1 vs the
    // mitred fallback's 1.016e-1). Opt-in until those passes are
    // adapted; the mitred fallback keeps today's serve.
    if std::env::var_os("POPKERN_SADDLE_RING").is_none_or(|v| v != "1") {
        return (
            None,
            Some("containment: POPKERN_SADDLE_RING not set (downstream repair passes misfire on the ring substrate)".to_string()),
        );
    }
    let Some(boundary) = n_gon_boundary else {
        return (None, None);
    };
    match crate::construct::setback_vertex_topology::try_saddle_boundary_ring(
        brep,
        vertex,
        boundary,
        fillet_surfaces_pre,
        tolerance,
    ) {
        Ok(face) => (Some(face), None),
        Err(reason) => {
            if std::env::var_os("PK_CARVE_TRACE").is_some() {
                eprintln!("  [saddle-ring] declined at {vertex:?}: {reason}");
            }
            (None, Some(reason))
        }
    }
}

/// Structured fact for the diagnostics pipeline — one per saddle
/// vertex, covering applied and declined outcomes (a decline means
/// the mitred fallback ring is serving with its known deviation set,
/// PK-236).
pub(super) fn saddle_ring_diagnostic(
    classification: VertexBlendClassification,
    use_rational_gregory: bool,
    n_sides: usize,
    ring_face: Option<FaceHandle>,
    decline: Option<&String>,
) -> Option<ApexConstructionDiagnostic> {
    if !matches!(classification, VertexBlendClassification::OppositeSide) || !use_rational_gregory {
        return None;
    }
    let stage = match (ring_face, decline) {
        (Some(_), _) => "saddle_boundary_ring.applied".to_string(),
        (None, Some(reason)) => format!("saddle_boundary_ring.declined: {reason}"),
        (None, None) => "saddle_boundary_ring.not_attempted".to_string(),
    };
    Some(ApexConstructionDiagnostic::ExactBoundaryCap(
        ExactBoundaryCapDiagnostic {
            stage,
            code: "saddle_boundary_ring.fact".to_string(),
            detail_code: None,
            sector_name: format!("saddle_ring_n{n_sides}"),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: None,
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    ))
}
