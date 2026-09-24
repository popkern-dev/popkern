//! The boundary-ring rung of `insert_vertex_blend` (PK-234 saddle
//! class, PK-239 convex generalization). Extracted from
//! `blend_insert.rs` for the structural file budget; see the call
//! site for the rung ordering contract.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::boundary::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};
use crate::vertex_classifier::VertexBlendClassification;

/// Gate + attempt for the boundary-driven ring carve. Returns
/// `(ring_eligible, ring_face, decline_reason)`.
///
/// Eligibility: the PK-234 saddle class (`OppositeSide`), plus the
/// PK-239 convex classes with no stronger serving route —
/// `SameSideConvex` where the setback surgery was attempted and
/// declined (the conic flat-rho class), where the sphere-ownership
/// prediction is measurably false (uniform cylinder strips on
/// mutually OBLIQUE axes — no rolling-ball sphere exists there,
/// PK-232's recorded 4.9e-2 fit decline), or where the boundary was
/// sphere-fit substituted (the substituted curves are the contract
/// the cap interpolates, and the mitred ring is not those curves).
/// Orthogonal sphere-owned corners never attempt setback, keep a
/// true ownership claim, and never substitute, so the rolling-ball
/// serve is unreachable by this gate by construction.
///
/// Convex corners keep the unpushed `n_gon_boundary` for
/// every downstream consumer — the ring attempt gets its OWN
/// boundary with station pushback, because the mitred carve ends
/// the strips short of the closest-to-vertex stations (the default
/// stations sit in trim-overshoot territory past the carved rails;
/// measured on the conic cube corner: junction x = -0.0156 vs rail
/// end x = 0.1). Only a SERVED ring substitutes the pushed boundary
/// into `n_gon_boundary`; a decline leaves dispatch bit-identical
/// to the mitred path. Saddles already computed `n_gon_boundary`
/// with pushback (their bridges are invalid without it) and use it
/// directly.
#[allow(clippy::too_many_arguments)]
pub(super) fn boundary_ring_rung(
    brep: &mut BRep,
    vertex: VertexHandle,
    classification: VertexBlendClassification,
    use_rational_gregory: bool,
    setback_declined: bool,
    sphere_claim_false: bool,
    n_gon_boundary: &mut Option<BlendBoundary>,
    boundary_options: crate::blend::BlendBoundaryOptions,
    local_indices_pre: &[usize],
    fillet_surfaces_pre: &[FilletSurfaceResult],
    n_sides: usize,
    tolerance: &Tolerance,
) -> (bool, Option<FaceHandle>, Option<String>) {
    let ring_eligible = use_rational_gregory
        && (matches!(classification, VertexBlendClassification::OppositeSide)
            || (matches!(classification, VertexBlendClassification::SameSideConvex)
                && (setback_declined
                    || sphere_claim_false
                    || n_gon_boundary
                        .as_ref()
                        .is_some_and(|b| b.sphere_fit_substituted))));
    let convex_ring_boundary: Option<BlendBoundary> =
        if ring_eligible && !matches!(classification, VertexBlendClassification::OppositeSide) {
            crate::blend::compute_blend_boundary_with_options(
                vertex,
                local_indices_pre,
                fillet_surfaces_pre,
                brep,
                tolerance,
                crate::blend::BlendBoundaryOptions {
                    ring_station_pushback: true,
                    ring_sense_from_solid: true,
                    ..boundary_options
                },
            )
            .ok()
            .filter(|b| b.segments.len() == 2 * n_sides)
        } else {
            None
        };
    let ring_attempt_boundary = convex_ring_boundary.as_ref().or(n_gon_boundary.as_ref());
    let (ring_face, decline) = try_boundary_ring_carve_gated(
        brep,
        vertex,
        ring_eligible,
        ring_attempt_boundary,
        fillet_surfaces_pre,
        tolerance,
    );
    if ring_face.is_some() && convex_ring_boundary.is_some() {
        // The ring IS the pushed boundary; downstream interpolation
        // must see the curves the carve committed.
        *n_gon_boundary = convex_ring_boundary;
    }
    (ring_eligible, ring_face, decline)
}

/// Attempt the transactional carve behind the kill-switch. On
/// decline the BRep is untouched and the reason is recorded for
/// [`boundary_ring_diagnostic`].
fn try_boundary_ring_carve_gated(
    brep: &mut BRep,
    vertex: VertexHandle,
    ring_eligible: bool,
    n_gon_boundary: Option<&BlendBoundary>,
    fillet_surfaces_pre: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> (Option<FaceHandle>, Option<String>) {
    if !ring_eligible {
        return (None, None);
    }
    // PK-234 stage 4: the ring carves machine-exact (4.7e-12), the
    // support-seam relocation no-ops on it, and with the blend loop
    // normalized to start on an arc the sub-patches bond exactly —
    // the saddle fixture serves validate_all CLEAN (zero check-9
    // edges) where the mitred fallback carries 14 failing edges at
    // worst 1.016e-1.
    let Some(boundary) = n_gon_boundary else {
        return (None, None);
    };
    match crate::construct::setback_vertex_topology::try_boundary_ring_carve(
        brep,
        vertex,
        boundary,
        fillet_surfaces_pre,
        tolerance,
    ) {
        Ok(face) => (Some(face), None),
        Err(reason) => {
            if popkern_core::trace::trace_enabled("corner-carve") {
                eprintln!("  [boundary-ring] declined at {vertex:?}: {reason}");
            }
            (None, Some(reason))
        }
    }
}

/// Structured fact for the diagnostics pipeline — one per saddle
/// vertex, covering applied and declined outcomes (a decline means
/// the mitred fallback ring is serving with its known deviation set,
/// PK-236).
pub(super) fn boundary_ring_diagnostic(
    ring_eligible: bool,
    n_sides: usize,
    ring_face: Option<FaceHandle>,
    decline: Option<&String>,
) -> Option<ApexConstructionDiagnostic> {
    if !ring_eligible {
        return None;
    }
    let stage = match (ring_face, decline) {
        (Some(_), _) => "boundary_ring.applied".to_string(),
        (None, Some(reason)) => format!("boundary_ring.declined: {reason}"),
        (None, None) => "boundary_ring.not_attempted".to_string(),
    };
    Some(ApexConstructionDiagnostic::ExactBoundaryCap(
        ExactBoundaryCapDiagnostic {
            stage,
            code: "boundary_ring.fact".to_string(),
            detail_code: None,
            sector_name: format!("boundary_ring_n{n_sides}"),
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
