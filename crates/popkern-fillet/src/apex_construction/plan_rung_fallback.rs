//! Phase 6c-DirectApex.1.f — Setback and Open rung fallback.
//!
//! When the Direct rung fails (sphere-fit error on one or more
//! fillets), the planner attempts the Setback rung: lenient
//! per-fillet termination computation that omits failures.
//! If <3 fillets succeed, falls through to Open rung (no cap).
//!
//! See ADR §4.3 for rung definitions.

use popkern_core::topology::{BRep, VertexHandle};

use crate::blend::boundary::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::FilletSpec;

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan::ApexPlan;
use super::plan_terminations::compute_per_fillet_terminations_lenient;
use super::plan_topology::{build_target_topology, decide_cap_kind};
use super::target_topology::{TargetCapFace, TargetTopology};
use super::telemetry::{ConstructionRung, PlanTelemetry};

/// Attempt the Setback rung: tolerate per-fillet sphere-fit failures
/// and build a cap from whatever fillets succeeded.  Requires at
/// least 3 successful fillets to form a cap; otherwise returns
/// `direct_err` as a sentinel (caller falls through to Open rung).
pub(crate) fn plan_setback_rung(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundary,
    fillet_specs: &[FilletSpec],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    direct_err: &ApexPlanError,
) -> Result<ApexPlan, ApexPlanError> {
    let terminations =
        compute_per_fillet_terminations_lenient(brep, apex, boundary, fillet_surfaces, tolerance)?;
    if terminations.len() < 3 {
        // Not enough successful fillets to build a cap; signal to
        // caller to fall through to Open rung.
        return Err(ApexPlanError::SphereCylinderNoIntersection {
            fillet_idx: 0, // sentinel; real failure recorded in direct_err
        });
    }
    let mut telemetry = PlanTelemetry::new();
    let cap_kind = decide_cap_kind(&terminations, fillet_specs, tolerance, &mut telemetry)?;
    let target_topology = build_target_topology(&cap_kind, &terminations);
    telemetry.rung = ConstructionRung::Setback {
        // No specific setback distance — Setback here means
        // "omit failed fillets from the cap".  A future enhancement
        // can compute a real setback distance per fillet.
        setback_distance: 0.0,
    };
    telemetry.fallback_reason = Some(format!(
        "Direct rung failed ({direct_err}); Setback rung succeeded by omitting failed fillets"
    ));
    Ok(ApexPlan {
        apex_vertex: apex,
        cap_kind,
        target_topology,
        telemetry,
    })
}

/// Build an Open-rung [`ApexPlan`]: empty topology, no cap face.
/// The fallback when neither Direct nor Setback can produce a
/// valid cap.  Builder's view of this plan: do nothing — leave
/// the apex region as the legacy rolling-ball termination.
pub(crate) fn plan_open_rung(
    apex: VertexHandle,
    n_fil: usize,
    direct_err: &ApexPlanError,
) -> ApexPlan {
    let cap_kind = CapKind::Open {
        n_fillets_attempted: n_fil,
    };
    let target_topology = TargetTopology {
        apex_vertices: Vec::new(),
        apex_edges: Vec::new(),
        cap_face: TargetCapFace::Open,
        face_boundary_updates: Vec::new(),
    };
    let mut telemetry = PlanTelemetry::new();
    telemetry.rung = ConstructionRung::Open;
    telemetry.fallback_reason = Some(format!(
        "Direct rung failed ({direct_err}); Setback rung also failed; \
         Open rung produces no apex cap.  {n_fil} fillets terminate at \
         their legacy rolling-ball footprints."
    ));
    telemetry.cap_kind_reason = "Open — no cap construction".to_string();
    ApexPlan {
        apex_vertex: apex,
        cap_kind,
        target_topology,
        telemetry,
    }
}
