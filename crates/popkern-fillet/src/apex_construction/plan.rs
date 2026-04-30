//! Phase 6c-A2 apex planner — orchestration.
//!
//! `plan_apex_resolution` is the **pure** entry point that produces
//! an [`ApexPlan`] from a BRep + apex vertex + fillet specs.
//! It does NOT mutate the BRep — callers can preview, inspect
//! telemetry, or cancel.
//!
//! Algorithm (per ADR §6.2):
//! 1. Validate inputs ([`plan_validate`](super::plan_validate)).
//! 2. Compute boundary segments
//!    ([`plan_terminations::compute_apex_boundary`](super::plan_terminations::compute_apex_boundary)).
//! 3. Try the [`ConstructionRung`] ladder Direct → Setback → Open
//!    ([`plan_with_rung_ladder`]):
//!    * Direct: full per-fillet sphere-fit + Sphere/Gregory cap-kind
//!      decision + target topology assembly.
//!    * Setback ([`plan_rung_fallback::plan_setback_rung`](super::plan_rung_fallback::plan_setback_rung)):
//!      lenient per-fillet sphere-fit, ≥3 successful fillets required.
//!    * Open ([`plan_rung_fallback::plan_open_rung`](super::plan_rung_fallback::plan_open_rung)):
//!      empty apex topology, fillets terminate at legacy rolling-ball
//!      footprints.
//! 4. Telemetry records the rung used + any fallback reason.

use popkern_core::topology::{BRep, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::FilletSpec;

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan_rung_fallback::{plan_open_rung, plan_setback_rung};
use super::plan_terminations::{compute_apex_boundary, compute_per_fillet_terminations};
use super::plan_topology::{build_target_topology, decide_cap_kind};
use super::plan_validate::validate_inputs;
use super::target_topology::TargetTopology;
use super::telemetry::{ConstructionRung, PlanTelemetry};

/// A plan for resolving an apex vertex's fillet topology.
///
/// Produced by [`plan_apex_resolution`]; consumed by
/// `construct_apex_topology` (when A2.2 lands).  Introspectable for
/// UI preview.
///
/// See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md` §3.1.
#[derive(Debug, Clone)]
pub struct ApexPlan {
    /// The apex vertex this plan resolves.
    pub apex_vertex: VertexHandle,
    /// Cap kind decision (Sphere or Gregory) with computed geometry.
    pub cap_kind: CapKind,
    /// Final topology spec — vertices, edges, faces with their geometry.
    pub target_topology: TargetTopology,
    /// Decision telemetry — explains why this plan was chosen, which
    /// rung was used, and exposes diagnostic numerics.  UI-friendly.
    pub telemetry: PlanTelemetry,
}

impl ApexPlan {
    /// Number of fillets meeting at the apex.
    pub fn n_fillets(&self) -> usize {
        self.cap_kind.n_fillets()
    }
}

/// Compute an [`ApexPlan`] for resolving an apex vertex's fillet
/// topology.  Pure — does not mutate the BRep.
///
/// # Inputs
///
/// * `brep` — the BRep being filleted (read-only).
/// * `apex` — the vertex where multiple fillets meet.
/// * `fillet_specs` — the fillet specs (one per incident fillet edge),
///   in the same order as the cap-loop walk would visit them.
/// * `fillet_surfaces` — fillet surface data (cylinder/torus surfaces
///   + trim curves), one per fillet, computed earlier in the
///     pipeline by `geometry_compute`.
/// * `tolerance` — geometric tolerance for cap-kind dispatch and
///   coincidence checks.
///
/// # Output
///
/// `Ok(ApexPlan)` on success; `Err(ApexPlanError)` describing the
/// failure mode for user-actionable feedback.
///
/// # Errors
///
/// See [`ApexPlanError`] for the full failure taxonomy.
pub fn plan_apex_resolution(
    brep: &BRep,
    apex: VertexHandle,
    fillet_specs: &[FilletSpec],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<ApexPlan, ApexPlanError> {
    // Step 1: validate inputs (A2.1.b).  Validation errors don't
    // benefit from rung fallback — invalid inputs are user errors.
    validate_inputs(brep, apex, fillet_specs)?;

    // Step 2: compute boundary segments (A2.1.c).  Same — boundary
    // failure is structural, not a sphere-fit failure.
    let boundary = compute_apex_boundary(brep, apex, fillet_surfaces, tolerance)?;

    // Steps 3-5 with rung fallback: try Direct → Setback → Open.
    plan_with_rung_ladder(
        brep,
        apex,
        &boundary,
        fillet_specs,
        fillet_surfaces,
        tolerance,
    )
}

/// Try [`ConstructionRung::Direct`] first, then Setback, then Open.
/// Records the successful rung (and any fallback reason) in
/// [`PlanTelemetry`].
fn plan_with_rung_ladder(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundaryRef,
    fillet_specs: &[FilletSpec],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<ApexPlan, ApexPlanError> {
    // --- Rung 1: Direct ---
    match compute_per_fillet_terminations(brep, apex, boundary, fillet_surfaces, tolerance) {
        Ok(terminations) => {
            let mut telemetry = PlanTelemetry::new();
            let cap_kind = decide_cap_kind(&terminations, fillet_specs, tolerance, &mut telemetry)?;
            let target_topology = build_target_topology(&cap_kind, &terminations);
            telemetry.rung = ConstructionRung::Direct;
            Ok(ApexPlan {
                apex_vertex: apex,
                cap_kind,
                target_topology,
                telemetry,
            })
        }
        Err(direct_err) => {
            // --- Rung 2: Setback ---
            //
            // Setback strategy: tolerate individual per-fillet
            // sphere-fit failures.  Successful fillets contribute
            // apex_vertices/edges; failed fillets are skipped.
            // Requires ≥3 successful fillets.
            if let Ok(plan) = plan_setback_rung(
                brep,
                apex,
                boundary,
                fillet_specs,
                fillet_surfaces,
                tolerance,
                &direct_err,
            ) {
                return Ok(plan);
            }

            // --- Rung 3: Open ---
            //
            // No cap construction.  Plan has empty apex topology;
            // builder leaves apex as raw rolling-ball terminations.
            Ok(plan_open_rung(apex, fillet_surfaces.len(), &direct_err))
        }
    }
}

// Type alias to keep the rung-ladder signature concise without
// re-importing BlendBoundary at the top of plan.rs (it's only
// referenced inside this private fn).
type BlendBoundaryRef = crate::blend::boundary::BlendBoundary;
