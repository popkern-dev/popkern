//! Phase 6c-A2 apex planner.
//!
//! `plan_apex_resolution` is the **pure** entry point that produces
//! an [`ApexPlan`] from a BRep + apex vertex + fillet specs.
//! It does NOT mutate the BRep — callers can preview, inspect
//! telemetry, or cancel.
//!
//! Algorithm (per ADR §6.2):
//! 1. Validate inputs (radii positive/finite, ≥3 fillets at apex).
//! 2. Compute boundary segments (via `compute_blend_boundary`).
//! 3. Per fillet: compute T_a[k], T_b[k] via `sphere_fit_termination`.
//! 4. Decide cap kind: Sphere if all radii equal AND bridges collapse,
//!    else Gregory.
//! 5. Build [`TargetTopology`] (apex_vertices, apex_edges, cap_face,
//!    face_boundary_updates).
//! 6. Build [`PlanTelemetry`] (cap_kind_reason, diagnostic_values).
//! 7. Return [`ApexPlan`].
//!
//! Failure paths fall back through the [`ConstructionRung`] ladder:
//! Direct → Setback → Open.  Each fallback emits a
//! [`BlendWarning::ContinuityDowngraded`] (legacy pattern).
//!
//! # Status
//!
//! Stub implementation (A2.1.a) — returns
//! `ApexPlanError::InvalidInput("planner not yet implemented")`.
//! A2.1.b through A2.1.f flesh out each step.

use popkern_core::topology::{BRep, VertexHandle};

use crate::types::FilletSpec;

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::target_topology::TargetTopology;
use super::telemetry::PlanTelemetry;

/// A plan for resolving an apex vertex's fillet topology.
///
/// Produced by [`plan_apex_resolution`]; consumed by
/// [`construct_apex_topology`](super::super::apex_construction::construct_apex_topology)
/// (when A2.2 lands).  Introspectable for UI preview.
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
///
/// # Status
///
/// **Stub**: A2.1.a scaffolds the type API.  Actual implementation
/// lands in A2.1.b through A2.1.f.
pub fn plan_apex_resolution(
    _brep: &BRep,
    apex: VertexHandle,
    fillet_specs: &[FilletSpec],
    _tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<ApexPlan, ApexPlanError> {
    // A2.1.b will replace this with the real validation pass.
    if fillet_specs.len() < 3 {
        return Err(ApexPlanError::InsufficientFillets {
            apex,
            n: fillet_specs.len(),
        });
    }
    Err(ApexPlanError::InvalidInput {
        field: "planner",
        value: "stub — A2.1.b through A2.1.f not yet implemented".into(),
    })
}
