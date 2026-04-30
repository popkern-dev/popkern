//! Phase 6c-A2.2 — direct apex construction via Euler ops.
//!
//! [`construct_apex`] is the **mutating** counterpart to
//! [`plan_apex_resolution`](super::plan::plan_apex_resolution): it
//! consumes a precomputed [`ApexPlan`] and applies Euler operators
//! to realise the plan's [`TargetTopology`] in the BRep.
//!
//! Per ADR §6.3 the construction is structured around three
//! snapshot levels:
//! 1. **Whole-plan snapshot** (taken at entry): on any builder
//!    failure, the BRep is restored from this snapshot before the
//!    error is returned, so the input BRep is never left in a
//!    half-mutated state.
//! 2. **Per-step checkpoints** (optional, off by default; helpful
//!    for tests): rollback to the last good state.
//! 3. **Op-level pre-condition checks**: every Euler op is
//!    validated before mutating, with typed errors on rejection.
//!
//! Build order (mirrors the four `TargetTopology` field groups):
//! * Apex vertices — find-or-create at each `apex_vertices[i]`
//!   point.
//! * Apex edges — `mef`-create with the curve assigned directly
//!   (no separate wire step).
//! * Face boundary updates — truncate non-cap face boundaries
//!   that previously ended at rolling-ball footprints.
//! * Cap face — assign Sphere or Gregory surface from the plan.
//!
//! A2.2.a (this slice) lands the public API surface and the
//! Open-rung no-op handler.  A2.2.b through A2.2.f flesh out
//! the per-stage logic.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::cancellation::CancellationToken;

use super::cap_kind::CapKind;
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::target_topology::TargetCapFace;
use super::telemetry::PlanTelemetry;

/// What `construct_apex` actually built in the BRep.
///
/// All handles refer to entities in the BRep that was passed to
/// [`construct_apex`].  For [`CapKind::Open`] plans the result is
/// trivially empty (the Open rung is a no-op for the builder).
///
/// Returned by both [`construct_apex`] and `resolve_apex` (A2.3).
#[derive(Debug, Clone)]
pub struct ApexBuildResult {
    /// New cap face(s) created by the builder.  Length 1 for
    /// Sphere, 1 for Gregory (the N-sided cap face — sub-patch
    /// decomposition happens in A2.4), 0 for Open.
    pub created_faces: Vec<FaceHandle>,
    /// Apex vertices in the realised topology, indexed parallel
    /// to `plan.target_topology.apex_vertices`.  Mix of newly
    /// created (via `mev`) and pre-existing vertices (within
    /// tolerance of the plan's target points).
    pub apex_vertices: Vec<VertexHandle>,
    /// New edges in the apex region, indexed parallel to
    /// `plan.target_topology.apex_edges`.  Carry the curve3d
    /// from the plan as their stored geometry.
    pub apex_edges: Vec<EdgeHandle>,
    /// Telemetry passed through from the plan (the builder
    /// doesn't add new diagnostic values; A2.5 may extend this
    /// with construction-side numerics).
    pub telemetry: PlanTelemetry,
}

/// Apply a precomputed [`ApexPlan`] to the BRep.
///
/// The function is **transactional**: on success the BRep
/// reflects the plan's target topology; on failure the BRep is
/// restored to its pre-call state and the error is returned.
///
/// # Inputs
///
/// * `brep` — the BRep being filleted.  Mutated in place on
///   success.
/// * `plan` — a previously-computed [`ApexPlan`] (typically from
///   [`plan_apex_resolution`](super::plan::plan_apex_resolution)).
/// * `tolerance` — geometric tolerance for `find-or-create`
///   vertex matching and pre-condition checks.
/// * `cancel` — optional cooperative cancellation token.  Polled
///   at each major build stage; on cancellation the BRep is
///   rolled back and [`ApexBuildError::Cancelled`] is returned.
///
/// # Output
///
/// `Ok(ApexBuildResult)` on success, with handles to the
/// realised topology.  `Err(ApexBuildError)` on failure (with
/// the BRep restored).
///
/// # Errors
///
/// See [`ApexBuildError`] for the failure taxonomy:
/// - [`ApexBuildError::EulerOpFailed`] for Euler-op
///   pre-condition violations.
/// - [`ApexBuildError::ValidationFailed`] for post-build
///   `validate()` errors.
/// - [`ApexBuildError::PlanMismatch`] when the plan refers to
///   BRep entities that no longer exist or have drifted.
/// - [`ApexBuildError::Cancelled`] when the supplied token is
///   tripped.
pub fn construct_apex(
    brep: &mut BRep,
    plan: &ApexPlan,
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexBuildResult, ApexBuildError> {
    // Whole-plan snapshot — on any failure below, restore.
    let snapshot = brep.clone();
    match construct_apex_inner(brep, plan, tolerance, cancel) {
        Ok(result) => Ok(result),
        Err(err) => {
            *brep = snapshot;
            Err(err)
        }
    }
}

/// Inner builder that performs the actual mutations.  Wrapped by
/// [`construct_apex`] so any error triggers a snapshot restore.
fn construct_apex_inner(
    brep: &mut BRep,
    plan: &ApexPlan,
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexBuildResult, ApexBuildError> {
    check_cancelled(cancel)?;

    // Open rung: nothing to construct.  The plan's
    // target_topology is empty and the cap face is
    // TargetCapFace::Open.  Return an empty result — fillets
    // remain at their legacy rolling-ball terminations as the
    // caller wired them.
    if plan.cap_kind.is_open() {
        debug_assert!(plan.target_topology.apex_vertices.is_empty());
        debug_assert!(plan.target_topology.apex_edges.is_empty());
        debug_assert!(matches!(plan.target_topology.cap_face, TargetCapFace::Open));
        return Ok(ApexBuildResult {
            created_faces: Vec::new(),
            apex_vertices: Vec::new(),
            apex_edges: Vec::new(),
            telemetry: plan.telemetry.clone(),
        });
    }

    // Sphere and Gregory paths come in A2.2.b through A2.2.f.
    // Until those slices land, callers driving non-Open plans
    // see a typed "unimplemented" error rather than a panic.
    let _ = (brep, tolerance);
    let cap_kind_label = match &plan.cap_kind {
        CapKind::Sphere { .. } => "Sphere",
        CapKind::Gregory { .. } => "Gregory",
        CapKind::Open { .. } => unreachable!("Open handled above"),
    };
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "construct_apex: {cap_kind_label} cap construction is not yet \
             implemented (A2.2.b–A2.2.f); only Open rung is currently supported"
        ),
    })
}

/// Polls the cancellation token (if any) and returns
/// [`ApexBuildError::Cancelled`] when tripped.
fn check_cancelled(cancel: Option<&CancellationToken>) -> Result<(), ApexBuildError> {
    if cancel.is_some_and(|c| c.is_cancelled()) {
        Err(ApexBuildError::Cancelled)
    } else {
        Ok(())
    }
}
