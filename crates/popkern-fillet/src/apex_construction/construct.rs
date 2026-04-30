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

use popkern_core::geometry::Point3;
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

    // Stage 1 (A2.2.b) — apex vertex resolution.
    //
    // For each plan apex vertex point, find or create a BRep
    // vertex at that position.  This slice implements the
    // "find" branch only: when no nearby BRep vertex exists,
    // the build returns a typed PlanMismatch.  A2.2.b+ extends
    // with edge-subdivision-based "create" once the wiring
    // context is available.
    let apex_vertex_handles = resolve_apex_vertices(brep, plan, tolerance)?;
    check_cancelled(cancel)?;

    // Stage 2 (A2.2.c) — apex edge resolution.
    //
    // For each plan apex edge (origin_idx, dest_idx, curve, role),
    // find or create a BRep edge between the two resolved
    // vertices.  This slice implements the "find" branch: an
    // edge whose endpoints match (in either direction) is
    // accepted.  The "create via mef" branch needs loop context
    // (which face boundary the new edge splits) and lands
    // alongside A2.2.d (face boundary updates).
    let _apex_edge_handles = resolve_apex_edges(brep, plan, &apex_vertex_handles)?;
    check_cancelled(cancel)?;

    // Stages 3-4 (A2.2.d–A2.2.f) — face boundary updates and
    // cap face surface assignment.  Until those slices land,
    // callers driving non-Open plans see a typed "in progress"
    // error rather than a panic.
    let cap_kind_label = match &plan.cap_kind {
        CapKind::Sphere { .. } => "Sphere",
        CapKind::Gregory { .. } => "Gregory",
        CapKind::Open { .. } => unreachable!("Open handled above"),
    };
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "construct_apex: {cap_kind_label} apex-vertex + apex-edge resolution \
             succeeded; face boundary updates (A2.2.d/e) and cap face assembly \
             (A2.2.f) are not yet implemented"
        ),
    })
}

/// Locate a BRep vertex within `tolerance.linear()` of `point`.
/// Returns the closest match if any qualifies.  O(N) over the
/// brep's vertex slotmap — acceptable for the apex use case
/// where N is small (≤16 in realistic corners).
fn find_vertex_within(brep: &BRep, point: Point3, tol: f64) -> Option<VertexHandle> {
    let mut best: Option<(VertexHandle, f64)> = None;
    for (vh, vertex) in &brep.vertices {
        let d2 = {
            let dx = vertex.position.x - point.x;
            let dy = vertex.position.y - point.y;
            let dz = vertex.position.z - point.z;
            dx * dx + dy * dy + dz * dz
        };
        if d2 <= tol * tol {
            match best {
                Some((_, b)) if b <= d2 => {}
                _ => best = Some((vh, d2)),
            }
        }
    }
    best.map(|(vh, _)| vh)
}

/// Resolve each plan apex vertex to a BRep vertex handle.
///
/// **A2.2.b — find-only.**  For each
/// `plan.target_topology.apex_vertices[i]: Point3`, search the
/// BRep for an existing vertex within `tolerance.linear()` and
/// return its handle.  When no nearby vertex exists, return
/// [`ApexBuildError::PlanMismatch`] identifying which plan
/// position couldn't be resolved.
///
/// Returned vector is index-aligned with
/// `plan.target_topology.apex_vertices`.
///
/// Future extension (A2.2.b's create branch): when no nearby
/// vertex exists, locate the BRep edge whose curve passes
/// through `point` and split it via edge subdivision; return
/// the handle of the newly-created vertex.  That requires a
/// stable edge-subdivision Euler op that the kernel doesn't
/// expose yet.
fn resolve_apex_vertices(
    brep: &BRep,
    plan: &ApexPlan,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<VertexHandle>, ApexBuildError> {
    let tol = tolerance.linear();
    let mut out: Vec<VertexHandle> = Vec::with_capacity(plan.target_topology.apex_vertices.len());
    for (i, &point) in plan.target_topology.apex_vertices.iter().enumerate() {
        match find_vertex_within(brep, point, tol) {
            Some(vh) => out.push(vh),
            None => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "apex_vertex[{i}] at {point:?} has no BRep vertex within {tol:.3e} m \
                         (A2.2.b's find-or-create extension not yet implemented)"
                    ),
                });
            }
        }
    }
    Ok(out)
}

/// Find an existing BRep edge whose endpoints match the given
/// vertex pair (in either direction).  Iterates the edge
/// slotmap; O(M) where M is total edges (small for the apex
/// use case).
fn find_edge_between(
    brep: &BRep,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
) -> Option<EdgeHandle> {
    for (eh, edge) in &brep.edges {
        let primary = edge.half_edge;
        let Some(primary_he) = brep.half_edges.get(primary) else {
            continue;
        };
        let twin = primary_he.twin;
        let Some(twin_he) = brep.half_edges.get(twin) else {
            continue;
        };
        let endpoints = (primary_he.origin, twin_he.origin);
        if endpoints == (v_origin, v_dest) || endpoints == (v_dest, v_origin) {
            return Some(eh);
        }
    }
    None
}

/// Resolve each plan apex edge to a BRep edge handle.
///
/// **A2.2.c — find-only.**  For each
/// `plan.target_topology.apex_edges[i]`, look up the BRep
/// vertices at `apex_vertex_handles[origin_idx]` and
/// `apex_vertex_handles[dest_idx]` and search the BRep for an
/// existing edge between them (in either direction).  Returns
/// one edge handle per plan edge.  When no such edge exists,
/// returns [`ApexBuildError::PlanMismatch`] identifying the
/// plan edge that couldn't be resolved.
///
/// Future extension (A2.2.c+'s create branch): when the edge
/// doesn't exist, identify the loop both endpoints share and
/// `mef`-create the edge with the plan's `curve` attached
/// directly.  This requires the loop context computed by stage
/// A2.2.d (face boundary updates), so the create branch lands
/// alongside it.
fn resolve_apex_edges(
    brep: &BRep,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
) -> Result<Vec<EdgeHandle>, ApexBuildError> {
    let mut out: Vec<EdgeHandle> = Vec::with_capacity(plan.target_topology.apex_edges.len());
    for (i, target_edge) in plan.target_topology.apex_edges.iter().enumerate() {
        let v_origin = apex_vertex_handles
            .get(target_edge.origin_vertex_idx)
            .copied()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_edge[{i}] references origin_vertex_idx={} but plan has only {} apex_vertices",
                    target_edge.origin_vertex_idx,
                    apex_vertex_handles.len(),
                ),
            })?;
        let v_dest = apex_vertex_handles
            .get(target_edge.dest_vertex_idx)
            .copied()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_edge[{i}] references dest_vertex_idx={} but plan has only {} apex_vertices",
                    target_edge.dest_vertex_idx,
                    apex_vertex_handles.len(),
                ),
            })?;
        match find_edge_between(brep, v_origin, v_dest) {
            Some(eh) => out.push(eh),
            None => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "apex_edge[{i}] (vertex_idx {} → {}) has no BRep edge between the resolved endpoints \
                         (A2.2.c's find-or-create extension not yet implemented)",
                        target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                    ),
                });
            }
        }
    }
    Ok(out)
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
