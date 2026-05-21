//! Legacy DirectApex mutating apex construction via Euler ops.
//!
//! The entrypoints here wrap the build in a whole-plan snapshot: on any
//! construction failure, the input BRep is restored before the typed error is
//! returned. Stage-specific implementation lives in `construct_*` modules.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::cancellation::CancellationToken;

use super::construct_cleanup::dissolve_v_bridges_until_valence_le_1;
use super::construct_edges::resolve_apex_edges;
pub use super::construct_result::ApexBuildResult;
use super::construct_surface::assign_cap_surface;
use super::construct_vertices::resolve_apex_vertices;
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::target_topology::TargetCapFace;

#[cfg(test)]
pub(crate) use super::construct_edges::attach_curve_to_edge;
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
    construct_apex_with_surfaces(brep, plan, &[], tolerance, cancel)
}

/// Variant of [`construct_apex`] that accepts the per-fillet
/// surface results.  These are required for DirectApex.4 Gregory cap
/// surface assembly (the bridge to legacy `three_quad_apply`
/// which needs `FilletSurfaceResult` to drive the K-P promotion
/// ladder).  Sphere caps don't need them and can call
/// [`construct_apex`] directly with an empty slice.
pub fn construct_apex_with_surfaces(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexBuildResult, ApexBuildError> {
    // Whole-plan snapshot — on any failure below, restore.
    let snapshot = brep.clone();
    match construct_apex_inner(brep, plan, fillet_surfaces, tolerance, cancel) {
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
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexBuildResult, ApexBuildError> {
    check_cancelled(cancel)?;

    // Open rung: nothing to construct.
    if plan.cap_kind.is_open() {
        debug_assert!(plan.target_topology.apex_vertices.is_empty());
        debug_assert!(plan.target_topology.apex_edges.is_empty());
        debug_assert!(matches!(plan.target_topology.cap_face, TargetCapFace::Open));
        return Ok(ApexBuildResult {
            created_faces: Vec::new(),
            apex_vertices: Vec::new(),
            apex_edges: Vec::new(),
            corner_patches: Vec::new(),
            exact_boundary_sectors: Vec::new(),
            exact_boundary_prepared_curve_count: 0,
            exact_boundary_inner_spoke_count: 0,
            exact_boundary_ribbon_count: 0,
            exact_boundary_solved_surface_count: 0,
            exact_boundary_pcurve_pair_count: 0,
            exact_boundary_performance: Default::default(),
            exact_boundary_report: crate::types::ExactBoundaryCapReport::not_applicable(
                "direct_apex",
            ),
            diagnostics: Vec::new(),
            telemetry: plan.telemetry.clone(),
        });
    }
    let apex_position = brep
        .vertices
        .get(plan.apex_vertex)
        .map(|v| v.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "apex vertex {:?} is not present in the BRep",
                plan.apex_vertex
            ),
        })?;

    // Perp-anchor relocation (§28.C): move anchors to plan targets before snap.
    super::perp_anchor_relocate::relocate_perp_anchors_to_plan_targets(brep, plan, tolerance)
        .map_err(|reason| ApexBuildError::PlanMismatch {
            reason: format!("§28.C relocate_perp_anchors_to_plan_targets: {reason}"),
        })?;
    check_cancelled(cancel)?;

    // Apex vertex resolution (DirectApex.2.b).
    //
    // For each plan apex vertex point, find or create a BRep
    // vertex at that position.  This slice implements the
    // "find" branch only: when no nearby BRep vertex exists,
    // the build returns a typed PlanMismatch.  DirectApex.2.b+ extends
    // with edge-subdivision-based "create" once the wiring
    // context is available.
    let apex_vertex_handles = resolve_apex_vertices(brep, plan, fillet_surfaces, tolerance)?;
    check_cancelled(cancel)?;

    check_cancelled(cancel)?;

    // Apex edge resolution (DirectApex.2.c + DirectApex.2.d).
    //
    // For each plan apex edge (origin_idx, dest_idx, curve, role),
    // find or create a BRep edge between the two resolved
    // vertices.  Find: matched in either direction.  Create: when
    // both endpoints lie on the same face's outer loop, `mef`
    // splits that face and the plan's curve is attached directly
    // to the new edge.  When no shared face exists (typical for
    // post-Phase-C BRep state where the rolling-ball footprints
    // sit on plane-face loops only), the create attempt errors
    // and the snapshot rolls back.
    let apex_edge_handles = resolve_apex_edges(brep, plan, &apex_vertex_handles, fillet_surfaces)?;
    check_cancelled(cancel)?;

    // Filter Some-edges for downstream consumers that index by edge,
    // not by plan-edge-idx.
    let valid_apex_edges: Vec<EdgeHandle> = apex_edge_handles.iter().filter_map(|e| *e).collect();

    // Apex teardown (DirectApex.2.e).
    //
    // After DirectApex.2.d's mef sequence each plane face is split into a
    // V-corner sub-face and a remainder.  V still has its
    // pre-mef stub edges (perp-edge splits from Phase A) which now
    // separate adjacent corner sub-faces.  Mirror legacy Phase
    // 2a/2b/3:
    //
    //   * `kef_cross_face_v_edges` walks V's incident edges and
    //     `kef`s each one separating distinct loops, merging the
    //     corner sub-faces.  Returns the surviving merged-face
    //     handle — the cap face.
    //   * `kev_v_bridge` removes V plus the trailing same-loop
    //     bridge in one Euler step (legacy Phase 2b).
    //   * `remove_orphaned_vertex` is the no-bridge safety net.
    //
    // The legacy helpers (`pub(crate)` in
    // `crate::construct::blend_insert`) handle the iterative kef
    // / outgoing repair / safety cap that DirectApex needs verbatim.
    if crate::direct_apex_trace::enabled() {
        let v = plan.apex_vertex;
        let mut incidences: Vec<(EdgeHandle, FaceHandle, FaceHandle, VertexHandle)> = Vec::new();
        for (eh, edge) in &brep.edges {
            let primary = edge.half_edge;
            let Some(he) = brep.half_edges.get(primary) else {
                continue;
            };
            let twin = brep.half_edges.get(he.twin);
            let touches = he.origin == v || twin.is_some_and(|t| t.origin == v);
            if !touches {
                continue;
            }
            let other = if he.origin == v {
                twin.map(|t| t.origin).unwrap_or(v)
            } else {
                he.origin
            };
            let face_a_h = brep.loops.get(he.loop_).map(|l| l.face);
            let face_b_h = twin.and_then(|t| brep.loops.get(t.loop_)).map(|l| l.face);
            let (Some(face_a), Some(face_b)) = (face_a_h, face_b_h) else {
                continue;
            };
            incidences.push((eh, face_a, face_b, other));
        }
        let v_pos = brep.vertices.get(v).map(|vv| vv.position);
        eprintln!(
            "[direct-apex-trace] pre-apex-teardown V={v:?}@{:?} has {} incident edge(s):",
            v_pos,
            incidences.len()
        );
        for (eh, fa, fb, other) in &incidences {
            let other_pos = brep.vertices.get(*other).map(|vv| vv.position);
            eprintln!(
                "[direct-apex-trace]   edge={eh:?} -> {other:?}@{:?} faces=({fa:?}, {fb:?})",
                other_pos
            );
        }
    }
    let kef_blend_face =
        crate::construct::blend_insert::kef_cross_face_v_edges(brep, plan.apex_vertex).map_err(
            |e| ApexBuildError::PlanMismatch {
                reason: format!("kef_cross_face_v_edges: {e}"),
            },
        )?;
    // §28.H Session 3.5 — multi-bridge V dissolve.  After Session
    // 3.6's stitches connect T_a/T_b to the perp-anchor mesh,
    // the cap face's outer loop is continuous and the dissolve
    // can safely remove V's residual same-loop bridges.
    dissolve_v_bridges_until_valence_le_1(brep, plan.apex_vertex);
    crate::construct::blend_insert::kev_v_bridge(brep, plan.apex_vertex).map_err(|e| {
        ApexBuildError::PlanMismatch {
            reason: format!("kev_v_bridge: {e}"),
        }
    })?;
    crate::construct::blend_insert::remove_orphaned_vertex(brep, plan.apex_vertex);

    // `kef_cross_face_v_edges` returns whichever face survived the
    // last kef.  When apex_edges are mef'd on multiple host faces
    // (Gregory cap = 2N edges on cyl + plane mix), the surviving
    // face may not be the one bounded by all 2N apex_edges; it may
    // be one of the host faces whose interior was unrelated to V.
    // Robust: find the face whose outer loop contains all 2N
    // apex_edge handles.  This is the actual cap face regardless
    // of which half-side the kef happened to keep.
    // Apex-region interior merge: targeted consolidation that only
    // kefs edges where one face is a cap-side accumulator (>=2
    // apex_edges), avoiding over-merge into host faces.
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] post-apex-teardown (kef_cross/kev_v): kef_blend_face={kef_blend_face:?}; \
             {} faces, {} edges",
            brep.faces.len(),
            brep.edges.len()
        );
        if let Some(face) = brep.faces.get(kef_blend_face)
            && let Some(loop_data) = brep.loops.get(face.outer_loop)
        {
            let start = loop_data.half_edge;
            let mut cursor = start;
            let mut iters = 0;
            let max_iters = brep.half_edges.len() + 1;
            eprintln!("[direct-apex-trace] kef_blend_face={kef_blend_face:?} outer_loop:");
            while let Some(he) = brep.half_edges.get(cursor) {
                let twin_face = brep
                    .half_edges
                    .get(he.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| l.face);
                let p = brep.vertices.get(he.origin).map(|v| v.position);
                eprintln!(
                    "[direct-apex-trace]   [{iters}] edge={:?} v={:?}@{:?} twin_face={twin_face:?}",
                    he.edge, he.origin, p
                );
                cursor = he.next;
                iters += 1;
                if cursor == start || iters >= max_iters {
                    break;
                }
            }
        }
        eprintln!("[direct-apex-trace] apex_edge face memberships:");
        for (i, &eh) in valid_apex_edges.iter().enumerate() {
            let mut faces: Vec<FaceHandle> = Vec::new();
            for (fh, face) in &brep.faces {
                let Some(loop_data) = brep.loops.get(face.outer_loop) else {
                    continue;
                };
                let start = loop_data.half_edge;
                let mut cursor = start;
                let mut iters = 0;
                let max_iters = brep.half_edges.len() + 1;
                let mut found = false;
                while let Some(he) = brep.half_edges.get(cursor) {
                    if he.edge == eh {
                        found = true;
                        break;
                    }
                    cursor = he.next;
                    iters += 1;
                    if cursor == start || iters >= max_iters {
                        break;
                    }
                }
                if found {
                    faces.push(fh);
                }
            }
            eprintln!("[direct-apex-trace]   apex_edge[{i}]={eh:?} → faces={faces:?}");
        }
    }
    super::gregory_cap::merge_cap_side_accumulators(brep, &valid_apex_edges).map_err(|e| {
        ApexBuildError::PlanMismatch {
            reason: format!("merge_cap_side_accumulators: {e}"),
        }
    })?;
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] post-merge_cap_side_accumulators: {} faces, {} edges; apex_edge face memberships:",
            brep.faces.len(),
            brep.edges.len()
        );
        for (i, &eh) in valid_apex_edges.iter().enumerate() {
            let mut faces: Vec<FaceHandle> = Vec::new();
            for (fh, face) in &brep.faces {
                let Some(loop_data) = brep.loops.get(face.outer_loop) else {
                    continue;
                };
                let start = loop_data.half_edge;
                let mut cursor = start;
                let mut iters = 0;
                let max_iters = brep.half_edges.len() + 1;
                let mut found = false;
                while let Some(he) = brep.half_edges.get(cursor) {
                    if he.edge == eh {
                        found = true;
                        break;
                    }
                    cursor = he.next;
                    iters += 1;
                    if cursor == start || iters >= max_iters {
                        break;
                    }
                }
                if found {
                    faces.push(fh);
                }
            }
            eprintln!("[direct-apex-trace]   apex_edge[{i}]={eh:?} → faces={faces:?}");
        }
    }
    let cap_face = super::gregory_cap::find_face_bounded_by_all_apex_edges(brep, &valid_apex_edges)
        .unwrap_or(kef_blend_face);

    if crate::direct_apex_trace::enabled()
        && let Some(face) = brep.faces.get(cap_face)
        && let Some(loop_data) = brep.loops.get(face.outer_loop)
    {
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut edges: Vec<(EdgeHandle, VertexHandle, FaceHandle)> = Vec::new();
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            let twin_face = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face)
                .unwrap_or(cap_face);
            edges.push((he.edge, he.origin, twin_face));
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        eprintln!(
            "[direct-apex-trace] cap_face={cap_face:?} outer_loop has {} edges:",
            edges.len()
        );
        for (i, (eh, vh, tf)) in edges.iter().enumerate() {
            let p = brep.vertices[*vh].position;
            eprintln!(
                "[direct-apex-trace]   [{i}] edge={eh:?} v={vh:?}@({:.3},{:.3},{:.3}) twin_face={tf:?}",
                p.x, p.y, p.z
            );
        }
    }
    check_cancelled(cancel)?;

    // Cap face surface assignment (DirectApex.2.f / DirectApex.4).
    //
    // Sphere caps land natively (DirectApex.2.f).  Gregory caps need 2N
    // sub-patch decomposition (DirectApex.4): bridge into legacy
    // `three_quad_apply::apply_three_quad_decomposition_to_blend`
    // by constructing a `BlendBoundary` from the post-stitch
    // BRep state + plan data + fillet_surfaces, then call the
    // K-P promotion ladder.
    let corner_patches = assign_cap_surface(
        brep,
        cap_face,
        plan,
        apex_position,
        &apex_vertex_handles,
        &apex_edge_handles,
        fillet_surfaces,
        tolerance,
    )?;
    check_cancelled(cancel)?;

    Ok(ApexBuildResult {
        created_faces: vec![cap_face],
        apex_vertices: apex_vertex_handles,
        apex_edges: apex_edge_handles,
        corner_patches,
        exact_boundary_sectors: Vec::new(),
        exact_boundary_prepared_curve_count: 0,
        exact_boundary_inner_spoke_count: 0,
        exact_boundary_ribbon_count: 0,
        exact_boundary_solved_surface_count: 0,
        exact_boundary_pcurve_pair_count: 0,
        exact_boundary_performance: Default::default(),
        exact_boundary_report: crate::types::ExactBoundaryCapReport::not_applicable("direct_apex"),
        diagnostics: Vec::new(),
        telemetry: plan.telemetry.clone(),
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
