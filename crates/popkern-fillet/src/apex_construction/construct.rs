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

use popkern_core::geometry::{CurveGeometry, CurveOnEdge, NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::cancellation::CancellationToken;
use crate::construct::blend_topo::find_halfedge_from_vertex_in_face;

use super::cap_kind::CapKind;
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::target_topology::{TargetCapFace, TargetEdge};
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
    construct_apex_with_surfaces(brep, plan, &[], tolerance, cancel)
}

/// Variant of [`construct_apex`] that accepts the per-fillet
/// surface results.  These are required for A2.4 Gregory cap
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

    // §27 commercial-grade design (Adopt + Create): we no longer
    // run the stitch-and-weld surgery here — that was a symptomatic
    // fix for the wrong-mef problem.  Stage 2 below now distinguishes
    // CapVsCyl (Adopt existing fillet-end arc) from CapVsPlane
    // (Create new chord via mef on V-incident face).  Re-enable the
    // stitch only as a last-resort fallback if both Adopt and
    // Create fail for cyl-side edges.
    //
    // §A2.4 option 2 surgery is now SKIPPED in the §27 Adopt path.
    // Set `POPKERN_A2_LEGACY_STITCH=1` to re-enable for diagnostics.
    if std::env::var("POPKERN_A2_LEGACY_STITCH").is_ok_and(|v| v == "1") {
        let stitched = super::cylinder_face_coherence::stitch_apex_edges_into_cyl_loops(
            brep,
            plan,
            &apex_vertex_handles,
        )
        .map_err(|reason| ApexBuildError::PlanMismatch {
            reason: format!("A2.4 option 2 stitch surgery: {reason}"),
        })?;
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!("[a2-trace] §A2.4 stitched {stitched} apex_edge(s) into cyl loops");
        }
    }
    check_cancelled(cancel)?;

    // Stage 2 (A2.2.c + A2.2.d) — apex edge resolution.
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

    // §27 Adopt+Create: skip interior vertex capture; stage 3's
    // V-spoke kef sequence produces a clean 2N-edge cap face when
    // stage 2 used Adopt for CapVsCyl.  Diagnostic re-enable via
    // POPKERN_A2_INTERIOR_ELIMINATION=1.
    let interior_vertex_set: std::collections::HashSet<VertexHandle> =
        if std::env::var("POPKERN_A2_INTERIOR_ELIMINATION").is_ok_and(|v| v == "1") {
            super::gregory_cap::compute_interior_vertex_set(
                brep,
                plan.apex_vertex,
                &apex_edge_handles,
            )
        } else {
            std::collections::HashSet::new()
        };

    // Stage 3 (A2.2.e) — apex teardown.
    //
    // After A2.2.d's mef sequence each plane face is split into a
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
    // / outgoing repair / safety cap that A2 needs verbatim.
    let kef_blend_face =
        crate::construct::blend_insert::kef_cross_face_v_edges(brep, plan.apex_vertex).map_err(
            |e| ApexBuildError::PlanMismatch {
                reason: format!("kef_cross_face_v_edges: {e}"),
            },
        )?;
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
    // Stage 3.5 — apex-region interior elimination ("carve out the
    // cap region").  Stage 3's `kef_cross_face_v_edges` only kefs
    // V's incident edges; Phase A perp-walks may have created perp-
    // anchor vertices INSIDE the cap perimeter that aren't V-
    // incident.  This pass kefs every cross-loop edge whose both
    // endpoints are in the interior set captured before stage 3.
    //
    // No-op when the post-stage-3 BRep is already coherent
    // (cube_uniform Sphere caps with no perp-anchor remnants).
    super::gregory_cap::eliminate_interior_edges(brep, &interior_vertex_set, &apex_edge_handles)
        .map_err(|e| ApexBuildError::PlanMismatch {
            reason: format!("eliminate_interior_edges: {e}"),
        })?;
    let cap_face =
        super::gregory_cap::find_face_bounded_by_all_apex_edges(brep, &apex_edge_handles)
            .unwrap_or(kef_blend_face);
    check_cancelled(cancel)?;

    // Stage 4 (A2.2.f / A2.4) — cap face surface assignment.
    //
    // Sphere caps land natively (A2.2.f).  Gregory caps need 2N
    // sub-patch decomposition (A2.4): bridge into legacy
    // `three_quad_apply::apply_three_quad_decomposition_to_blend`
    // by constructing a `BlendBoundary` from the post-stitch
    // BRep state + plan data + fillet_surfaces, then call the
    // K-P promotion ladder.
    assign_cap_surface(
        brep,
        cap_face,
        plan,
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
        telemetry: plan.telemetry.clone(),
    })
}

/// Assign the planned cap surface to the merged face.
///
/// * **Sphere** (A2.2.f) — maps directly to
///   [`SurfaceGeometry::Sphere`] (analytic, no NURBS conversion).
/// * **Gregory** (A2.4) — bridges into legacy
///   `three_quad_apply::apply_three_quad_decomposition_to_blend`
///   by constructing a [`BlendBoundary`] from the post-stitch
///   BRep state + plan + fillet_surfaces.  The K-P promotion
///   ladder handles the actual surface authoring; this function
///   just builds the input shape and assigns the resulting
///   sub-patches' surfaces to their (poked) sub-faces.
/// * **Open** is handled by [`construct_apex_inner`]'s early
///   return and is unreachable here.
#[allow(clippy::too_many_arguments)]
fn assign_cap_surface(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
    apex_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    _tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<(), ApexBuildError> {
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    match &plan.cap_kind {
        CapKind::Sphere { center, radius, .. } => {
            brep.faces[cap_face].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Sphere {
                    center: *center,
                    radius: *radius,
                },
                orientation: true,
            });
            Ok(())
        }
        CapKind::Gregory { .. } => super::gregory_cap::assemble_gregory_cap_surface(
            brep,
            cap_face,
            plan,
            apex_vertex_handles,
            apex_edge_handles,
            fillet_surfaces,
        ),
        CapKind::Open { .. } => unreachable!("Open handled in construct_apex_inner"),
    }
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

/// Resolve each plan apex vertex to a BRep vertex handle, finding
/// or creating as needed.
///
/// **A2.2.b — find + create.**  For each
/// `plan.target_topology.apex_vertices[i]: Point3`:
/// 1. Search the BRep for an existing vertex within
///    `tolerance.linear()`; if found, use it.
/// 2. Else, scan all BRep edges for one whose chord passes through
///    the target position and call
///    [`BRep::split_edge`](popkern_core::topology::BRep::split_edge)
///    to insert a new vertex at that position.
/// 3. If no edge passes through the target either, return a typed
///    [`ApexBuildError::PlanMismatch`] — the snapshot rolls back.
///
/// Returned vector is index-aligned with
/// `plan.target_topology.apex_vertices`.
///
/// **Chord projection caveat.**
/// [`crate::construct::helpers::find_edge_and_param_at`] uses
/// straight-line projection between edge endpoints — accurate for
/// `CurveGeometry::Line` edges (cube trim curves on planar faces
/// are Line variants) but approximate for arc / spline edges.  For
/// the cube fixtures landing today this is sufficient; richer
/// fixtures (curved supports, tori) may need a curve-aware
/// `closest_point_on_curve` projection later.
fn resolve_apex_vertices(
    brep: &mut BRep,
    plan: &ApexPlan,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<VertexHandle>, ApexBuildError> {
    let tol = tolerance.linear();
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let mut out: Vec<VertexHandle> = Vec::with_capacity(plan.target_topology.apex_vertices.len());
    for (i, &point) in plan.target_topology.apex_vertices.iter().enumerate() {
        if let Some(vh) = find_vertex_within(brep, point, tol) {
            if trace {
                eprintln!(
                    "[a2-trace] apex_vertex[{i}] at {point:?} → FOUND existing {vh:?} \
                     (within {tol:.3e})"
                );
            }
            out.push(vh);
            continue;
        }
        // Create branch: split the edge passing through `point`.
        let Some((edge, t)) = crate::construct::helpers::find_edge_and_param_at(brep, point) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_vertex[{i}] at {point:?}: no existing BRep vertex within \
                     {tol:.3e} m AND no edge whose chord passes through {point:?} \
                     (A2.2.b's create branch needs an edge to subdivide)"
                ),
            });
        };
        // Capture the split-edge's adjacent face handles before
        // mutation so the trace can show which faces inherit the
        // new vertex.
        let pre_faces = if trace {
            edge_adjacent_faces(brep, edge)
        } else {
            (None, None)
        };
        let split = brep
            .split_edge(edge, t)
            .map_err(|e| ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_vertex[{i}] at {point:?}: split_edge on edge {edge:?} at \
                     t={t:.6e} failed: {e}"
                ),
            })?;
        if trace {
            let loops = vertex_loop_face_handles(brep, split.vertex);
            eprintln!(
                "[a2-trace] apex_vertex[{i}] at {point:?} → CREATED {:?} via \
                 split_edge(edge={edge:?}, t={t:.6e}); pre-split faces=({:?}, {:?}); \
                 post-split outer-loop face memberships={loops:?}",
                split.vertex, pre_faces.0, pre_faces.1
            );
        }
        out.push(split.vertex);
    }
    Ok(out)
}

/// Helper for the `POPKERN_A2_TRACE` diagnostic: return the
/// `(face_a, face_b)` pair adjacent to a manifold edge.  Returns
/// `(None, None)` when handles are stale.
fn edge_adjacent_faces(brep: &BRep, edge: EdgeHandle) -> (Option<FaceHandle>, Option<FaceHandle>) {
    let Some(edata) = brep.edges.get(edge) else {
        return (None, None);
    };
    let primary = edata.half_edge;
    let Some(primary_he) = brep.half_edges.get(primary) else {
        return (None, None);
    };
    let face_a = brep.loops.get(primary_he.loop_).map(|l| l.face);
    let face_b = brep
        .half_edges
        .get(primary_he.twin)
        .and_then(|t| brep.loops.get(t.loop_))
        .map(|l| l.face);
    (face_a, face_b)
}

/// Helper for the `POPKERN_A2_TRACE` diagnostic: return every
/// face whose outer loop contains a half-edge originating from
/// `vertex`.
fn vertex_loop_face_handles(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut out: Vec<FaceHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        if let Some(loop_data) = brep.loops.get(he.loop_)
            && !out.contains(&loop_data.face)
        {
            out.push(loop_data.face);
        }
    }
    out
}

/// §27 Adopt path: scan for a pre-existing fillet-end arc on
/// `fillet_surfaces[fillet_idx]`'s cylinder face that connects
/// `v_origin` and `v_dest`.  Returns the BRep edge handle if such
/// an arc exists, `None` otherwise.
///
/// Phase B's fillet surface authoring writes each fillet's
/// transverse termination arc as a discrete BRep edge between
/// two perp-anchors.  When A2 inserts apex_vertices via
/// `split_edge` on the same perp-anchor positions, the existing
/// arc edge ends up between two apex_vertex_handles, ready to be
/// adopted as a CapVsCyl apex_edge with no further Euler op.
///
/// Match criterion: edge endpoints equal (v_origin, v_dest) in
/// either direction AND at least one adjacent face's surface is
/// a `Cylinder` matching the fillet's `(origin, axis, radius)`
/// within tolerance (1e-4 relative on radius, 1e-3 on axis dot,
/// 1e-3 m on origin distance to fillet axis).
fn adopt_fillet_end_arc(
    brep: &BRep,
    fillet_idx: usize,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Option<EdgeHandle> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let fillet = fillet_surfaces.get(fillet_idx)?;
    let SurfaceGeometry::Cylinder {
        origin: f_origin,
        axis: f_axis,
        radius: f_radius,
    } = &fillet.surface
    else {
        // Non-cylindrical fillet (Torus etc.) — adopt path doesn't
        // apply; caller falls back to find_edge_between or mef.
        return None;
    };
    let f_axis_unit = f_axis.normalized()?;
    for (eh, edge) in &brep.edges {
        let primary = edge.half_edge;
        let he = brep.half_edges.get(primary)?;
        let twin_he = brep.half_edges.get(he.twin)?;
        let endpoints = (he.origin, twin_he.origin);
        if endpoints != (v_origin, v_dest) && endpoints != (v_dest, v_origin) {
            continue;
        }
        let face_a = brep.loops.get(he.loop_).map(|l| l.face);
        let face_b = brep.loops.get(twin_he.loop_).map(|l| l.face);
        for fh_opt in [face_a, face_b] {
            let Some(fh) = fh_opt else { continue };
            let Some(face) = brep.faces.get(fh) else {
                continue;
            };
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let SurfaceGeometry::Cylinder {
                origin: face_origin,
                axis: face_axis,
                radius: face_radius,
            } = &sof.geometry
            else {
                continue;
            };
            // Match (origin on axis, axis direction, radius).
            if (face_radius - f_radius).abs() > f_radius.abs() * 1e-4 + 1e-9 {
                continue;
            }
            let Some(face_axis_unit) = face_axis.normalized() else {
                continue;
            };
            if face_axis_unit.dot(&f_axis_unit).abs() < 1.0 - 1e-3 {
                continue;
            }
            // Origin lies on the fillet's axis line.
            let d = popkern_core::geometry::Vec3::new(
                face_origin.x - f_origin.x,
                face_origin.y - f_origin.y,
                face_origin.z - f_origin.z,
            );
            let along = d.dot(&f_axis_unit);
            let perp = popkern_core::geometry::Vec3::new(
                d.x - along * f_axis_unit.x,
                d.y - along * f_axis_unit.y,
                d.z - along * f_axis_unit.z,
            );
            if perp.length() > 1e-3 {
                continue;
            }
            return Some(eh);
        }
    }
    None
}

/// Collect every face whose outer loop contains `vertex` as a
/// half-edge origin.  Used by the §27 V-incident chord preference
/// in `resolve_apex_edges` for CapVsPlane edges.
fn faces_with_vertex_on_outer_loop(
    brep: &BRep,
    vertex: VertexHandle,
) -> std::collections::HashSet<FaceHandle> {
    let mut out: std::collections::HashSet<FaceHandle> = std::collections::HashSet::new();
    for (fh, face) in &brep.faces {
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if he.origin == vertex {
                out.insert(fh);
                break;
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    out
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

/// Resolve each plan apex edge to a BRep edge handle, finding or
/// creating as needed.
///
/// **A2.2.c (find) + A2.2.d (create).**  For each
/// `plan.target_topology.apex_edges[i]`:
/// 1. Look up the resolved endpoint handles via
///    `apex_vertex_handles[origin_vertex_idx / dest_vertex_idx]`.
/// 2. **Find:** scan `brep.edges` for an existing edge between
///    those handles (in either direction).  If found, return it.
/// 3. **Create:** locate a face whose outer loop contains both
///    endpoints, call `mef` to split that face along the new
///    edge, and store the plan's `curve` directly on the new
///    edge as a [`CurveGeometry::Nurbs`].
///
/// The create branch errors out (rolling back via
/// [`construct_apex`]'s snapshot) when no face has both endpoints
/// on its outer loop.  This is the typical mid-pipeline state
/// today: rolling-ball footprint vertices live on plane-face
/// outer loops, but the planner picks one anchor per
/// geometrically-coincident pair and the picked anchors for two
/// adjacent apex vertices may sit on different plane faces — no
/// shared loop, no `mef` site.  Resolving that requires the
/// stage A2.2.e/f handoff contract; for now the dispatch falls
/// back to legacy with a typed warning.
fn resolve_apex_edges(
    brep: &mut BRep,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Result<Vec<EdgeHandle>, ApexBuildError> {
    use super::target_topology::TargetEdgeRole;
    let mut out: Vec<EdgeHandle> = Vec::with_capacity(plan.target_topology.apex_edges.len());
    // Track the cap-side accumulator face — the face that
    // successive mef calls should split, so each new apex_edge
    // attaches to the previously-mef'd cap-side sub-face.  Without
    // this preference, `find_shared_face_for_mef` picks any
    // face containing both endpoints, which on Gregory caps with
    // mixed cyl/plane host faces scatters apex_edges across
    // multiple faces and leaves stage 3's V-spoke kef sequence
    // unable to fully merge them.
    //
    // Initialised empty; populated with the new face produced by
    // each mef.  `find_shared_face_for_mef_preferring` consults
    // this set when multiple faces qualify.
    let mut cap_side_faces: std::collections::HashSet<FaceHandle> =
        std::collections::HashSet::new();
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
        // §27 Adopt path for CapVsCyl: scan existing BRep edges
        // between v_origin and v_dest, prefer one whose adjacent face
        // is the fillet's cylinder.  This re-uses Phase B's already-
        // authored fillet-end arc rather than mef'ing a new edge.
        if let TargetEdgeRole::CapVsCyl { fillet_idx } = target_edge.role
            && let Some(eh) =
                adopt_fillet_end_arc(brep, fillet_idx, v_origin, v_dest, fillet_surfaces)
        {
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!(
                    "[a2-trace] apex_edge[{i}] (vidx {} → {}) → ADOPTED existing \
                     fillet-end arc {eh:?} on fillet[{fillet_idx}]",
                    target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                );
            }
            out.push(eh);
            continue;
        }
        if let Some(eh) = find_edge_between(brep, v_origin, v_dest) {
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!(
                    "[a2-trace] apex_edge[{i}] (vidx {} → {}) → FOUND existing edge \
                     {eh:?} between {v_origin:?} and {v_dest:?}",
                    target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                );
            }
            out.push(eh);
            continue;
        }
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            let origin_faces = vertex_loop_face_handles(brep, v_origin);
            let dest_faces = vertex_loop_face_handles(brep, v_dest);
            let shared: Vec<_> = origin_faces
                .iter()
                .filter(|f| dest_faces.contains(f))
                .copied()
                .collect();
            eprintln!(
                "[a2-trace] apex_edge[{i}] (vidx {} → {}) NEED to create; \
                 v_origin={v_origin:?} on faces {origin_faces:?}; \
                 v_dest={v_dest:?} on faces {dest_faces:?}; \
                 SHARED faces = {shared:?}",
                target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
            );
        }
        // §27 Create path for CapVsPlane: prefer a V-incident face
        // (the V-corner sub-region of the plane between two adjacent
        // fillets).  Without this preference, mef can land on any
        // face containing both endpoints — including cylinder faces
        // (after stitch surgery) or far-side plane sub-faces.  The
        // V-incident preference ensures the new chord edge is on the
        // V-corner triangle, so stage 3's V-spoke kef merges the
        // corner sub-faces into the cap face naturally.
        let v_incident_faces = if matches!(target_edge.role, TargetEdgeRole::CapVsPlane { .. }) {
            faces_with_vertex_on_outer_loop(brep, plan.apex_vertex)
        } else {
            std::collections::HashSet::new()
        };
        let prefer_set: std::collections::HashSet<FaceHandle> = if !v_incident_faces.is_empty() {
            v_incident_faces
        } else {
            cap_side_faces.clone()
        };
        let (edge, mef_face_pair) =
            create_apex_edge_preferring(brep, v_origin, v_dest, target_edge, Some(&prefer_set))
                .map_err(|reason| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "apex_edge[{i}] (vertex_idx {} → {}) create via mef failed: {reason}",
                        target_edge.origin_vertex_idx, target_edge.dest_vertex_idx,
                    ),
                })?;
        // Record the new sub-face produced by this mef.  Either
        // half-edge of `edge` lives on a face: the "old" face
        // (which already existed pre-mef) and the "new" face
        // (created by mef).  We track both — preferred-face
        // matching just needs *any* face the previous mef created.
        cap_side_faces.insert(mef_face_pair.0);
        cap_side_faces.insert(mef_face_pair.1);
        out.push(edge);
    }
    Ok(out)
}

/// Locate a face whose outer loop contains both `v_origin` and
/// `v_dest` as half-edge origins.  Returns `(face, he_origin,
/// he_dest)` where the two half-edges originate from the matching
/// vertices on that face's outer loop.  O(F · L) where F is total
/// faces and L the average loop length — small for the apex use
/// case.
/// Search every face for one whose outer loop contains both
/// `v_origin` and `v_dest` as half-edge origins.  When `prefer` is
/// supplied, faces in the preference set are returned over equally
/// valid alternatives — used by Gregory cap construction to keep
/// successive `mef`-created apex_edges contiguous on a single
/// cap-side accumulator face.  Returns `(face, he_origin, he_dest)`
/// on success; `None` when no face qualifies.
fn find_shared_face_for_mef_preferring(
    brep: &BRep,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    prefer: Option<&std::collections::HashSet<FaceHandle>>,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let mut preferred: Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> = None;
    let mut fallback: Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> = None;
    for (fh, face) in &brep.faces {
        let outer_loop = face.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut he_origin: Option<HalfEdgeHandle> = None;
        let mut he_dest: Option<HalfEdgeHandle> = None;
        let mut cursor = start;
        let max_iters = brep.half_edges.len() + 1;
        let mut iters = 0;
        while let Some(he) = brep.half_edges.get(cursor) {
            if he.origin == v_origin && he_origin.is_none() {
                he_origin = Some(cursor);
            }
            if he.origin == v_dest && he_dest.is_none() {
                he_dest = Some(cursor);
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        if let (Some(ho), Some(hd)) = (he_origin, he_dest) {
            if prefer.is_some_and(|s| s.contains(&fh)) {
                preferred = Some((fh, ho, hd));
                break;
            }
            if fallback.is_none() {
                fallback = Some((fh, ho, hd));
            }
        }
    }
    preferred.or(fallback)
}

/// Create one apex edge by splitting a shared face via `mef` and
/// attaching the plan's curve directly.  Returns the new edge
/// handle and the `(old_face, new_face)` mef pair so callers can
/// chain successive mef sites onto a contiguous cap-side
/// accumulator via the `prefer` argument.  The caller is
/// responsible for rolling back (typically via the whole-plan
/// snapshot in [`construct_apex`]).
fn create_apex_edge_preferring(
    brep: &mut BRep,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    target_edge: &TargetEdge,
    prefer: Option<&std::collections::HashSet<FaceHandle>>,
) -> Result<(EdgeHandle, (FaceHandle, FaceHandle)), String> {
    let (face, he_origin, he_dest) =
        find_shared_face_for_mef_preferring(brep, v_origin, v_dest, prefer).ok_or_else(|| {
            format!(
                "no face has both v_origin {v_origin:?} and v_dest {v_dest:?} on its outer loop"
            )
        })?;
    // Defensive — find_halfedge_from_vertex_in_face guarantees
    // the same property find_shared_face_for_mef just verified,
    // but cheap to re-check and surfaces a clearer error if the
    // BRep is in an unexpected state.
    let he_origin_chk = find_halfedge_from_vertex_in_face(brep, face, v_origin)
        .map_err(|e| format!("origin half-edge lookup: {e}"))?;
    let he_dest_chk = find_halfedge_from_vertex_in_face(brep, face, v_dest)
        .map_err(|e| format!("dest half-edge lookup: {e}"))?;
    debug_assert_eq!(he_origin_chk, he_origin);
    debug_assert_eq!(he_dest_chk, he_dest);
    let mef_result = brep
        .mef(he_origin, he_dest)
        .map_err(|e| format!("mef on face {face:?}: {e:?}"))?;
    attach_curve_to_edge(brep, mef_result.edge, &target_edge.curve);
    Ok((mef_result.edge, (face, mef_result.face)))
}

/// Wrap the plan's `NurbsCurve` in a [`CurveOnEdge`] and attach
/// it to the edge.  Parameter range is read from the curve's
/// own knot vector — the planner builds the curve with endpoints
/// at the apex vertex positions, so the natural knot range is
/// the correct `t_start`/`t_end`.
///
/// `pub(crate)` for white-box test access; production callers go
/// through [`create_apex_edge`].
pub(crate) fn attach_curve_to_edge(brep: &mut BRep, edge: EdgeHandle, curve: &NurbsCurve) {
    let p = curve.degree;
    let n = curve.control_points.len() - 1;
    let t_start = curve.knots[p];
    let t_end = curve.knots[n + 1];
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve.clone()),
        t_start,
        t_end,
    });
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
