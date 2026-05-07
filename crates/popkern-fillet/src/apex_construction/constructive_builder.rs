//! Phase 6c-A2 §28.H Session 4 — constructive apex builder.
//!
//! Per ADR `2026-05-02-phase-6c-a2-constructive-apex-builder.md`. Replaces
//! the legacy carve-and-cleanup `construct_apex` with a deterministic
//! Euler-op sequence that takes an [`ApexPlan`] and authors the cap
//! topology directly: no carve, no merge accumulator, no dissolve cleanup.
//!
//! # Construction Sequence
//!
//! Per ADR §3 + §14.8:
//!
//! - **Apex vertex authoring** — author T-vertices. Sphere uniform: tight-match existing
//!   trim-trim corners. Gregory variable: `mev`-spurs from V into host plane
//!   sub-face interiors at the planner-given T-positions (§14.8).
//! - **Face-edge construction** — `mef` FaceEdge bridges across V-incident host plane
//!   sub-faces (Gregory only).
//! - **Cap-loop construction** — `kef` V's preexisting cross-face edge
//!   segments; merges V-corner fragments into one cap face. The last
//!   residual edge becomes a same-face bridge consumed by apex removal.
//! - **Fillet-end edge construction** — `mef` FilletEnd arc edges on the cap
//!   face between T_a[k] and T_b[k] (§14.8). Gregory only.
//! - **Apex removal** — `kev V` via the residual same-face bridge plus
//!   `dissolve_same_loop_bridge` for any MEV spurs that became same-face
//!   bridges after cap-loop construction.
//! - **Cap surface assignment** — assign cap surface (Sphere directly; Gregory delegates
//!   to legacy `assemble_gregory_cap_surface`).
//!
//! # Atomicity
//!
//! All mutations happen inside [`construct_apex_inner_constructive`], wrapped
//! by [`construct_apex_constructively`]'s `BRep::clone()` snapshot. On any
//! error the BRep is restored byte-equivalently before the error is returned.
//!
//! # Status
//!
//! The builder now routes supported cap plans through these semantic topology
//! operations. Unsupported plans still decline through structured errors and
//! snapshot rollback.
#![allow(dead_code)]

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::cancellation::CancellationToken;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};

use super::apex_vertex_authoring::{author_gregory_apex_vertices, author_sphere_apex_vertices};
use super::apex_vertex_removal::remove_apex_vertex;
use super::cap_kind::CapKind;
use super::cap_loop_construction::{
    collect_required_fillet_ends, construct_cap_loop_from_preexisting_cross_face_edges,
    create_remaining_fillet_end_edges,
};
use super::cap_surface_assignment::{
    assign_gregory_cap_surface, assign_sphere_cap_surface, gregory_sectors_from_handles,
};
use super::construct::ApexBuildResult;
use super::constructive_geometry::cylinder_arc_through_edge_endpoints;
use super::error::ApexBuildError;
use super::face_edge_construction::{
    construct_sphere_cap_face_from_v_incident_edges, create_gregory_face_edges,
};
use super::plan::ApexPlan;
use super::target_topology::TargetCapFace;

/// Apply an [`ApexPlan`] to the BRep using the constructive topology
/// algorithm. Drop-in replacement for [`super::construct::construct_apex_with_surfaces`]
/// for supported apex plans.
///
/// Snapshot-rollback wrapper around [`construct_apex_inner_constructive`].
pub fn construct_apex_constructively(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexBuildResult, ApexBuildError> {
    let snapshot = brep.clone();
    match construct_apex_inner_constructive(brep, plan, fillet_surfaces, tolerance, cancel) {
        Ok(result) => {
            if let Err(reason) = validate_constructive_topology_postconditions(brep) {
                *brep = snapshot;
                Err(ApexBuildError::PlanMismatch {
                    reason: format!("constructive topology postcondition failed: {reason}"),
                })
            } else {
                Ok(result)
            }
        }
        Err(err) => {
            *brep = snapshot;
            Err(err)
        }
    }
}

fn validate_constructive_topology_postconditions(brep: &BRep) -> Result<(), String> {
    brep.validate_euler_poincare()
        .and_then(|_| brep.validate_half_edge_pairing())
        .and_then(|_| brep.validate_half_edge_cycles())
        .and_then(|_| brep.validate_vertex_consistency())
        .and_then(|_| brep.validate_manifold_edges())
        .and_then(|_| brep.validate_vertex_valence())
        .map_err(|e| format!("{e:?}"))
}

/// Inner builder; mutations happen here. Wrapped by
/// [`construct_apex_constructively`] for atomicity.
fn construct_apex_inner_constructive(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexBuildResult, ApexBuildError> {
    check_cancelled(cancel)?;

    // Open rung — no topology to construct.
    if plan.cap_kind.is_open() {
        debug_assert!(plan.target_topology.apex_vertices.is_empty());
        debug_assert!(plan.target_topology.apex_edges.is_empty());
        debug_assert!(matches!(plan.target_topology.cap_face, TargetCapFace::Open));
        return Ok(ApexBuildResult {
            created_faces: Vec::new(),
            apex_vertices: Vec::new(),
            apex_edges: Vec::new(),
            corner_patches: Vec::new(),
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

    match &plan.cap_kind {
        CapKind::Sphere {
            center: cap_center,
            radius: cap_radius,
            ..
        } => {
            // Sphere path: adopt nearby apex-neighborhood anchors to the
            // planned target positions, then tight-match T-vertices, remove
            // V-incident edges, remove V, and assign the Sphere surface.
            // No Gregory FaceEdge MEFs and no cap-loop FilletEnd MEFs.
            super::perp_anchor_relocate::relocate_perp_anchors_to_plan_targets(
                brep, plan, tolerance,
            )
            .map_err(|reason| ApexBuildError::PlanMismatch {
                reason: format!("sphere apex anchor adoption: {reason}"),
            })?;
            let apex_vertex_handles =
                author_sphere_apex_vertices(brep, plan, fillet_surfaces, tolerance)?;
            check_cancelled(cancel)?;
            let cap_face = construct_sphere_cap_face_from_v_incident_edges(brep, plan)?;
            check_cancelled(cancel)?;
            remove_apex_vertex(brep, plan, cap_face, false)?;
            check_cancelled(cancel)?;
            assign_sphere_cap_surface(brep, cap_face, *cap_center, *cap_radius)?;
            let corner_patches = vec![CornerPatchInfo {
                vertex: plan.apex_vertex,
                n_sides: plan.cap_kind.n_fillets(),
                method: CornerMethod::RollingBallSphere,
                achieved_continuity: ContinuityLevel::G1,
                classification: crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
                face: cap_face,
                bridge_fallbacks: Vec::new(),
                overlap_remediated: false,
                relaxed_continuity_coplanarity: None,
                cap_sphere: Some((*cap_center, *cap_radius)),
            }];

            let n = plan.cap_kind.n_fillets();
            Ok(ApexBuildResult {
                created_faces: vec![cap_face],
                apex_vertices: apex_vertex_handles,
                // Sphere apex_edges are all CapVsCyl by §28.D — represented
                // as `None` because they're not separately mef'd; they're
                // existing trim-trim segments on the cap face's outer loop.
                apex_edges: vec![None; n],
                corner_patches,
                telemetry: plan.telemetry.clone(),
            })
        }
        CapKind::Gregory { .. } => {
            // Gregory path: MEV T-vertices into host plane sub-faces (§14.8)
            // -> MEF FaceEdges -> construct the cap loop from preexisting
            // cross-face edges -> MEF FilletEnds on the cap face -> remove V
            // -> delegate to the Gregory surface assembler.
            //
            // On any construction error the snapshot wrapper rolls back and
            // the caller falls through to legacy A2.

            // Capture V's original cross-face incidences before anchor MEVs.
            // Cap-loop construction removes only these edges; anchor spurs
            // become cross-face later and must survive into fillet-end
            // construction and apex removal.
            let preexisting_cross_face_edges: Vec<popkern_core::topology::EdgeHandle> = brep
                .half_edges
                .iter()
                .filter(|(_, h)| h.origin == plan.apex_vertex)
                .map(|(_, h)| h.edge)
                .collect();

            let apex_vertex_handles =
                author_gregory_apex_vertices(brep, plan, fillet_surfaces, tolerance)?;
            check_cancelled(cancel)?;
            let face_edge_handles = create_gregory_face_edges(brep, plan, &apex_vertex_handles)?;
            check_cancelled(cancel)?;
            let (mut cap_face, mut fillet_end_handles) =
                construct_cap_loop_from_preexisting_cross_face_edges(
                    brep,
                    plan,
                    &preexisting_cross_face_edges,
                    &apex_vertex_handles,
                    fillet_surfaces,
                )?;
            check_cancelled(cancel)?;
            let repeated_apex_normalized = create_remaining_fillet_end_edges(
                brep,
                plan,
                &mut cap_face,
                &apex_vertex_handles,
                &face_edge_handles,
                fillet_surfaces,
                &mut fillet_end_handles,
            )?;
            let fillet_end_handles = collect_required_fillet_ends(fillet_end_handles)?;
            check_cancelled(cancel)?;
            remove_apex_vertex(brep, plan, cap_face, repeated_apex_normalized)?;
            check_cancelled(cancel)?;
            let corner_patches = assign_gregory_cap_surface(
                brep,
                cap_face,
                plan,
                apex_position,
                &apex_vertex_handles,
                &face_edge_handles,
                &fillet_end_handles,
                &gregory_sectors_from_handles(plan, &fillet_end_handles, &face_edge_handles)?,
                fillet_surfaces,
                tolerance,
            )?;
            if repeated_apex_normalized
                && std::env::var_os("POPKERN_A2_ENABLE_TOPOLOGY_EDGE_CURVE_CONTAINMENT").is_some()
            {
                complete_cap_patch_edge_geometry_from_topology(brep, &corner_patches);
            }

            let n = plan.cap_kind.n_fillets();
            // Build apex_edges return value: 2N entries, alternating
            // FilletEnd (even) / FaceEdge (odd), each Some(handle).
            let mut apex_edges: Vec<Option<popkern_core::topology::EdgeHandle>> =
                Vec::with_capacity(2 * n);
            for k in 0..n {
                apex_edges.push(Some(fillet_end_handles[k]));
                apex_edges.push(Some(face_edge_handles[k]));
            }
            let mut created_faces = vec![cap_face];
            for patch in &corner_patches {
                if !created_faces.contains(&patch.face) {
                    created_faces.push(patch.face);
                }
            }
            Ok(ApexBuildResult {
                created_faces,
                apex_vertices: apex_vertex_handles,
                apex_edges,
                corner_patches,
                telemetry: plan.telemetry.clone(),
            })
        }
        CapKind::Open { .. } => unreachable!("Open handled above"),
    }
}

// ---------------------------------------------------------------------------
// Cap edge completion
// ---------------------------------------------------------------------------

/// Diagnostic containment for repeated-apex topology experiments.
///
/// This remains opt-in through
/// `POPKERN_A2_ENABLE_TOPOLOGY_EDGE_CURVE_CONTAINMENT`: accepted native A2
/// construction must author topology, 3D curves, pcurves, and Gregory boundary
/// sectors from the same physical FilletEnd/FaceEdge ownership at construction
/// time.
pub(crate) fn complete_cap_patch_edge_geometry_from_topology(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
) {
    let cap_faces = corner_patches
        .iter()
        .map(|patch| patch.face)
        .collect::<std::collections::BTreeSet<_>>();
    let edge_handles = brep
        .half_edges
        .iter()
        .filter_map(|(_, he)| {
            let face = brep.loops.get(he.loop_)?.face;
            cap_faces.contains(&face).then_some(he.edge)
        })
        .collect::<std::collections::BTreeSet<_>>();

    for edge in edge_handles {
        complete_edge_geometry_from_topological_endpoints(brep, edge);
    }
}

fn complete_edge_geometry_from_topological_endpoints(brep: &mut BRep, edge: EdgeHandle) {
    let Some(edge_data) = brep.edges.get(edge) else {
        return;
    };
    let he = edge_data.half_edge;
    let Some(he_data) = brep.half_edges.get(he) else {
        return;
    };
    let twin = he_data.twin;
    let (Some(a), Some(b)) = (
        brep.vertices.get(he_data.origin),
        brep.half_edges
            .get(twin)
            .and_then(|twin_data| brep.vertices.get(twin_data.origin)),
    ) else {
        return;
    };
    let needs_curve = match edge_data.curve.as_ref() {
        Some(curve_on_edge) => {
            let p_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
            let p_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
            let forward_cost = a.position.distance_to(&p_start) + b.position.distance_to(&p_end);
            let reversed_cost = a.position.distance_to(&p_end) + b.position.distance_to(&p_start);
            forward_cost.min(reversed_cost) > 1.0e-5
        }
        None => true,
    };
    if needs_curve {
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/edge-curve] completing cap edge geometry on {edge:?} from topology"
            );
        }
        brep.edges[edge].curve = curve_on_support_or_line(brep, he, twin, a.position, b.position);
    }
    if let Some(curve) = brep.edges[edge].curve.clone() {
        write_projected_pcurve_from_curve(brep, he, &curve);
        write_projected_pcurve_from_curve(brep, twin, &curve);
    }
}

fn curve_on_support_or_line(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    twin: popkern_core::topology::HalfEdgeHandle,
    start: popkern_core::geometry::Point3,
    end: popkern_core::geometry::Point3,
) -> Option<CurveOnEdge> {
    for half_edge in [he, twin] {
        let face = brep.loops[brep.half_edges[half_edge].loop_].face;
        let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        if let Some((geometry, t_start, t_end)) =
            cylinder_arc_through_edge_endpoints(&surface.geometry, start, end)
        {
            return Some(CurveOnEdge {
                geometry,
                t_start,
                t_end,
            });
        }
    }
    CurveOnEdge::line_segment(start, end, 1.0e-9)
}

fn write_projected_pcurve_from_curve(
    brep: &mut BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    curve: &CurveOnEdge,
) {
    let Some(he_data) = brep.half_edges.get(he) else {
        return;
    };
    let face = brep.loops[he_data.loop_].face;
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    else {
        return;
    };
    let Some(origin) = brep.vertices.get(he_data.origin).map(|v| v.position) else {
        return;
    };
    let p_start = curve.geometry.point_at(curve.t_start);
    let p_end = curve.geometry.point_at(curve.t_end);
    let forward = origin.distance_to(&p_start) <= origin.distance_to(&p_end);
    let Some(pcurve) = fit_projected_pcurve(&surface, curve, forward) else {
        return;
    };
    brep.half_edges[he].pcurve = Some(pcurve.into());
}

fn fit_projected_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveOnEdge,
    forward: bool,
) -> Option<NurbsCurve> {
    let n_samples = 12usize;
    let mut points = Vec::with_capacity(n_samples);
    let mut params = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let oriented_s = if forward { s } else { 1.0 - s };
        let t = curve.t_start + oriented_s * (curve.t_end - curve.t_start);
        let point = curve.geometry.point_at(t);
        let (u, v) = crate::construct::helpers::project_to_uv(surface, &point);
        points.push(popkern_core::geometry::Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    popkern_core::geometry::fitting::fit_curve(&points, &params, degree, n_cps).ok()
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn check_cancelled(cancel: Option<&CancellationToken>) -> Result<(), ApexBuildError> {
    if let Some(token) = cancel
        && token.is_cancelled()
    {
        return Err(ApexBuildError::Cancelled);
    }
    Ok(())
}

#[cfg(test)]
#[path = "constructive_builder_tests.rs"]
mod tests;
