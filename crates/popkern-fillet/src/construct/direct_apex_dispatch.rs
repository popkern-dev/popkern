//! Phase 6c-direct-apex dispatch.b — direct-apex apex-construction dispatch helper.
//!
//! Bridges the legacy `apply_blend_geometry` apex-vertex-blend
//! loop to the direct-apex planner + builder
//! ([`crate::apex_construction::plan_apex_resolution`] +
//! [`crate::apex_construction::construct_apex`]).
//!
//! Called from `apply_blend_geometry` only when
//! `FilletOptions::use_direct_apex_construction == true`.  Returns
//! `Ok(DirectApexDispatchOutcome { built_faces, corner_patches })` on success (caller
//! skips the legacy `insert_vertex_blend` path, merges `built_faces` into the
//! post-loop protected-face set, and threads patch reports to diagnostics).
//! Returns
//! `Err(reason)` when direct-apex can't finish; caller emits
//! [`BlendWarning::DirectApexFellBackToLegacy`](crate::error::BlendWarning::DirectApexFellBackToLegacy)
//! and runs legacy as a backup.

#[path = "direct_apex_acceptance.rs"]
mod direct_apex_acceptance;

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::EdgeGeometry;

use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::construct::construct_apex_with_surfaces;
use crate::apex_construction::constructive_builder::construct_apex_constructively;
use crate::apex_construction::plan::plan_apex_resolution_lenient;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{CornerPatchInfo, FilletRadius, FilletSpec};

/// What `try_direct_apex_construction` produced on success.
///
/// `built_faces` lists the cap face(s) `construct_apex` authored
/// (`ApexBuildResult::created_faces`).  The caller threads this
/// list into the protected-face set passed to `sliver_absorb` /
/// `cap_insert` / `late_wire_jobs` so those post-loop stages
/// don't re-author or absorb the direct-apex-built cap.
///
/// `corner_patches` carries any cap patch reports produced by
/// the apex builder so continuity diagnostics and validation see
/// native direct-apex cap faces.
#[derive(Debug, Clone, Default)]
pub(crate) struct DirectApexDispatchOutcome {
    pub built_faces: Vec<FaceHandle>,
    pub corner_patches: Vec<CornerPatchInfo>,
    pub exact_boundary_sector_count: usize,
    pub exact_boundary_prepared_curve_count: usize,
    pub exact_boundary_ribbon_count: usize,
    pub diagnostics: Vec<crate::types::ApexConstructionDiagnostic>,
}

/// Try the direct-apex apex-construction path for one vertex.
///
/// On success, the BRep has been mutated (apex region built by
/// `construct_apex`) and `Ok(DirectApexDispatchOutcome { built_faces })`
/// is returned.  On any direct-apex error (plan or build), the BRep is
/// rolled back by `construct_apex`'s own snapshot logic and
/// `Err(reason)` is returned describing why direct-apex couldn't finish.
pub(super) fn try_direct_apex_construction(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    _descendants: &HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> Result<DirectApexDispatchOutcome, String> {
    // Edge handles in `edge_geometries` are the ORIGINAL pre-
    // pipeline handles.  Phase C kef's the fillet edges so
    // those handles are no longer live in the BRep — but the
    // planner's lenient entry point (`plan_apex_resolution_
    // lenient`) skips the edge-existence and edge-incidence
    // checks that would otherwise fire.  The planner reads
    // geometry from `fillet_surfaces`, not from `spec.edge`,
    // so the stale handle is harmless metadata.
    let specs: Vec<FilletSpec> = edge_indices
        .iter()
        .map(|&i| {
            let eg = &edge_geometries[i];
            let radius = if (eg.radius_at_start - eg.radius_at_end).abs() < 1e-12 {
                FilletRadius::Constant(eg.radius_at_start)
            } else {
                FilletRadius::Variable {
                    start: eg.radius_at_start,
                    end: eg.radius_at_end,
                }
            };
            FilletSpec {
                edge: eg.analysis.edge,
                radius,
            }
        })
        .collect();
    let surfaces: Vec<FilletSurfaceResult> = edge_indices
        .iter()
        .map(|&i| {
            let eg = &edge_geometries[i];
            // EdgeGeometry's continuity_downgrade is
            // Option<(ContinuityLevel, String)>; the surface
            // result's is Option<String> (just the reason).
            // Map by extracting the second element.
            let downgrade_reason = eg.continuity_downgrade.as_ref().map(|(_, s)| s.clone());
            FilletSurfaceResult {
                surface: eg.surface.clone(),
                trim_curve_1: eg.trim_3d_1.clone(),
                trim_curve_2: eg.trim_3d_2.clone(),
                continuity_downgrade: downgrade_reason,
            }
        })
        .collect();

    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] try_direct_apex_construction: vertex={vertex:?}; \
             pre-call BRep has {} faces, {} edges, {} vertices",
            brep.faces.len(),
            brep.edges.len(),
            brep.vertices.len(),
        );
        for (fh, face) in &brep.faces {
            use popkern_core::geometry::surface::SurfaceGeometry;
            let kind = face
                .surface
                .as_ref()
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Plane { origin, normal, .. } => {
                        format!("Plane(origin={origin:?}, normal={normal:?})")
                    }
                    SurfaceGeometry::Cylinder {
                        origin,
                        axis,
                        radius,
                        ..
                    } => format!("Cylinder(origin={origin:?}, axis={axis:?}, r={radius:.3e})"),
                    SurfaceGeometry::Sphere { center, radius } => {
                        format!("Sphere(center={center:?}, r={radius:.3e})")
                    }
                    SurfaceGeometry::Torus { .. } => "Torus".to_string(),
                    SurfaceGeometry::Cone { .. } => "Cone".to_string(),
                    SurfaceGeometry::Nurbs(_) => "Nurbs".to_string(),
                    SurfaceGeometry::Gregory(_) => "Gregory".to_string(),
                })
                .unwrap_or_else(|| "None".to_string());
            eprintln!("[direct-apex-trace]   {fh:?}: {kind}");
        }
        // Dump every cylinder face's outer-loop walk so we can
        // see which edges bound which cylinder, and whether two
        // cylinder face handles share a seam edge.
        for (fh, face) in &brep.faces {
            let is_cyl = face.surface.as_ref().is_some_and(|s| {
                matches!(
                    s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                )
            });
            if !is_cyl {
                continue;
            }
            let outer_loop = face.outer_loop;
            let Some(loop_data) = brep.loops.get(outer_loop) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut cursor = start;
            let mut chain: Vec<String> = Vec::new();
            for _ in 0..(brep.half_edges.len() + 1) {
                let Some(he) = brep.half_edges.get(cursor) else {
                    break;
                };
                let origin_pos = brep.vertices[he.origin].position;
                chain.push(format!(
                    "v={:?}@({:.3},{:.3},{:.3})",
                    he.origin, origin_pos.x, origin_pos.y, origin_pos.z
                ));
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
            eprintln!(
                "[direct-apex-trace]   {fh:?} outer loop ({} half-edges): {}",
                chain.len(),
                chain.join(" -> ")
            );
        }
    }
    let plan = plan_apex_resolution_lenient(brep, vertex, &specs, &surfaces, tolerance, None)
        .map_err(|plan_err| format!("plan_apex_resolution: {plan_err}"))?;

    // §direct-apex cap-surface option 2 diagnostic: detect Phase-C cylinder-face
    // incoherence (loops with off-cyl vertices, per ADR §25).
    // Today this is purely informational under
    // `POPKERN_DIRECT_APEX_TRACE=1` — it fires on both Sphere fixtures
    // (where direct-apex still succeeds because the planner finds shared
    // faces despite off-cyl loop vertices) and Gregory fixtures
    // (where direct-apex will fail downstream).  Future slices will use
    // this detection to drive the actual `kef`/`split_edge`
    // surgery that rebuilds coherent cylinder loops.
    if crate::direct_apex_trace::enabled() {
        let incoherent =
            crate::apex_construction::cylinder_face_coherence::find_off_cyl_loop_vertices(brep);
        if !incoherent.is_empty() {
            eprintln!(
                "[direct-apex-trace] cylinder-face coherence: {} face(s) with off-cyl vertices — {}",
                incoherent.len(),
                crate::apex_construction::cylinder_face_coherence::describe_incoherence(
                    &incoherent
                ),
            );
            for entry in &incoherent {
                for (vh, pos, dist) in &entry.off_cyl_vertices {
                    eprintln!(
                        "[direct-apex-trace]   {:?} (r={:.4}) v={vh:?} pos=({:.4},{:.4},{:.4}) radial_dist={:.4}",
                        entry.face, entry.radius, pos.x, pos.y, pos.z, dist,
                    );
                }
            }
        }
    }

    // Stage 2/3/4/5 dispatch: Sphere AND Gregory both try the constructive
    // builder first; on PlanMismatch (e.g. pyramid Sphere with sphere-fit
    // T-positions that don't coincide with existing trim-trim corners,
    // any Gregory edge case the constructive path doesn't yet handle),
    // fall through to legacy direct-apex (§28.C-relocate path) so we don't regress
    // fixtures that were succeeding under legacy. Open is unchanged
    // (legacy short-circuits it as no-op).
    let build_result = match plan.cap_kind {
        CapKind::Sphere { .. } | CapKind::Gregory { .. } => {
            try_constructive_with_legacy_fallthrough(brep, &plan, &surfaces, tolerance)?
        }
        CapKind::Open { .. } => {
            construct_apex_with_surfaces(brep, &plan, &surfaces, tolerance, None)
                .map_err(|build_err| format!("construct_apex: {build_err}"))?
        }
    };
    Ok(DirectApexDispatchOutcome {
        built_faces: build_result.created_faces,
        corner_patches: build_result.corner_patches,
        exact_boundary_sector_count: build_result.exact_boundary_sectors.len(),
        exact_boundary_prepared_curve_count: build_result.exact_boundary_prepared_curve_count,
        exact_boundary_ribbon_count: build_result.exact_boundary_ribbon_count,
        diagnostics: build_result.diagnostics,
    })
}

/// Try the constructive builder for one apex; on PlanMismatch, fall through
/// to the legacy `construct_apex_with_surfaces`. Either result is returned
/// as `Ok(ApexBuildResult)`. If both fail, the legacy error is propagated
/// upward via `Err(String)`.
///
/// The constructive builder's snapshot-rollback ensures that on its own
/// failure, the BRep is byte-equivalent to the pre-call state — so the
/// subsequent legacy call sees the same input as if the constructive path
/// were never attempted.
fn try_constructive_with_legacy_fallthrough(
    brep: &mut BRep,
    plan: &crate::apex_construction::plan::ApexPlan,
    surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<crate::apex_construction::construct::ApexBuildResult, String> {
    let snapshot = brep.clone();
    let execution_options =
        crate::apex_construction::constructive_builder_options::ConstructiveApexExecutionOptions::from_environment_adapter();
    match construct_apex_constructively(brep, plan, surfaces, tolerance, None, execution_options) {
        Ok(mut result) => {
            if matches!(plan.cap_kind, CapKind::Gregory { .. }) {
                let constructive_err =
                    if !direct_apex_acceptance::native_gregory_has_expected_cylinder_seams(
                        brep,
                        &result,
                        plan.n_fillets(),
                    ) {
                        Some(format!(
                            "native Gregory cap has fewer than {} cap-to-cylinder seams; \
                             treating as incomplete direct-apex topology",
                            plan.n_fillets()
                        ))
                    } else {
                        direct_apex_acceptance::native_gregory_path_a_error(brep)
                            .map(|reason| {
                                format!(
                                    "native Gregory cap failed structural acceptance: {reason}; \
                                 treating as incomplete direct-apex topology"
                                )
                            })
                            .or_else(|| {
                                direct_apex_acceptance::native_gregory_cap_outer_g0_error(
                                    brep, &result, 5.0e-3,
                                )
                                .map(|reason| {
                                    format!(
                                        "native Gregory cap failed geometric acceptance: {reason}; \
                                         treating as incomplete direct-apex topology"
                                    )
                                })
                            })
                    };
                if let Some(constructive_err) = constructive_err {
                    if crate::direct_apex_trace::enabled() {
                        eprintln!(
                            "[direct-apex-trace] constructive fell through to legacy DirectApex: {constructive_err}"
                        );
                    }
                    *brep = snapshot;
                    return construct_apex_with_surfaces(brep, plan, surfaces, tolerance, None)
                        .map_err(|build_err| {
                            format!(
                                "construct_apex (after constructive fail: {constructive_err}): \
                                 {build_err}"
                            )
                        });
                }
                if let Some(diagnostic) =
                    direct_apex_acceptance::native_gregory_cap_outer_g1_diagnostic(
                        brep, &result, 5.0e-2,
                    )
                {
                    result.diagnostics.push(diagnostic);
                }
            }
            Ok(result)
        }
        Err(constructive_err) => {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace] constructive fell through to legacy DirectApex: {constructive_err}"
                );
            }
            construct_apex_with_surfaces(brep, plan, surfaces, tolerance, None).map_err(
                |build_err| {
                    format!(
                        "construct_apex (after constructive fail: {constructive_err}): \
                         {build_err}"
                    )
                },
            )
        }
    }
}
