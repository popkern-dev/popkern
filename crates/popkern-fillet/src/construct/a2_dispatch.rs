//! Phase 6c-A2.3.b — A2 apex-construction dispatch helper.
//!
//! Bridges the legacy `apply_blend_geometry` apex-vertex-blend
//! loop to the A2 planner + builder
//! ([`crate::apex_construction::plan_apex_resolution`] +
//! [`crate::apex_construction::construct_apex`]).
//!
//! Called from `apply_blend_geometry` only when
//! `FilletOptions::use_a2_apex_construction == true`.  Returns
//! `Ok(A2DispatchOutcome { built_faces })` on success (caller
//! skips the legacy `insert_vertex_blend` path and merges
//! `built_faces` into the post-loop protected-face set).  Returns
//! `Err(reason)` when A2 can't finish; caller emits
//! [`BlendWarning::A2FellBackToLegacy`](crate::error::BlendWarning::A2FellBackToLegacy)
//! and runs legacy as a backup.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::EdgeGeometry;

use crate::apex_construction::construct_apex;
use crate::apex_construction::plan::plan_apex_resolution_lenient;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{FilletRadius, FilletSpec};

/// What `try_a2_apex_construction` produced on success.
///
/// `built_faces` lists the cap face(s) `construct_apex` authored
/// (`ApexBuildResult::created_faces`).  The caller threads this
/// list into the protected-face set passed to `sliver_absorb` /
/// `cap_insert` / `late_wire_jobs` so those post-loop stages
/// don't re-author or absorb the A2-built cap.
#[derive(Debug, Clone, Default)]
pub(crate) struct A2DispatchOutcome {
    pub built_faces: Vec<FaceHandle>,
}

/// Try the A2 apex-construction path for one vertex.
///
/// On success, the BRep has been mutated (apex region built by
/// `construct_apex`) and `Ok(A2DispatchOutcome { built_faces })`
/// is returned.  On any A2 error (plan or build), the BRep is
/// rolled back by `construct_apex`'s own snapshot logic and
/// `Err(reason)` is returned describing why A2 couldn't finish.
pub(super) fn try_a2_apex_construction(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    _descendants: &HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> Result<A2DispatchOutcome, String> {
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

    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        eprintln!(
            "[a2-trace] try_a2_apex_construction: vertex={vertex:?}; \
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
            eprintln!("[a2-trace]   {fh:?}: {kind}");
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
                "[a2-trace]   {fh:?} outer loop ({} half-edges): {}",
                chain.len(),
                chain.join(" -> ")
            );
        }
    }
    let plan = plan_apex_resolution_lenient(brep, vertex, &specs, &surfaces, tolerance, None)
        .map_err(|plan_err| format!("plan_apex_resolution: {plan_err}"))?;

    // §A2.4 option 2 diagnostic: detect Phase-C cylinder-face
    // incoherence (loops with off-cyl vertices, per ADR §25).
    // Today this is purely informational under
    // `POPKERN_A2_TRACE=1` — it fires on both Sphere fixtures
    // (where A2 still succeeds because the planner finds shared
    // faces despite off-cyl loop vertices) and Gregory fixtures
    // (where A2 will fail downstream).  Future slices will use
    // this detection to drive the actual `kef`/`split_edge`
    // surgery that rebuilds coherent cylinder loops.
    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        let incoherent =
            crate::apex_construction::cylinder_face_coherence::find_off_cyl_loop_vertices(brep);
        if !incoherent.is_empty() {
            eprintln!(
                "[a2-trace] cylinder-face coherence: {} face(s) with off-cyl vertices — {}",
                incoherent.len(),
                crate::apex_construction::cylinder_face_coherence::describe_incoherence(
                    &incoherent
                ),
            );
            for entry in &incoherent {
                for (vh, pos, dist) in &entry.off_cyl_vertices {
                    eprintln!(
                        "[a2-trace]   {:?} (r={:.4}) v={vh:?} pos=({:.4},{:.4},{:.4}) radial_dist={:.4}",
                        entry.face, entry.radius, pos.x, pos.y, pos.z, dist,
                    );
                }
            }
        }
    }

    let build_result = construct_apex(brep, &plan, tolerance, None)
        .map_err(|build_err| format!("construct_apex: {build_err}"))?;
    // TODO(A2.4): when construct_apex returns substantive cap
    // geometry for Gregory plans, bridge `ApexBuildResult` into the
    // `corner_patches` list so downstream stages (continuity
    // diagnostics, `validate_all`) see this apex's patch info.
    Ok(A2DispatchOutcome {
        built_faces: build_result.created_faces,
    })
}
