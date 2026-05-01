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

    let plan = plan_apex_resolution_lenient(brep, vertex, &specs, &surfaces, tolerance, None)
        .map_err(|plan_err| format!("plan_apex_resolution: {plan_err}"))?;
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
