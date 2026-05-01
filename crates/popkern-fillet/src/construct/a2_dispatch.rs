//! Phase 6c-A2.3.b — A2 apex-construction dispatch helper.
//!
//! Bridges the legacy `apply_blend_geometry` apex-vertex-blend
//! loop to the A2 planner + builder
//! ([`crate::apex_construction::plan_apex_resolution`] +
//! [`crate::apex_construction::construct_apex`]).
//!
//! Called from `apply_blend_geometry` only when
//! `FilletOptions::use_a2_apex_construction == true`.  Returns
//! `None` on A2 success (the caller skips the legacy
//! `insert_vertex_blend` path), or `Some(reason)` with a
//! human-readable error string when the A2 path could not
//! finish (the caller emits
//! [`BlendWarning::A2FellBackToLegacy`](crate::error::BlendWarning::A2FellBackToLegacy)
//! and runs the legacy path as a backup).

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use super::EdgeGeometry;

use crate::apex_construction::construct_apex;
use crate::apex_construction::plan::plan_apex_resolution_lenient;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{FilletRadius, FilletSpec};

/// Try the A2 apex-construction path for one vertex.
///
/// On success, the BRep has been mutated (apex region built by
/// `construct_apex`) and `None` is returned.  On any A2 error
/// (plan or build), the BRep is rolled back by `construct_apex`'s
/// own snapshot logic and `Some(reason)` is returned describing
/// why A2 couldn't finish.
///
/// **Currently never returns `None` for non-Open-rung apices** —
/// `construct_apex` errors on Sphere / Gregory cap construction
/// (A2.2.d/e/f are not yet implemented), so the dispatch always
/// hits the fallback path.  Open-rung plans (degenerate cases
/// where both Direct and Setback failed) succeed but produce an
/// empty `ApexBuildResult`, which leaves the apex region as the
/// pre-call BRep state — typically not what the caller wants for
/// production fillets but acceptable for the parallel-path
/// transition period.
pub(super) fn try_a2_apex_construction(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    _descendants: &HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> Option<String> {
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

    let plan = match plan_apex_resolution_lenient(brep, vertex, &specs, &surfaces, tolerance, None)
    {
        Ok(plan) => plan,
        Err(plan_err) => {
            return Some(format!("plan_apex_resolution: {plan_err}"));
        }
    };
    match construct_apex(brep, &plan, tolerance, None) {
        Ok(_build_result) => {
            // A2 success — caller skips legacy.  TODO(A2.2.d/e/f):
            // when construct_apex returns substantive cap
            // geometry, bridge ApexBuildResult into the
            // corner_patches list so downstream stages see this
            // apex's patch info.
            None
        }
        Err(build_err) => Some(format!("construct_apex: {build_err}")),
    }
}
