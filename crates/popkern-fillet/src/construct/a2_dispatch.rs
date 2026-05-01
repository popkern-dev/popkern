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

use crate::apex_construction::{construct_apex, plan_apex_resolution};
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
    descendants: &HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> Option<String> {
    // Phase A's perpendicular-edge splits replaced each
    // original fillet edge with multiple descendants.  For the
    // planner's edge-incidence validation to pass we have to
    // pass a current edge handle that's incident to the apex
    // vertex.  Walk the descendants map and pick the one
    // descendant edge whose primary or twin half-edge has
    // `vertex` as `origin`.
    let resolve_current_edge = |orig: EdgeHandle| -> Option<EdgeHandle> {
        let descs = descendants.get(&orig)?;
        for &d in descs {
            let edge = brep.edges.get(d)?;
            let primary = brep.half_edges.get(edge.half_edge)?;
            let twin = brep.half_edges.get(primary.twin)?;
            if primary.origin == vertex || twin.origin == vertex {
                return Some(d);
            }
        }
        None
    };

    let mut specs: Vec<FilletSpec> = Vec::with_capacity(edge_indices.len());
    for (k, &i) in edge_indices.iter().enumerate() {
        let eg = &edge_geometries[i];
        let radius = if (eg.radius_at_start - eg.radius_at_end).abs() < 1e-12 {
            FilletRadius::Constant(eg.radius_at_start)
        } else {
            FilletRadius::Variable {
                start: eg.radius_at_start,
                end: eg.radius_at_end,
            }
        };
        let edge_now = match resolve_current_edge(eg.analysis.edge) {
            Some(e) => e,
            None => {
                return Some(format!(
                    "post-Phase-A descendant of original edge \
                     {:?} (apex_edges[{k}]) not incident to vertex {vertex:?}",
                    eg.analysis.edge
                ));
            }
        };
        specs.push(FilletSpec {
            edge: edge_now,
            radius,
        });
    }
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

    let plan = match plan_apex_resolution(brep, vertex, &specs, &surfaces, tolerance, None) {
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
