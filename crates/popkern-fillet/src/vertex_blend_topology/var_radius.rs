//! Phase 6c-A continuation §4.2 — variable-radius apex-region
//! topology rebuild.
//!
//! The uniform-radius rebuild in `super::{detect,plan,execute,collapse,rebuild}`
//! produces a single-Sphere cap face with a 3-edge canonical
//! cap-tangent-arc loop.  Variable-radius corners cannot share
//! one sphere across the whole apex (per-fillet spheres are
//! nested, not coincident), so the canonical form is N
//! per-fillet sphere caps, each with a 3-edge outer loop:
//!
//! 1. Cap-tangent arc on `sphere_k ∩ cylinder_k` (perfect G1
//!    against fillet `k`'s cylinder).
//! 2. Inner-spoke edge to cap `k+1` (cap-vs-cap, G0).
//! 3. Inner-spoke edge to cap `k-1` (cap-vs-cap, G0).
//!
//! Today the variable-radius caps end up with *5*-edge loops
//! whose outer twins are 3 cylinders + 2 caps — i.e., each cap
//! shares boundary with cylinder `k-1`, `k`, and `k+1`, but its
//! sphere is tangent only to cylinder `k`.  The two extra
//! cap-vs-cylinder edges are responsible for the residual G1
//! mismatch on `cube_mixed` and `cube_extreme`.
//!
//! Surgery overview (per pair of adjacent caps `k, k+1`):
//!
//! * Identify the legacy "footprint" vertex shared between caps
//!   k, k+1 and cylinders k, k+1.
//! * Splice the cap-side outer loops so the two extra cap-vs-
//!   cylinder edges (cap_k bordering cylinder_{k+1}, and
//!   cap_{k+1} bordering cylinder_k) become cap-vs-cap edges
//!   between cap_k and cap_{k+1}.
//! * The cylinder faces lose 1 edge from their outer loops; their
//!   outer-loop continuity is patched with a `kemr` move-edge
//!   operation (or equivalent kernel primitive).
//!
//! This module is the planning + execution scaffold for that
//! surgery.  Currently in Phase 6c-A.5 *detection only* — the
//! mutation primitives are stubbed pending careful BRep-state
//! testing.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

#[path = "var_radius_classify.rs"]
mod var_radius_classify;
#[path = "var_radius_geometry.rs"]
mod var_radius_geometry;
#[path = "var_radius_legacy.rs"]
mod var_radius_legacy;
#[path = "var_radius_path_b.rs"]
mod var_radius_path_b;
#[path = "var_radius_path_b_vertex_move.rs"]
mod var_radius_path_b_vertex_move;
#[path = "var_radius_vertex_split.rs"]
mod var_radius_vertex_split;

pub(crate) use var_radius_classify::{TwinKind, classify_cap_loops};
#[cfg(test)]
pub(crate) use var_radius_legacy::dump_neighcyl_connectivity;
use var_radius_path_b::{apply_path_b_plan, plan_path_b};
pub(crate) use var_radius_vertex_split::split_apex_footprint_vertex;

/// Errors from `reroute_neighbour_cylinder_edges`.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RerouteError {
    /// No `OwnCylinder` edge found for some cap — surgery can't
    /// anchor to identify which radius this cap "owns".
    NoOwnCylinder { cap_idx: usize },
    /// Couldn't find the destination cap for a NeighbourCylinder
    /// edge — the twin cylinder doesn't match any other cap's
    /// primary radius.
    NoDestinationCap { cap_idx: usize },
    /// Loop linkage check failed after the move — would corrupt
    /// the BRep. Snapshot-rolled back.
    LoopBroken,
    /// `validate()` failed after the surgery. Snapshot-rolled back.
    ValidationFailed,
    /// Not yet implemented — the design is documented in the
    /// function body but the actual mutation primitives need
    /// careful BRep-state testing.
    NotImplemented,
}

/// Attempt the variable-radius cap-region surgery: re-route each
/// cap face's `NeighbourCylinder` edges into cap-vs-cap edges
/// between the cap and its appropriate neighbour cap.
///
/// For each `(cap_k, NeighbourCylinder_edge `e`)`:
/// 1. Identify the twin cylinder face `cyl_m` and the destination
///    cap `cap_n` whose primary cylinder is `cyl_m`.
/// 2. Locate the cap-side half-edge in `cap_n`'s outer loop that
///    has the SAME endpoints as `e` and a Cylinder twin pointing
///    to `cyl_m` — call this the "mirror" edge `e'`. The two
///    edges `e` and `e'` are geometrically identical (same
///    cylinder, same two vertices) but topologically distinct:
///    `e` shares cap_k with cyl_m, `e'` shares cap_n with cyl_m.
///    The fix: kill `e'` (the redundant duplicate), and re-parent
///    `e`'s twin from cyl_m to cap_n.
///
/// On any sub-step failure, snapshot-rollback restores `brep`.
///
/// Returns the number of edges successfully re-routed.
#[allow(dead_code)]
pub(crate) fn reroute_neighbour_cylinder_edges(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(popkern_core::geometry::Point3, f64)],
) -> Result<usize, RerouteError> {
    let snapshot = brep.clone();
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    // Step 2: iteratively classify + apply.  Each surgery
    // mutates the BRep, so we re-classify after each iteration
    // to find the next NeighbourCylinder edge. If a surgery
    // fails (rolled back via snapshot inside apply_path_b_plan),
    // we mark its target as tried-and-failed and skip in
    // subsequent iterations to avoid infinite loops.
    let mut succeeded = 0usize;
    let trace = std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1");
    let mut tried_failed: std::collections::BTreeSet<HalfEdgeHandle> =
        std::collections::BTreeSet::new();
    for iter in 0..16 {
        // Find one NeighbourCylinder target not previously tried.
        let classification = classify_cap_loops(brep, cap_faces, per_fillet_caps);
        let mut target: Option<HalfEdgeHandle> = None;
        for c in &classification {
            for entry in &c.entries {
                if entry.twin_kind == TwinKind::NeighbourCylinder
                    && !tried_failed.contains(&entry.half_edge)
                {
                    target = Some(entry.half_edge);
                    break;
                }
            }
            if target.is_some() {
                break;
            }
        }
        let Some(he_cap) = target else {
            if trace {
                eprintln!("[reroute] iter {iter}: no more NeighCyl targets");
            }
            break;
        };
        let cap_k = brep.loops[brep.half_edges[he_cap].loop_].face;
        let plan = plan_path_b(brep, he_cap, cap_k, cap_faces, per_fillet_caps);
        let Some(plan) = plan else {
            if trace {
                eprintln!("[reroute] iter {iter}: plan_path_b None for he={he_cap:?}");
            }
            tried_failed.insert(he_cap);
            continue;
        };
        match apply_path_b_plan(brep, &plan) {
            Ok(()) => {
                if trace {
                    eprintln!("[reroute] iter {iter}: applied OK for he={he_cap:?}");
                }
                succeeded += 1;
            }
            Err(e) => {
                if trace {
                    eprintln!("[reroute] iter {iter}: apply failed for he={he_cap:?}: {e:?}");
                }
                tried_failed.insert(he_cap);
            }
        }
    }
    let _ = cap_set;

    // Phase 6c-D Stage 2 fallback: if Path B made zero successful
    // surgeries (typical pyramid case — `plan_path_b` returns None
    // because per-fillet sphere centres don't lie on neighbour
    // cylinder axes), run vertex-snap-only as a generalised
    // surrogate. Snap projects each cap-loop vertex onto its
    // incident sphere∩cylinder / sphere∩plane circle, restoring
    // analytic alignment between cap surface and adjacent fillet
    // cylinder.
    if succeeded == 0 {
        let snap_result =
            super::var_radius_snap::apply_vertex_snap_for_caps(brep, cap_faces, per_fillet_caps);
        if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
            eprintln!(
                "[reroute] Path B → 0 surgeries; vertex-snap fallback moved {} verts, updated {} edge curves",
                snap_result.vertices_moved, snap_result.edges_with_curves_updated
            );
        }
        if snap_result.vertices_moved > 0 {
            // Phase 6c-D Stage 3a (arc-rebuild) lives in
            // `var_radius_arc_rebuild` but is currently disabled —
            // the analytic Ellipse 3-D curve writes correctly, but
            // the sphere-side pcurve refit fails to capture the
            // great-circle UV trajectory on a TILTED-axis cylinder
            // (sphere u/v are world-aligned, cylinder axis isn't).
            // The existing `set_cut_edge_geometry` Linear sphere
            // pcurve relies on axis-aligned cylinders (cube body-
            // diagonal); pyramid's tilted axes don't satisfy that.
            // Closing this needs either a non-Linear sphere pcurve
            // (degree-3 NURBS through dense samples, with proper
            // endpoint snapping for pole noise — see Phase 6c-A's
            // `pole_radial_eps` machinery) OR a tilted-axis-aware
            // sphere reparameterisation. Tracked as Stage 3b.
            let _arc_rebuild_disabled =
                super::var_radius_arc_rebuild::rebuild_cap_tangent_arc_curves;
            if brep.validate().is_err() {
                *brep = snapshot;
                return Err(RerouteError::ValidationFailed);
            }
            return Ok(snap_result.vertices_moved);
        }
    }

    // Final validate; if anything went wrong, roll back.
    if succeeded > 0 && brep.validate().is_err() {
        *brep = snapshot;
        return Err(RerouteError::ValidationFailed);
    }

    Ok(succeeded)
}

#[cfg(test)]
#[path = "var_radius_tests_basic.rs"]
mod tests_basic;
#[cfg(test)]
#[path = "var_radius_tests_classify.rs"]
mod tests_classify;
#[cfg(test)]
#[path = "var_radius_tests_diagnostics.rs"]
mod tests_diagnostics;
#[cfg(test)]
#[path = "var_radius_tests_g1.rs"]
mod tests_g1;
#[cfg(test)]
#[path = "var_radius_tests_reroute.rs"]
mod tests_reroute;
#[cfg(test)]
#[path = "var_radius_tests_support.rs"]
mod var_radius_tests_support;
