//! In-place fillet radius modification (§14.7).
//!
//! [`modify_fillet_radius`] is the no-record path for changing a
//! fillet's radius after the fact — given just the BRep and the
//! fillet face's handle, it composes algorithmic unfillet (Round
//! 5b/5c) with a fresh fillet at the new radius. The record-based
//! path lives on [`crate::FilletOperationRecord::with_uniform_radius`]
//! and should be preferred when the caller has the original record
//! (it skips the unfillet step).
//!
//! # §14.7 acceptance criteria
//!
//! - **Identical to remove+recreate**: yes — that's literally what
//!   the implementation does (algorithmic unfillet → fillet at new
//!   radius). Any topology/geometry difference would be a bug in
//!   one of the underlying primitives.
//! - **≥2x faster than remove+recreate**: not yet. Today the
//!   function pays the full cost of unfillet + fillet. A future
//!   round will add a true in-place geometry update that only
//!   touches the fillet surface, contact edges, and corner
//!   vertices without rebuilding topology.
//!
//! # When to use which path
//!
//! | Caller has...               | Use                                              |
//! |-----------------------------|--------------------------------------------------|
//! | `FilletOperationRecord`     | [`crate::FilletOperationRecord::with_uniform_radius`] |
//! | only the BRep + face handle | [`modify_fillet_radius`]                         |

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::error::BlendError;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::unfillet::{compute_unfillet_plan, unfillet_face};

/// Modify the radius of a fillet face by composing algorithmic
/// unfillet with a fresh fillet at `new_radius`.
///
/// This is the "no record" §14.7 path — when the caller has only a
/// BRep and a fillet face handle (no [`crate::FilletOperationRecord`]),
/// they can still change the fillet's radius in one call.
///
/// `new_radius` accepts any value convertible into [`FilletRadius`]:
///   - A bare `f64` becomes [`FilletRadius::Constant`] (the common case).
///   - [`FilletRadius::Variable { start, end }`] applies a linearly-
///     varying radius along the edge.
///   - [`FilletRadius::ByPoints`] applies a per-parameter radius profile.
///
/// The constant case takes the in-place fast path (Round 5d-2/5d-3).
/// Variable and `ByPoints` cases go through algorithmic unfillet
/// followed by a fresh fillet at the new spec — slower but
/// supported transparently.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if `fillet_face` isn't
///   recognised as a fillet, or if the restored edge can't be
///   located in the unfilleted BRep.
/// - Any error from the underlying [`unfillet_face`] /
///   [`crate::fillet_edges`] calls.
pub fn modify_fillet_radius(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: impl Into<FilletRadius>,
    tolerance: Tolerance,
) -> Result<BRep, BlendError> {
    let radius = new_radius.into();
    if is_closed_annulus(brep, fillet_face) {
        return modify_closed_fillet_radius(brep, fillet_face, radius, tolerance);
    }
    modify_open_fillet_radius(brep, fillet_face, radius, tolerance)
}

/// Open (plane-plane) fillet path — tries the in-place geometry
/// update first (Round 5d-3 for constant, Round 5h for variable).
/// Falls back to remove+recreate via `compute_unfillet_plan` +
/// `unfillet_face` + `fillet_edges` when the in-place path can't
/// handle the spec (e.g. `ByPoints`) or fails for any reason.
fn modify_open_fillet_radius(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: FilletRadius,
    tolerance: Tolerance,
) -> Result<BRep, BlendError> {
    // Fast path: in-place geometry update.
    match &new_radius {
        FilletRadius::Constant(r) => {
            if let Ok(out) = crate::modify_inplace::modify_open_in_place(brep, fillet_face, *r) {
                return Ok(out);
            }
        }
        FilletRadius::Variable { start, end } => {
            if let Ok(out) = crate::modify_inplace::modify_open_in_place_variable(
                brep,
                fillet_face,
                *start,
                *end,
            ) {
                return Ok(out);
            }
        }
        FilletRadius::ByPoints(_) => {
            // ByPoints carries arbitrary parameter→radius pairs; for
            // a 4-corner plane-plane fillet only the start (parameter
            // 0) and end (parameter 1) values matter. Fall through
            // to remove+recreate which honours the full profile.
        }
    }

    // Fallback: remove + recreate.
    let plan = compute_unfillet_plan(brep, fillet_face)?;
    let p_start = plan.restored_edge_start;
    let p_end = plan.restored_edge_end;

    let unfilleted = unfillet_face(brep, fillet_face)?;

    let restored_edge = find_edge_with_endpoints(&unfilleted, &p_start, &p_end, &tolerance)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: format!("{fillet_face:?}"),
            reason: "could not locate restored sharp edge in the unfilleted BRep".to_string(),
        })?;

    let specs = [FilletSpec {
        edge: restored_edge,
        radius: new_radius,
    }];
    let result = crate::fillet_edges(
        &unfilleted,
        &specs,
        &FilletOptions::default(),
        tolerance,
        None,
    )?;
    Ok(result.brep)
}

/// Closed-fillet (annulus) path — tries the in-place geometry
/// update first (Round 5d-2 — true §14.7 ≥2× speedup) when the
/// radius is constant, and falls back to remove+recreate for any
/// variable-radius spec or for constant cases that don't fit the
/// fast path.
fn modify_closed_fillet_radius(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: FilletRadius,
    tolerance: Tolerance,
) -> Result<BRep, BlendError> {
    // Fast path: in-place geometry update (constant only).
    if let FilletRadius::Constant(r) = &new_radius
        && let Ok(out) = crate::modify_inplace::modify_closed_in_place(brep, fillet_face, *r)
    {
        return Ok(out);
    }

    // Fallback: remove + recreate via algorithmic unfillet.
    let unfilleted = crate::unfillet_closed::unfillet_closed_annulus(brep, fillet_face)?;
    let restored_edge = find_closed_plane_cylinder_edge(&unfilleted).ok_or_else(|| {
        BlendError::InvalidParameter {
            parameter: "fillet_face".to_string(),
            value: format!("{fillet_face:?}"),
            reason: "could not locate restored closed rim edge after unfillet".to_string(),
        }
    })?;
    let specs = [FilletSpec {
        edge: restored_edge,
        radius: new_radius,
    }];
    let result = crate::fillet_edges(
        &unfilleted,
        &specs,
        &FilletOptions::default(),
        tolerance,
        None,
    )?;
    Ok(result.brep)
}

/// True when `face` is a closed-fillet annulus (one outer + one
/// inner loop, both single-element self-closed). Mirrors the
/// detector in [`crate::unfillet`] but inlined here so this module
/// has no cross-module visibility dependency.
fn is_closed_annulus(brep: &BRep, face: FaceHandle) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    if face_data.inner_loops.len() != 1 {
        return false;
    }
    let single = |lh: popkern_core::topology::LoopHandle| -> bool {
        let Some(l) = brep.loops.get(lh) else {
            return false;
        };
        let start = l.half_edge;
        let Some(he) = brep.half_edges.get(start) else {
            return false;
        };
        he.next == start && he.prev == start
    };
    single(face_data.outer_loop) && single(face_data.inner_loops[0])
}

/// Find an edge whose two endpoint vertices are at `p_start` and
/// `p_end` (in either order) within `tolerance.linear`.
fn find_edge_with_endpoints(
    brep: &BRep,
    p_start: &Point3,
    p_end: &Point3,
    tolerance: &Tolerance,
) -> Option<EdgeHandle> {
    let tol = tolerance.linear().max(1e-9);
    let close = |a: &Point3, b: &Point3| -> bool {
        let dx = a.x - b.x;
        let dy = a.y - b.y;
        let dz = a.z - b.z;
        (dx * dx + dy * dy + dz * dz).sqrt() < tol
    };
    for (eh, edata) in brep.edges.iter() {
        let he_h = edata.half_edge;
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let v0_pos = brep.vertices[he.origin].position;
        let v1_pos = brep.vertices[twin.origin].position;
        if (close(&v0_pos, p_start) && close(&v1_pos, p_end))
            || (close(&v0_pos, p_end) && close(&v1_pos, p_start))
        {
            return Some(eh);
        }
    }
    None
}

/// Find a closed (v_start == v_end) edge whose two adjacent faces
/// are a Plane and a Cylinder — the topological signature of the
/// rim edge restored by [`crate::unfillet_closed::unfillet_closed_annulus`].
fn find_closed_plane_cylinder_edge(brep: &BRep) -> Option<EdgeHandle> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    for (eh, edata) in brep.edges.iter() {
        let he_h = edata.half_edge;
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if he.origin != twin.origin {
            continue;
        }
        let face_a = brep.loops.get(he.loop_)?.face;
        let face_b = brep.loops.get(twin.loop_)?.face;
        let surf_a = brep.faces.get(face_a)?.surface.as_ref()?;
        let surf_b = brep.faces.get(face_b)?.surface.as_ref()?;
        let is_plane_cyl = matches!(
            (&surf_a.geometry, &surf_b.geometry),
            (
                SurfaceGeometry::Plane { .. },
                SurfaceGeometry::Cylinder { .. },
            ) | (
                SurfaceGeometry::Cylinder { .. },
                SurfaceGeometry::Plane { .. },
            )
        );
        if is_plane_cyl {
            return Some(eh);
        }
    }
    None
}

#[cfg(test)]
#[path = "modify_tests.rs"]
mod tests;
