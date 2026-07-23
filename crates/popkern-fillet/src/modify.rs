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
    // PK-143 stale-cap guard: an in-place update on a strip whose
    // neighbouring cap surface already disagrees with its boundary
    // would compound the lie silently — decline toward the record
    // re-dispatch instead.
    let neighbours = crate::identify::blend_neighbours(brep, fillet_face);
    let stale = crate::stale_cap::detect_stale_caps(brep, &tolerance);
    if let Some(s) = stale.iter().find(|s| neighbours.contains(&s.face)) {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".into(),
            value: format!("{fillet_face:?}"),
            reason: format!(
                "neighbouring cap {:?} is STALE (boundary deviates {:.3e} from its surface) —                  modify via modify_fillet_radius_via_record, which re-dispatches the corner",
                s.face, s.max_deviation
            ),
        });
    }
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
        FilletRadius::ByPoints(points) => {
            // Round 5m: ByPoints in-place. The corners are
            // positioned by the radii at parameters 0 and 1, but
            // the fillet surface is built as a NURBS sweep driven
            // by a piecewise-linear closure over the full points
            // list — intermediate radii shape the surface.
            if let Ok(out) =
                crate::modify_inplace::modify_open_in_place_bypoints(brep, fillet_face, points)
            {
                return Ok(out);
            }
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
/// update first (Round 5d-2 for constant; Round 5p extended to
/// Variable/ByPoints by averaging) and falls back to remove+
/// recreate on failure.
fn modify_closed_fillet_radius(
    brep: &BRep,
    fillet_face: FaceHandle,
    new_radius: FilletRadius,
    tolerance: Tolerance,
) -> Result<BRep, BlendError> {
    // Round 5s: `ByPoints` closed-fillet modify takes the true
    // around-the-rim NURBS-sweep fast path, interpreting each
    // `(t, r)` pair's `t` as the normalized angular position
    // around the cylinder axis (t = 0 at the seam). The resulting
    // fillet surface is a tensor-product NURBS with the local
    // minor radius following the piecewise-linear profile, and
    // both contact loops become closed NURBS curves with
    // actually-varying effective radii. See
    // `modify_inplace::modify_closed_in_place_bypoints`.
    //
    // `Variable { start, end }` takes the same sweep path by
    // synthesizing a symmetric 3-point profile:
    // `[(0.0, start), (0.5, end), (1.0, start)]`. The canonical
    // interpretation for a closed rim is *"`start` at the seam
    // (θ = 0), linearly interpolated to `end` at the opposite
    // side of the rim (θ = π), then linearly back to `start`"*.
    // The profile is continuous at the seam (t = 0 and t = 1
    // both map to `start`), which is required for the closed
    // NURBS-sweep builder to produce a well-formed surface.
    if let FilletRadius::ByPoints(pts) = &new_radius
        && !pts.is_empty()
        && let Ok(out) =
            crate::modify_inplace::modify_closed_in_place_bypoints(brep, fillet_face, pts)
    {
        return Ok(out);
    }
    if let FilletRadius::Variable { start, end } = &new_radius {
        let synthetic: Vec<(f64, f64)> = vec![(0.0, *start), (0.5, *end), (1.0, *start)];
        if let Ok(out) =
            crate::modify_inplace::modify_closed_in_place_bypoints(brep, fillet_face, &synthetic)
        {
            return Ok(out);
        }
    }

    // Fast path: averaged in-place geometry update (Constant
    // only; Variable + ByPoints have already tried the true
    // sweep path above, and this averaging branch is reached
    // only if the sweep path failed for a structural reason).
    let averaged = match &new_radius {
        FilletRadius::Constant(r) => *r,
        FilletRadius::Variable { start, end } => (*start + *end) * 0.5,
        FilletRadius::ByPoints(pts) if !pts.is_empty() => {
            pts.iter().map(|(_, r)| *r).sum::<f64>() / pts.len() as f64
        }
        FilletRadius::ByPoints(_) => 0.0,
    };
    if averaged > 0.0
        && let Ok(out) = crate::modify_inplace::modify_closed_in_place(brep, fillet_face, averaged)
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

/// PK-143: record-based radius modify — the cap-family re-dispatch
/// path. Strips whose ends meet a corner cap (sphere / Gregory /
/// setback families) cannot be modified by the local in-place or
/// unfillet-plan paths: the cap end has no restored-edge endpoint
/// ("found 1 valid restored-edge endpoints; need exactly 2"), and a
/// corner's radii are COUPLED — changing one strip re-dispatches the
/// whole corner construction. With the operation record the exact
/// commercial semantics are a deterministic REPLAY: the target
/// strip's spec is rewritten to the new radius and every spec is
/// re-applied to the recorded input, so the corner (whatever cap
/// family it dispatches to at the new radius) is rebuilt by the
/// same machinery that built it.
///
/// The target spec is identified by proximity: the spec whose
/// recorded edge midpoint lies nearest the strip face's boundary
/// (within 2× the larger radius). Returns the full [`BlendResult`]
/// of the replay.
pub fn modify_fillet_radius_via_record(
    record: &crate::history::FilletOperationRecord,
    fillet_face: FaceHandle,
    new_radius: impl Into<FilletRadius>,
) -> Result<crate::BlendResult, BlendError> {
    use popkern_core::geometry::curve::CurveEval;
    let new_radius = new_radius.into();
    // Strip boundary sample points in the RESULT body.
    let face =
        record
            .result
            .faces
            .get(fillet_face)
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "fillet_face".into(),
                value: format!("{fillet_face:?}"),
                reason: "face not present in the recorded result".into(),
            })?;
    let mut samples = Vec::new();
    if let Some(l) = record.result.loops.get(face.outer_loop) {
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..record.result.half_edges.len() + 1 {
            let Some(he) = record.result.half_edges.get(cur) else {
                break;
            };
            if let Some(v) = record.result.vertices.get(he.origin) {
                samples.push(v.position);
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
    if samples.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".into(),
            value: format!("{fillet_face:?}"),
            reason: "strip face has no boundary vertices".into(),
        });
    }
    // Nearest recorded spec by edge-midpoint distance to the strip.
    let mut best: Option<(usize, f64)> = None;
    for (i, spec) in record.specs.iter().enumerate() {
        let Some(edge) = record.input.edges.get(spec.edge) else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        // Sample the spec edge along its span (the strip's stored
        // boundary vertices sit at the ENDS, up to half the edge
        // length from the midpoint).
        let mut d = f64::INFINITY;
        for k in 0..=8 {
            let t = curve.t_start + (curve.t_end - curve.t_start) * f64::from(k) / 8.0;
            let q = curve.geometry.point_at(t);
            for p in &samples {
                d = d.min(p.distance_to(&q));
            }
        }
        if best.map(|(_, bd)| d < bd).unwrap_or(true) {
            best = Some((i, d));
        }
    }
    let spec_radius_scale = |r: &FilletRadius| -> f64 {
        match r {
            FilletRadius::Constant(c) => *c,
            FilletRadius::Variable { start, end } => start.max(*end),
            _ => 1.0,
        }
    };
    let (idx, d) = best.ok_or_else(|| BlendError::InvalidParameter {
        parameter: "record".into(),
        value: "specs".into(),
        reason: "no recorded spec has a live edge curve".into(),
    })?;
    let gate = 2.0 * spec_radius_scale(&record.specs[idx].radius).max(1e-6);
    if d > gate {
        return Err(BlendError::InvalidParameter {
            parameter: "fillet_face".into(),
            value: format!("{fillet_face:?}"),
            reason: format!(
                "no recorded spec edge lies within {gate:.3e} of the strip (nearest {d:.3e}) — \
                 the face does not belong to this record"
            ),
        });
    }
    let mut specs = record.specs.clone();
    specs[idx].radius = new_radius;
    crate::fillet_edges(
        &record.input,
        &specs,
        &record.options,
        record.tolerance,
        None,
    )
}
