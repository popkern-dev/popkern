//! Thicken operation: convert a sheet body into a solid.
//!
//! Takes a sheet body (a BRep with open boundary edges) and thickens it
//! by offsetting in both directions and connecting the free boundaries
//! with side wall faces.

use crate::ShellWarning;
use crate::assemble::assemble_shell;
use crate::error::ShellError;
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Thicken a body into a solid wall body.
///
/// Offsets the input by `thickness / 2` in both material directions and
/// assembles the two skins: a CLOSED body becomes a two-shell hollow
/// body; an OPEN sheet (free boundary half-edges) becomes a single
/// closed solid, its free rims closed with wall faces adopting the rim
/// topology (design §6 open-sheet row, PK-324). Non-planar free rims
/// decline with the edge named — ruled free-rim walls are Stage-3
/// work (PK-325).
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, the wall does not
/// fit (`OffsetTooLarge`), a free rim needs a non-planar wall, or the
/// assembled result fails structural validation.
pub fn thicken(
    brep: &BRep,
    thickness: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    thicken_impl(
        brep,
        thickness,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Honest per-run stage count of [`thicken_impl`] for progress totals
/// (PK-276): one skin ladder per side, assembly, pcurve recompute.
pub(crate) const THICKEN_STAGE_COUNT: usize = 2 * crate::skin_ladder::SKIN_LADDER_STAGE_COUNT + 2;

/// Shared implementation for [`thicken`], with PK-276 interactive
/// control threaded through both directional skin ladders and the
/// assembly stage.
pub(crate) fn thicken_impl(
    brep: &BRep,
    thickness: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<ShellResult, ShellError> {
    use std::time::Instant;

    // Normalize the open-boundary representation: constructors register
    // sheet rims as null-sentinel twins (`boundary_edges`), while the
    // skin pipeline works in the loopless-twin form (both half-edges
    // present; the free one carries a dangling loop). Materialize the
    // missing twins on a working copy so one rim representation flows
    // through classify → recompute → assembly → walls.
    let normalized;
    let brep = if brep
        .half_edges
        .iter()
        .any(|(_, he)| !brep.half_edges.contains_key(he.twin))
    {
        let mut owned = brep.clone();
        materialize_boundary_twins(&mut owned)?;
        normalized = owned;
        &normalized
    } else {
        brep
    };

    let half = thickness / 2.0;
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);

    // Offset outward to create the outer skin. The ladder's collapse
    // gate reports the caller's full-thickness units (scale 2.0).
    let (outer, outer_map, offset_out, joins_out) = crate::skin_ladder::offset_skin_ladder(
        brep,
        half,
        &[],
        ShellDirection::Outward,
        &[],
        false,
        options.allow_partial,
        options.join_style,
        &stage_tol,
        2.0,
        [
            "thicken.outward.classify_faces",
            "thicken.outward.offset_surfaces",
            "thicken.outward.skin",
            "thicken.outward.recompute",
            "thicken.outward.apply_offset_geometry",
        ],
        tracker,
    )?;

    // Offset inward to create the inner skin.
    let (mut inner, inner_map, offset_in, joins_in) = crate::skin_ladder::offset_skin_ladder(
        brep,
        half,
        &[],
        ShellDirection::Inward,
        &[],
        false,
        options.allow_partial,
        options.join_style,
        &stage_tol,
        2.0,
        [
            "thicken.inward.classify_faces",
            "thicken.inward.offset_surfaces",
            "thicken.inward.skin",
            "thicken.inward.recompute",
            "thicken.inward.apply_offset_geometry",
        ],
        tracker,
    )?;

    // Mirror the inner skin's material side honestly: winding and
    // orientation flags flip together (a flag-only flip ships the
    // PK-244 winding-dishonest mirror class).
    inner.reverse_body_orientation();

    // Build a combined skin map from outer to inner.
    let combined_skin_map = crate::skin_ladder::combine_skin_maps(&outer_map, &inner_map);

    // Assemble: outer + inner; free rims (open sheets) are closed with
    // wall faces inside assembly.
    if tracker.cancelled() {
        return Err(ShellError::Cancelled);
    }
    let t0 = Instant::now();
    let mut result = assemble_shell(&outer, &inner, &[], &combined_skin_map, tolerance)?;
    tracker.finish_stage("thicken.assemble", t0, result.faces.len());

    // Assembly copies topology without pcurves — re-author them all
    // (trim domains, loop winding, downstream tessellation all read
    // them; serving a pcurve-less body meshed every periodic face as
    // an untrimmed grid).
    if tracker.cancelled() {
        return Err(ShellError::Cancelled);
    }
    let t0 = Instant::now();
    result
        .recompute_all_pcurves_winding_repaired()
        .map_err(|e| ShellError::Topology(e.to_string()))?;
    tracker.finish_stage("thicken.pcurve_recompute", t0, result.faces.len());

    // Orientation flags are DERIVED state (a function of finalized
    // loop winding + surface geometry): reconcile them at the serving
    // boundary. Join faces on charts whose handedness opposes their
    // construction-time material measurement (the cone-rim torus
    // strips) get their honest value here (PK-325).
    result.reconcile_orientation_flags();

    // Serve gate: never hand back a body a validator would reject.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    let mut warnings = offset_out.warnings;
    warnings.extend(offset_in.warnings);
    warnings.extend(joins_out.warnings);
    warnings.extend(joins_in.warnings);
    for face in offset_out
        .auto_removed
        .iter()
        .chain(&offset_in.auto_removed)
    {
        warnings.push(ShellWarning::FaceAutoRemoved {
            face: *face,
            reason: "offset failed".into(),
        });
    }

    Ok(ShellResult {
        brep: result,
        warnings,
    })
}

/// Materialize the missing twin half-edges of null-sentinel boundary
/// rims (the constructor representation of an open sheet) into the
/// loopless-twin form: each rim gains a real half-edge running against
/// its survivor, with `next`/`prev` chained around the hole ring so
/// vertex fan walks stay closed. The dangling `loop_` marks it free.
fn materialize_boundary_twins(brep: &mut BRep) -> Result<(), ShellError> {
    use popkern_core::topology::{HalfEdge, HalfEdgeHandle, VertexHandle};
    use std::collections::HashMap;

    let survivors: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| !brep.half_edges.contains_key(he.twin))
        .map(|(h, _)| h)
        .collect();

    // Create the twins first (next/prev wired below).
    let mut new_by_origin: HashMap<VertexHandle, HalfEdgeHandle> = HashMap::new();
    let mut created: Vec<HalfEdgeHandle> = Vec::new();
    for s in &survivors {
        let next = brep.half_edges[*s].next;
        let dest = brep
            .half_edges
            .get(next)
            .map(|n| n.origin)
            .ok_or_else(|| ShellError::Topology("boundary survivor has no next".into()))?;
        let edge = brep.half_edges[*s].edge;
        let new = brep.half_edges.insert(HalfEdge {
            origin: dest,
            twin: *s,
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge,
            loop_: Default::default(),
            pcurve: None,
        });
        brep.half_edges[*s].twin = new;
        if new_by_origin.insert(dest, new).is_some() {
            return Err(ShellError::Topology(format!(
                "non-manifold sheet rim at vertex {dest:?} (two outgoing boundary half-edges)"
            )));
        }
        created.push(new);
    }

    // Chain each hole ring: next(b) = the boundary half-edge starting
    // where b ends (its survivor's origin).
    for b in created {
        let survivor = brep.half_edges[b].twin;
        let end_vertex = brep.half_edges[survivor].origin;
        let next_b = *new_by_origin.get(&end_vertex).ok_or_else(|| {
            ShellError::Topology(format!(
                "sheet rim ring broken at vertex {end_vertex:?} — no continuing boundary edge"
            ))
        })?;
        brep.half_edges[b].next = next_b;
        brep.half_edges[next_b].prev = b;
    }

    Ok(())
}

#[cfg(test)]
#[path = "thicken_tests.rs"]
mod tests;
