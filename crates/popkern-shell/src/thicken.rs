//! Thicken operation: convert a sheet body into a solid.
//!
//! Takes a sheet body (a BRep with open boundary edges) and thickens it
//! by offsetting in both directions and connecting the free boundaries
//! with side wall faces.

use crate::ShellWarning;
use crate::assemble::assemble_shell;
use crate::classify::classify_faces;
use crate::edge_recompute::recompute_edges;
use crate::error::ShellError;
use crate::inner_skin::{apply_offset_geometry, build_inner_skin};
use crate::skin_map::SkinMap;
use crate::surface_offset::compute_offset_surfaces;
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use crate::vertex_recompute::recompute_vertices;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Thicken a closed body into a hollow wall solid.
///
/// Offsets the input by `thickness / 2` in both material directions and
/// assembles the two skins into a two-shell hollow body.
///
/// Open sheets (free boundary half-edges) DECLINE with
/// [`ShellError::OpenBoundary`]: closing a free rim needs the side-wall
/// join machinery, which is not built yet (temporary containment — the
/// historical doc claimed rim walls were generated; they never were).
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, has free boundary
/// half-edges, the wall does not fit (`OffsetTooLarge`), or the
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
/// (PK-276): one skin ladder per side plus assembly.
pub(crate) const THICKEN_STAGE_COUNT: usize = 2 * SKIN_STAGE_COUNT + 1;
const SKIN_STAGE_COUNT: usize = 5;

/// One directional skin ladder (classify → offset surfaces → skin →
/// vertex/edge recompute → apply), shared by the outward and inward
/// halves of [`thicken`]. `label_*` are the per-side stage labels.
#[allow(clippy::too_many_arguments)]
fn thicken_skin(
    brep: &BRep,
    half: f64,
    direction: ShellDirection,
    allow_partial: bool,
    stage_tol: &Tolerance,
    labels: [&'static str; SKIN_STAGE_COUNT],
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<(BRep, SkinMap, crate::surface_offset::OffsetResult), ShellError> {
    use std::time::Instant;
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(ShellError::Cancelled);
            }
        };
    }
    stage_gate!();
    let t0 = Instant::now();
    let classification = classify_faces(brep, half, &[], direction, &[])?;
    tracker.finish_stage(labels[0], t0, brep.faces.len());

    stage_gate!();
    let t0 = Instant::now();
    let offset = compute_offset_surfaces(brep, &classification, allow_partial, stage_tol)?;
    tracker.finish_stage(labels[1], t0, offset.surfaces.len());

    stage_gate!();
    let t0 = Instant::now();
    let (mut skin, skin_map) = build_inner_skin(brep, &[])?;
    tracker.finish_stage(labels[2], t0, skin.faces.len());

    stage_gate!();
    let t0 = Instant::now();
    recompute_vertices(&mut skin, &offset.surfaces, &skin_map, brep, stage_tol)?;
    recompute_edges(&mut skin, &offset.surfaces, &skin_map, brep, stage_tol)?;
    // Collapse gate (part of the recompute stage): a half-thickness past
    // the local feature size inverts this skin. The gate runs in
    // half-thickness units; the decline reports the caller's full
    // thickness so a GUI can bound its slider directly.
    crate::collapse_check::check_skin_collapse(brep, &skin, &skin_map, half, stage_tol).map_err(
        |e| match e {
            ShellError::OffsetTooLarge {
                requested,
                max_safe,
                binding,
            } => ShellError::OffsetTooLarge {
                requested: requested * 2.0,
                max_safe: max_safe * 2.0,
                binding,
            },
            other => other,
        },
    )?;
    tracker.finish_stage(labels[3], t0, skin.edges.len());

    stage_gate!();
    let t0 = Instant::now();
    apply_offset_geometry(&mut skin, &offset.surfaces, &skin_map)?;
    tracker.finish_stage(labels[4], t0, skin.faces.len());

    Ok((skin, skin_map, offset))
}

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

    // Open-boundary precondition (temporary containment, recorded in
    // docs/architecture/solid-offset-shell-design.md §6): a free rim
    // needs side walls the join stage does not build yet. Without this
    // gate an open sheet used to serve as TWO disconnected skins with an
    // unclosed rim while this function's docs claimed otherwise.
    let free_half_edges = brep
        .half_edges
        .iter()
        .filter(|(_, he)| !brep.loops.contains_key(he.loop_))
        .count();
    if free_half_edges > 0 {
        return Err(ShellError::OpenBoundary { free_half_edges });
    }

    let half = thickness / 2.0;
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);

    // Offset outward to create the outer skin.
    let (outer, outer_map, offset_out) = thicken_skin(
        brep,
        half,
        ShellDirection::Outward,
        options.allow_partial,
        &stage_tol,
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
    let (mut inner, inner_map, offset_in) = thicken_skin(
        brep,
        half,
        ShellDirection::Inward,
        options.allow_partial,
        &stage_tol,
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
    let combined_skin_map = build_combined_skin_map(&outer_map, &inner_map);

    // Assemble: outer + inner, no removed faces (side walls at free edges).
    if tracker.cancelled() {
        return Err(ShellError::Cancelled);
    }
    let t0 = Instant::now();
    let result = assemble_shell(&outer, &inner, &[], &combined_skin_map, tolerance)?;
    tracker.finish_stage("thicken.assemble", t0, result.faces.len());

    // Serve gate: never hand back a body a validator would reject.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    let mut warnings = offset_out.warnings;
    warnings.extend(offset_in.warnings);
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

/// Build a combined skin map mapping outer skin handles to inner skin handles.
///
/// Both `outer_map` and `inner_map` map from the original BRep to their
/// respective skins. We need a map from outer skin to inner skin.
fn build_combined_skin_map(outer_map: &SkinMap, inner_map: &SkinMap) -> SkinMap {
    let reverse_outer_vertices = outer_map.reverse_vertex_map();
    let reverse_outer_edges = outer_map.reverse_edge_map();
    let reverse_outer_he = outer_map
        .half_edges
        .iter()
        .map(|(&k, &v)| (v, k))
        .collect::<std::collections::HashMap<_, _>>();
    let reverse_outer_faces = outer_map.reverse_face_map();

    let mut vertices = std::collections::HashMap::new();
    let mut edges = std::collections::HashMap::new();
    let mut half_edges = std::collections::HashMap::new();
    let mut faces = std::collections::HashMap::new();

    for (&outer_skin_vh, &orig_vh) in &reverse_outer_vertices {
        if let Some(&inner_skin_vh) = inner_map.vertices.get(&orig_vh) {
            vertices.insert(outer_skin_vh, inner_skin_vh);
        }
    }
    for (&outer_skin_eh, &orig_eh) in &reverse_outer_edges {
        if let Some(&inner_skin_eh) = inner_map.edges.get(&orig_eh) {
            edges.insert(outer_skin_eh, inner_skin_eh);
        }
    }
    for (&outer_skin_heh, &orig_heh) in &reverse_outer_he {
        if let Some(&inner_skin_heh) = inner_map.half_edges.get(&orig_heh) {
            half_edges.insert(outer_skin_heh, inner_skin_heh);
        }
    }
    for (&outer_skin_fh, &orig_fh) in &reverse_outer_faces {
        if let Some(&inner_skin_fh) = inner_map.faces.get(&orig_fh) {
            faces.insert(outer_skin_fh, inner_skin_fh);
        }
    }

    SkinMap {
        vertices,
        edges,
        half_edges,
        faces,
    }
}

#[cfg(test)]
#[path = "thicken_tests.rs"]
mod tests;
