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
pub(crate) const THICKEN_STAGE_COUNT: usize = 2 * crate::skin_ladder::SKIN_LADDER_STAGE_COUNT + 1;

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

    // Offset outward to create the outer skin. The ladder's collapse
    // gate reports the caller's full-thickness units (scale 2.0).
    let (outer, outer_map, offset_out) = crate::skin_ladder::offset_skin_ladder(
        brep,
        half,
        &[],
        ShellDirection::Outward,
        &[],
        false,
        options.allow_partial,
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
    let (mut inner, inner_map, offset_in) = crate::skin_ladder::offset_skin_ladder(
        brep,
        half,
        &[],
        ShellDirection::Inward,
        &[],
        false,
        options.allow_partial,
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

#[cfg(test)]
#[path = "thicken_tests.rs"]
mod tests;
