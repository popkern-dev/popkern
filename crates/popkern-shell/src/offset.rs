//! Offset solid operations.
//!
//! Offsets all faces of a solid by a uniform or per-face distance,
//! producing a new solid that is larger (positive offset) or smaller
//! (negative offset) than the original.

use crate::ShellWarning;
use crate::classify::classify_faces;
use crate::error::ShellError;
use crate::inner_skin::{apply_offset_geometry, build_inner_skin};
use crate::surface_offset::compute_offset_surfaces;
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Honest per-run stage count of [`offset_impl`] for progress totals
/// (PK-276).
pub(crate) const OFFSET_STAGE_COUNT: usize = 8;

/// Offset all faces of a solid by a uniform distance.
///
/// A positive `distance` expands the solid outward; a negative distance
/// shrinks it inward.  No faces are removed and no side walls are built.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces or all faces fail.
pub fn offset_solid(
    brep: &BRep,
    distance: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    offset_solid_impl(
        brep,
        distance,
        &[],
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )
    .map_err(|e| {
        let direction = if distance >= 0.0 {
            crate::types::ShellDirection::Outward
        } else {
            crate::types::ShellDirection::Inward
        };
        crate::max_safe::upgrade_trim_decline(
            brep,
            direction,
            distance.abs(),
            tolerance,
            options.trim_policy,
            e,
        )
    })
}

/// Offset all faces with per-face variable distances.
///
/// Faces not listed in `face_distances` use `default_distance`.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces or all faces fail.
pub fn offset_solid_variable(
    brep: &BRep,
    default_distance: f64,
    face_distances: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    offset_solid_impl(
        brep,
        default_distance,
        face_distances,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Shared offset ladder for [`offset_solid`] and
/// [`offset_solid_variable`] (identical stage sequence; the variable
/// entry adds per-face distances).
///
/// `tracker` threads PK-276 interactive control: a cancellation poll
/// before every stage (`ShellError::Cancelled`, model untouched) and
/// one progress event per stage.
pub(crate) fn offset_solid_impl(
    brep: &BRep,
    default_distance: f64,
    face_distances: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<ShellResult, ShellError> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(ShellError::Cancelled);
            }
        };
    }
    // Directions are MATERIAL directions (classify resolves the
    // component material sense): positive distance grows the solid —
    // every face moves out of the material — and negative shrinks it.
    // (The historical mapping compensated for an inverted flag chain
    // that moved "Inward" skins outward; both lies died together.)
    let direction = if default_distance >= 0.0 {
        ShellDirection::Outward
    } else {
        ShellDirection::Inward
    };
    let abs_default = default_distance.abs();
    let abs_faces: Vec<(FaceHandle, f64)> = face_distances
        .iter()
        .map(|&(fh, d)| (fh, d.abs()))
        .collect();

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let classification = classify_faces(brep, abs_default, &abs_faces, direction, &[])?;
    tracker.finish_stage("offset.classify_faces", t0, brep.faces.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);
    let offset_result =
        compute_offset_surfaces(brep, &classification, options.allow_partial, &stage_tol)?;
    tracker.finish_stage("offset.offset_surfaces", t0, offset_result.surfaces.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let (mut result, skin_map) = build_inner_skin(brep, &[])?;
    tracker.finish_stage("offset.skin", t0, result.faces.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    // Face-collapse resolution (PK-326): features thinner than the
    // offset dissolve before the join plan sees the skin.
    let collapse_warnings = crate::collapse_resolve::resolve_face_collapse(
        brep,
        &mut result,
        &skin_map,
        &offset_result.surfaces,
        &stage_tol,
    )?;
    let join_plan = crate::join_classification::classify_joins(
        &result,
        &skin_map,
        &offset_result.surfaces,
        &stage_tol,
        options.join_style,
        false,
    )?;
    crate::vertex_recompute::recompute_vertices_with_gap_corners(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
        &join_plan.gap_corners,
    )?;
    tracker.finish_stage("offset.vertex_recompute", t0, result.vertices.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let undone = crate::edge_recompute::recompute_edges_with_skips(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
        &join_plan.arc_skip_edges(),
    )?;
    crate::join_classification::intersection_gap_postcheck(&join_plan, &undone)?;
    // Collapse gate (part of the recompute stage): a shrink past the
    // local feature size inverts the skin. Tangential zero-length
    // collapses decline here; INVERTED chords (transversal collision)
    // defer their decline to the post-assembly trim stage (PK-327).
    let trim_routing = crate::collapse_check::check_skin_collapse(
        brep,
        &result,
        &skin_map,
        default_distance,
        &tolerance,
    )?;
    tracker.finish_stage("offset.edge_recompute", t0, result.edges.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    apply_offset_geometry(&mut result, &offset_result.surfaces, &skin_map)?;
    // Join construction (PK-325): rolling-ball strips + sphere corner
    // windows, via Euler operators only.
    let arc = crate::arc_join::construct_arc_joins(
        &mut result,
        &join_plan,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    tracker.finish_stage("offset.apply_offset_geometry", t0, result.faces.len());

    // The clone carries the ORIGINAL pcurves. Most exact offsets are
    // parametrization-preserving (planes, cylinders, spheres, tori) so
    // those stay valid, but the exact cone arm SHIFTS the axial
    // parameter (apex translation) and the inherited pcurves/trim
    // windows go geometrically stale. The recompute's survival rule
    // keeps every still-accurate pcurve and re-authors only the stale
    // ones.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    result
        .recompute_all_pcurves_winding_repaired()
        .map_err(|e| ShellError::Topology(e.to_string()))?;
    tracker.finish_stage("offset.pcurve_recompute", t0, result.faces.len());

    // Orientation flags are DERIVED state (a function of finalized
    // loop winding + surface geometry): reconcile them at the serving
    // boundary. Join faces on charts whose handedness opposes their
    // construction-time material measurement (the cone-rim torus
    // strips) get their honest value here (PK-325).
    result.reconcile_orientation_flags();

    // Global self-intersection trim stage (PK-327): scan the
    // assembled skin for transversal collisions; trim-and-re-close
    // under TrimPolicy::Trim, decline with region reports under
    // Decline. The chord gate's deferred decline lands here when the
    // trim does not apply.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let mut provenance = crate::global_si_trim::TrimProvenance {
        default_distance: abs_default,
        ..Default::default()
    };
    for (outer_fh, skin_fh) in &skin_map.faces {
        if result.faces.contains_key(*skin_fh) {
            let d = classification
                .offset_faces
                .iter()
                .find(|(fh, _, _)| fh == outer_fh)
                .map(|(_, d, _)| d.abs())
                .unwrap_or(abs_default);
            provenance.expected.insert(*skin_fh, d);
        }
    }
    provenance.free_form_offsets = crate::global_si_trim::collect_free_form_pairs(
        brep,
        &offset_result.surfaces,
        |fh| {
            skin_map
                .faces
                .get(&fh)
                .copied()
                .filter(|s| result.faces.contains_key(*s))
        },
        |fh| {
            classification
                .offset_faces
                .iter()
                .find(|(f, _, _)| *f == fh)
                .map(|(_, d, _)| *d)
                .unwrap_or(abs_default)
        },
    );
    let trim_warnings = crate::global_si_trim::run_trim_stage(
        &mut result,
        brep,
        &provenance,
        options.trim_policy,
        trim_routing.deferred,
        trim_routing.inverted_edges.len(),
        &tolerance,
    )
    .map_err(|e| {
        crate::max_safe::upgrade_trim_decline(
            brep,
            direction,
            abs_default,
            tolerance,
            options.trim_policy,
            e,
        )
    })?;
    tracker.finish_stage("offset.global_si_trim", t0, result.faces.len());

    // Serve gate: never hand back a body a validator would reject.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    let mut warnings = offset_result.warnings;
    warnings.extend(trim_warnings);
    warnings.extend(collapse_warnings);
    warnings.extend(join_plan.warnings.iter().cloned());
    crate::arc_join_apex::retire_served_apex_warnings(&mut warnings, &arc.served_apexes);
    for face in &offset_result.auto_removed {
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
#[path = "offset_tests.rs"]
mod tests;
