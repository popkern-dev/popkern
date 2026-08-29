#![warn(missing_docs)]
//! popkern-shell: Shell and offset solid operations for the popkern BRep kernel.
//!
//! Provides the [`shell`] and [`shell_variable`] functions for creating
//! thin-walled solids, plus [`offset_solid`] and [`thicken`] operations.

mod arc_join;
mod arc_join_apex;
mod arc_join_corner;
mod arc_join_corner_mixed;
mod arc_join_geometry;
mod arc_join_transversal;
mod arc_join_verify;
pub mod assemble;
mod body_distance;
pub mod classify;
mod classify_sense;
pub mod collapse_check;
mod collapse_resolve;
mod curved_recompute;
mod curved_recompute_similarity;
mod edge_curve_dispatch;
mod edge_reauthor;
pub mod edge_recompute;
pub mod error;
mod global_si_trim;
pub mod inner_skin;
pub mod interactive;
mod join_classification;
pub mod max_safe;
mod max_safe_motion;
mod midplane;
pub mod offset;
mod parallel_dispatch;
pub mod preview;
mod shell_naming;
pub mod side_walls;
mod skin_ladder;
pub mod skin_map;
mod spindle_pinch;
pub mod surface_offset;
pub mod thicken;
pub mod types;
pub mod vertex_recompute;

pub use error::{ShellError, ShellWarning};
pub use interactive::{
    offset_solid_variable_with_control, offset_solid_with_control, shell_variable_with_control,
    shell_with_control, thicken_with_control,
};
pub use max_safe::{OffsetConstraint, SafeOffsetReport, max_safe_offset};
pub use preview::{PreviewMesh, shell_preview};
pub use types::*;

pub(crate) use shell_naming::{ShellProvenance, assign_shell_names, join_source_name};

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Create a thin-walled solid by offsetting faces of a BRep body.
///
/// Removes the specified faces and offsets the remaining faces inward (or
/// outward, depending on `direction`) by `thickness`, producing a hollow
/// shell solid.  When faces are removed, side wall faces are built to close
/// the gap between the outer and inner skins.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, all faces fail to
/// offset, or the offset distance is too large.
#[allow(clippy::too_many_arguments)]
pub fn shell(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<ShellResult, ShellError> {
    let mut provenance = ShellProvenance::default();
    let mut result = shell_impl(
        brep,
        thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
        tracking.is_some().then_some(&mut provenance),
    )
    .map_err(|e| {
        max_safe::upgrade_trim_decline(
            brep,
            direction,
            thickness,
            tolerance,
            options.trim_policy,
            e,
        )
    })?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, brep, &provenance);
    }
    Ok(result)
}

/// Create a shell with per-face variable thickness.
///
/// Like [`shell`], but the `face_thicknesses` parameter provides
/// per-face thickness overrides. Faces not listed use `default_thickness`.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, all faces fail to
/// offset, or the offset distance is too large.
#[allow(clippy::too_many_arguments)]
pub fn shell_variable(
    brep: &BRep,
    default_thickness: f64,
    face_thicknesses: &[(FaceHandle, f64)],
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<ShellResult, ShellError> {
    let mut provenance = ShellProvenance::default();
    let mut result = shell_impl(
        brep,
        default_thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
        tracking.is_some().then_some(&mut provenance),
    )
    .map_err(|e| {
        max_safe::upgrade_trim_decline(
            brep,
            direction,
            default_thickness,
            tolerance,
            options.trim_policy,
            e,
        )
    })?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, brep, &provenance);
    }
    Ok(result)
}

/// Honest per-run stage count of [`shell_impl`] for progress totals
/// (PK-276).
pub(crate) const SHELL_STAGE_COUNT: usize = 8;

/// Shared implementation for [`shell`] and [`shell_variable`].
///
/// `tracker` threads PK-276 interactive control: a cancellation poll
/// before every stage (`ShellError::Cancelled`, model untouched — the
/// working skins are private) and one progress event per stage.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shell_impl(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
    naming_sink: Option<&mut ShellProvenance>,
) -> Result<ShellResult, ShellError> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(ShellError::Cancelled);
            }
        };
    }
    // Midplane moves BOTH skins; it runs its own two-ladder pipeline
    // (PK-330) instead of the single-skin ladder below.
    if matches!(direction, ShellDirection::Midplane) {
        return midplane::midplane_shell_impl(
            brep,
            thickness,
            removed_faces,
            face_thicknesses,
            options,
            tolerance,
            tracker,
            naming_sink,
        );
    }

    // Stage 1: classify faces.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let classification =
        classify::classify_faces(brep, thickness, face_thicknesses, direction, removed_faces)?;
    tracker.finish_stage("shell.classify_faces", t0, brep.faces.len());

    // Stage 2: compute offset surfaces.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);
    let offset_result = surface_offset::compute_offset_surfaces(
        brep,
        &classification,
        options.allow_partial,
        &stage_tol,
    )?;
    tracker.finish_stage("shell.offset_surfaces", t0, offset_result.surfaces.len());

    // Merge auto-removed faces into the removed list.
    let mut all_removed = classification.removed_faces.clone();
    all_removed.extend_from_slice(&offset_result.auto_removed);

    // Stage 3: build inner skin.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let (mut inner, skin_map) = inner_skin::build_inner_skin(brep, &all_removed)?;
    tracker.finish_stage("shell.inner_skin", t0, inner.faces.len());

    // Stage 4: recompute inner skin geometry. Face-collapse
    // resolution (PK-326) runs FIRST: features thinner than the
    // offset dissolve out of the skin, so the join plan and the
    // recomputes see the surviving topology only.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let collapse_warnings = collapse_resolve::resolve_face_collapse(
        brep,
        &mut inner,
        &skin_map,
        &offset_result.surfaces,
        &stage_tol,
    )?;
    // PK-340: opening-rim gap corners serve on the classic shell path
    // (the connector is the strip × removed-original-surface
    // transversal; the side walls grow arc sectors).
    let join_plan = join_classification::classify_joins(
        &inner,
        &skin_map,
        &offset_result.surfaces,
        &stage_tol,
        options.join_style,
        true,
    )?;
    vertex_recompute::recompute_vertices_with_gap_corners(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
        &join_plan.gap_corners,
    )?;
    tracker.finish_stage("shell.vertex_recompute", t0, inner.vertices.len());
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let undone = edge_recompute::recompute_edges_with_skips(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
        &join_plan.arc_skip_edges(),
    )?;
    join_classification::intersection_gap_postcheck(&join_plan, &undone)?;
    // Collapse gate (part of the recompute stage, O(E) and sub-stage
    // cheap): a thickness past the local feature size inverts skin
    // edges; the fold would serve structurally-clean wrong geometry.
    // Inverted chords defer their decline to the trim stage (PK-327).
    let trim_routing =
        collapse_check::check_skin_collapse(brep, &inner, &skin_map, thickness, &tolerance)?;
    tracker.finish_stage("shell.edge_recompute", t0, inner.edges.len());
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    inner_skin::apply_offset_geometry(&mut inner, &offset_result.surfaces, &skin_map)?;
    // Join construction (PK-325): rolling-ball strips + sphere corner
    // windows on the recomputed skin, via Euler operators only.
    let arc = arc_join::construct_arc_joins(
        &mut inner,
        &join_plan,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    let join_records = arc.records;
    tracker.finish_stage("shell.apply_offset_geometry", t0, inner.faces.len());

    // Mirror the inner skin's material side honestly: winding and
    // orientation flags flip together, so the void skin keeps the
    // interior-left invariant (a flag-only flip ships the PK-244
    // winding-dishonest mirror class).
    inner.reverse_body_orientation();

    // Stage 5: assemble outer + inner + side walls.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let (mut result, assembly_maps) =
        assemble::assemble_shell_with_maps(brep, &inner, &all_removed, &skin_map, tolerance)?;
    tracker.finish_stage("shell.assemble", t0, result.faces.len());

    // Exact face provenance for stable naming (PK-328): handles are
    // final from assembly on; trim-stage deletions are checked at
    // assignment time.
    if let Some(sink) = naming_sink {
        for (&src, &res) in &assembly_maps.outer_faces {
            sink.modified.push((res, src));
        }
        for (&src, &skin_fh) in &skin_map.faces {
            if let Some(&res) = assembly_maps.inner_faces.get(&skin_fh) {
                sink.offset.push((res, src));
            }
        }
        for record in &join_records {
            if let Some(&res) = assembly_maps.inner_faces.get(&record.face) {
                sink.joins.push((res, join_source_name(brep, record.kind)));
            }
        }
        sink.sort();
    }

    // Stage 6: recompute pcurves for all faces.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    result
        .recompute_all_pcurves_winding_repaired()
        .map_err(|e| ShellError::Topology(e.to_string()))?;
    tracker.finish_stage("shell.pcurve_recompute", t0, result.faces.len());

    // Orientation flags are DERIVED state (a function of finalized
    // loop winding + surface geometry): reconcile them at the serving
    // boundary. Join faces on charts whose handedness opposes their
    // construction-time material measurement (the cone-rim torus
    // strips) get their honest value here (PK-325).
    result.reconcile_orientation_flags();

    // Global self-intersection trim stage (PK-327): the outer skin IS
    // the source boundary (expected distance 0, trivially realized);
    // the moved inner skin expects the shell thickness; walls are
    // exempt.
    let mut provenance = global_si_trim::TrimProvenance {
        default_distance: thickness.abs(),
        ..Default::default()
    };
    for new_fh in assembly_maps.outer_faces.values() {
        provenance.expected.insert(*new_fh, 0.0);
    }
    let mapped: std::collections::BTreeSet<_> = assembly_maps
        .outer_faces
        .values()
        .chain(assembly_maps.inner_faces.values())
        .copied()
        .collect();
    for fh in result.faces.keys() {
        if !mapped.contains(&fh) {
            provenance.always_keep.insert(fh);
        }
    }
    provenance.free_form_offsets = global_si_trim::collect_free_form_pairs(
        brep,
        &offset_result.surfaces,
        |fh| {
            skin_map
                .faces
                .get(&fh)
                .and_then(|s| assembly_maps.inner_faces.get(s))
                .copied()
        },
        |_| thickness.abs(),
    );
    let trim_warnings = global_si_trim::run_trim_stage(
        &mut result,
        brep,
        &provenance,
        options.trim_policy,
        trim_routing.deferred,
        trim_routing.inverted_edges.len(),
        &tolerance,
    )
    .map_err(|e| {
        max_safe::upgrade_trim_decline(
            brep,
            direction,
            thickness,
            tolerance,
            options.trim_policy,
            e,
        )
    })?;

    // Serve gate: never hand back a body a validator would reject. The
    // structural collection is the affordable in-pipeline form (PK-273);
    // tests and acceptance batteries additionally hold results to
    // validate_all.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        if popkern_core::trace::trace_enabled("shell-serve") {
            for e in &validation_errors {
                eprintln!("[serve-trace] {e}");
            }
            for (eh, edge) in &result.edges {
                let ends = {
                    let he = &result.half_edges[edge.half_edge];
                    let tw = &result.half_edges[he.twin];
                    (
                        result.vertices[he.origin].position,
                        result.vertices[tw.origin].position,
                    )
                };
                eprintln!(
                    "[serve-trace] {eh:?} curve={:?} ends=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    edge.curve.as_ref().map(|c| format!("{:?}", c.geometry)
                        .chars()
                        .take(90)
                        .collect::<String>()),
                    ends.0.x,
                    ends.0.y,
                    ends.0.z,
                    ends.1.x,
                    ends.1.y,
                    ends.1.z
                );
            }
        }
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    // Collect warnings.
    let mut warnings = offset_result.warnings;
    warnings.extend(trim_warnings);
    warnings.extend(collapse_warnings);
    warnings.extend(join_plan.warnings.iter().cloned());
    arc_join_apex::retire_served_apex_warnings(&mut warnings, &arc.served_apexes);
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
#[path = "shell_entry_tests.rs"]
mod tests;
