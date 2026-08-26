//! Midplane shells: both skins move — outer out of the material by half
//! the thickness, inner into it by half — so the wall straddles the
//! reference surface. Openings (removed faces) get rim walls spanning
//! the two MOVED skins; the opening rim stays flush with each removed
//! face's original surface (the rim constraints in vertex/edge
//! recompute use the unclassified faces' original geometry).
//!
//! Until PK-321 this direction silently served an inward shell at half
//! thickness (the outer skin never moved), then declined as named
//! containment; PK-330 built this construction.

use crate::error::ShellError;
use crate::skin_ladder::{SKIN_LADDER_STAGE_COUNT, combine_skin_maps, offset_skin_ladder};
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Stages a midplane shell emits (PK-276 progress totals): two skin
/// ladders, assembly, pcurve recompute.
pub(crate) const MIDPLANE_STAGE_COUNT: usize = 2 * SKIN_LADDER_STAGE_COUNT + 2;

/// Build a midplane shell: outer skin at +thickness/2 (material-outward),
/// inner skin at −thickness/2, side walls at removed-face openings.
#[allow(clippy::too_many_arguments)]
pub(crate) fn midplane_shell_impl(
    brep: &BRep,
    thickness: f64,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
    join_sink: Option<&mut Vec<crate::JoinFaceData>>,
) -> Result<ShellResult, ShellError> {
    use std::time::Instant;
    let half = thickness / 2.0;
    let half_faces: Vec<(FaceHandle, f64)> = face_thicknesses
        .iter()
        .map(|&(fh, t)| (fh, t / 2.0))
        .collect();
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);

    // Outer skin: material-outward by half. Removed faces stay PRESENT
    // (unclassified) so assembly can walk their loops for rim discovery;
    // the copy filter drops them and the walls adopt their rims.
    let (outer_skin, outer_map, offset_out, joins_out, routing_out) = offset_skin_ladder(
        brep,
        half,
        &half_faces,
        ShellDirection::Outward,
        removed_faces,
        false,
        options.allow_partial,
        options.join_style,
        true,
        &stage_tol,
        2.0,
        [
            "shell.midplane.outward.classify_faces",
            "shell.midplane.outward.offset_surfaces",
            "shell.midplane.outward.skin",
            "shell.midplane.outward.recompute",
            "shell.midplane.outward.apply_offset_geometry",
        ],
        tracker,
    )?;

    // Faces whose OUTWARD offset failed under allow_partial become
    // openings for the inner ladder too, so both skins agree on the
    // removed set. A face failing only the INWARD ladder stays present
    // in the outer skin without a counterpart — the serve gate declines
    // that mismatch rather than serving it.
    let mut all_removed: Vec<FaceHandle> = removed_faces.to_vec();
    all_removed.extend_from_slice(&offset_out.auto_removed);

    // Inner skin: material-inward by half, removed faces deleted (the
    // classic void-skin shape whose free rim half-edges the walls adopt).
    let (mut inner_skin, inner_map, offset_in, joins_in, mut routing_in) = offset_skin_ladder(
        brep,
        half,
        &half_faces,
        ShellDirection::Inward,
        &all_removed,
        true,
        options.allow_partial,
        options.join_style,
        true,
        &stage_tol,
        2.0,
        [
            "shell.midplane.inward.classify_faces",
            "shell.midplane.inward.offset_surfaces",
            "shell.midplane.inward.skin",
            "shell.midplane.inward.recompute",
            "shell.midplane.inward.apply_offset_geometry",
        ],
        tracker,
    )?;
    all_removed.extend_from_slice(&offset_in.auto_removed);

    // Join-face naming data (PK-325 row 474): captured on the skins
    // before assembly re-homes every handle.
    if let Some(sink) = join_sink {
        sink.extend(crate::collect_join_face_data(
            &outer_skin,
            brep,
            &joins_out.records,
        ));
        sink.extend(crate::collect_join_face_data(
            &inner_skin,
            brep,
            &joins_in.records,
        ));
    }

    // Mirror the inner skin honestly (winding + flag together, PK-244).
    inner_skin.reverse_body_orientation();

    // Assembly keys on the OUTER SKIN's handle space.
    let combined = combine_skin_maps(&outer_map, &inner_map);
    let mapped_removed: Vec<FaceHandle> = all_removed
        .iter()
        .filter_map(|fh| outer_map.faces.get(fh).copied())
        .collect();

    if tracker.cancelled() {
        return Err(ShellError::Cancelled);
    }
    let t0 = Instant::now();
    let (mut result, assembly_maps) = crate::assemble::assemble_shell_with_maps(
        &outer_skin,
        &inner_skin,
        &mapped_removed,
        &combined,
        tolerance,
    )?;
    tracker.finish_stage("shell.midplane.assemble", t0, result.faces.len());

    if tracker.cancelled() {
        return Err(ShellError::Cancelled);
    }
    let t0 = Instant::now();
    result
        .recompute_all_pcurves_winding_repaired()
        .map_err(|e| ShellError::Topology(e.to_string()))?;
    tracker.finish_stage("shell.midplane.pcurve_recompute", t0, result.faces.len());

    // Orientation flags are DERIVED state (a function of finalized
    // loop winding + surface geometry): reconcile them at the serving
    // boundary. Join faces on charts whose handedness opposes their
    // construction-time material measurement (the cone-rim torus
    // strips) get their honest value here (PK-325).
    result.reconcile_orientation_flags();

    // Global self-intersection trim stage (PK-327): both skins moved
    // by half the thickness; walls and corner sectors are exempt.
    let mut routing = routing_out;
    routing
        .inverted_edges
        .append(&mut routing_in.inverted_edges);
    if routing.deferred.is_none() {
        routing.deferred = routing_in.deferred.take();
    }
    let mut provenance = crate::global_si_trim::TrimProvenance {
        default_distance: half,
        ..Default::default()
    };
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
    for (offset, ladder_map, assembly_faces) in [
        (&offset_out, &outer_map, &assembly_maps.outer_faces),
        (&offset_in, &inner_map, &assembly_maps.inner_faces),
    ] {
        provenance
            .free_form_offsets
            .extend(crate::global_si_trim::collect_free_form_pairs(
                brep,
                &offset.surfaces,
                |fh| {
                    ladder_map
                        .faces
                        .get(&fh)
                        .and_then(|s| assembly_faces.get(s))
                        .copied()
                },
                |fh| {
                    half_faces
                        .iter()
                        .find(|(f, _)| *f == fh)
                        .map(|(_, d)| *d)
                        .unwrap_or(half)
                },
            ));
    }
    let trim_warnings = crate::global_si_trim::run_trim_stage(
        &mut result,
        brep,
        &provenance,
        options.trim_policy,
        routing.deferred,
        routing.inverted_edges.len(),
        &tolerance,
    )?;

    // Serve gate: never hand back a body a validator would reject.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    let mut warnings = offset_out.warnings;
    warnings.extend(trim_warnings);
    warnings.extend(offset_in.warnings);
    warnings.extend(joins_out.warnings);
    warnings.extend(joins_in.warnings);
    for face in offset_out
        .auto_removed
        .iter()
        .chain(&offset_in.auto_removed)
    {
        warnings.push(crate::ShellWarning::FaceAutoRemoved {
            face: *face,
            reason: "offset failed".into(),
        });
    }

    Ok(ShellResult {
        brep: result,
        warnings,
    })
}
