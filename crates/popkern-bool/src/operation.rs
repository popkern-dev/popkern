//! Pipeline orchestrator: connects all 5 stages of the boolean pipeline.
//!
//! The [`execute_boolean`] function drives intersection, split, classify,
//! select, and merge in sequence, handling the no-intersection shortcut
//! and accumulating diagnostics on failure.

#[path = "operation_helpers.rs"]
mod helpers;

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::classify::classify_fragments;
use crate::error::{BoolDiagnostics, BoolError, BoolOperation};
use crate::intersect::find_intersections;
use crate::merge::graph_merge::graph_merge;
use crate::select::{Operand, select_faces};
use crate::split::split_along_curves;

use self::helpers::{handle_empty_selection, handle_no_intersection};

/// Execute a boolean operation on two BRep solids.
///
/// Runs the full five-stage pipeline: intersect, split, classify, select,
/// and merge. When no intersection curves are found, takes a shortcut based
/// on the spatial relationship between the two solids.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
///
/// # Errors
///
/// Returns [`BoolError`] describing the failure, wrapping partial diagnostics
/// from whichever stage failed.
pub(crate) fn execute_boolean(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    // Stage 1: Intersect — compute face-face intersection curves.
    let intersection_result =
        find_intersections(&a, &b, tolerance).map_err(|diag| diag.into_error())?;

    // No-intersection shortcut.
    if intersection_result.records.is_empty() {
        return handle_no_intersection(
            op,
            a,
            b,
            &intersection_result.aabbs_a,
            &intersection_result.aabbs_b,
            tolerance,
        );
    }

    // Coincident-body shortcut: if all face overlaps cover at least half
    // the faces of either body, treat as full containment.
    if !intersection_result.coincidence.face_overlaps.is_empty() {
        let n_a = a.faces.len();
        let n_b = b.faces.len();
        let n_overlaps = intersection_result.coincidence.face_overlaps.len();
        if n_overlaps >= n_a.min(n_b) / 2 {
            return handle_no_intersection(
                op,
                a,
                b,
                &intersection_result.aabbs_a,
                &intersection_result.aabbs_b,
                tolerance,
            );
        }
    }

    // Clone originals for graph merge (graph merge reads originals, not split).
    let original_a = a.clone();
    let original_b = b.clone();

    // Stage 2: Split — trim faces along intersection curves.
    let split_result = split_along_curves(a, b, &intersection_result.records, tolerance)
        .map_err(|diag| diag.into_error())?;

    // Stage 3: Classify — label each fragment as inside/outside/on-boundary.
    let classifications_a = classify_fragments(
        &split_result.a,
        &split_result.b,
        &intersection_result.aabbs_b,
        &split_result.fragments_a.intersection_edges,
        tolerance,
    )?;

    let classifications_b = classify_fragments(
        &split_result.b,
        &split_result.a,
        &intersection_result.aabbs_a,
        &split_result.fragments_b.intersection_edges,
        tolerance,
    )?;

    // Build the combined classification list for Stage 4.
    let mut all_classifications = Vec::new();
    for (face, class) in &classifications_a {
        all_classifications.push((*face, Operand::A, *class));
    }
    for (face, class) in &classifications_b {
        all_classifications.push((*face, Operand::B, *class));
    }

    // Stage 4: Select — choose which fragments belong to the result.
    let selection = select_faces(op, &all_classifications);

    // If selection is empty, the split/classify pipeline classified all fragments
    // the same way (e.g. all Inside). Infer containment and return the
    // appropriate operand directly instead of producing a 0-face result.
    if selection.keep_a.is_empty() && selection.keep_b.is_empty() {
        let result =
            handle_empty_selection(op, split_result, &classifications_a, &classifications_b);
        log_validation_issues(&result, op);
        return Ok(result);
    }

    // Stage 5: Cellular graph merge — order-independent fragment selection.
    // Falls back to standard graph-merge if cellular fails or produces
    // invalid topology.
    let result = crate::cellular::boolean::cellular_graph_merge(
        &original_a,
        &original_b,
        &intersection_result,
        op,
        tolerance,
    )
    .ok()
    .filter(|r| {
        let valid = r.validate().is_ok();
        #[cfg(debug_assertions)]
        if !valid {
            eprintln!("[boolean] cellular assembly invalid, falling back to graph-merge");
        }
        valid
    })
    .map(Ok)
    .unwrap_or_else(|| {
        // Fallback: standard graph-merge
        graph_merge(
            &original_a,
            &original_b,
            &intersection_result,
            op,
            tolerance,
        )
    })
    .map_err(|diag| diag.into_error())?;

    // Recompute pcurves on all faces so that the tessellation pipeline
    // can build accurate trim loops. Boolean merge may leave seam edges
    // and original-operand edges without pcurves.
    let mut result = result;
    let _ = result.recompute_all_pcurves();

    // Log validation issues (don't fail — diagnostic only for now).
    log_validation_issues(&result, op);

    // Post-processing: remove redundant edges if the result topology is valid.
    if result.validate().is_ok() {
        remove_redundant_edges_safe(&mut result, tolerance);
    }

    // Post-processing: face simplification and zero-thickness sheet removal
    // are disabled — the current implementations are too aggressive and
    // remove valid faces from boolean results. They compare all face pairs
    // globally, but should only compare ADJACENT face pairs.
    // fix_zero_thickness_sheets(&mut result, tolerance);
    // crate::postprocess::face_simplify::apply_face_simplification(&mut result, tolerance);

    // Post-processing: collapse micro-edges (position merge).
    let micro_edges = crate::imprint::micro_edge::find_micro_edges(&result, tolerance.linear);
    if !micro_edges.is_empty() {
        crate::imprint::micro_edge::collapse_micro_edge_positions(&mut result, &micro_edges);
    }

    // Post-processing: propagate tolerances to intersection entities (§10.3).
    // Intersection-created vertices/edges may have larger tolerance than
    // original entities due to numerical approximation.
    propagate_output_tolerances(&mut result, tolerance);

    // Post-processing diagnostics: sliver detection and watertightness.
    log_post_processing_diagnostics(&result, tolerance);

    Ok(result)
}

/// Log any validation issues on a boolean result BRep.
///
/// Runs `validate_all()` and prints a summary to stderr. Does not fail the
/// operation — this is purely diagnostic to surface topological problems early.
/// Catches panics from validation on corrupted topology.
fn log_validation_issues(result: &BRep, op: BoolOperation) {
    let op_name = match op {
        BoolOperation::Union => "union",
        BoolOperation::Intersection => "intersection",
        BoolOperation::Subtraction => "subtraction",
    };

    // Validation may panic on corrupted topology (invalid SlotMap keys).
    // Catch panics so diagnostics don't crash the pipeline.
    let validation_result =
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| result.validate_all()));

    match validation_result {
        Ok(errors) if !errors.is_empty() => {
            eprintln!(
                "[boolean] validation: {} issues in {} result",
                errors.len(),
                op_name,
            );
        }
        Err(_) => {
            eprintln!(
                "[boolean] validation: panicked on {} result (corrupted topology)",
                op_name,
            );
        }
        _ => {}
    }
}

/// Propagate tolerances to output entities (spec §10.3).
///
/// Vertices and edges that were created during the boolean (not present
/// in either input) get an explicit tolerance reflecting the accuracy
/// of intersection computation. Original entities keep `None` (inherit
/// from BRep default).
fn propagate_output_tolerances(result: &mut BRep, tolerance: &Tolerance) {
    // Intersection-created vertices have tolerance = linear tolerance
    // (the accuracy of the Newton solver). Original vertices keep None.
    // Since we can't distinguish them here (the graph info is gone),
    // set tolerance on vertices that have None and edges that have None
    // only if they don't match the BRep's default tolerance. This is
    // a no-op for most cases — the real benefit comes when the boolean
    // is called with a tolerance different from the BRep's default.
    let bool_tol = tolerance.linear;

    // Set per-entity tolerance on all entities to the boolean tolerance
    // so downstream operations know the achieved accuracy.
    for (_, vertex) in &mut result.vertices {
        if vertex.tolerance.is_none() {
            vertex.tolerance = Some(bool_tol);
        }
    }
    for (_, edge) in &mut result.edges {
        if edge.tolerance.is_none() {
            edge.tolerance = Some(bool_tol);
        }
    }
    for (_, face) in &mut result.faces {
        if face.tolerance.is_none() {
            face.tolerance = Some(bool_tol);
        }
    }
}

/// Fix zero-thickness sheet faces: two faces sharing the same surface
/// with the same orientation, which creates a non-manifold internal sheet.
/// Remove one of the duplicate faces.
#[allow(dead_code)]
fn fix_zero_thickness_sheets(result: &mut BRep, tolerance: &Tolerance) {
    use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;

    let face_handles: Vec<_> = result.faces.keys().collect();
    let mut to_remove = Vec::new();

    for i in 0..face_handles.len() {
        if to_remove.contains(&face_handles[i]) {
            continue;
        }
        for j in (i + 1)..face_handles.len() {
            if to_remove.contains(&face_handles[j]) {
                continue;
            }

            let surf_a = match result
                .faces
                .get(face_handles[i])
                .and_then(|f| f.surface.as_ref())
            {
                Some(s) => s,
                None => continue,
            };
            let surf_b = match result
                .faces
                .get(face_handles[j])
                .and_then(|f| f.surface.as_ref())
            {
                Some(s) => s,
                None => continue,
            };

            // Same surface, same orientation → zero-thickness sheet
            if surf_a.orientation == surf_b.orientation
                && are_surfaces_coincident_analytic(&surf_a.geometry, &surf_b.geometry, tolerance)
            {
                to_remove.push(face_handles[j]);
            }
        }
    }

    // Remove duplicate faces (simple removal — topology may become invalid,
    // but this is a best-effort cleanup).
    for fh in &to_remove {
        result.faces.remove(*fh);
    }

    if !to_remove.is_empty() {
        #[cfg(debug_assertions)]
        eprintln!(
            "[boolean] removed {} zero-thickness sheet face(s)",
            to_remove.len()
        );
    }
}

/// Safely remove redundant edges, validating after each removal.
///
/// Finds redundant edges and removes them one at a time via KEF.
/// After each removal, validates the topology. If validation fails,
/// stops immediately (leaving remaining redundant edges in place).
fn remove_redundant_edges_safe(result: &mut BRep, tolerance: &Tolerance) {
    use crate::postprocess::redundant_edge::find_redundant_edges;

    let report = find_redundant_edges(result, tolerance);

    if report.redundant_edges.is_empty() {
        return;
    }

    let mut removed = 0;
    for eh in &report.redundant_edges {
        let he = match result.edges.get(*eh) {
            Some(e) => e.half_edge,
            None => continue,
        };

        // Snapshot vertex/edge/face counts before KEF
        let pre_v = result.vertices.len();
        let pre_e = result.edges.len();
        let pre_f = result.faces.len();

        if result.kef(he).is_ok() {
            // Validate after removal — if topology is corrupted, stop
            if result.validate().is_err() {
                // Revert is impossible, but at least stop making it worse
                #[cfg(debug_assertions)]
                eprintln!(
                    "[boolean] KEF corrupted topology after removing edge, stopping. \
                     Before: V={pre_v} E={pre_e} F={pre_f}"
                );
                break;
            }
            removed += 1;
        }
    }

    if removed > 0 {
        #[cfg(debug_assertions)]
        eprintln!("[boolean] post-processing: removed {removed} redundant edge(s)");
    }
}

/// Log post-processing diagnostics: sliver detection and watertightness.
///
/// Runs sliver analysis on each face and watertightness validation on all
/// edges. Results are logged to stderr (diagnostic only, does not fail).
fn log_post_processing_diagnostics(result: &BRep, tolerance: &Tolerance) {
    use crate::clip::face_vertex_positions;
    use crate::imprint::shared_geometry::validate_watertightness;
    use crate::imprint::sliver::{SliverThresholds, analyze_face_polygon};

    let sliver_thresholds = SliverThresholds::from_tolerance(tolerance.linear * 5.0);

    // Sliver detection
    let mut sliver_count = 0;
    for fh in result.faces.keys() {
        let positions = face_vertex_positions(result, fh);
        if positions.len() >= 3 {
            let report = analyze_face_polygon(&positions, &sliver_thresholds);
            if report.is_sliver {
                sliver_count += 1;
            }
        }
    }
    if sliver_count > 0 {
        #[cfg(debug_assertions)]
        eprintln!("[boolean] post-processing: {sliver_count} sliver face(s) detected");
    }

    // Watertightness validation
    let wt_report = validate_watertightness(result, tolerance.linear * 10.0, 3);
    if !wt_report.is_watertight() {
        #[cfg(debug_assertions)]
        eprintln!(
            "[boolean] post-processing: {} watertightness gap(s), max deviation: {:.2e}",
            wt_report.gaps.len(),
            wt_report.max_gap()
        );
    }
}

impl BoolDiagnostics {
    /// Extract the error from diagnostics, or produce a generic error if none.
    fn into_error(self) -> BoolError {
        self.error.unwrap_or_else(|| BoolError::MergeFailed {
            stage: format!("{:?}", self.failed_stage),
            reason: "unknown error in boolean pipeline".to_string(),
        })
    }
}

#[cfg(test)]
#[path = "operation_tests.rs"]
mod tests;
