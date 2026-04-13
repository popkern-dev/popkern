//! Pipeline orchestrator: connects all 5 stages of the boolean pipeline.
//!
//! The [`execute_boolean`] function drives intersection, split, classify,
//! select, and merge in sequence, handling the no-intersection shortcut
//! and accumulating diagnostics on failure.

#[path = "operation_helpers.rs"]
pub(crate) mod helpers;

use std::time::Instant;

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::cancellation::CancellationToken;
use crate::classify::classify_fragments;
use crate::error::{BoolDiagnostics, BoolError, BoolOperation, PhaseMetric, Stage};
use crate::intersect::find_intersections;
use crate::merge::graph_merge::graph_merge;
use crate::progress::{ProgressEvent, ProgressReporter};
use crate::select::{Operand, select_faces};
use crate::split::split_along_curves;

const TOTAL_PHASES: usize = 6;

/// Notify the reporter (if any) that a stage just completed.
///
/// Takes `&mut Option<...>` so the caller can re-borrow `progress`
/// across multiple calls without moving it.
fn report_phase(
    reporter: &mut Option<&mut dyn ProgressReporter>,
    stage: Stage,
    phase_index: usize,
    duration: std::time::Duration,
    entity_count: usize,
) {
    if let Some(r) = reporter.as_deref_mut() {
        r.on_phase(ProgressEvent {
            stage,
            phase_index,
            total_phases: TOTAL_PHASES,
            duration,
            entity_count,
        });
    }
}

/// Check `cancel` and short-circuit with `BoolError::Cancelled` if set.
///
/// `Option<&CancellationToken>` keeps existing call sites that pass
/// `None` zero-cost — there is no atomic load when the caller did not
/// supply a token.
#[inline]
fn check_cancel(cancel: Option<&CancellationToken>) -> Result<(), BoolError> {
    if let Some(t) = cancel
        && t.is_cancelled()
    {
        return Err(BoolError::Cancelled);
    }
    Ok(())
}

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
    execute_boolean_full(op, a, b, tolerance, &mut Vec::new(), None, None)
}

/// Execute a boolean and populate `metrics` with per-stage timing samples.
///
/// Backwards-compatible wrapper used by `boolean_with_diagnostics`.
pub(crate) fn execute_boolean_with_metrics(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
) -> Result<BRep, BoolError> {
    execute_boolean_full(op, a, b, tolerance, metrics, None, None)
}

/// Execute a boolean with full pipeline knobs: telemetry, cancellation,
/// and progress reporting.
///
/// Appends a [`PhaseMetric`] after each stage that runs to completion,
/// polls `cancel` between stages, and notifies `progress` (if any)
/// after each stage with a [`ProgressEvent`]. Returns
/// [`BoolError::Cancelled`] as soon as a cancellation is observed;
/// the input model is unchanged.
pub(crate) fn execute_boolean_full(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
    cancel: Option<&CancellationToken>,
    mut progress: Option<&mut dyn ProgressReporter>,
) -> Result<BRep, BoolError> {
    check_cancel(cancel)?;

    // Stage 1: Intersect — compute face-face intersection curves.
    let t_intersect = Instant::now();
    let intersection_result =
        find_intersections(&a, &b, tolerance).map_err(|diag| diag.into_error())?;
    let elapsed_intersect = t_intersect.elapsed();
    let count_intersect = intersection_result.records.len();
    metrics.push(PhaseMetric {
        stage: Stage::Intersection,
        duration: elapsed_intersect,
        entity_count: count_intersect,
    });
    report_phase(
        &mut progress,
        Stage::Intersection,
        0,
        elapsed_intersect,
        count_intersect,
    );
    check_cancel(cancel)?;

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

    // §6.4: deferred ON audit. When the coincidence map contains
    // partial-overlap pairs, log the count so callers (and debug
    // builds) know how many ambiguous-on cases exist. The audit
    // itself is read-only and runs after assembly via the
    // public `classify_on_deferred::audit_ambiguous_on_classifications`
    // helper for callers that want to inspect or override the
    // immediate ON resolution.
    #[cfg(debug_assertions)]
    {
        let ambiguous =
            crate::classify_on_deferred::find_ambiguous_pairs(&intersection_result.coincidence);
        if !ambiguous.is_empty() {
            eprintln!(
                "[boolean] {} ambiguous-on face pair(s) detected (partial overlaps)",
                ambiguous.len()
            );
        }
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
    let t_split = Instant::now();
    let split_result = split_along_curves(a, b, &intersection_result.records, tolerance)
        .map_err(|diag| diag.into_error())?;
    let split_fragment_count = split_result.a.faces.len() + split_result.b.faces.len();
    let elapsed_split = t_split.elapsed();
    metrics.push(PhaseMetric {
        stage: Stage::Split,
        duration: elapsed_split,
        entity_count: split_fragment_count,
    });
    report_phase(
        &mut progress,
        Stage::Split,
        1,
        elapsed_split,
        split_fragment_count,
    );
    check_cancel(cancel)?;

    // Stage 3: Classify — label each fragment as inside/outside/on-boundary.
    let t_classify = Instant::now();
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
    let elapsed_classify = t_classify.elapsed();
    let count_classify = classifications_a.len() + classifications_b.len();
    metrics.push(PhaseMetric {
        stage: Stage::Classification,
        duration: elapsed_classify,
        entity_count: count_classify,
    });
    report_phase(
        &mut progress,
        Stage::Classification,
        2,
        elapsed_classify,
        count_classify,
    );
    check_cancel(cancel)?;

    // Build the combined classification list for Stage 4.
    let mut all_classifications = Vec::new();
    for (face, class) in &classifications_a {
        all_classifications.push((*face, Operand::A, *class));
    }
    for (face, class) in &classifications_b {
        all_classifications.push((*face, Operand::B, *class));
    }

    // Stage 4: Select — choose which fragments belong to the result.
    let t_select = Instant::now();
    let selection = select_faces(op, &all_classifications);
    let elapsed_select = t_select.elapsed();
    let count_select = selection.keep_a.len() + selection.keep_b.len();
    metrics.push(PhaseMetric {
        stage: Stage::Selection,
        duration: elapsed_select,
        entity_count: count_select,
    });
    report_phase(
        &mut progress,
        Stage::Selection,
        3,
        elapsed_select,
        count_select,
    );
    check_cancel(cancel)?;

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
    let t_merge = Instant::now();
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
    let elapsed_merge = t_merge.elapsed();
    let count_merge = result.faces.len();
    metrics.push(PhaseMetric {
        stage: Stage::Merge,
        duration: elapsed_merge,
        entity_count: count_merge,
    });
    report_phase(&mut progress, Stage::Merge, 4, elapsed_merge, count_merge);

    // Recompute pcurves on all faces so that the tessellation pipeline
    // can build accurate trim loops. Boolean merge may leave seam edges
    // and original-operand edges without pcurves.
    let mut result = result;
    let _ = result.recompute_all_pcurves();

    // §12.5 / §7.3: gated post-validation. Run the boolean-specific
    // post-check (manifold, Euler, half-edge cycles, watertightness).
    // Surface any errors as a structured failure rather than silently
    // returning a corrupted result. The merge-phase cut-behavior fix
    // (commit history) made this consistently green; if a future
    // change regresses, fail loudly.
    let t_validate = Instant::now();
    let pc_report = crate::validate::post_check::post_check(&result, tolerance, false);
    let pc_errors: Vec<String> = pc_report
        .findings
        .iter()
        .filter(|f| {
            matches!(
                f.severity,
                crate::error::Severity::Error | crate::error::Severity::Critical
            )
        })
        .map(|f| f.message.clone())
        .collect();
    log_validation_issues(&result, op);
    let elapsed_validate = t_validate.elapsed();
    let count_validate = result.faces.len();
    metrics.push(PhaseMetric {
        stage: Stage::Validation,
        duration: elapsed_validate,
        entity_count: count_validate,
    });
    report_phase(
        &mut progress,
        Stage::Validation,
        5,
        elapsed_validate,
        count_validate,
    );
    if !pc_errors.is_empty() {
        return Err(BoolError::ValidationFailed(pc_errors.join("; ")));
    }

    // Post-processing: prune stranded half-edges left behind by
    // broken merge loops. The merge phase occasionally tags
    // half-edges with a loop handle they're not reachable from;
    // this pruning step removes them so the visible topology
    // matches the Euler-Poincaré accounting.
    let _ = crate::postprocess::loop_repair::repair_stranded_loops(&mut result);

    // Post-processing: remove redundant edges if the result topology is valid.
    if result.validate().is_ok() {
        remove_redundant_edges_safe(&mut result, tolerance);
    }

    // Post-processing: remove sliver faces (§8.3.1) if the topology
    // is still valid. The pass merges each sliver into a non-sliver
    // neighbor via KEF, validating after every removal and bailing
    // out on any corruption rather than compounding damage.
    if result.validate().is_ok() {
        let _ = crate::postprocess::sliver_removal::remove_slivers(&mut result, tolerance);
    }

    // Post-processing: face simplification (§11.2). Merges adjacent
    // co-surfacial faces via KEF, validating after each merge with
    // rollback on regression. The previous implementation was
    // disabled because of a stale-edge-handle bug; the iterate-and-
    // recompute pattern fixes that.
    if result.validate().is_ok() {
        crate::postprocess::face_simplify::apply_face_simplification(&mut result, tolerance);
    }
    // Note: fix_zero_thickness_sheets remains disabled — its global
    // face-pair comparison is unrelated to this fix and still risky.
    // fix_zero_thickness_sheets(&mut result, tolerance);

    // Post-processing: collapse micro-edges (position merge).
    let micro_edges = crate::imprint::micro_edge::find_micro_edges(&result, tolerance.linear());
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
    let bool_tol = tolerance.linear();

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

    let sliver_thresholds = SliverThresholds::from_tolerance(tolerance.linear() * 5.0);

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
    let wt_report = validate_watertightness(result, tolerance.linear() * 10.0, 3);
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
    pub(crate) fn into_error(self) -> BoolError {
        self.error.unwrap_or_else(|| BoolError::MergeFailed {
            stage: format!("{:?}", self.failed_stage),
            reason: "unknown error in boolean pipeline".to_string(),
        })
    }
}

#[cfg(test)]
#[path = "operation_tests.rs"]
mod tests;
