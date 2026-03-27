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

    // Stage 5: Graph Merge — build intersection graph and assemble result.
    let result = graph_merge(
        &original_a,
        &original_b,
        &intersection_result,
        op,
        tolerance,
    )
    .map_err(|diag| diag.into_error())?;

    // Recompute pcurves on all faces so that the tessellation pipeline
    // can build accurate trim loops. Boolean merge may leave seam edges
    // and original-operand edges without pcurves.
    let mut result = result;
    let _ = result.recompute_all_pcurves();

    // Log validation issues (don't fail — diagnostic only for now).
    log_validation_issues(&result, op);

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
