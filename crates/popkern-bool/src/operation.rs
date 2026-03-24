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
use crate::merge::merge_faces;
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

    // Stage 5: Merge — stitch selected fragments into a valid BRep.
    let mut result =
        merge_faces(&split_result, &selection, tolerance).map_err(|diag| diag.into_error())?;

    // For union: strip inner loops from result faces. Inner loops were created
    // by the split stage to represent intersection regions, but in a union those
    // regions are covered by the other operand's faces.
    if op == BoolOperation::Union {
        strip_inner_loops(&mut result);
    }

    // Log validation issues (don't fail — diagnostic only for now).
    log_validation_issues(&result, op);

    Ok(result)
}

/// Strip all inner loops from faces in a BRep, removing the associated
/// topology (loops, half-edges, edges, vertices) that become orphaned.
///
/// Used for union results where inner loops represent intersection regions
/// that are already covered by the other operand's faces.
fn strip_inner_loops(brep: &mut BRep) {
    use std::collections::HashSet;

    // Collect all inner loop handles and clear them from faces.
    let mut inner_loops_to_remove: Vec<popkern_core::topology::LoopHandle> = Vec::new();
    for (_fh, face) in &mut brep.faces {
        inner_loops_to_remove.append(&mut face.inner_loops);
    }

    if inner_loops_to_remove.is_empty() {
        return;
    }

    // Collect half-edges belonging to inner loops.
    let mut he_to_remove = HashSet::new();
    for &lh in &inner_loops_to_remove {
        if let Some(loop_) = brep.loops.get(lh) {
            let start = loop_.half_edge;
            if start != popkern_core::topology::HalfEdgeHandle::default() {
                let mut current = start;
                loop {
                    he_to_remove.insert(current);
                    current = brep.half_edges[current].next;
                    if current == start {
                        break;
                    }
                }
            }
        }
    }

    // Remove inner loops.
    for lh in &inner_loops_to_remove {
        brep.loops.remove(*lh);
    }

    // Collect edges referenced only by removed half-edges.
    let mut edges_to_remove = HashSet::new();
    for &heh in &he_to_remove {
        if let Some(he) = brep.half_edges.get(heh) {
            edges_to_remove.insert(he.edge);
        }
    }
    // Keep edges that are still referenced by surviving half-edges.
    for (heh, he) in &brep.half_edges {
        if !he_to_remove.contains(&heh) {
            edges_to_remove.remove(&he.edge);
        }
    }

    // Collect vertices referenced only by removed half-edges.
    let mut vertices_to_remove = HashSet::new();
    for &heh in &he_to_remove {
        if let Some(he) = brep.half_edges.get(heh) {
            vertices_to_remove.insert(he.origin);
        }
    }
    // Keep vertices that are still referenced by surviving half-edges.
    for (heh, he) in &brep.half_edges {
        if !he_to_remove.contains(&heh) {
            vertices_to_remove.remove(&he.origin);
        }
    }

    // Remove the topology entities.
    for heh in &he_to_remove {
        brep.half_edges.remove(*heh);
    }
    for eh in &edges_to_remove {
        brep.edges.remove(*eh);
    }
    for vh in &vertices_to_remove {
        brep.vertices.remove(*vh);
    }
}

/// Log any validation issues on a boolean result BRep.
///
/// Runs `validate_all()` and prints a summary to stderr. Does not fail the
/// operation — this is purely diagnostic to surface topological problems early.
fn log_validation_issues(result: &BRep, op: BoolOperation) {
    let validation_errors = result.validate_all();
    if !validation_errors.is_empty() {
        eprintln!(
            "[boolean] validation: {} issues in {} result",
            validation_errors.len(),
            match op {
                BoolOperation::Union => "union",
                BoolOperation::Intersection => "intersection",
                BoolOperation::Subtraction => "subtraction",
            }
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
