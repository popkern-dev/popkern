//! Pipeline orchestrator: connects all 5 stages of the boolean pipeline.
//!
//! The [`execute_boolean`] function drives intersection, split, classify,
//! select, and merge in sequence, handling the no-intersection shortcut
//! and accumulating diagnostics on failure.

#[path = "operation_helpers.rs"]
pub(crate) mod helpers;
#[path = "operation_pcurve_coherence.rs"]
mod helpers_pcurve_coherence;
#[path = "operation_no_intersection.rs"]
pub(crate) mod no_intersection;

#[path = "operation_contact.rs"]
mod contact;
#[path = "operation_containment.rs"]
mod containment;
#[path = "operation_postprocess.rs"]
mod postprocess;
#[path = "operation_postprocess_gate.rs"]
mod postprocess_gate;
#[path = "operation_serve.rs"]
mod serve;

use postprocess::log_validation_issues;

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::cancellation::CancellationToken;
use crate::classify::classify_fragments;
use crate::error::{BoolDiagnostics, BoolError, BoolOperation, PhaseMetric, Stage};
use crate::intersect::find_intersections;
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
    duration: Option<std::time::Duration>,
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

use self::helpers::{handle_empty_selection, verify_whole_body_containment};
use self::no_intersection::handle_no_intersection;

#[path = "operation_sense_rescue.rs"]
mod sense_rescue;
pub(crate) use sense_rescue::execute_boolean_full;

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

/// The staged boolean pipeline body (intersect → split → classify →
/// select → merge → validate). Runs the operands exactly as given;
/// the mixed-sense rescue orchestration lives in
/// [`sense_rescue::execute_boolean_full`].
pub(super) fn execute_boolean_pipeline(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
    cancel: Option<&CancellationToken>,
    mut progress: Option<&mut dyn ProgressReporter>,
) -> Result<BRep, BoolError> {
    check_cancel(cancel)?;
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [exec] execute_boolean {op:?} ({} vs {} faces)",
            a.faces.len(),
            b.faces.len()
        );
    }

    // Stage 1: Intersect — compute face-face intersection curves.
    let t_intersect = popkern_core::time::Stopwatch::start();
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
    if popkern_core::trace::trace_enabled("boolean") {
        helpers::trace_records(
            &intersection_result.records,
            intersection_result.coincidence.face_overlaps.len(),
        );
    }
    if intersection_result.records.is_empty() {
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!("  [exec] -> no-records containment path");
        }
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

    // Coincident-body containment shortcut — verified geometrically
    // (see `operation_containment.rs`).
    if containment::coincident_containment_applies(&a, &b, &intersection_result, tolerance) {
        return handle_no_intersection(
            op,
            a,
            b,
            &intersection_result.aabbs_a,
            &intersection_result.aabbs_b,
            tolerance,
        );
    }

    // Bodies that only touch (opposite contact faces, point records).
    if let Some(result) = contact::touching_only(op, &a, &b, &intersection_result, tolerance) {
        return Ok(result);
    }

    // Clone originals for graph merge (graph merge reads originals, not split).
    let original_a = a.clone();
    let original_b = b.clone();

    // Stage 2: Split — trim faces along intersection curves.
    let t_split = popkern_core::time::Stopwatch::start();
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
    // The classified-against body is the SPLIT operand, so its face
    // AABBs must be computed from the split body too: the pre-split
    // operand boxes collide by handle with the split faces and pair a
    // face with another face's box (latent until the per-hit face-box
    // gate cross-referenced them — a contained-cubes intersection
    // rejected every real hit and served empty, PK-464).
    let t_classify = popkern_core::time::Stopwatch::start();
    let split_aabbs_a = crate::intersect::compute_face_aabbs(&split_result.a);
    let split_aabbs_b = crate::intersect::compute_face_aabbs(&split_result.b);
    let classifications_a = classify_fragments(
        &split_result.a,
        &split_result.b,
        &split_aabbs_b,
        &split_result.fragments_a.intersection_edges,
        tolerance,
    )?;

    let classifications_b = classify_fragments(
        &split_result.b,
        &split_result.a,
        &split_aabbs_a,
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
    let t_select = popkern_core::time::Stopwatch::start();
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
    //
    // NOT when coincident face overlaps exist: the face-level
    // classifier's samples sit in the other body's boundary plane
    // there, its parity verdicts are coin flips, and an all-drop
    // selection is noise, not containment (measured: a dovetail-rail
    // cutter whose void floor is flush with the plate underside
    // emptied the whole tip piece). The cellular stage below carries
    // the coincident-aware offset-pair classification — let it run.
    if selection.keep_a.is_empty()
        && selection.keep_b.is_empty()
        && intersection_result.coincidence.face_overlaps.is_empty()
    {
        // The inferred containment must be VERIFIED against the
        // original operands before it is served: the classifications
        // come from one seed per propagation component, and an
        // under-split operand fuses both sides of the trim into one
        // component whose single seed wipes everything (measured: a
        // rim-straddling key subtract served an EMPTY body, PK-173).
        // A refuted claim falls through to the cellular stage, which
        // classifies per fragment.
        let all_a_inside = !classifications_a.is_empty()
            && classifications_a
                .values()
                .all(|c| *c == crate::error::Classification::Inside);
        let all_b_inside = !classifications_b.is_empty()
            && classifications_b
                .values()
                .all(|c| *c == crate::error::Classification::Inside);
        let claim_verified = if all_a_inside {
            verify_whole_body_containment(&original_a, &original_b)
        } else if all_b_inside {
            verify_whole_body_containment(&original_b, &original_a)
        } else {
            false
        };
        if claim_verified {
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!("  [exec] -> empty-selection containment path (verified)");
            }
            let result =
                handle_empty_selection(op, split_result, &classifications_a, &classifications_b);
            log_validation_issues(&result, op);
            return Ok(result);
        }
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!("  [exec] empty-selection containment REFUTED — continuing to cellular");
        }
    }

    serve::select_and_serve(
        &original_a,
        &original_b,
        &intersection_result,
        op,
        tolerance,
        metrics,
        &mut progress,
    )
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

/// First face with surface geometry whose loop half-edges are missing
/// a pcurve — the cellular acceptance gate's completeness scan
/// (PK-463). `None` = trim inventory complete.
fn first_missing_pcurve_face(r: &BRep) -> Option<popkern_core::topology::FaceHandle> {
    for (fh, fd) in &r.faces {
        if fd.surface.is_none() {
            continue;
        }
        for lh in std::iter::once(fd.outer_loop).chain(fd.inner_loops.iter().copied()) {
            let Some(ld) = r.loops.get(lh) else { continue };
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..r.half_edges.len() + 1 {
                let Some(h) = r.half_edges.get(he) else { break };
                if h.pcurve.is_none() {
                    return Some(fh);
                }
                he = h.next;
                if he == start {
                    break;
                }
            }
        }
    }
    None
}

#[cfg(test)]
#[path = "operation_tests.rs"]
mod tests;
