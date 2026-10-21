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
#[path = "operation_contact_rim.rs"]
mod contact_rim;
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
use crate::error::{BoolDiagnostics, BoolError, BoolOperation, PhaseMetric, Stage};
use crate::intersect::find_intersections;
use crate::progress::{ProgressEvent, ProgressReporter};

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
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    execute_boolean_full(op, a, b, tolerance, &mut Vec::new(), None, None)
}

/// Execute a boolean and populate `metrics` with per-stage timing samples.
///
/// Backwards-compatible wrapper used by `boolean_with_diagnostics`.
pub(crate) fn execute_boolean_with_metrics(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
) -> Result<BRep, BoolError> {
    execute_boolean_full(op, a, b, tolerance, metrics, None, None)
}

/// The staged boolean pipeline body (intersect → split → classify →
/// select → merge → validate). Runs the operands exactly as given,
/// borrowed: only the no-intersection shortcut, which combines them,
/// takes copies (PK-815: the late cut of a 255-face disc copied it twice
/// per boolean); the mixed-sense rescue orchestration lives in
/// [`sense_rescue::execute_boolean_full`].
pub(super) fn execute_boolean_pipeline(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
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
        find_intersections(a, b, tolerance).map_err(|diag| diag.into_error())?;
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
            a.clone(),
            b.clone(),
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
    if containment::coincident_containment_applies(a, b, &intersection_result, tolerance) {
        return handle_no_intersection(
            op,
            a.clone(),
            b.clone(),
            &intersection_result.aabbs_a,
            &intersection_result.aabbs_b,
            tolerance,
        );
    }

    // Bodies that only touch (opposite contact faces, point records).
    if let Some(result) = contact::touching_only(op, a, b, &intersection_result, tolerance) {
        return Ok(result);
    }

    // The cellular stage serves from the operands and the
    // intersection. A face-level split, classification and selection
    // used to run first, and fed nothing but an empty-selection
    // containment shortcut the cellular stage answers by itself: with
    // them skipped the workspace suite (10101 tests) passed unchanged,
    // and a 64-slot disc's late cut took 8.0 ms instead of 10.2 (the
    // first 3.9 instead of 5.0; PK-814).
    serve::select_and_serve(
        a,
        b,
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
