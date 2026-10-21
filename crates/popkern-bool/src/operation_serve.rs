//! Stage 5 and the served-body tail of the boolean pipeline (split
//! from `operation.rs` under the file-size ratchet): cellular
//! selection with its graph-merge fallback, pcurve recompute,
//! post-validation and post-processing — and the CLOSURE postcondition
//! on the served body (PK-555).
//!
//! Closure: when both operands are closed (no free half-edge, i.e. no
//! half-edge that is its own twin), the served body must be closed
//! too. A free half-edge that survives post-processing is an open
//! seam in the mesh (measured: the angled tapped hole's entry-chamfer
//! cut, served by the cellular pipeline once the tapped wall's
//! pcurves were exact, kept 5 free half-edges — 215 open mesh edges
//! and a split chamfer cone — while graph-merge serves the same cut
//! closed). A cellular body failing closure is discarded for the
//! graph-merge body (recorded as a fallback); a graph-merge body
//! failing it declines the boolean.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::error::{BoolError, BoolOperation, PhaseMetric, Stage};
use crate::intersect::IntersectionResult;
use crate::merge::graph_merge::graph_merge;
use crate::progress::ProgressReporter;

use super::postprocess::{
    log_post_processing_diagnostics, log_validation_issues, run_output_postprocess,
};

/// Half-edges that are their own twin (free boundary half-edges).
pub(super) fn free_half_edges(brep: &BRep) -> usize {
    brep.half_edges
        .iter()
        .filter(|(handle, he)| he.twin == *handle)
        .count()
}

/// Select fragments (cellular first, graph-merge on fallback) and
/// serve the result through the post-merge tail and the closure
/// postcondition (see the module docs).
pub(super) fn select_and_serve(
    original_a: &BRep,
    original_b: &BRep,
    intersection_result: &IntersectionResult,
    op: BoolOperation,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
    progress: &mut Option<&mut dyn ProgressReporter>,
) -> Result<BRep, BoolError> {
    // Stage 5: Cellular graph merge — order-independent fragment selection.
    // Falls back to standard graph-merge if cellular fails or produces
    // invalid topology. Every fallback firing is recorded structurally
    // and said once on stderr (`merge_fallback`) — the cellular error
    // was previously discarded with `.ok()` on all builds, and a whole
    // campaign's worth of pcurve archaeology traced a defect to
    // graph-merge before discovering graph-merge was even running.
    let t_merge = popkern_core::time::Stopwatch::start();
    // POPKERN_CELLULAR=0 forces the graph-merge path — diagnostic
    // kill-switch (PK-239 precedent).
    let cellular_enabled = std::env::var_os("POPKERN_CELLULAR").is_none_or(|v| v != "0");
    let cellular_result = if cellular_enabled {
        // Local scope (PK-814): far faces carried, the rest merged; any
        // scoped failure or invalid result reruns unscoped.
        let merge = |scopes: Option<&crate::scope::OperandScopes>| {
            crate::cellular::boolean::cellular_graph_merge(
                original_a,
                original_b,
                intersection_result,
                op,
                tolerance,
                scopes,
            )
        };
        let attempt = |scopes: &crate::scope::OperandScopes| match merge(Some(scopes)) {
            Ok((r, carried)) if r.validate().is_ok() => Some(Ok((r, carried))),
            outcome => {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [scope] scoped run declined (holes {}): {:?}",
                        scopes.a.as_ref().map_or(0, |s| s.far_holes.len()),
                        outcome.err().map(|d| d.into_error().to_string())
                    );
                }
                None
            }
        };
        // Hole scope, then face scope, then unscoped.
        let scopes = crate::scope::serving_scopes(
            original_a,
            original_b,
            intersection_result,
            op,
            tolerance,
        );
        let scoped = (scopes.a.is_some() || scopes.b.is_some())
            .then(|| attempt(&scopes).or_else(|| attempt(&scopes.faces_only()?)))
            .flatten();
        // A scoped result arrives validated (its validity decided the
        // fallback); the arm below does not validate it twice.
        let prevalidated = scoped.is_some();
        match scoped.unwrap_or_else(|| merge(None)) {
            Ok((r, _))
                if let Some(why) = crate::scope::compare_violation(
                    (original_a, original_b),
                    intersection_result,
                    op,
                    tolerance,
                    &r,
                ) =>
            {
                return Err(BoolError::MergeFailed {
                    stage: "ScopeCompare".to_string(),
                    reason: format!("local scope premise violated: {why}"),
                });
            }
            Ok((r, carried)) => match if prevalidated { Ok(()) } else { r.validate() } {
                // Pcurve completeness is part of the acceptance
                // (PK-463): the cellular assembly leaves record-use
                // pcurves to the recompute, and a recompute skip
                // (projection failure) is silent — the served body
                // meshed torn on the chamfered rod while graph-merge
                // served it watertight. Attempt the fill HERE and
                // decline to graph-merge when any face still misses
                // one; the accepted result's later shared recompute
                // is then a no-op.
                Ok(()) => {
                    let mut r = r;
                    let _ = r.recompute_pcurves_where_missing();
                    match super::first_missing_pcurve_face(&r) {
                        None => Some((r, carried)),
                        Some(fh) => {
                            crate::merge_fallback::record_merge_fallback(
                                op,
                                crate::merge_fallback::MergeFallbackReason::CellularInvalid(
                                    format!("face {fh:?} pcurve-incomplete after recompute"),
                                ),
                            );
                            None
                        }
                    }
                }
                Err(e) => {
                    super::helpers_pcurve_coherence::trace_cellular_nonmanifold(&r);
                    crate::merge_fallback::record_merge_fallback(
                        op,
                        crate::merge_fallback::MergeFallbackReason::CellularInvalid(e.to_string()),
                    );
                    None
                }
            },
            Err(diag) => {
                crate::merge_fallback::record_merge_fallback(
                    op,
                    crate::merge_fallback::MergeFallbackReason::CellularError(
                        diag.into_error().to_string(),
                    ),
                );
                None
            }
        }
    } else {
        crate::merge_fallback::record_merge_fallback(
            op,
            crate::merge_fallback::MergeFallbackReason::CellularDisabled,
        );
        None
    };
    let closed_operands = free_half_edges(original_a) == 0 && free_half_edges(original_b) == 0;
    let graph = || {
        graph_merge(original_a, original_b, intersection_result, op, tolerance)
            .map_err(|diag| diag.into_error())
    };
    let from_cellular = cellular_result.is_some();
    // A cellular result arrives validated in full (its acceptance);
    // a graph-merge one does not.
    let (result, carried, validated) = match cellular_result {
        Some((r, carried)) => (r, carried, true),
        None => (graph()?, crate::serve_region::Carried::default(), false),
    };
    let merged = (result, t_merge.elapsed());
    let served = serve(
        merged,
        (&carried, validated),
        op,
        tolerance,
        original_a,
        original_b,
        metrics,
        progress,
    )?;
    let open = free_half_edges(&served);
    if !closed_operands || open == 0 {
        return Ok(served);
    }
    if from_cellular {
        crate::merge_fallback::record_merge_fallback(
            op,
            crate::merge_fallback::MergeFallbackReason::CellularInvalid(format!(
                "served body carries {open} free half-edge(s) after post-processing"
            )),
        );
        let t_retry = popkern_core::time::Stopwatch::start();
        let retried = graph()?;
        let merged = (retried, t_retry.elapsed());
        let served = serve(
            merged,
            (&crate::serve_region::Carried::default(), false),
            op,
            tolerance,
            original_a,
            original_b,
            metrics,
            progress,
        )?;
        let open = free_half_edges(&served);
        if open == 0 {
            return Ok(served);
        }
        return Err(unclosed(open));
    }
    Err(unclosed(open))
}

fn unclosed(open: usize) -> BoolError {
    BoolError::ValidationFailed(format!(
        "result of two closed operands carries {open} free half-edge(s) after \
         post-processing — declining rather than serving an open body"
    ))
}

/// The post-merge tail on a selected body: pcurve recompute,
/// post-validation, post-processing and stub repair.
#[allow(clippy::too_many_arguments)]
fn serve(
    (result, elapsed_merge): (BRep, Option<std::time::Duration>),
    (carried, validated): (&crate::serve_region::Carried, bool),
    op: BoolOperation,
    tolerance: &Tolerance,
    original_a: &BRep,
    original_b: &BRep,
    metrics: &mut Vec<PhaseMetric>,
    progress: &mut Option<&mut dyn ProgressReporter>,
) -> Result<BRep, BoolError> {
    let count_merge = result.faces.len();
    metrics.push(PhaseMetric {
        stage: Stage::Merge,
        duration: elapsed_merge,
        entity_count: count_merge,
    });
    super::report_phase(progress, Stage::Merge, 4, elapsed_merge, count_merge);

    // Recompute pcurves so that the tessellation pipeline can build
    // accurate trim loops. Boolean merge may leave seam edges and
    // original-operand edges without pcurves — only faces with a gap
    // in their pcurve inventory are recomputed (PK-274): a face the
    // merge carried over intact keeps its authoritative pcurves.
    // Pcurve edits keep a validity verdict: `validate()` reads none.
    let mut gate = if validated {
        super::postprocess_gate::ValidityGate::known(true)
    } else {
        super::postprocess_gate::ValidityGate::new()
    };
    let mut result = result;
    if let Err(e) = result.recompute_pcurves_where_missing()
        && popkern_core::trace::trace_enabled("boolean")
    {
        eprintln!("[boolean] pcurve recompute declined: {e}");
    }
    super::helpers_pcurve_coherence::normalize_output_loop_charts(&mut result);

    let region = crate::serve_region::ServeRegion::of(&result, carried);
    gate.edited(super::helpers::pre_check_retire_and_empty_gate(
        &mut result,
        op,
        (original_a, original_b),
        tolerance,
        &region,
    )?);

    // §12.5 / §7.3: gated post-validation. Run the boolean-specific
    // post-check (manifold, Euler, half-edge cycles, watertightness).
    // Surface any errors as a structured failure rather than silently
    // returning a corrupted result. The merge-phase cut-behavior fix
    // (commit history) made this consistently green; if a future
    // change regresses, fail loudly.
    let t_validate = popkern_core::time::Stopwatch::start();
    let pc_errors: Vec<String> = super::helpers::validate_with_junction_reconcile_retry(
        &mut result,
        tolerance,
        &carried.edges,
        &mut gate,
    );
    let elapsed_validate = t_validate.elapsed();
    let count_validate = result.faces.len();
    metrics.push(PhaseMetric {
        stage: Stage::Validation,
        duration: elapsed_validate,
        entity_count: count_validate,
    });
    super::report_phase(
        progress,
        Stage::Validation,
        5,
        elapsed_validate,
        count_validate,
    );
    if !pc_errors.is_empty() {
        super::helpers::trace_invalid_result(&result);
        return Err(BoolError::ValidationFailed(pc_errors.join("; ")));
    }
    run_output_postprocess(&mut result, tolerance, gate, region, carried);
    super::postprocess::decline_on_unrepairable_stubs(&mut result, tolerance)?;
    crate::postprocess::curveless_edge::repair_or_decline_curveless_edges(&mut result, tolerance)?;

    // Pcurve COHERENCE sweep for doubly-periodic (torus) faces —
    // mechanics and rationale in `operation_pcurve_coherence.rs`.
    super::helpers_pcurve_coherence::recompute_incoherent_periodic_pcurves(&mut result);
    // The served body's faces must bound regions of their charts
    // (PK-836, PK-862): coincident cooriented cylinder faces once served
    // with loops winding round the period (validate() passes them, the
    // body measured 16587 mm^3 of nothing). Checked here, after the
    // output repairs re-register pcurves; a face that still winds gets
    // its pcurves re-derived once, else the boolean declines.
    super::helpers_pcurve_coherence::gate_loop_period_winding(&mut result)
        .map_err(BoolError::ValidationFailed)?;

    // Fitted trims record their measured agreement on the edge (PK-826);
    // edges carried over unsplit keep their operands' (scoped: the
    // whole-body measurement made a late slot cut 25x the first).
    let raised = result.record_measured_edge_tolerances_where(|e| !carried.edges.contains_key(e));
    // An edge whose trims miss its curve past the model's ceiling fails
    // `validate()`: decline instead of serving it (PK-911: a mixed-sense
    // disc ∩ inward-boss served four record edges whose pcurves ran
    // against their curves, 4.0 off). The decline is what lets the
    // mixed-sense rescue retry with the inward operand reversed.
    if let Some((edge, measured)) = raised
        .iter()
        .filter_map(|&e| Some((e, result.edges.get(e)?.tolerance?)))
        .find(|&(_, t)| t > result.max_imprecision)
    {
        return Err(BoolError::ValidationFailed(format!(
            "edge {edge:?}: its curve and trims disagree by {measured:.3e}, past the model's \
             max imprecision {:.1e}; declining rather than serving it",
            result.max_imprecision
        )));
    }

    // Debug/the boolean trace topic diagnostics on the SERVED result. This runs
    // after run_output_postprocess deliberately (PK-272): the earlier
    // placement fired before reconcile_orientation_flags and reported
    // transient winding findings that the served result doesn't have.
    log_validation_issues(&result, op);

    // Post-processing diagnostics: sliver detection and watertightness.
    log_post_processing_diagnostics(&result, tolerance);

    Ok(result)
}
