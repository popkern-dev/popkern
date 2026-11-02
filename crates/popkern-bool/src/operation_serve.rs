//! Stage 5 and the served-body tail of the boolean pipeline (split
//! from `operation.rs` under the file-size ratchet): cellular
//! selection (the complex primary, the cellular classifier its
//! fallback, the complex again as the rescue), pcurve recompute,
//! post-validation and post-processing — and the CLOSURE postcondition
//! on the served body (PK-555).
//!
//! Closure: when both operands are closed (no free half-edge, i.e. no
//! half-edge that is its own twin), the served body must be closed
//! too. A free half-edge that survives post-processing is an open
//! seam in the mesh (measured: the angled tapped hole's entry-chamfer
//! cut once kept 5 free half-edges — 215 open mesh edges and a split
//! chamfer cone). A body failing closure declines the boolean.
//!
//! Graph-merge selection was the last rung until plan Stage 7d
//! retired it: over two consecutive whole-suite runs it served no
//! boolean outside its own kill switch (`MergeFallback` lines still
//! mark the cellular declines it used to catch).

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::error::{BoolError, BoolOperation, PhaseMetric, Stage};
use crate::intersect::IntersectionResult;
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
#[allow(clippy::too_many_arguments)]
pub(super) fn select_and_serve(
    original_a: &BRep,
    original_b: &BRep,
    intersection_result: &IntersectionResult,
    op: BoolOperation,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
    progress: &mut Option<&mut dyn ProgressReporter>,
    policy: crate::options::NonManifoldPolicy,
) -> Result<BRep, BoolError> {
    // Stage 5: Cellular graph merge — order-independent fragment selection.
    // Every cellular decline is recorded structurally and said once on
    // stderr (`merge_fallback`) — the cellular error was once discarded
    // with `.ok()` on all builds, and a whole campaign's worth of
    // pcurve archaeology traced a defect to graph-merge before
    // discovering graph-merge was even running.
    let t_merge = popkern_core::time::Stopwatch::start();
    // Why the cellular stage declined: the boolean's error when no
    // later rung serves.
    let mut cellular_error: Option<BoolError> = None;
    let cellular_result = {
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
                crate::complex::serve::primary_source(),
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
                // meshed torn on the chamfered rod. Attempt the fill
                // HERE and decline when any face still misses one; the
                // accepted result's later shared recompute is then a
                // no-op.
                Ok(()) => {
                    let mut r = r;
                    let _ = r.recompute_pcurves_where_missing();
                    match super::first_missing_pcurve_face(&r) {
                        None => Some((r, carried)),
                        Some(fh) => {
                            let why = format!("face {fh:?} pcurve-incomplete after recompute");
                            cellular_error = Some(BoolError::ValidationFailed(why.clone()));
                            crate::merge_fallback::record_merge_fallback(
                                op,
                                crate::merge_fallback::MergeFallbackReason::CellularInvalid(why),
                            );
                            None
                        }
                    }
                }
                Err(e) => {
                    super::helpers_pcurve_coherence::trace_cellular_nonmanifold(&r);
                    cellular_error = Some(BoolError::ValidationFailed(format!(
                        "cellular result invalid: {e}"
                    )));
                    crate::merge_fallback::record_merge_fallback(
                        op,
                        crate::merge_fallback::MergeFallbackReason::CellularInvalid(e.to_string()),
                    );
                    None
                }
            },
            Err(diag) => {
                let e = diag.into_error();
                crate::merge_fallback::record_merge_fallback(
                    op,
                    crate::merge_fallback::MergeFallbackReason::CellularError(e.to_string()),
                );
                cellular_error = Some(e);
                None
            }
        }
    };
    // The cell complex (PK-976): where the cellular classifier could
    // not serve, the complex classifies the same graph and walk by
    // region membership. A non-manifold verdict declines with the
    // contacts (D1), or under the split policy the complex splits it
    // into manifold pieces.
    let mut nonmanifold: Option<BoolError> = None;
    let cellular_result = match cellular_result {
        None if crate::complex::serve::complex_rescue_enabled() => {
            match crate::complex::serve::rescue(
                original_a,
                original_b,
                intersection_result,
                op,
                tolerance,
            ) {
                Ok(served) => served,
                Err(error) => {
                    nonmanifold = Some(error);
                    // Under the split policy the complex splits its own
                    // verdict into manifold pieces (PK-987).
                    (policy == crate::options::NonManifoldPolicy::SplitToManifold)
                        .then(|| {
                            crate::complex::serve::rescue_split(
                                original_a,
                                original_b,
                                intersection_result,
                                op,
                                tolerance,
                            )
                        })
                        .and_then(|r| r.ok().flatten())
                }
            }
        }
        other => other,
    };
    let closed_operands = free_half_edges(original_a) == 0 && free_half_edges(original_b) == 0;
    // No rung served: the non-manifold verdict, else why the cellular
    // stage declined.
    let Some((result, carried)) = cellular_result else {
        return Err(nonmanifold
            .or(cellular_error)
            .unwrap_or_else(|| BoolError::MergeFailed {
                stage: "Cellular".to_string(),
                reason: "no selection served".to_string(),
            }));
    };
    // A cellular result arrives validated in full (its acceptance).
    let merged = (result, t_merge.elapsed());
    let served = serve(
        merged,
        (&carried, true),
        op,
        tolerance,
        original_a,
        original_b,
        metrics,
        progress,
    );
    let served = match (served, nonmanifold) {
        // The complex split of a non-manifold verdict (the split
        // policy): it serves once every contact is recorded and
        // `validate_all` accepts it.
        (Ok(body), Some(nm)) => {
            // The caller declines any served body that touches itself
            // (`execute_boolean_full`): only a contact-free body serves
            // as it is.
            let touches = || {
                !crate::postprocess::zero_thickness::zero_thickness_contacts(&body, tolerance)
                    .is_empty()
                    || !crate::postprocess::shell_contact::shell_contacts(&body, tolerance)
                        .is_empty()
            };
            if body.validate_all().is_empty() && !touches() {
                return Ok(body);
            }
            return match policy {
                crate::options::NonManifoldPolicy::Decline => Err(nm),
                crate::options::NonManifoldPolicy::SplitToManifold => {
                    crate::postprocess::contact_record::accept_split(body, tolerance).map_err(|e| {
                        crate::complex_shadow::record_complex_rescue(
                            op,
                            crate::complex_shadow::ComplexRescueOutcome::SplitRejected(
                                e.to_string(),
                            ),
                        );
                        nm
                    })
                }
            };
        }
        (Err(_), Some(nm)) => return Err(nm),
        (served, None) => served?,
    };
    let open = free_half_edges(&served);
    if !closed_operands || open == 0 {
        return Ok(served);
    }
    crate::merge_fallback::record_merge_fallback(
        op,
        crate::merge_fallback::MergeFallbackReason::CellularInvalid(format!(
            "served body carries {open} free half-edge(s) after post-processing"
        )),
    );
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
    // A certified empty result (the gate above declines the rest) is
    // served as a bare `BRep::new`, like every other empty path — not
    // as the builder's solid and shell with nothing in them (a
    // line-contact intersection served one shell of zero faces,
    // PK-966).
    if result.faces.is_empty() {
        return Ok(BRep::new(*tolerance));
    }

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

    // The output repairs can remove the last faces (a zero-height
    // sliver ring is all a line contact keeps): the same certificate
    // decides whether that emptiness is the answer.
    if result.faces.is_empty() {
        super::empty::certify_empty_result(op, original_a, original_b)?;
        return Ok(BRep::new(*tolerance));
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
