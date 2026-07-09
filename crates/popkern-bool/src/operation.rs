//! Pipeline orchestrator: connects all 5 stages of the boolean pipeline.
//!
//! The [`execute_boolean`] function drives intersection, split, classify,
//! select, and merge in sequence, handling the no-intersection shortcut
//! and accumulating diagnostics on failure.

#[path = "operation_helpers.rs"]
pub(crate) mod helpers;

#[path = "operation_postprocess.rs"]
mod postprocess;

use postprocess::{
    log_post_processing_diagnostics, log_validation_issues, propagate_output_tolerances,
    remove_redundant_edges_safe,
};

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
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        eprintln!(
            "  [exec] execute_boolean {op:?} ({} vs {} faces)",
            a.faces.len(),
            b.faces.len()
        );
    }

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
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        eprintln!(
            "  [exec] {} records, {} overlaps",
            intersection_result.records.len(),
            intersection_result.coincidence.face_overlaps.len()
        );
    }
    if intersection_result.records.is_empty() {
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
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

    // Coincident-body shortcut. The old gate was a FACE-COUNT
    // heuristic — "overlaps ≥ half the smaller body's faces ⇒ full
    // containment" — and it silently swallowed a flush-sided gusset
    // whose 2 coplanar caps tripped the count on a 5-face tool: the
    // union returned the bare arm, watertight and χ-clean, minus
    // 44 cm³ of material. Containment is now VERIFIED geometrically:
    // every vertex of the candidate body must classify inside-or-on
    // the other. Bodies that merely SHARE some coplanar faces fall
    // through to the full pipeline.
    if !intersection_result.coincidence.face_overlaps.is_empty() {
        let n_a = a.faces.len();
        let n_b = b.faces.len();
        let n_overlaps = intersection_result.coincidence.face_overlaps.len();
        if n_overlaps >= n_a.min(n_b) / 2 {
            // Probe FACE INTERIOR samples, not vertices: corner
            // vertices sit on several faces of the other body at
            // once and the ray votes degrade (two identical cubes
            // misclassified their own corners as Outside). A face
            // interior point of a genuinely contained body is
            // Inside-or-On the other; a poking face has interior
            // points strictly Outside.
            // Each probe is offset INWARD from its face sample (the
            // Newell winding of the outer loop gives the outward
            // side): an on-boundary sample of a coincident face has
            // ill-defined ray parity, but its inward offset is
            // strictly interior to the candidate — Inside the other
            // body iff containment holds. Any face that cannot be
            // probed declines the shortcut (conservative: the full
            // pipeline runs).
            let contained = |inner: &popkern_core::topology::BRep,
                             outer: &popkern_core::topology::BRep,
                             outer_aabbs: &[(
                popkern_core::topology::FaceHandle,
                popkern_core::geometry::intersection::surface_subdivision::Aabb,
            )]|
             -> bool {
                let eps = tolerance.linear() * 50.0;
                for (fh, face) in inner.faces.iter() {
                    let sample = crate::classify::pick_sample_point(inner, fh);
                    // Outward normal from the outer loop's Newell
                    // vector (loops wind CCW seen from outside).
                    let Some(loop_data) = inner.loops.get(face.outer_loop) else {
                        return false;
                    };
                    let start = loop_data.half_edge;
                    if !inner.half_edges.contains_key(start) {
                        return false;
                    }
                    let mut pts = Vec::new();
                    let mut he = start;
                    loop {
                        pts.push(inner.vertices[inner.half_edges[he].origin].position);
                        he = inner.half_edges[he].next;
                        if he == start {
                            break;
                        }
                    }
                    if pts.len() < 3 {
                        return false;
                    }
                    let mut n = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
                    for i in 0..pts.len() {
                        let p = pts[i];
                        let q = pts[(i + 1) % pts.len()];
                        n.x += (p.y - q.y) * (p.z + q.z);
                        n.y += (p.z - q.z) * (p.x + q.x);
                        n.z += (p.x - q.x) * (p.y + q.y);
                    }
                    let len = n.length();
                    if len < 1.0e-12 {
                        return false;
                    }
                    // Measured convention: this kernel's loop Newell
                    // vector points INWARD (a unit cube's bottom
                    // face yielded +z), so the interior probe adds.
                    let probe = sample + n * (eps / len);
                    let verdict = crate::classify::point_in_solid_final(
                        &probe,
                        outer,
                        outer_aabbs,
                        tolerance,
                    );
                    if std::env::var_os("PK_CONTAIN_TRACE").is_some() {
                        eprintln!(
                            "[contain] {fh:?} probe ({:.4},{:.4},{:.4}) -> {verdict:?}",
                            probe.x, probe.y, probe.z
                        );
                    }
                    match verdict {
                        Ok(crate::error::PointClassification::Inside) => {}
                        Ok(_) => return false,
                        Err(_) => return false,
                    }
                }
                true
            };
            let a_in_b = n_a <= n_b && contained(&a, &b, &intersection_result.aabbs_b);
            let b_in_a = !a_in_b && n_b <= n_a && contained(&b, &a, &intersection_result.aabbs_a);
            if a_in_b || b_in_a {
                return handle_no_intersection(
                    op,
                    a,
                    b,
                    &intersection_result.aabbs_a,
                    &intersection_result.aabbs_b,
                    tolerance,
                );
            } else if std::env::var_os("PK_BOOL_TRACE").is_some() {
                eprintln!(
                    "[boolean] coincident shortcut declined: {n_overlaps} overlaps but \
                     containment not verified — full pipeline"
                );
            }
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
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!("  [exec] -> empty-selection containment path");
        }
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
        // Explicit fallback boundary: record why the stronger
        // (order-independent cellular) path didn't apply. Debug-only
        // logging hid this in release — a whole campaign's worth of
        // pcurve archaeology traced a defect to graph-merge before
        // discovering graph-merge was even running.
        match r.validate() {
            Ok(()) => true,
            Err(e) => {
                if cfg!(debug_assertions) || std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "[boolean] cellular assembly invalid ({e}), falling back to graph-merge"
                    );
                }
                if std::env::var_os("PK_CELL_TRACE").is_some() {
                    // Dump every edge's use count + position for
                    // non-manifold diagnosis.
                    for (eh, edge) in r.edges.iter() {
                        let mut uses = 0;
                        for (_, he) in r.half_edges.iter() {
                            if he.edge == eh && r.loops.contains_key(he.loop_) {
                                uses += 1;
                            }
                        }
                        if uses != 2 {
                            use popkern_core::geometry::curve::CurveEval;
                            let pos = edge
                                .curve
                                .as_ref()
                                .map(|c| {
                                    let p = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                                    format!("({:.2},{:.2},{:.2})", p.x, p.y, p.z)
                                })
                                .unwrap_or_else(|| {
                                    let p =
                                        r.vertices[r.half_edges[edge.half_edge].origin].position;
                                    format!("v({:.2},{:.2},{:.2})", p.x, p.y, p.z)
                                });
                            let mut endpoints = String::new();
                            {
                                use std::fmt::Write as _;
                                let he0 = edge.half_edge;
                                if r.half_edges.contains_key(he0) {
                                    let a = r.vertices[r.half_edges[he0].origin].position;
                                    let tw = r.half_edges[he0].twin;
                                    if r.half_edges.contains_key(tw) {
                                        let b = r.vertices[r.half_edges[tw].origin].position;
                                        let _ = write!(
                                            endpoints,
                                            " span ({:.1},{:.1},{:.1})->({:.1},{:.1},{:.1})",
                                            a.x, a.y, a.z, b.x, b.y, b.z
                                        );
                                    }
                                }
                            }
                            eprintln!("  [nm] {eh:?} uses={uses} at {pos}{endpoints}");
                        }
                    }
                }
                false
            }
        }
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
