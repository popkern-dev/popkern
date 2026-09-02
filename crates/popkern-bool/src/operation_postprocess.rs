//! Post-processing helpers for the boolean pipeline: validation
//! logging, tolerance propagation, redundant-edge removal, and
//! sliver/watertightness diagnostics.
//!
//! Split out of `operation.rs` to keep the parent under the
//! 500-line structural cap.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::error::BoolOperation;

/// Log any validation issues on a boolean result BRep.
///
/// Runs `validate_all()` and prints a summary to stderr. Does not fail the
/// operation — this is purely diagnostic to surface topological problems early.
/// Catches panics from validation on corrupted topology.
///
/// Debug/the `boolean` trace topic-only (PK-272): `validate_all()` includes the
/// O(n²) surface-surface self-intersection sweep, and running it on
/// every release boolean purely to print stderr was measurable
/// overhead the caller can't act on. The correctness gate for boolean
/// results is `validate::post_check` in `operation.rs`, which runs
/// unconditionally either way. On the main pipeline this is called
/// after `run_output_postprocess` (so it describes the served result,
/// not a transient pre-reconcile state); the no-intersection paths
/// serve their result directly and log it as-is.
pub(crate) fn log_validation_issues(result: &BRep, op: BoolOperation) {
    if !(cfg!(debug_assertions) || popkern_core::trace::trace_enabled("boolean")) {
        return;
    }
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
            const MAX_LISTED: usize = 8;
            for error in errors.iter().take(MAX_LISTED) {
                eprintln!("[boolean] validation:   {error}");
            }
            if errors.len() > MAX_LISTED {
                eprintln!(
                    "[boolean] validation:   ... {} more suppressed",
                    errors.len() - MAX_LISTED,
                );
            }
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
pub(crate) fn propagate_output_tolerances(result: &mut BRep, tolerance: &Tolerance) {
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

/// Safely remove redundant edges, validating after each removal.
///
/// Finds redundant edges and removes them one at a time via KEF.
/// After each removal, validates the topology. If validation fails,
/// stops immediately (leaving remaining redundant edges in place).
pub(crate) fn remove_redundant_edges_safe(result: &mut BRep, tolerance: &Tolerance) {
    use crate::postprocess::redundant_edge::find_redundant_edges;

    let report = find_redundant_edges(result, tolerance);

    if report.redundant_edges.is_empty() {
        dissolve_redundant_chains(result, tolerance);
        return;
    }

    let mut removed = 0;
    for eh in &report.redundant_edges {
        // Re-check against the CURRENT topology: the report is
        // computed once, and earlier KEFs in this loop change both
        // which edges two faces share (a stale single-share verdict
        // dissolves one edge of a multi-edge boundary and leaves the
        // survivor self-twinned) and how far a merged face spans the
        // chart (a merge that was span-safe at report time can wrap
        // a period by the time its edge comes up).
        if !crate::postprocess::redundant_edge::edge_redundant_now(result, *eh, tolerance) {
            continue;
        }
        if !crate::clip::chart_cell::edge_dissolve_keeps_cell_resolvable(result, *eh, tolerance) {
            continue;
        }
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
                if cfg!(debug_assertions) {
                    eprintln!(
                        "[boolean] KEF corrupted topology after removing edge, stopping. \
                         Before: V={pre_v} E={pre_e} F={pre_f}"
                    );
                }
                break;
            }
            removed += 1;
        }
    }

    if removed > 0 && cfg!(debug_assertions) {
        eprintln!("[boolean] post-processing: removed {removed} redundant edge(s)");
    }
    dissolve_redundant_chains(result, tolerance);
}

/// Merge co-surfacial face pairs whose shared boundary is a multi-edge
/// PATH (PK-379): staged per pair, kept only when the body still
/// validates. Single-edge pairs were handled above; pairs whose
/// shared edges are not one path stay split (a ring would need an
/// inner loop, not a merge).
fn dissolve_redundant_chains(result: &mut BRep, tolerance: &Tolerance) {
    use crate::postprocess::redundant_chain::{dissolve_chain, shared_coplanar_chain};
    let mut merged = 0usize;
    loop {
        let mut progressed = false;
        let edges: Vec<_> = result.edges.keys().collect();
        for eh in edges {
            let Some(chain) = shared_coplanar_chain(result, eh, tolerance) else {
                continue;
            };
            if chain.len() < 2 {
                continue;
            }
            if !crate::clip::chart_cell::edge_dissolve_keeps_cell_resolvable(
                result, chain[0], tolerance,
            ) {
                continue;
            }
            let mut staged = result.clone();
            if dissolve_chain(&mut staged, &chain).is_ok() && staged.validate().is_ok() {
                *result = staged;
                merged += 1;
                progressed = true;
                break;
            }
        }
        if !progressed {
            break;
        }
    }
    if merged > 0 && cfg!(debug_assertions) {
        eprintln!("[boolean] post-processing: dissolved {merged} redundant shared chain(s)");
    }
}

/// Log post-processing diagnostics: sliver detection and watertightness.
///
/// Runs sliver analysis on each face and watertightness validation on all
/// edges. Results are logged to stderr (diagnostic only, does not fail).
///
/// Debug/the `boolean` trace topic-only (PK-272): before the gate, release
/// builds computed the full sliver and watertightness sweeps and then
/// discarded them — the eprintlns were `#[cfg(debug_assertions)]` but
/// the computation was not. Now nothing is computed unless the output
/// can be observed.
pub(crate) fn log_post_processing_diagnostics(result: &BRep, tolerance: &Tolerance) {
    use crate::clip::face_vertex_positions;
    use crate::imprint::shared_geometry::validate_watertightness;
    use crate::imprint::sliver::{SliverThresholds, analyze_face_polygon};

    if !(cfg!(debug_assertions) || popkern_core::trace::trace_enabled("boolean")) {
        return;
    }

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
        eprintln!("[boolean] post-processing: {sliver_count} sliver face(s) detected");
    }

    // Watertightness validation
    let wt_report = validate_watertightness(result, tolerance.linear() * 10.0, 3);
    if !wt_report.is_watertight() {
        eprintln!(
            "[boolean] post-processing: {} watertightness gap(s), max deviation: {:.2e}",
            wt_report.gaps.len(),
            wt_report.max_gap()
        );
    }
}

/// the boolean trace topic-gated face/edge census line for the post-processing
/// stages — sliver removal silently deleted a cylinder wall once
/// (PK-150); the per-stage census makes any face-count change in the
/// post pipeline visible.
pub(crate) fn trace_face_census(label: &str, brep: &BRep) {
    if popkern_core::trace::trace_enabled("boolean") {
        let self_twins = brep
            .half_edges
            .iter()
            .filter(|(heh, he)| he.twin == *heh)
            .count();
        eprintln!(
            "  [post] {label}: {} faces {} edges {} self-twinned",
            brep.faces.len(),
            brep.edges.len(),
            self_twins
        );
    }
}

/// Output post-processing pipeline shared by every boolean opcode:
/// stranded-loop pruning, redundant-edge removal, sliver removal,
/// face simplification, micro-edge collapse, and tolerance
/// propagation — each stage census-traced under the boolean trace topic
/// (sliver removal silently deleted a cylinder wall once, PK-150).
pub(crate) fn run_output_postprocess(
    result: &mut BRep,
    op: crate::BoolOperation,
    tolerance: &Tolerance,
) {
    // Post-processing: prune stranded half-edges left behind by
    // broken merge loops. The merge phase occasionally tags
    // half-edges with a loop handle they're not reachable from;
    // this pruning step removes them so the visible topology
    // matches the Euler-Poincaré accounting.
    trace_face_census("pre-postprocess", result);
    let _ = crate::postprocess::loop_repair::repair_stranded_loops(result);
    trace_face_census("post-stranded", result);

    // Post-processing: remove redundant edges if the result topology is valid.
    if result.validate().is_ok() {
        remove_redundant_edges_safe(result, tolerance);
    }
    trace_face_census("post-redundant", result);

    // Post-processing: remove sliver faces (§8.3.1) if the topology
    // is still valid. The pass merges each sliver into a non-sliver
    // neighbor via KEF, validating after every removal and bailing
    // out on any corruption rather than compounding damage.
    if result.validate().is_ok() {
        let _ = crate::postprocess::sliver_removal::remove_slivers(result, tolerance);
    }
    trace_face_census("post-sliver", result);

    // Post-processing: face simplification (§11.2). Merges adjacent
    // co-surfacial faces via KEF, validating after each merge with
    // rollback on regression. The previous implementation was
    // disabled because of a stale-edge-handle bug; the iterate-and-
    // recompute pattern fixes that.
    if result.validate().is_ok() {
        crate::postprocess::face_simplify::apply_face_simplification(result, tolerance);
    }
    trace_face_census("post-simplify", result);
    // Note: fix_zero_thickness_sheets remains disabled — its global
    // face-pair comparison is unrelated to this fix and still risky.
    // fix_zero_thickness_sheets(result, tolerance);

    // Post-processing: dissolve redundant collinear valence-2
    // vertices (flush-weld residue). Segment semantics leak into
    // dress-up selection otherwise — a chord-length filter sees half
    // of a physical corner line (the L-union fillet class). Each
    // dissolution is a transactional kve; pcurves cleared by the pass
    // are re-derived immediately below.
    if result.validate().is_ok() {
        let dissolved = crate::postprocess::collinear_vertex::dissolve_collinear_valence2_vertices(
            result, tolerance,
        );
        if dissolved > 0 {
            // PK-274: each dissolution is a transactional kve that
            // cleared the pcurves it touched — refill exactly those.
            let _ = result.recompute_pcurves_where_missing();
            if popkern_core::trace::trace_enabled("boolean") || cfg!(debug_assertions) {
                eprintln!(
                    "[boolean] post-processing: dissolved {dissolved} collinear valence-2 vertex(es)"
                );
            }
        }
    }
    trace_face_census("post-collinear-vertex", result);

    // Post-processing: collapse micro-edges (position merge).
    let micro_edges = crate::imprint::micro_edge::find_micro_edges(result, tolerance.linear());
    if !micro_edges.is_empty() {
        crate::imprint::micro_edge::collapse_micro_edge_positions(result, &micro_edges);
    }

    // Post-processing: retire degenerate spur edges (PK-293). Runs
    // AFTER the micro-edge position pass so spurs that pass squashed
    // to zero length are caught along with walk-emitted stubs; each
    // removal is a transactional kev.
    let spur_report =
        crate::postprocess::spur_edge::collapse_degenerate_spur_edges(result, tolerance);
    if spur_report.spurs_removed > 0
        && (popkern_core::trace::trace_enabled("boolean") || cfg!(debug_assertions))
    {
        eprintln!(
            "[boolean] post-processing: retired {} degenerate spur edge(s)",
            spur_report.spurs_removed
        );
    }
    trace_face_census("post-spur", result);

    // Slit-seam dissolution (PK-295): a fragment of a periodic face
    // that no longer wraps the period can carry the original seam as
    // a chart-co-located slit, which pinches the UV ring and forces
    // the tessellator's grid fallback (57 open mesh edges on the
    // rim-filleted-cylinder subtraction's seam-side bulge).
    let slit_report = crate::postprocess::slit_seam::dissolve_slit_seams(result, tolerance);
    if slit_report.dissolved > 0
        && (popkern_core::trace::trace_enabled("boolean") || cfg!(debug_assertions))
    {
        eprintln!(
            "[boolean] post-processing: dissolved {} slit seam edge(s)",
            slit_report.dissolved
        );
    }
    trace_face_census("post-slit-seam", result);

    // Weld the two masked halves of a collapsed near-tangential gore
    // into one twin pair (PK-309) — postprocess/parallel_weld.rs.
    if result.validate().is_ok() {
        let weld = crate::postprocess::parallel_weld::weld_parallel_self_twins(result, tolerance);
        if weld.welded > 0
            && (cfg!(debug_assertions) || popkern_core::trace::trace_enabled("boolean"))
        {
            eprintln!(
                "[boolean] post-processing: welded {} collapsed-gore rim pair(s) ({} declined)",
                weld.welded, weld.declined
            );
        }
    }
    trace_face_census("post-parallel-weld", result);

    // Exact rims (PK-417): intersection edges that trace a circle are
    // re-spelled as the analytic Arc they are, with co-parameterised
    // pcurves on both faces, so a bore's rims select as `circular()`
    // and keep the analytic path downstream. Cut rims only for now:
    // the G2 closed-fillet ladder does not yet take a union's
    // junction rim as an arc (PK-419 carries the measured evidence).
    let rims = if op == crate::BoolOperation::Subtraction {
        crate::postprocess::exact_rims::promote_exact_rims(result, tolerance)
    } else {
        crate::postprocess::exact_rims::ExactRimReport::default()
    };
    if rims.promoted + rims.declined > 0
        && (cfg!(debug_assertions) || popkern_core::trace::trace_enabled("boolean"))
    {
        eprintln!(
            "[boolean] post-processing: {} circular rim(s) promoted to exact arcs ({} declined)",
            rims.promoted, rims.declined
        );
    }
    trace_face_census("post-exact-rims", result);

    // Post-processing: propagate tolerances to intersection entities (§10.3).
    // Intersection-created vertices/edges may have larger tolerance than
    // original entities due to numerical approximation.
    propagate_output_tolerances(result, tolerance);

    // PK-244 serving-boundary orientation-flag reconcile: the flag is
    // derived state (loop winding × surface geometry), and boolean
    // merge authors stamped tool-body faces (bore walls, countersink
    // cones, slot planes) with constant flags that the corpus audit
    // measured dishonest. Runs last, after every pcurve author above.
    let flips = result.reconcile_orientation_flags();
    if flips > 0 && popkern_core::trace::trace_enabled("boolean") {
        eprintln!("[boolean] post-processing: reconciled {flips} orientation flag(s)");
    }
}
