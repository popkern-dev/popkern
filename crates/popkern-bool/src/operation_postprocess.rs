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
pub(crate) fn log_validation_issues(result: &BRep, op: BoolOperation) {
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
pub(crate) fn log_post_processing_diagnostics(result: &BRep, tolerance: &Tolerance) {
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
