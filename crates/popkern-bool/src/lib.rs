#![warn(missing_docs)]
#![allow(clippy::result_large_err)]

//! popkern-bool: Boolean set operations on BRep solids.
//!
//! Provides union, intersection, and subtraction operations on
//! [`BRep`](popkern_core::topology::BRep) solid models.
//!
//! # Pipeline
//!
//! Each boolean operation proceeds through five stages:
//!
//! 1. **Intersect** — compute face-face intersection curves.
//! 2. **Split** — trim faces along those curves.
//! 3. **Classify** — label each fragment as inside, outside, or on-boundary.
//! 4. **Select** — choose which fragments belong to the result.
//! 5. **Merge** — stitch fragments into a valid closed solid.
//!
//! # Example (placeholder — not yet implemented)
//!
//! ```rust,ignore
//! use popkern_bool::{BoolOperation, boolean};
//! use popkern_core::tolerance::Tolerance;
//!
//! let result = boolean(BoolOperation::Union, solid_a, solid_b, Tolerance::default())?;
//! ```

pub mod adaptive_precision;
mod bool_api_extras;
pub mod cancellation;
pub mod cellular;
pub mod chart_walk_fallback;
pub(crate) mod classify;
pub mod classify_on_deferred;
pub mod classify_on_face;
pub mod classify_propagation;
pub(crate) mod classify_uv_loops;
pub mod classify_winding;
pub(crate) mod clip;
mod clip_boundary_edges;
mod clip_endpoint_far;
mod clip_endpoint_refine;
mod clip_far_junction;
mod clip_triple_junction;
mod closed_source_window;
pub mod coincidence;
pub mod csg;
mod curve_projection;
mod entity_lineage;
mod entity_naming;
mod entity_naming_edges;
pub mod error;
mod error_remediation;
pub mod history;
pub mod imprint;
pub mod instancing;
pub(crate) mod intersect;
pub(crate) mod intersect_edge_face;
pub(crate) mod intersect_edge_face_coincident;
pub(crate) mod intersect_residue_filter;
pub(crate) mod intersect_residue_split;
pub(crate) mod merge;
pub mod merge_fallback;
mod pcurve_box;
pub(crate) mod periodic_chart_parity;
pub(crate) mod periodic_chart_parity_torus;
pub mod propagation_report;
pub(crate) mod record_side_witness;

/// PK-478 E1 diagnostic export (see `merge::closed_side_census`):
/// closed-edge witness census for the sense-aware-merge exploration.
#[doc(hidden)]
pub use merge::graph_build::closed_side_census::{ClosedSideCensusRow, closed_side_census};
pub mod multi_body;
pub(crate) mod operation;
mod operation_postprocess_trace;
pub mod ops;
pub mod parallel;
pub mod perturbation;
pub mod postprocess;
pub mod progress;
pub mod provenance;
pub(crate) mod select;
pub mod sheet;
pub(crate) mod split;
pub(crate) mod split_inner_loop;
pub mod tangent_analysis;
pub mod tolerance;
pub mod validate;

pub use bool_api_extras::{boolean_strict, boolean_with_provenance, boolean_with_remediation};
pub use classify::{point_in_solid, point_in_solid_final};
pub use csg::{CsgNode, csg_eval};

pub use error::{
    BoolDiagnostics, BoolError, BoolOperation, Classification, PointClassification, Stage,
};
pub use ops::interference::{InterferenceMode, InterferenceReport, interference_check};
pub use postprocess::quality::{QualityReport, QualityThresholds, quality_report};
pub use postprocess::sliver_removal::{SliverRemovalReport, remove_slivers};
pub use progress::{CaptureReporter, NullReporter, ProgressEvent, ProgressReporter};
pub use sheet::{
    SheetBooleanResult, SheetOperationKind, SideOfSheet, boolean_with_sheets,
    classify_sheet_operation, is_closed_solid, is_sheet_body, side_of_sheet,
};

#[cfg(test)]
mod error_and_tracking_tests;
#[cfg(test)]
mod error_remediation_tests;
#[cfg(test)]
mod integration_tests;
#[cfg(test)]
mod integration_tests_curved;
#[cfg(test)]
mod integration_tests_curved_regression;
#[cfg(test)]
mod integration_tests_hole_debug;
#[cfg(test)]
mod integration_tests_multicurve;
#[cfg(test)]
mod translate_helper_tests;

use popkern_core::history::TrackingContext;
use popkern_core::{tolerance::Tolerance, topology::BRep};

/// Perform a boolean operation on two solids, consuming both inputs.
///
/// Returns the result solid on success, or a [`BoolError`] describing the failure.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tol` — geometric tolerance for all distance and coincidence tests.
pub fn boolean(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    // §2.3: reject sheet body inputs with a structured error rather
    // than letting volume-based classification produce undefined
    // results. Future work will add boundary-relative semantics; for
    // now this is the only safe response.
    let a_is_sheet = sheet::is_sheet_body(&a);
    let b_is_sheet = sheet::is_sheet_body(&b);
    if a_is_sheet || b_is_sheet {
        return Err(BoolError::OpenShellNotSupported {
            a: a_is_sheet,
            b: b_is_sheet,
        });
    }

    // Collect input face data before the boolean consumes the inputs.
    let mut input_a_faces = if tracking.is_some() {
        entity_naming::collect_face_data(&a)
    } else {
        Vec::new()
    };
    let mut input_b_faces = if tracking.is_some() {
        entity_naming::collect_face_data(&b)
    } else {
        Vec::new()
    };
    // Edges and vertices too (PK-386): a dress-up after a boolean is
    // recorded by stable name, so the survivors must keep theirs.
    // Where operand b's entries start in the combined edge / vertex
    // lists (their names come from b's pool).
    let (mut input_edges, mut input_vertices, b_edges_from, b_vertices_from) = if tracking.is_some()
    {
        let mut edges = entity_naming_edges::collect_edge_data(&a);
        let b_edges_from = edges.len();
        edges.extend(entity_naming_edges::collect_edge_data(&b));
        let mut vertices = entity_naming_edges::collect_vertex_data(&a);
        let b_vertices_from = vertices.len();
        vertices.extend(entity_naming_edges::collect_vertex_data(&b));
        (edges, vertices, b_edges_from, b_vertices_from)
    } else {
        (Vec::new(), Vec::new(), 0, 0)
    };

    let mut result = operation::execute_boolean(op, a.clone(), b.clone(), &tol).or_else(|e| {
        // Perturbation fallback: retry with slightly perturbed geometry
        // if the error is the kind that perturbation can fix.
        if matches!(
            e,
            BoolError::CoincidentFaces { .. }
                | BoolError::IntersectionFailed(_)
                | BoolError::ClassificationFailed { .. }
        ) {
            perturbation::boolean_with_perturbation(
                op,
                &a,
                &b,
                tol,
                &perturbation::PerturbationConfig::default(),
            )
        } else {
            Err(e)
        }
    })?;
    if let Some(ctx) = tracking {
        // Naming scheme 3 (PK-731): the result pool starts as the
        // target operand's (a shared clone), as every other operation's
        // result carries its input's pool: the target's names keep
        // their ids, so importing them is the identity. Scheme 2 built
        // the pool by importing every operand name, ancestors first,
        // which re-interned the target's whole history on every
        // boolean (49k names at the 64th slot of a disc).
        let seeded = ctx.naming_scheme >= 3;
        if seeded {
            result.name_pool = a.name_pool.clone();
        }
        if ctx.naming_scheme >= 2 {
            // Naming scheme 2 (PK-677): each operand name comes into the
            // result pool by value, ancestors first, so a tool-side
            // name's ancestor ids are this pool's (scheme 1 interned the
            // tool's ids as they were).
            let pool = &mut result.name_pool;
            for (name, ..) in input_a_faces.iter_mut().filter(|_| !seeded) {
                entity_naming::import_name(name, &a, pool);
            }
            for (name, ..) in &mut input_b_faces {
                entity_naming::import_name(name, &b, pool);
            }
            for (i, (name, ..)) in input_edges.iter_mut().enumerate() {
                if i >= b_edges_from || !seeded {
                    entity_naming::import_name(name, if i < b_edges_from { &a } else { &b }, pool);
                }
            }
            for (i, (name, _)) in input_vertices.iter_mut().enumerate() {
                if i >= b_vertices_from || !seeded {
                    let source = if i < b_vertices_from { &a } else { &b };
                    entity_naming::import_name(name, source, pool);
                }
            }
        }
        entity_naming::assign_boolean_names(&mut result, &input_a_faces, &input_b_faces, ctx, &tol);
        // Lineage beside the names (PK-633): which input faces each
        // output face is a piece of, so labels follow trimmed faces.
        entity_lineage::record_face_lineage(&mut result, &a, &b, &tol);
        entity_naming_edges::assign_boolean_edge_vertex_names(
            &mut result,
            &input_edges,
            &input_vertices,
            ctx,
            &tol,
        );
        entity_lineage::record_edge_lineage(&mut result, &a, &b, &tol);
        popkern_core::naming_lineage::prune_unreadable_lineage(&mut result);
    }
    Ok(result)
}

/// Concatenate two solids into ONE multi-solid [`BRep`] with NO
/// boolean semantics: every shell of `b` is copied beside `a`'s.
/// This is the whole-document view (an export or render of several
/// bodies at once) — no classification runs, so kissing or even
/// overlapping bodies concatenate verbatim. For material semantics
/// use [`boolean`].
///
/// # Errors
///
/// Returns [`BoolError`] if the structural copy fails.
pub fn concatenate_solids(a: BRep, b: BRep) -> Result<BRep, BoolError> {
    operation::helpers::combine_breps(a, b)
}

/// Perform a boolean operation on two solids, borrowing both inputs.
///
/// Clones `a` and `b` internally before operating. Prefer [`boolean`] when
/// ownership can be transferred to avoid the clone cost.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tol` — geometric tolerance for all distance and coincidence tests.
/// * `tracking` — optional tracking context for stable naming.
pub fn boolean_ref(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(op, a.clone(), b.clone(), tol, tracking)
}

/// Compute the union of two solids (A ∪ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Union, a, b, tol, tracking)`.
pub fn union(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Union, a, b, tol, tracking)
}

/// Compute the union of two solids (A ∪ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Union, a, b, tol, tracking)`.
pub fn union_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Union, a, b, tol, tracking)
}

/// Compute the intersection of two solids (A ∩ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Intersection, a, b, tol, tracking)`.
pub fn intersection(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Intersection, a, b, tol, tracking)
}

/// Compute the intersection of two solids (A ∩ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Intersection, a, b, tol, tracking)`.
pub fn intersection_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Intersection, a, b, tol, tracking)
}

/// Subtract `b` from `a` (A − B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Subtraction, a, b, tol, tracking)`.
pub fn subtraction(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Subtraction, a, b, tol, tracking)
}

/// Subtract `b` from `a` (A − B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Subtraction, a, b, tol, tracking)`.
pub fn subtraction_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Subtraction, a, b, tol, tracking)
}

/// Perform a boolean operation and return per-phase telemetry alongside the result.
///
/// Identical to [`boolean`] but populates a [`BoolDiagnostics`] with
/// [`PhaseMetric`](crate::error::PhaseMetric) samples for each pipeline stage
/// (intersection, split, classification, selection, merge, validation). Use
/// this entry point for performance regression detection and tooling that
/// needs to attribute time to specific phases.
///
/// # Arguments
///
/// * `op` — the boolean set operation.
/// * `a`, `b` — operand solids (cloned).
/// * `tol` — geometric tolerance.
pub fn boolean_with_diagnostics(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
) -> Result<(BRep, BoolDiagnostics), BoolError> {
    let mut metrics = Vec::new();
    let result =
        operation::execute_boolean_with_metrics(op, a.clone(), b.clone(), &tol, &mut metrics)?;
    let diagnostics = BoolDiagnostics {
        phase_metrics: metrics,
        ..BoolDiagnostics::default()
    };
    Ok((result, diagnostics))
}

/// Perform a boolean with cooperative cancellation support (§15) and
/// the structural "< 1 s clean termination" guarantee.
///
/// Behaves identically to [`boolean`] when the token never fires and
/// returns [`BoolError::Cancelled`] as soon as a cancellation request
/// is observed. The input model is unchanged.
///
/// The pipeline runs on a worker thread over private clones (PK-276,
/// the same detach-on-cancel layer as
/// `popkern_blend::fillet_edges_cancellable`): when the token fires,
/// this function returns immediately — the latency does not depend on
/// any stage's internal poll granularity. The abandoned worker
/// observes the token at its next stage-boundary poll, discards its
/// clones, and exits; completion wins a tie so a finished result is
/// never thrown away.
///
/// Pair this with [`cancellation::CancellationToken`] held by a UI
/// progress dialog or signal handler.
pub fn boolean_cancellable(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    cancel: &cancellation::CancellationToken,
) -> Result<BRep, BoolError> {
    let (worker_a, worker_b, worker_token) = (a.clone(), b.clone(), cancel.clone());
    match popkern_core::cancellation::run_detach_on_cancel(cancel, move || {
        let mut metrics = Vec::new();
        operation::execute_boolean_full(
            op,
            worker_a,
            worker_b,
            &tol,
            &mut metrics,
            Some(&worker_token),
            None,
        )
    }) {
        popkern_core::cancellation::DetachOutcome::Completed(result) => result,
        popkern_core::cancellation::DetachOutcome::Detached => Err(BoolError::Cancelled),
    }
}

/// Perform a boolean and report progress through a [`progress::ProgressReporter`] (§15).
///
/// Behaves identically to [`boolean`] but invokes `reporter.on_phase`
/// after every pipeline stage with a [`progress::ProgressEvent`].
/// Combine with [`boolean_cancellable`] in your own wrapper if you
/// need both progress notifications and cancellation.
pub fn boolean_with_progress(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    reporter: &mut dyn progress::ProgressReporter,
) -> Result<BRep, BoolError> {
    let mut metrics = Vec::new();
    operation::execute_boolean_full(
        op,
        a.clone(),
        b.clone(),
        &tol,
        &mut metrics,
        None,
        Some(reporter),
    )
}

/// Wedge-matching trace gate (fragment-walk successor diagnostics).
pub(crate) fn trace_enabled_wedge() -> bool {
    popkern_core::trace::trace_enabled("wedge-match")
}
