#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.
//!
//! Provides [`fillet_edges`] and [`chamfer_edges`] for applying fillets and
//! chamfers to BRep edges, plus convenience wrappers [`fillet_edge`] and
//! [`chamfer_edge`] for single-edge operations.

pub mod analysis;
pub mod blend;
pub mod cancellation;
pub mod center_curve;
pub mod center_curve_helpers;
pub mod chamfer_curve;
pub mod construct;
pub mod error;
pub mod fillet_surface;
pub(crate) mod geometry_compute;
pub mod history;
pub mod identify;
pub(crate) mod identify_helpers;
pub mod modify;
pub(crate) mod modify_inplace;
pub(crate) mod modify_inplace_geometry;
pub(crate) mod modify_inplace_sweep;
pub(crate) mod naming;
pub mod offset;
pub mod overlap;
pub mod precheck;
pub mod preview;
pub mod progress;
pub mod quality;
pub mod remediation;
pub mod trim;
pub mod types;
pub mod unfillet;
pub(crate) mod unfillet_closed;
pub(crate) mod unfillet_closed_helpers;
pub(crate) mod unfillet_helpers;
pub(crate) mod unfillet_surgery;
pub mod validate;

pub use cancellation::CancellationToken;
pub use error::{
    BlendError, BlendFailure, BlendWarning, Remediation, RemediationAction, RemediationContext,
};
pub use history::{
    ChamferOperationRecord, FilletOperationRecord, chamfer_recorded, fillet_recorded,
};
pub use identify::{
    BlendInfo, BlendKind, IdentificationSource, blend_neighbours, blend_radius,
    blend_support_faces, identify_all_blends, identify_blend,
};
pub use modify::modify_fillet_radius;
pub use precheck::{
    BlendQuality, ComputeQualityInput, FilletFeasibility, FilletIssue, FilletRadiusBounds,
    feasible_radius_bounds, fillet_quality, precheck_fillet,
};
pub use preview::{PreviewEdgeGeometry, PreviewGeometry, preview_chamfer, preview_fillet};
pub use progress::{CaptureReporter, FilletPhase, NullReporter, ProgressEvent, ProgressReporter};
pub use types::*;
pub use unfillet::{UnfilletPlan, compute_unfillet_plan, unfillet_face, unfillet_via_record};

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Fillet the specified edges of a BRep solid.
///
/// Analyzes each edge, computes the rolling-ball center curve and fillet
/// surface, derives trim curves, and applies the topology modifications.
///
/// # Errors
///
/// Returns [`BlendError::NoValidEdges`] if analysis fails for all edges, or
/// [`BlendError::PartialFailure`] if some edges fail (with a partial result).
pub fn fillet_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(brep, specs, options, tolerance, tracking, None, None)
}

/// Fillet the specified edges, checking a [`CancellationToken`]
/// between pipeline phases.
///
/// If the token reports `is_cancelled()` between any two phases,
/// the operation returns [`BlendError::Cancelled`] and the input
/// BRep is guaranteed unchanged. Otherwise behaves identically to
/// [`fillet_edges`].
///
/// # Errors
///
/// Same as [`fillet_edges`], plus [`BlendError::Cancelled`] when
/// cancellation was observed.
pub fn fillet_edges_cancellable(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: &cancellation::CancellationToken,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(
        brep,
        specs,
        options,
        tolerance,
        tracking,
        Some(cancel),
        None,
    )
}

/// Fillet the specified edges, emitting progress events between
/// pipeline phases via a [`ProgressReporter`].
///
/// The reporter is called once per pipeline phase in
/// [`FilletPhase`]'s natural order. A slow reporter adds directly
/// to the fillet's wall time, so implementations should be cheap.
///
/// # Errors
///
/// Same as [`fillet_edges`].
pub fn fillet_edges_with_progress(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    reporter: &mut dyn progress::ProgressReporter,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(
        brep,
        specs,
        options,
        tolerance,
        tracking,
        None,
        Some(reporter),
    )
}

/// The real body of [`fillet_edges`]. The public entry points all
/// dispatch here with their respective combinations of cancellation
/// token and progress reporter. Keeping the pipeline in one place
/// means cancellation / progress wiring cannot drift between the
/// variants.
fn fillet_edges_internal(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    mut reporter: Option<&mut dyn progress::ProgressReporter>,
) -> Result<BlendResult, BlendError> {
    use std::time::Instant;

    macro_rules! emit {
        ($phase:expr, $start:expr, $count:expr) => {
            if let Some(r) = reporter.as_deref_mut() {
                r.on_phase(progress::ProgressEvent {
                    phase: $phase,
                    duration: $start.elapsed(),
                    entity_count: $count,
                });
            }
        };
    }

    macro_rules! check_cancel {
        () => {
            if let Some(tok) = cancel
                && tok.is_cancelled()
            {
                return Err(BlendError::Cancelled);
            }
        };
    }

    // PHASE 0: Validate.
    let t0 = Instant::now();
    validate::validate_fillet_specs(brep, specs)?;
    emit!(progress::FilletPhase::Validate, t0, specs.len());
    check_cancel!();

    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    // PHASE 1: Analyze edges.
    let t1 = Instant::now();
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;
    emit!(progress::FilletPhase::Analyze, t1, analyses.len());
    check_cancel!();

    // PHASE 2: Compute geometry for each edge.
    let t2 = Instant::now();
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        // Per-edge cancellation point — lets the user kill a long
        // NURBS surface intersection mid-stride.
        check_cancel!();
        match compute_fillet_geometry(&analyses[i], &spec.radius, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }
    emit!(
        progress::FilletPhase::ComputeGeometry,
        t2,
        edge_geometries.len()
    );
    check_cancel!();

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // PHASE 3: Detect overlaps (warnings only — don't block the operation).
    let t3 = Instant::now();
    let mut overlap_warnings = Vec::new();
    let (overlap_trims, _face_trim_map) = build_overlap_inputs(&edge_geometries, brep);

    let overlaps = overlap::detect_overlaps(&overlap_trims, &tolerance);
    for ov in &overlaps {
        overlap_warnings.push(BlendWarning::FilletOverlap {
            fillet_a: ov.fillet_a,
            fillet_b: ov.fillet_b,
        });
    }
    emit!(
        progress::FilletPhase::DetectOverlaps,
        t3,
        overlap_warnings.len()
    );
    check_cancel!();

    // Face disappearance is checked geometrically using the actual fillet
    // radius and face edge lengths, not the unbounded trim curves (which
    // extend beyond the face boundary and produce false positives).

    if edge_geometries.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // PHASE 4: Apply topology.
    let t4 = Instant::now();
    let mut apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &options.vertex_mode)
            .map_err(|edge_failures| BlendError::PartialFailure {
            failed: edge_failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: edge_failures,
        })?;
    emit!(
        progress::FilletPhase::ApplyTopology,
        t4,
        edge_geometries.len()
    );
    check_cancel!();

    // PHASE 5: Ensure pcurve coverage and correctness.
    let t5 = Instant::now();
    let _ = apply_result.brep.recompute_all_pcurves();
    apply_result.brep.fix_stale_pcurves();
    emit!(
        progress::FilletPhase::RepairPcurves,
        t5,
        apply_result.brep.half_edges.len()
    );

    // Merge overlap warnings with topology warnings.
    let mut all_warnings = apply_result.warnings;
    for warn in overlap_warnings {
        let edge_h = specs
            .first()
            .map(|s| s.edge)
            .unwrap_or(edge_geometries[0].analysis.edge);
        all_warnings.push((edge_h, warn));
    }

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: all_warnings,
    })
}

/// Chamfer the specified edges of a BRep solid.
///
/// Analyzes each edge, computes offset curves on each adjacent face, builds
/// a ruled surface between them, and applies topology modifications.
///
/// # Errors
///
/// Returns [`BlendError::NoValidEdges`] if analysis fails for all edges, or
/// [`BlendError::PartialFailure`] if some edges fail (with a partial result).
pub fn chamfer_edges(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    // §10.1 — reject invalid parameters at the API boundary.
    validate::validate_chamfer_specs(brep, specs)?;

    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology. Chamfers always use mitred vertex mode.
    let mut apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &VertexMode::Mitred)
            .map_err(|edge_failures| BlendError::PartialFailure {
                failed: edge_failures.len(),
                total: specs.len(),
                partial_result: Box::new(brep.clone()),
                failures: edge_failures,
            })?;

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: apply_result.warnings,
    })
}

/// Fillet a set of edges and return the result bundled with
/// post-operation quality metrics ([`BlendQuality`]).
///
/// This is the commercial-grade entry point for UI callers: the
/// quality record drives any "fillet health" indicator in a
/// property panel, and the raw metrics (`R_min`, `patch_count`,
/// `surface_area`) feed into manufacturing-readiness reports.
///
/// The operation is otherwise identical to [`fillet_edges`] — same
/// parameter validation, same error taxonomy, same topology
/// modifications. On success the returned [`BlendQuality`] is
/// computed from the new-face delta between the input BRep and the
/// output BRep using [`precheck::fillet_quality`].
///
/// # Errors
///
/// Same as [`fillet_edges`].
pub fn fillet_edges_with_quality(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<(BlendResult, BlendQuality), BlendError> {
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();
    let result = fillet_edges(brep, specs, options, tolerance, tracking)?;
    // Use the spec's largest requested radius as the scoring
    // denominator — for constant radius this is the request, for
    // variable specs it's the peak.
    let requested_radius = specs
        .iter()
        .map(|s| geometry_compute::largest_requested_radius(&s.radius))
        .fold(0.0_f64, f64::max);
    let quality = precheck::fillet_quality(&precheck::ComputeQualityInput {
        post_brep: &result.brep,
        pre_faces: &pre_faces,
        requested_radius,
    });
    Ok((result, quality))
}

/// Convenience: fillet a single edge with a constant radius.
///
/// Equivalent to calling [`fillet_edges`] with a single [`FilletSpec`] and
/// default options.
///
/// # Errors
///
/// Returns a [`BlendError`] if the edge cannot be filleted.
pub fn fillet_edge(
    brep: &BRep,
    edge: EdgeHandle,
    radius: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    fillet_edges(
        brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions::default(),
        tolerance,
        tracking,
    )
}

/// Convenience: chamfer a single edge with equal distances.
///
/// Equivalent to calling [`chamfer_edges`] with a single [`ChamferSpec`] and
/// default options.
///
/// # Errors
///
/// Returns a [`BlendError`] if the edge cannot be chamfered.
pub fn chamfer_edge(
    brep: &BRep,
    edge: EdgeHandle,
    distance: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    chamfer_edges(
        brep,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(distance),
        }],
        &ChamferOptions::default(),
        tolerance,
        tracking,
    )
}

// Internal geometry computation is in `geometry_compute` module.
use geometry_compute::{build_overlap_inputs, compute_chamfer_geometry, compute_fillet_geometry};
use naming::assign_fillet_names;

#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
