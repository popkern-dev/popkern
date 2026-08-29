//! Convenience wrappers over [`super::fillet_edges`] and
//! [`super::chamfer_edges`] — single-edge helpers and the
//! quality-scored pipeline.

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::error::BlendError;
use crate::precheck::{self, BlendQuality};
use crate::types::{
    BlendResult, ChamferKind, ChamferOptions, ChamferSpec, FilletOptions, FilletRadius, FilletSpec,
};
use crate::{chamfer_edges, fillet_edges, geometry_compute_feasibility};

/// Fillet the specified edges and compute a post-operation quality
/// score. See [`super::fillet_edges`] for failure modes and
/// [`BlendQuality`] for the metric definition.
pub fn fillet_edges_with_quality(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<(BlendResult, BlendQuality), BlendError> {
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();
    let result = fillet_edges(brep, specs, options, tolerance, tracking)?;
    let requested_radius = specs
        .iter()
        .map(|s| geometry_compute_feasibility::largest_requested_radius(&s.radius))
        .fold(0.0_f64, f64::max);
    let quality = precheck::fillet_quality(&precheck::ComputeQualityInput {
        post_brep: &result.brep,
        pre_faces: &pre_faces,
        requested_radius,
    });
    Ok((result, quality))
}

/// Convenience: fillet a single edge with a constant radius.
/// Equivalent to [`super::fillet_edges`] with one [`FilletSpec`]
/// + default options. Returns a [`BlendError`] if the edge cannot be filleted.
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
/// Equivalent to [`super::chamfer_edges`] with one [`ChamferSpec`]
/// + default options. Returns a [`BlendError`] if the edge cannot be chamfered.
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

/// Cancellable variant of [`fillet_edges`] — returns
/// [`BlendError::Cancelled`] within one second of the token firing
/// (§17: clean termination < 1 s, model unchanged).
///
/// The pipeline runs on a worker thread. When the token fires, this
/// function returns immediately; the worker observes the same token
/// at its next cooperative poll (between phases and between corner
/// constructions — measured worst gap ≈ one G2 corner, a few
/// seconds), then discards its private clone and exits. The caller's
/// BRep is never touched by a cancelled run, and no partial result
/// escapes.
///
/// `tracking`: the context is moved into the worker for the run. On
/// normal completion it is moved back; on cancellation the caller's
/// context is left freshly re-created for the same operation id —
/// consistent with the [`TrackingContext`] per-operation contract
/// (a cancelled operation records nothing).
pub fn fillet_edges_cancellable(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: &crate::cancellation::CancellationToken,
) -> Result<BlendResult, BlendError> {
    if cancel.is_cancelled() {
        return Err(BlendError::Cancelled);
    }
    // Move the caller's tracking context into the worker; a fresh
    // same-operation context takes its place until (unless) the
    // worker hands the recorded one back.
    let mut tracking = tracking;
    let worker_tracking: Option<TrackingContext> = tracking.as_deref_mut().map(|ctx| {
        let operation = ctx.operation;
        std::mem::replace(ctx, TrackingContext::new(operation))
    });
    let worker_brep = brep.clone();
    let worker_specs = specs.to_vec();
    let worker_options = options.clone();
    let worker_token = cancel.clone();
    // Detach-on-cancel via the shared kernel helper (PK-276): the
    // worker runs over private clones; completion wins a tie so
    // finished work is never thrown away; a worker panic is re-raised
    // here. On detach the worker sees the token at its next
    // cooperative poll and discards its clones.
    match popkern_core::cancellation::run_detach_on_cancel(cancel, move || {
        let mut worker_tracking = worker_tracking;
        let result = crate::fillet_edges_internal(
            &worker_brep,
            &worker_specs,
            &worker_options,
            tolerance,
            worker_tracking.as_mut(),
            Some(&worker_token),
            None,
        );
        (result, worker_tracking)
    }) {
        popkern_core::cancellation::DetachOutcome::Completed((result, recorded)) => {
            if let (Some(caller_ctx), Some(recorded)) = (tracking.take(), recorded) {
                *caller_ctx = recorded;
            }
            result
        }
        popkern_core::cancellation::DetachOutcome::Detached => Err(BlendError::Cancelled),
    }
}

/// Progress-reporting variant of [`fillet_edges`]. See [`ProgressReporter`].
#[rustfmt::skip]
pub fn fillet_edges_with_progress(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    reporter: &mut dyn crate::progress::ProgressReporter,
) -> Result<BlendResult, BlendError> {
    crate::fillet_edges_internal(brep, specs, options, tolerance, tracking, None, Some(reporter))
}
