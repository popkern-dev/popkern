//! PK-276 interactive-control variants of the smooth-mesh entry
//! points (extracted from `smooth_mesh.rs` for the file-size cap).
//! Contract: with a cancel token the construction runs detached on a
//! worker thread over a private clone of the cage; without one it runs
//! inline. See `docs/architecture/interactive-operation-control.md`.

use popkern_core::tolerance::Tolerance;

use crate::error::PartError;
use crate::part::Part;
use popkern_core::geometry::polyhedral_spline::ControlNet;

use crate::smooth_mesh::{
    ATTEMPT_STAGE_COUNT, SmoothMeshOptions, SmoothMeshReport, smooth_solid_impl,
};

/// [`smooth_solid_from_control_net`] with PK-276 interactive control
/// (cooperative cancellation + per-stage progress).
///
/// With a cancel token the construction runs detached on a worker
/// thread over a private clone of the cage: the call returns
/// [`PartError::Cancelled`] within milliseconds of the token firing,
/// and the abandoned worker exits at its next stage boundary (worst
/// gap: one opaque stage — the spline-complex build or the O(n²)
/// self-intersection screen, multi-second on cage-scale bodies).
/// Without a token the pipeline runs inline. A call with a never-fired
/// token and a null reporter matches the plain entry point's result.
pub fn smooth_solid_from_control_net_with_control(
    net: &ControlNet,
    tolerance: Tolerance,
    options: SmoothMeshOptions,
    control: popkern_core::progress::OperationControl<'_>,
) -> Result<(Part, SmoothMeshReport), PartError> {
    use popkern_core::cancellation::{DetachOutcome, run_detach_on_cancel_with_events};
    use popkern_core::progress::StageTracker;
    let total_stages = (options.max_refinement_steps + 1) * ATTEMPT_STAGE_COUNT;
    let popkern_core::progress::OperationControl {
        cancel,
        mut reporter,
    } = control;
    match cancel {
        None => smooth_solid_impl(
            net,
            tolerance,
            options,
            &mut StageTracker::new(None, reporter, total_stages),
        ),
        Some(token) => {
            let (worker_net, worker_token) = (net.clone(), token.clone());
            match run_detach_on_cancel_with_events(
                token,
                move |mut sink| {
                    let mut tracker =
                        StageTracker::new(Some(&worker_token), Some(&mut sink), total_stages);
                    smooth_solid_impl(&worker_net, tolerance, options, &mut tracker)
                },
                |event| {
                    if let Some(r) = reporter.as_mut() {
                        r.on_stage(event);
                    }
                },
            ) {
                DetachOutcome::Completed(result) => result,
                DetachOutcome::Detached => Err(PartError::Cancelled),
            }
        }
    }
}

/// [`smooth_solid_from_obj`] with PK-276 interactive control: the
/// (cheap) OBJ parse runs on the calling thread; the construction runs
/// under the same contract as
/// [`smooth_solid_from_control_net_with_control`].
pub fn smooth_solid_from_obj_with_control(
    obj_text: &str,
    tolerance: Tolerance,
    options: SmoothMeshOptions,
    control: popkern_core::progress::OperationControl<'_>,
) -> Result<(Part, SmoothMeshReport), PartError> {
    let net = ControlNet::from_obj_str(obj_text).map_err(|source| PartError::SmoothMesh {
        stage: "obj parse",
        detail: source.to_string(),
    })?;
    smooth_solid_from_control_net_with_control(&net, tolerance, options, control)
}
