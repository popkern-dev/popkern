//! PK-276 interactive-control variants of the sweep-family entry
//! points: cooperative cancellation with the structural "< 1 s clean
//! termination" detach layer, plus per-stage progress.
//!
//! Contract (shared with `popkern-shell::interactive`): with a cancel
//! token the pipeline runs detached on a worker thread over private
//! clones; the worker polls the same token at every stage boundary
//! (bounding post-cancel burn to one stage) and its progress events
//! reach the caller's reporter through the event sink. Without a
//! token it runs inline. A call with a never-fired token and a null
//! reporter produces the same result as the plain entry point.
//!
//! Tracking handoff (the fillet contract): the caller's
//! [`TrackingContext`] moves into the worker; on completion the
//! recorded context moves back, on cancellation the caller holds a
//! fresh context for the same operation — a cancelled operation
//! records nothing.

use popkern_core::cancellation::{
    CancellationToken, DetachOutcome, EventSink, run_detach_on_cancel_with_events,
};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::progress::{OperationControl, ProgressEvent, StageTracker};
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::types::{
    HelixSpec, LoftOptions, LoftSection, PipeRadius, SweepOptions, SweepProfile, SweepResult,
};

type WorkerOutput = (Result<SweepResult, SweepError>, Option<TrackingContext>);

/// Run a sweep-family pipeline under an [`OperationControl`]: inline
/// when no token is present, detached (with sink-forwarded progress
/// and tracking handoff) when one is.
fn run_controlled<F, M, W>(
    control: OperationControl<'_>,
    tracking: Option<&mut TrackingContext>,
    total_stages: usize,
    inline: F,
    make_worker: M,
) -> Result<SweepResult, SweepError>
where
    F: FnOnce(
        Option<&mut TrackingContext>,
        &mut StageTracker<'_, '_>,
    ) -> Result<SweepResult, SweepError>,
    M: FnOnce(CancellationToken, Option<TrackingContext>) -> W,
    W: FnOnce(EventSink<ProgressEvent>) -> WorkerOutput + Send + 'static,
{
    let OperationControl {
        cancel,
        mut reporter,
    } = control;
    match cancel {
        None => inline(
            tracking,
            &mut StageTracker::new(None, reporter, total_stages),
        ),
        Some(token) => {
            let mut tracking = tracking;
            let worker_tracking: Option<TrackingContext> = tracking.as_mut().map(|ctx| {
                let operation = ctx.operation;
                let scheme = ctx.naming_scheme;
                std::mem::replace(
                    &mut **ctx,
                    TrackingContext::new(operation).with_naming_scheme(scheme),
                )
            });
            let worker = make_worker(token.clone(), worker_tracking);
            match run_detach_on_cancel_with_events(token, worker, |event| {
                if let Some(r) = reporter.as_mut() {
                    r.on_stage(event);
                }
            }) {
                DetachOutcome::Completed((result, recorded)) => {
                    if let (Some(ctx), Some(recorded)) = (tracking, recorded) {
                        *ctx = recorded;
                    }
                    result
                }
                DetachOutcome::Detached => Err(SweepError::Cancelled),
            }
        }
    }
}

/// [`crate::sweep`] with PK-276 interactive control.
pub fn sweep_with_control(
    profile: &SweepProfile,
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    run_controlled(
        control,
        tracking,
        crate::pipeline::SWEEP_STAGE_COUNT,
        |tracking, tracker| {
            crate::pipeline::sweep_impl(profile, path, options, tolerance, tracking, tracker)
        },
        |token, mut worker_tracking| {
            let (profile, path, options) = (profile.clone(), path.clone(), options.clone());
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::pipeline::SWEEP_STAGE_COUNT,
                );
                let result = crate::pipeline::sweep_impl(
                    &profile,
                    &path,
                    &options,
                    tolerance,
                    worker_tracking.as_mut(),
                    &mut tracker,
                );
                (result, worker_tracking)
            }
        },
    )
}

/// [`crate::loft`] with PK-276 interactive control.
pub fn loft_with_control(
    sections: &[LoftSection],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    guided_loft_with_control(sections, &[], options, tolerance, tracking, control)
}

/// [`crate::guided_loft`] with PK-276 interactive control (empty
/// `guides` = plain loft; the guide stage's event is then never
/// emitted).
pub fn guided_loft_with_control(
    sections: &[LoftSection],
    guides: &[NurbsCurve],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    run_controlled(
        control,
        tracking,
        crate::pipeline::LOFT_STAGE_COUNT,
        |tracking, tracker| {
            crate::pipeline::loft_impl(sections, guides, options, tolerance, tracking, tracker)
        },
        |token, mut worker_tracking| {
            let (sections, guides, options) = (sections.to_vec(), guides.to_vec(), options.clone());
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::pipeline::LOFT_STAGE_COUNT,
                );
                let result = crate::pipeline::loft_impl(
                    &sections,
                    &guides,
                    &options,
                    tolerance,
                    worker_tracking.as_mut(),
                    &mut tracker,
                );
                (result, worker_tracking)
            }
        },
    )
}

/// [`crate::pipe`] with PK-276 interactive control.
pub fn pipe_with_control(
    path: &NurbsCurve,
    radius: &PipeRadius,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    run_controlled(
        control,
        tracking,
        crate::pipeline::PIPE_STAGE_COUNT,
        |tracking, tracker| {
            crate::pipeline::pipe_impl(path, radius, options, tolerance, tracking, tracker)
        },
        |token, mut worker_tracking| {
            let (path, radius, options) = (path.clone(), radius.clone(), options.clone());
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::pipeline::PIPE_STAGE_COUNT,
                );
                let result = crate::pipeline::pipe_impl(
                    &path,
                    &radius,
                    &options,
                    tolerance,
                    worker_tracking.as_mut(),
                    &mut tracker,
                );
                (result, worker_tracking)
            }
        },
    )
}

/// [`crate::helical_sweep`] with PK-276 interactive control. The
/// (cheap, analytic) helix path builds on the calling thread; the
/// sweep itself runs under the same contract as
/// [`sweep_with_control`].
pub fn helical_sweep_with_control(
    profile: &SweepProfile,
    spec: &HelixSpec,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    let (path, mut warnings) = crate::helix::build_helix_path(spec)?;
    let mut result = sweep_with_control(profile, &path, options, tolerance, tracking, control)?;
    result.warnings.append(&mut warnings);
    Ok(result)
}

/// [`crate::guided_sweep`] with PK-276 interactive control.
pub fn guided_sweep_with_control(
    profile: &SweepProfile,
    path: &NurbsCurve,
    guides: &[NurbsCurve],
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    run_controlled(
        control,
        tracking,
        crate::pipeline::GUIDED_SWEEP_STAGE_COUNT,
        |tracking, tracker| {
            crate::pipeline::guided_sweep_impl(
                profile, path, guides, options, tolerance, tracking, tracker,
            )
        },
        |token, mut worker_tracking| {
            let (profile, path, guides, options) = (
                profile.clone(),
                path.clone(),
                guides.to_vec(),
                options.clone(),
            );
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::pipeline::GUIDED_SWEEP_STAGE_COUNT,
                );
                let result = crate::pipeline::guided_sweep_impl(
                    &profile,
                    &path,
                    &guides,
                    &options,
                    tolerance,
                    worker_tracking.as_mut(),
                    &mut tracker,
                );
                (result, worker_tracking)
            }
        },
    )
}

/// [`crate::swept_blend`] with PK-276 interactive control.
pub fn swept_blend_with_control(
    profiles: &[(f64, SweepProfile)],
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    run_controlled(
        control,
        tracking,
        crate::pipeline::SWEPT_BLEND_STAGE_COUNT,
        |tracking, tracker| {
            crate::pipeline::swept_blend_impl(profiles, path, options, tolerance, tracking, tracker)
        },
        |token, mut worker_tracking| {
            let (profiles, path, options) = (profiles.to_vec(), path.clone(), options.clone());
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::pipeline::SWEPT_BLEND_STAGE_COUNT,
                );
                let result = crate::pipeline::swept_blend_impl(
                    &profiles,
                    &path,
                    &options,
                    tolerance,
                    worker_tracking.as_mut(),
                    &mut tracker,
                );
                (result, worker_tracking)
            }
        },
    )
}

/// [`crate::revolve`] with PK-276 interactive control.
pub fn revolve_with_control(
    profile: &SweepProfile,
    axis_origin: &Point3,
    axis_direction: &Vec3,
    angle: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<SweepResult, SweepError> {
    run_controlled(
        control,
        tracking,
        crate::revolve::REVOLVE_STAGE_COUNT,
        |tracking, tracker| {
            crate::revolve::revolve_impl(
                profile,
                axis_origin,
                axis_direction,
                angle,
                tolerance,
                tracking,
                tracker,
            )
        },
        |token, mut worker_tracking| {
            let (profile, axis_origin, axis_direction) =
                (profile.clone(), *axis_origin, *axis_direction);
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::revolve::REVOLVE_STAGE_COUNT,
                );
                let result = crate::revolve::revolve_impl(
                    &profile,
                    &axis_origin,
                    &axis_direction,
                    angle,
                    tolerance,
                    worker_tracking.as_mut(),
                    &mut tracker,
                );
                (result, worker_tracking)
            }
        },
    )
}
