//! PK-276 interactive-control variants of the shell/offset/thicken
//! entry points: cooperative cancellation with the structural "< 1 s
//! clean termination" detach layer, plus per-stage progress.
//!
//! Every wrapper takes one extra
//! [`popkern_core::progress::OperationControl`] argument. With a
//! cancel token the pipeline runs detached on a worker thread over
//! private clones
//! ([`popkern_core::cancellation::run_detach_on_cancel_with_events`]);
//! the worker polls the same token at every stage boundary (bounding
//! post-cancel CPU burn to one stage) and its progress events reach
//! the caller's reporter through the event sink. Without a token the
//! pipeline runs inline. Either way, a call with a never-fired token
//! and a null reporter produces the same result as the plain entry
//! point, and a cancelled run leaves the caller's model untouched.

use popkern_core::cancellation::{
    CancellationToken, DetachOutcome, EventSink, run_detach_on_cancel_with_events,
};
use popkern_core::history::TrackingContext;
use popkern_core::progress::{OperationControl, ProgressEvent, StageTracker};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::ShellError;
use crate::types::{ShellDirection, ShellOptions, ShellResult};

/// Run a shell-crate pipeline under an [`OperationControl`]: inline
/// when no token is present, detached with sink-forwarded progress
/// when one is. `inline` borrows the caller's data; `make_worker`
/// receives the worker's token clone and returns the detached body
/// owning pre-cloned inputs (cloning happens only on the detach path).
fn run_controlled<T, F, M, W>(
    control: OperationControl<'_>,
    total_stages: usize,
    inline: F,
    make_worker: M,
) -> Result<T, ShellError>
where
    T: Send + 'static,
    F: FnOnce(&mut StageTracker<'_, '_>) -> Result<T, ShellError>,
    M: FnOnce(CancellationToken) -> W,
    W: FnOnce(EventSink<ProgressEvent>) -> Result<T, ShellError> + Send + 'static,
{
    let OperationControl {
        cancel,
        mut reporter,
    } = control;
    match cancel {
        None => inline(&mut StageTracker::new(None, reporter, total_stages)),
        Some(token) => {
            let worker = make_worker(token.clone());
            match run_detach_on_cancel_with_events(token, worker, |event| {
                if let Some(r) = reporter.as_mut() {
                    r.on_stage(event);
                }
            }) {
                DetachOutcome::Completed(result) => result,
                DetachOutcome::Detached => Err(ShellError::Cancelled),
            }
        }
    }
}

/// [`crate::shell`] with PK-276 interactive control. A cancelled run
/// leaves the caller's model untouched and records no names.
#[allow(clippy::too_many_arguments)]
pub fn shell_with_control(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<ShellResult, ShellError> {
    // Midplane runs its own two-ladder pipeline with a different honest
    // stage total (PK-276: the promised count must match the emitted
    // events).
    let stage_count = if matches!(direction, ShellDirection::Midplane) {
        crate::midplane::MIDPLANE_STAGE_COUNT
    } else {
        crate::SHELL_STAGE_COUNT
    };
    // Provenance rides the generic payload across the detach boundary
    // (PK-328): the worker owns its record and returns it with the
    // result, so the interactive path names exactly like the plain one.
    let track = tracking.is_some();
    let (mut result, provenance) = run_controlled(
        control,
        stage_count,
        |tracker| {
            let mut provenance = crate::ShellProvenance::default();
            let result = crate::shell_impl(
                brep,
                thickness,
                direction,
                removed_faces,
                face_thicknesses,
                options,
                tolerance,
                tracker,
                track.then_some(&mut provenance),
            )?;
            Ok((result, provenance))
        },
        |token| {
            let brep = brep.clone();
            let removed_faces = removed_faces.to_vec();
            let face_thicknesses = face_thicknesses.to_vec();
            let options = options.clone();
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(Some(&token), Some(&mut sink), stage_count);
                let mut provenance = crate::ShellProvenance::default();
                let result = crate::shell_impl(
                    &brep,
                    thickness,
                    direction,
                    &removed_faces,
                    &face_thicknesses,
                    &options,
                    tolerance,
                    &mut tracker,
                    track.then_some(&mut provenance),
                )?;
                Ok((result, provenance))
            }
        },
    )?;
    if let Some(ctx) = tracking {
        crate::assign_shell_names(&mut result.brep, ctx, brep, &provenance);
        let pre = popkern_blend::naming::PreOperationSnapshot::capture(brep);
        popkern_blend::naming::assign_output_edge_vertex_names(
            &mut result.brep,
            ctx,
            &pre,
            &tolerance,
        );
    }
    Ok(result)
}

/// [`crate::shell_variable`] with PK-276 interactive control.
#[allow(clippy::too_many_arguments)]
pub fn shell_variable_with_control(
    brep: &BRep,
    default_thickness: f64,
    face_thicknesses: &[(FaceHandle, f64)],
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    control: OperationControl<'_>,
) -> Result<ShellResult, ShellError> {
    // Same pipeline as `shell_with_control` — `shell_variable` and
    // `shell` differ only in argument spelling.
    shell_with_control(
        brep,
        default_thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
        tracking,
        control,
    )
}

/// [`crate::thicken`] with PK-276 interactive control.
pub fn thicken_with_control(
    brep: &BRep,
    thickness: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
    control: OperationControl<'_>,
) -> Result<ShellResult, ShellError> {
    run_controlled(
        control,
        crate::thicken::THICKEN_STAGE_COUNT,
        |tracker| crate::thicken::thicken_impl(brep, thickness, options, tolerance, tracker),
        |token| {
            let brep = brep.clone();
            let options = options.clone();
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::thicken::THICKEN_STAGE_COUNT,
                );
                crate::thicken::thicken_impl(&brep, thickness, &options, tolerance, &mut tracker)
            }
        },
    )
}

/// [`crate::offset_solid`] with PK-276 interactive control.
pub fn offset_solid_with_control(
    brep: &BRep,
    distance: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
    control: OperationControl<'_>,
) -> Result<ShellResult, ShellError> {
    offset_solid_variable_with_control(brep, distance, &[], options, tolerance, control)
}

/// [`crate::offset_solid_variable`] with PK-276 interactive control.
pub fn offset_solid_variable_with_control(
    brep: &BRep,
    default_distance: f64,
    face_distances: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    control: OperationControl<'_>,
) -> Result<ShellResult, ShellError> {
    run_controlled(
        control,
        crate::offset::OFFSET_STAGE_COUNT,
        |tracker| {
            crate::offset::offset_solid_impl(
                brep,
                default_distance,
                face_distances,
                options,
                tolerance,
                tracker,
            )
        },
        |token| {
            let brep = brep.clone();
            let face_distances = face_distances.to_vec();
            let options = options.clone();
            move |mut sink: EventSink<ProgressEvent>| {
                let mut tracker = StageTracker::new(
                    Some(&token),
                    Some(&mut sink),
                    crate::offset::OFFSET_STAGE_COUNT,
                );
                crate::offset::offset_solid_impl(
                    &brep,
                    default_distance,
                    &face_distances,
                    &options,
                    tolerance,
                    &mut tracker,
                )
            }
        },
    )
}
