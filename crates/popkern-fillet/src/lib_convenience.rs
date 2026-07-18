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
/// [`BlendError::Cancelled`] if the token fires between phases.
#[rustfmt::skip]
pub fn fillet_edges_cancellable(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: &crate::cancellation::CancellationToken,
) -> Result<BlendResult, BlendError> {
    crate::fillet_edges_internal(brep, specs, options, tolerance, tracking, Some(cancel), None)
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
