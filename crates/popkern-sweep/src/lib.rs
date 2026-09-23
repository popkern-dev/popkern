#![warn(missing_docs)]
//! popkern-sweep: Sweep and loft operations for the popkern BRep kernel.
//!
//! Provides sweep, loft, pipe, helical sweep, guided sweep, and swept blend
//! operations that produce BRep solids from profile curves and paths.

pub mod arc_length;
pub mod boolean;
pub mod error;
pub mod frame;
pub(crate) mod frame_budget;
pub mod helix;
pub mod interactive;
pub mod pipe;
pub(crate) mod pipeline;
pub mod profile;
pub mod revolve;
pub mod self_intersect;
pub mod surface;
pub mod topology;
pub mod types;

pub use boolean::{loft_add, loft_cut, pipe_add, pipe_cut, sweep_add, sweep_cut, sweep_intersect};
pub use error::{SweepError, SweepWarning};
pub use interactive::{
    guided_loft_with_control, guided_sweep_with_control, helical_sweep_with_control,
    loft_with_control, pipe_with_control, revolve_with_control, sweep_with_control,
    swept_blend_with_control,
};
pub use revolve::revolve;
pub use types::*;

// Re-exported for `lib_tests.rs` (historic unqualified uses).
#[cfg(test)]
pub(crate) use frame_budget::{MIN_FRAMES, compute_sample_count};
#[cfg(test)]
pub(crate) use pipeline::check_path_closed;

use popkern_core::geometry::NurbsCurve;
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;

/// Sweep a profile along a path, producing a BRep solid.
///
/// The profile is analyzed for tangent discontinuities, centered at the
/// origin, and swept along the path using orientation frames. The resulting
/// NURBS surfaces are assembled into a BRep with cap faces at each end
/// (for closed profiles) or a sheet body (for open profiles).
///
/// # Errors
///
/// Returns `SweepError` if the profile or path is invalid, frame computation
/// fails, or self-intersection is detected (unless allowed by options).
pub fn sweep(
    profile: &SweepProfile,
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    pipeline::sweep_impl(
        profile,
        path,
        options,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Loft through a sequence of cross-sections, producing a BRep solid.
///
/// Sections are made compatible, parameterized by chord length, and
/// interpolated to create a NURBS surface. Cap faces are added for
/// closed sections.
///
/// # Errors
///
/// Returns `SweepError` if fewer than 2 sections are provided or
/// section compatibility cannot be achieved.
pub fn loft(
    sections: &[LoftSection],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    pipeline::loft_impl(
        sections,
        &[],
        options,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Create a pipe by sweeping a circular cross-section along a path.
///
/// Supports both constant and variable radius pipes.
///
/// # Errors
///
/// Returns `SweepError` if the radius is invalid, path is degenerate,
/// or frame computation fails.
pub fn pipe(
    path: &NurbsCurve,
    radius: &PipeRadius,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    pipeline::pipe_impl(
        path,
        radius,
        options,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Sweep a profile along a helical path.
///
/// Builds a helical path curve from the specification, then performs a
/// standard sweep along that path.
///
/// # Errors
///
/// Returns `SweepError` if the helix specification is invalid or the
/// sweep fails.
pub fn helical_sweep(
    profile: &SweepProfile,
    spec: &HelixSpec,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    let (path, mut warnings) = helix::build_helix_path(spec)?;
    let mut result = sweep(profile, &path, options, tolerance, tracking)?;
    result.warnings.append(&mut warnings);
    Ok(result)
}

/// Sweep a profile along a path with guide rail constraints.
///
/// After generating the sweep surface, validates that it passes through
/// the guide rail curves within tolerance.
///
/// # Errors
///
/// Returns `SweepError` if the sweep fails or the guide rails deviate
/// beyond tolerance.
pub fn guided_sweep(
    profile: &SweepProfile,
    path: &NurbsCurve,
    guides: &[NurbsCurve],
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    pipeline::guided_sweep_impl(
        profile,
        path,
        guides,
        options,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Loft through sections with guide rail constraints.
///
/// # Errors
///
/// Returns `SweepError` if the loft fails or guide rails deviate.
pub fn guided_loft(
    sections: &[LoftSection],
    guides: &[NurbsCurve],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    pipeline::loft_impl(
        sections,
        guides,
        options,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Generate a swept blend from multiple profiles along a path.
///
/// Each profile is positioned at its specified arc-length fraction along
/// the path, then the positioned sections are lofted together.
///
/// # Errors
///
/// Returns `SweepError` if fewer than 2 profiles are provided or loft fails.
pub fn swept_blend(
    profiles: &[(f64, SweepProfile)],
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<SweepResult, SweepError> {
    pipeline::swept_blend_impl(
        profiles,
        path,
        options,
        tolerance,
        tracking,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
