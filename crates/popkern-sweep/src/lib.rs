#![warn(missing_docs)]
//! popkern-sweep: Sweep and loft operations for the popkern BRep kernel.
//!
//! Provides sweep, loft, pipe, helical sweep, guided sweep, and swept blend
//! operations that produce BRep solids from profile curves and paths.

pub mod arc_length;
pub mod boolean;
pub mod error;
pub mod frame;
pub mod helix;
pub mod pipe;
pub mod profile;
pub mod revolve;
pub mod self_intersect;
pub mod surface;
pub mod topology;
pub mod types;

pub use boolean::{loft_add, loft_cut, pipe_add, pipe_cut, sweep_add, sweep_cut, sweep_intersect};
pub use error::{SweepError, SweepWarning};
pub use revolve::revolve;
pub use types::*;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;

/// Default number of arc-length samples for path parameterization.
const ARC_LENGTH_SAMPLES: usize = 100;

/// Minimum number of frames for sweep operations.
const MIN_FRAMES: usize = 10;

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
    // 1. Analyze profile.
    let analyzed = profile::analyze_profile(profile, &tolerance)?;

    // 2. Build arc-length table.
    let table = arc_length::ArcLengthTable::build(path, ARC_LENGTH_SAMPLES);

    // 3. Compute frames.
    let n_samples = compute_sample_count(&table);
    let is_closed_path = check_path_closed(path, &tolerance);
    let (frames, mut warnings) = frame::compute_frames(
        path,
        &table,
        &options.frame,
        n_samples,
        options.twist,
        is_closed_path,
    )?;

    // 4. Generate surfaces.
    let surface_result =
        surface::generate_sweep_surface(&analyzed, &frames, &options.scale, &tolerance)?;

    // 5. Self-intersection check.
    for surf in &surface_result.surfaces {
        if let Some(t) = self_intersect::check_self_intersection(surf, &tolerance) {
            if options.allow_self_intersection {
                warnings.push(SweepWarning::SelfIntersectionAllowed { t });
            } else {
                return Err(SweepError::SelfIntersection { t });
            }
        }
    }

    // 6. Build topology.
    let brep = topology::side_faces::build_sweep_topology(
        &surface_result.surfaces,
        &surface_result.start_cap,
        &surface_result.end_cap,
        analyzed.is_closed,
        tolerance,
        tracking,
    )?;

    Ok(SweepResult { brep, warnings })
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
    let result = surface::generate_loft_surface(sections, options, &tolerance)?;

    // Check if the first section is closed.
    let is_closed = if !sections.is_empty() {
        let curve = &sections[0].curve;
        let (t_min, t_max) = curve.parameter_range();
        let start = curve.point_at(t_min);
        let end = curve.point_at(t_max);
        tolerance.points_coincident(&start, &end)
    } else {
        false
    };

    let brep = topology::side_faces::build_sweep_topology(
        &[result.surface],
        std::slice::from_ref(&result.start_cap),
        std::slice::from_ref(&result.end_cap),
        is_closed,
        tolerance,
        tracking,
    )?;

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
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
    let surface_result = pipe::generate_pipe_surface(path, radius, options, &tolerance)?;

    let brep = topology::side_faces::build_sweep_topology(
        &surface_result.surfaces,
        &surface_result.start_cap,
        &surface_result.end_cap,
        true, // pipes are always closed profiles
        tolerance,
        tracking,
    )?;

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
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
    // Generate the base sweep.
    let analyzed = profile::analyze_profile(profile, &tolerance)?;
    let table = arc_length::ArcLengthTable::build(path, ARC_LENGTH_SAMPLES);
    let n_samples = compute_sample_count(&table);
    let is_closed_path = check_path_closed(path, &tolerance);
    let (frames, warnings) = frame::compute_frames(
        path,
        &table,
        &options.frame,
        n_samples,
        options.twist,
        is_closed_path,
    )?;
    let surface_result =
        surface::generate_sweep_surface(&analyzed, &frames, &options.scale, &tolerance)?;

    // Validate guides against each surface.
    for surf in &surface_result.surfaces {
        surface::validate_guide_rails(surf, guides, &tolerance)?;
    }

    let brep = topology::side_faces::build_sweep_topology(
        &surface_result.surfaces,
        &surface_result.start_cap,
        &surface_result.end_cap,
        analyzed.is_closed,
        tolerance,
        tracking,
    )?;

    Ok(SweepResult { brep, warnings })
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
    let result = surface::generate_loft_surface(sections, options, &tolerance)?;

    // Validate guides against the lofted surface.
    surface::validate_guide_rails(&result.surface, guides, &tolerance)?;

    let is_closed = if !sections.is_empty() {
        let curve = &sections[0].curve;
        let (t_min, t_max) = curve.parameter_range();
        let start = curve.point_at(t_min);
        let end = curve.point_at(t_max);
        tolerance.points_coincident(&start, &end)
    } else {
        false
    };

    let brep = topology::side_faces::build_sweep_topology(
        &[result.surface],
        std::slice::from_ref(&result.start_cap),
        std::slice::from_ref(&result.end_cap),
        is_closed,
        tolerance,
        tracking,
    )?;

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
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
    let result =
        surface::swept_blend::generate_swept_blend_surface(profiles, path, options, &tolerance)?;

    let brep = topology::side_faces::build_sweep_topology(
        &[result.surface],
        std::slice::from_ref(&result.start_cap),
        std::slice::from_ref(&result.end_cap),
        false, // swept blend sections may not be closed
        tolerance,
        tracking,
    )?;

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
}

/// Compute sample count based on path arc length.
fn compute_sample_count(table: &arc_length::ArcLengthTable) -> usize {
    // Roughly 10 samples per unit length, minimum MIN_FRAMES.
    let count = (table.total_length * 10.0).ceil() as usize;
    count.clamp(MIN_FRAMES, 200)
}

/// Check if a path curve is closed (start point equals end point).
fn check_path_closed(path: &NurbsCurve, tolerance: &Tolerance) -> bool {
    let (t_min, t_max) = path.parameter_range();
    let start = path.point_at(t_min);
    let end = path.point_at(t_max);
    tolerance.points_coincident(&start, &end)
}

#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
