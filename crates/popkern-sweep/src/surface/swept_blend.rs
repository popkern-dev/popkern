//! Swept blend: multiple profiles interpolated along a path.
//!
//! A swept blend positions multiple cross-section profiles along a sweep
//! path using orientation frames, then lofts through the positioned sections
//! to create the final surface.

use popkern_core::geometry::NurbsCurve;
use popkern_core::tolerance::Tolerance;

use crate::arc_length::ArcLengthTable;
use crate::error::SweepError;
use crate::frame::compute_frames;
use crate::profile::analyze_profile;
use crate::surface::loft_surface::{LoftSurfaceResult, generate_loft_surface};
use crate::types::{LoftOptions, LoftSection, SweepOptions, SweepProfile};

/// Minimum number of frames to compute along the path.
const MIN_FRAMES: usize = 10;

/// Generate a swept blend surface from multiple profiles along a path.
///
/// Each profile is specified as a (t, profile) pair where `t` is the
/// arc-length fraction along the path (0.0 = start, 1.0 = end) and
/// `profile` is the cross-section at that position.
///
/// The profiles are positioned using orientation frames from the path,
/// then lofted together to produce the final surface.
///
/// # Errors
///
/// Returns an error if fewer than 2 profiles are provided, if any profile
/// is invalid, or if the loft fails.
pub fn generate_swept_blend_surface(
    profiles: &[(f64, SweepProfile)],
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: &Tolerance,
) -> Result<LoftSurfaceResult, SweepError> {
    if profiles.len() < 2 {
        return Err(SweepError::IncompatibleSections(format!(
            "swept blend needs at least 2 profiles, got {}",
            profiles.len()
        )));
    }

    // Validate all t values are in [0, 1] and sorted.
    for (i, &(t, _)) in profiles.iter().enumerate() {
        if !(0.0..=1.0).contains(&t) {
            return Err(SweepError::InvalidPath(format!(
                "profile {} has t={} outside [0, 1]",
                i, t
            )));
        }
        if i > 0 && t < profiles[i - 1].0 {
            return Err(SweepError::InvalidPath(
                "profile t values must be non-decreasing".into(),
            ));
        }
    }

    // Build arc-length table and frames.
    let table = ArcLengthTable::build(path, 100);
    let n_frames = MIN_FRAMES.max(profiles.len() * 5);
    let (frames, _warnings) =
        compute_frames(path, &table, &options.frame, n_frames, options.twist, false)?;

    // For each profile, find the frame at the given t value and transform
    // the profile into world coordinates.
    let mut sections = Vec::with_capacity(profiles.len());

    for &(t, ref profile) in profiles {
        let analyzed = analyze_profile(profile, tolerance)?;

        // Find the frame closest to arc-length fraction t.
        let frame_idx = find_nearest_frame(&frames.s_values, t, table.total_length);
        let frame = &frames.frames[frame_idx];

        // Transform the analyzed profile spans into world coordinates and
        // combine them into a single section curve.
        let section_curve = if analyzed.spans.is_empty() {
            return Err(SweepError::InvalidProfile("profile has no spans".into()));
        } else {
            // Use the first span as the section curve (main outer wire).
            transform_span_to_world(&analyzed.spans[0], frame)
        };

        sections.push(LoftSection {
            curve: section_curve,
            position: Some(frame.position),
        });
    }

    // Loft through the positioned sections.
    let loft_options = LoftOptions::default();
    generate_loft_surface(&sections, &loft_options, tolerance)
}

/// Find the frame index closest to a given arc-length fraction.
fn find_nearest_frame(s_values: &[f64], t_frac: f64, total_arc: f64) -> usize {
    let target_s = t_frac * total_arc;
    let mut best_idx = 0;
    let mut best_diff = f64::MAX;

    for (i, &s) in s_values.iter().enumerate() {
        let diff = (s - target_s).abs();
        if diff < best_diff {
            best_diff = diff;
            best_idx = i;
        }
    }

    best_idx
}

/// Transform a profile span from local frame coordinates to world coordinates.
fn transform_span_to_world(span: &NurbsCurve, frame: &crate::frame::PathFrame) -> NurbsCurve {
    let control_points = span
        .control_points
        .iter()
        .map(|cp| {
            let offset_n = frame.normal * cp.x;
            let offset_b = frame.binormal * cp.y;
            frame.position + offset_n + offset_b
        })
        .collect();

    NurbsCurve {
        degree: span.degree,
        knots: span.knots.clone(),
        control_points,
        weights: span.weights.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    /// A straight-line path along the Z axis.
    fn z_axis_path() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 10.0)],
            weights: vec![1.0, 1.0],
        }
    }

    /// A unit circle in the XY plane.
    fn circle_profile(radius: f64) -> SweepProfile {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let r = radius;
        SweepProfile::Wire(NurbsCurve {
            degree: 2,
            knots: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            control_points: vec![
                Point3::new(r, 0.0, 0.0),
                Point3::new(r, r, 0.0),
                Point3::new(0.0, r, 0.0),
                Point3::new(-r, r, 0.0),
                Point3::new(-r, 0.0, 0.0),
                Point3::new(-r, -r, 0.0),
                Point3::new(0.0, -r, 0.0),
                Point3::new(r, -r, 0.0),
                Point3::new(r, 0.0, 0.0),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        })
    }

    #[test]
    fn swept_blend_two_identical_circles() {
        let path = z_axis_path();
        let profiles = vec![(0.0, circle_profile(1.0)), (1.0, circle_profile(1.0))];
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_swept_blend_surface(&profiles, &path, &options, &tol);
        assert!(result.is_ok(), "swept blend failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        // Should produce a surface with sections at z=0 and z=10.
        use popkern_core::geometry::surface::SurfaceEval;
        let p_start = result.surface.point_at(0.0, 0.0);
        let p_end = result.surface.point_at(0.0, 1.0);
        // Start and end z should be near 0 and 10.
        assert!(
            p_start.z.abs() < 1.0,
            "expected z~0 at start, got {}",
            p_start.z
        );
        assert!(
            (p_end.z - 10.0).abs() < 1.0,
            "expected z~10 at end, got {}",
            p_end.z
        );
    }

    #[test]
    fn swept_blend_expanding_tube() {
        let path = z_axis_path();
        let profiles = vec![(0.0, circle_profile(1.0)), (1.0, circle_profile(3.0))];
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_swept_blend_surface(&profiles, &path, &options, &tol);
        assert!(result.is_ok(), "expanding tube failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        use popkern_core::geometry::surface::SurfaceEval;

        // At v=0, radius should be ~1.0.
        let p0 = result.surface.point_at(0.0, 0.0);
        let r0 = (p0.x * p0.x + p0.y * p0.y).sqrt();
        assert!(
            (r0 - 1.0).abs() < 0.5,
            "expected radius ~1.0 at start, got {}",
            r0
        );

        // At v=1, radius should be ~3.0.
        let p1 = result.surface.point_at(0.0, 1.0);
        let r1 = (p1.x * p1.x + p1.y * p1.y).sqrt();
        assert!(
            (r1 - 3.0).abs() < 0.5,
            "expected radius ~3.0 at end, got {}",
            r1
        );
    }

    #[test]
    fn swept_blend_too_few_profiles() {
        let path = z_axis_path();
        let profiles = vec![(0.5, circle_profile(1.0))];
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_swept_blend_surface(&profiles, &path, &options, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn swept_blend_invalid_t_range() {
        let path = z_axis_path();
        let profiles = vec![(0.0, circle_profile(1.0)), (1.5, circle_profile(1.0))];
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_swept_blend_surface(&profiles, &path, &options, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn swept_blend_unsorted_t_returns_error() {
        let path = z_axis_path();
        let profiles = vec![(0.8, circle_profile(1.0)), (0.2, circle_profile(1.0))];
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_swept_blend_surface(&profiles, &path, &options, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn find_nearest_frame_boundary() {
        let s_values = vec![0.0, 2.5, 5.0, 7.5, 10.0];
        assert_eq!(find_nearest_frame(&s_values, 0.0, 10.0), 0);
        assert_eq!(find_nearest_frame(&s_values, 1.0, 10.0), 4);
        assert_eq!(find_nearest_frame(&s_values, 0.5, 10.0), 2);
    }
}
