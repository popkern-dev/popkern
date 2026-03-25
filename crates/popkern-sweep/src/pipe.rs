//! Pipe (circular sweep) operations.
//!
//! Generates a pipe surface by sweeping a circular cross-section along a
//! NURBS path curve. Supports both constant and variable radius pipes.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;

use crate::arc_length::ArcLengthTable;
use crate::error::SweepError;
use crate::frame::compute_frames;
use crate::profile::analyze_profile;
use crate::surface::sweep_surface::{SweepSurfaceResult, generate_sweep_surface};
use crate::types::{PipeRadius, ScaleMethod, SweepOptions, SweepProfile};

/// Number of frames to sample along the path for pipe generation.
///
/// A moderate count balances surface fidelity against file size.
/// Too many frames produce oversized NURBS surfaces that some CAD
/// viewers struggle with; too few lose path curvature detail.
const PIPE_FRAME_COUNT: usize = 17;

/// Generate a pipe surface by sweeping a circle along a path.
///
/// For constant radius, sweeps an exact rational NURBS circle.
/// For variable radius, uses uniform scaling to vary the circle size
/// along the path.
pub fn generate_pipe_surface(
    path: &NurbsCurve,
    radius: &PipeRadius,
    options: &SweepOptions,
    tolerance: &Tolerance,
) -> Result<SweepSurfaceResult, SweepError> {
    let base_radius = match radius {
        PipeRadius::Constant(r) => {
            if *r <= 0.0 {
                return Err(SweepError::InvalidProfile(
                    "pipe radius must be positive".into(),
                ));
            }
            *r
        }
        PipeRadius::Variable(pts) => {
            if pts.is_empty() {
                return Err(SweepError::InvalidProfile(
                    "variable radius must have at least one point".into(),
                ));
            }
            let r = pts[0].1;
            if r <= 0.0 {
                return Err(SweepError::InvalidProfile(
                    "pipe radius must be positive".into(),
                ));
            }
            r
        }
    };

    let circle = build_circle_nurbs(base_radius);
    let profile = SweepProfile::Wire(circle);
    let analyzed = analyze_profile(&profile, tolerance)?;

    let table = ArcLengthTable::build(path, 100);
    let (frames, _warnings) = compute_frames(
        path,
        &table,
        &options.frame,
        PIPE_FRAME_COUNT,
        options.twist,
        false,
    )?;

    let scale = match radius {
        PipeRadius::Constant(_) => ScaleMethod::None,
        PipeRadius::Variable(pts) => {
            let base = pts[0].1;
            if base.abs() < 1e-30 {
                return Err(SweepError::InvalidProfile(
                    "base radius too small for variable pipe".into(),
                ));
            }
            ScaleMethod::Uniform(pts.iter().map(|&(t, r)| (t, r / base)).collect())
        }
    };

    generate_sweep_surface(&analyzed, &frames, &scale, tolerance)
}

/// Build a rational degree-2 NURBS circle in the XY plane.
///
/// Uses 9 control points (4 quarter arcs) with weights of 1/sqrt(2)
/// on the midpoints to produce an exact circle.
fn build_circle_nurbs(radius: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = radius;
    NurbsCurve {
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    /// A straight-line path along the Z axis from (0,0,0) to (0,0,5).
    fn z_axis_path() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn constant_radius_pipe_along_z() {
        let path = z_axis_path();
        let radius = PipeRadius::Constant(2.0);
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_pipe_surface(&path, &radius, &options, &tol);
        assert!(result.is_ok(), "pipe failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        assert!(!result.surfaces.is_empty());

        // Evaluate at (u=0, v=0.5) — should be near radius 2 from center at z~2.5.
        let surf = &result.surfaces[0];
        let p = surf.point_at(0.0, 0.5);
        let radial_dist = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (radial_dist - 2.0).abs() < 0.5,
            "expected radius ~2.0, got {}",
            radial_dist
        );
    }

    #[test]
    fn variable_radius_pipe() {
        let path = z_axis_path();
        let radius = PipeRadius::Variable(vec![(0.0, 1.0), (1.0, 3.0)]);
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_pipe_surface(&path, &radius, &options, &tol);
        assert!(result.is_ok(), "variable pipe failed: {:?}", result.err());

        let result = result.unwrap(); // ok in tests
        assert!(!result.surfaces.is_empty());

        // At v=0, radius should be ~1.0; at v=1.0, radius should be ~3.0.
        let surf = &result.surfaces[0];
        let p_start = surf.point_at(0.0, 0.0);
        let r_start = (p_start.x * p_start.x + p_start.y * p_start.y).sqrt();
        assert!(
            (r_start - 1.0).abs() < 0.5,
            "expected start radius ~1.0, got {}",
            r_start
        );

        let p_end = surf.point_at(0.0, 1.0);
        let r_end = (p_end.x * p_end.x + p_end.y * p_end.y).sqrt();
        assert!(
            (r_end - 3.0).abs() < 0.5,
            "expected end radius ~3.0, got {}",
            r_end
        );
    }

    #[test]
    fn zero_radius_returns_error() {
        let path = z_axis_path();
        let radius = PipeRadius::Constant(0.0);
        let options = SweepOptions::default();
        let tol = Tolerance::default();
        let result = generate_pipe_surface(&path, &radius, &options, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn negative_radius_returns_error() {
        let path = z_axis_path();
        let radius = PipeRadius::Constant(-1.0);
        let options = SweepOptions::default();
        let tol = Tolerance::default();
        let result = generate_pipe_surface(&path, &radius, &options, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn empty_variable_radius_returns_error() {
        let path = z_axis_path();
        let radius = PipeRadius::Variable(vec![]);
        let options = SweepOptions::default();
        let tol = Tolerance::default();
        let result = generate_pipe_surface(&path, &radius, &options, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn build_circle_has_correct_structure() {
        let circle = build_circle_nurbs(3.0);
        assert_eq!(circle.degree, 2);
        assert_eq!(circle.control_points.len(), 9);
        assert_eq!(circle.weights.len(), 9);
        assert_eq!(circle.knots.len(), 12);

        // First and last control points should coincide (closed).
        let first = &circle.control_points[0];
        let last = &circle.control_points[8];
        assert!((first.x - last.x).abs() < 1e-10);
        assert!((first.y - last.y).abs() < 1e-10);

        // Evaluate at t=0 should be at (3, 0, 0).
        let p = circle.point_at(0.0);
        assert!((p.x - 3.0).abs() < 0.01, "expected x=3.0, got {}", p.x);
        assert!(p.y.abs() < 0.01, "expected y=0.0, got {}", p.y);
    }

    #[test]
    fn pipe_generates_caps() {
        let path = z_axis_path();
        let radius = PipeRadius::Constant(1.0);
        let options = SweepOptions::default();
        let tol = Tolerance::default();

        let result = generate_pipe_surface(&path, &radius, &options, &tol);
        let result = result.unwrap(); // ok in tests
        assert!(!result.start_cap.is_empty());
        assert!(!result.end_cap.is_empty());
    }
}
