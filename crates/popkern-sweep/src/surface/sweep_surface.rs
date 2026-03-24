//! Sweep surface generation via cross-section skinning.
//!
//! For each smooth profile span, transforms the profile control points
//! by each frame in the sequence (position + normal/binormal orientation),
//! applies scaling, and skins the result into a NURBS surface.

use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::frame::{FrameSequence, PathFrame};
use crate::profile::AnalyzedProfile;
use crate::types::ScaleMethod;

/// Result of sweep surface generation.
#[derive(Debug)]
pub struct SweepSurfaceResult {
    /// One NURBS surface per smooth profile span.
    pub surfaces: Vec<NurbsSurface>,
    /// Start cap profile (transformed to path start).
    pub start_cap: Vec<NurbsCurve>,
    /// End cap profile (transformed to path end).
    pub end_cap: Vec<NurbsCurve>,
}

/// Generate sweep surface(s) from a profile and frame sequence.
///
/// Each profile span is transformed by every frame in the sequence,
/// producing a grid of control points that becomes a NURBS surface.
/// The u-direction follows the profile and the v-direction follows the path.
pub fn generate_sweep_surface(
    profile: &AnalyzedProfile,
    frames: &FrameSequence,
    scale: &ScaleMethod,
    _tolerance: &Tolerance,
) -> Result<SweepSurfaceResult, SweepError> {
    let n_frames = frames.frames.len();
    if n_frames < 2 {
        return Err(SweepError::InvalidPath(
            "need at least 2 frames for sweep".into(),
        ));
    }

    if profile.spans.is_empty() {
        return Err(SweepError::InvalidProfile("no profile spans".into()));
    }

    let total_arc = frames.s_values.last().copied().unwrap_or(1.0).max(1e-30);

    let mut surfaces = Vec::with_capacity(profile.spans.len());

    for span in &profile.spans {
        let surface = skin_span(span, frames, scale, total_arc)?;
        surfaces.push(surface);
    }

    let start_cap = transform_profile_at_frame(&profile.spans, &frames.frames[0], scale, 0.0);
    let end_cap =
        transform_profile_at_frame(&profile.spans, &frames.frames[n_frames - 1], scale, 1.0);

    Ok(SweepSurfaceResult {
        surfaces,
        start_cap,
        end_cap,
    })
}

/// Skin a single profile span into a NURBS surface using the frame sequence.
fn skin_span(
    span: &NurbsCurve,
    frames: &FrameSequence,
    scale: &ScaleMethod,
    total_arc: f64,
) -> Result<NurbsSurface, SweepError> {
    let num_u = span.control_points.len();
    let num_v = frames.frames.len();

    if num_u == 0 {
        return Err(SweepError::InvalidProfile(
            "span has no control points".into(),
        ));
    }

    // Build control point grid in row-major order: for each u (profile cp),
    // store all v (frame) values. This matches NurbsSurface layout where
    // cp(i_u, i_v) = control_points[i_u * num_v + i_v].
    let mut control_points = Vec::with_capacity(num_u * num_v);
    let mut weights = Vec::with_capacity(num_u * num_v);

    for (cp_idx, cp) in span.control_points.iter().enumerate() {
        for (frame_idx, frame) in frames.frames.iter().enumerate() {
            let s_frac = frames.s_values[frame_idx] / total_arc;
            let s_factor = evaluate_scale(scale, s_frac);
            let transformed = transform_point(cp, frame, s_factor);
            control_points.push(transformed);
            weights.push(span.weights[cp_idx]);
        }
    }

    let degree_u = span.degree;
    let knots_u = span.knots.clone();

    let degree_v = compute_v_degree(num_v);
    let knots_v = clamped_uniform_knots(num_v, degree_v);

    Ok(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights,
    })
}

/// Transform a profile control point by a frame with scaling.
///
/// Profile points are in local coordinates where x maps to normal
/// and y maps to binormal. The z component is ignored (profile is 2D).
fn transform_point(cp: &Point3, frame: &PathFrame, scale: f64) -> Point3 {
    let x = cp.x;
    let y = cp.y;
    let offset_n = frame.normal * (scale * x);
    let offset_b = frame.binormal * (scale * y);
    frame.position + offset_n + offset_b
}

/// Evaluate scale factor at a given arc-length fraction.
fn evaluate_scale(method: &ScaleMethod, s_frac: f64) -> f64 {
    match method {
        ScaleMethod::None => 1.0,
        ScaleMethod::Uniform(pairs) => interpolate_pairs(pairs, s_frac),
    }
}

/// Linearly interpolate a value from sorted (parameter, value) pairs.
fn interpolate_pairs(pairs: &[(f64, f64)], t: f64) -> f64 {
    if pairs.is_empty() {
        return 1.0;
    }
    if pairs.len() == 1 || t <= pairs[0].0 {
        return pairs[0].1;
    }
    if t >= pairs[pairs.len() - 1].0 {
        return pairs[pairs.len() - 1].1;
    }

    let idx = pairs.partition_point(|&(s, _)| s < t);
    if idx == 0 {
        return pairs[0].1;
    }

    let (s0, v0) = pairs[idx - 1];
    let (s1, v1) = pairs[idx];
    let denom = s1 - s0;
    if denom.abs() < 1e-30 {
        return v0;
    }
    let frac = (t - s0) / denom;
    v0 + frac * (v1 - v0)
}

/// Choose v-direction degree based on the number of frames.
fn compute_v_degree(num_v: usize) -> usize {
    if num_v <= 2 {
        1
    } else if num_v <= 4 {
        num_v - 1
    } else {
        3
    }
}

/// Generate clamped uniform knot vector for `n` control points with given degree.
fn clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1; // total number of knots
    let mut knots = vec![0.0; degree + 1];

    // Interior knots
    let n_interior = m.saturating_sub(2 * (degree + 1));
    for i in 1..=n_interior {
        knots.push(i as f64 / (n_interior + 1) as f64);
    }

    // Last (degree + 1) knots = 1.0
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

/// Transform profile spans at a single frame (for cap generation).
fn transform_profile_at_frame(
    spans: &[NurbsCurve],
    frame: &PathFrame,
    scale: &ScaleMethod,
    s_frac: f64,
) -> Vec<NurbsCurve> {
    let s_factor = evaluate_scale(scale, s_frac);
    spans
        .iter()
        .map(|span| {
            let control_points = span
                .control_points
                .iter()
                .map(|cp| transform_point(cp, frame, s_factor))
                .collect();
            NurbsCurve {
                degree: span.degree,
                knots: span.knots.clone(),
                control_points,
                weights: span.weights.clone(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::curve::CurveEval;

    /// Build a minimal frame sequence along the Z axis.
    fn z_axis_frames(n: usize, length: f64) -> FrameSequence {
        let mut frames = Vec::with_capacity(n);
        let mut s_values = Vec::with_capacity(n);

        for i in 0..n {
            let frac = i as f64 / (n - 1).max(1) as f64;
            let z = frac * length;
            frames.push(PathFrame {
                position: Point3::new(0.0, 0.0, z),
                tangent: Vec3::new(0.0, 0.0, 1.0),
                normal: Vec3::new(1.0, 0.0, 0.0),
                binormal: Vec3::new(0.0, 1.0, 0.0),
            });
            s_values.push(frac * length);
        }

        FrameSequence { frames, s_values }
    }

    /// A simple unit circle profile span (centered at origin).
    fn circle_span() -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve {
            degree: 2,
            knots: vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(1.0, -1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
        }
    }

    /// A simple line segment profile along the X axis.
    fn line_span() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn sweep_circle_along_z_axis() {
        let profile = AnalyzedProfile {
            spans: vec![circle_span()],
            is_closed: true,
            centroid: Point3::new(0.0, 0.0, 0.0),
        };
        let frames = z_axis_frames(10, 5.0);
        let tol = Tolerance::default();
        let result = generate_sweep_surface(&profile, &frames, &ScaleMethod::None, &tol);
        assert!(result.is_ok());
        let result = result.unwrap(); // ok in tests
        assert_eq!(result.surfaces.len(), 1);

        let surf = &result.surfaces[0];
        assert_eq!(surf.num_u, 9);
        assert_eq!(surf.num_v, 10);

        // Check that the midpoint of the sweep (v=0.5) at u=0 is at approximately (1,0,2.5).
        use popkern_core::geometry::surface::SurfaceEval;
        let p = surf.point_at(0.0, 0.5);
        assert!((p.x - 1.0).abs() < 0.2, "expected x≈1.0, got {}", p.x);
        assert!(p.y.abs() < 0.2, "expected y≈0.0, got {}", p.y);
        assert!((p.z - 2.5).abs() < 0.5, "expected z≈2.5, got {}", p.z);
    }

    #[test]
    fn sweep_line_along_z_axis() {
        let profile = AnalyzedProfile {
            spans: vec![line_span()],
            is_closed: false,
            centroid: Point3::new(0.0, 0.0, 0.0),
        };
        let frames = z_axis_frames(5, 3.0);
        let tol = Tolerance::default();
        let result = generate_sweep_surface(&profile, &frames, &ScaleMethod::None, &tol);
        assert!(result.is_ok());
        let result = result.unwrap(); // ok in tests

        let surf = &result.surfaces[0];
        assert_eq!(surf.num_u, 2);
        assert_eq!(surf.num_v, 5);

        // At (u=0, v=0) should be (-1, 0, 0).
        use popkern_core::geometry::surface::SurfaceEval;
        let p = surf.point_at(0.0, 0.0);
        assert!((p.x - (-1.0)).abs() < 0.1, "expected x≈-1.0, got {}", p.x);
        assert!(p.y.abs() < 0.1, "expected y≈0, got {}", p.y);
        assert!(p.z.abs() < 0.1, "expected z≈0, got {}", p.z);
    }

    #[test]
    fn sweep_with_uniform_scale() {
        let profile = AnalyzedProfile {
            spans: vec![line_span()],
            is_closed: false,
            centroid: Point3::new(0.0, 0.0, 0.0),
        };
        let frames = z_axis_frames(5, 3.0);
        let tol = Tolerance::default();
        let scale = ScaleMethod::Uniform(vec![(0.0, 2.0), (1.0, 2.0)]);
        let result = generate_sweep_surface(&profile, &frames, &scale, &tol);
        assert!(result.is_ok());
        let result = result.unwrap(); // ok in tests

        // Control points should be 2x further from the path center.
        let surf = &result.surfaces[0];
        // First control point of first frame should be at (-2, 0, 0).
        let first_cp = surf.control_points[0];
        assert!(
            (first_cp.x - (-2.0)).abs() < 0.01,
            "expected x≈-2.0 with scale=2, got {}",
            first_cp.x
        );
    }

    #[test]
    fn too_few_frames_returns_error() {
        let profile = AnalyzedProfile {
            spans: vec![line_span()],
            is_closed: false,
            centroid: Point3::new(0.0, 0.0, 0.0),
        };
        let frames = FrameSequence {
            frames: vec![PathFrame {
                position: Point3::new(0.0, 0.0, 0.0),
                tangent: Vec3::new(0.0, 0.0, 1.0),
                normal: Vec3::new(1.0, 0.0, 0.0),
                binormal: Vec3::new(0.0, 1.0, 0.0),
            }],
            s_values: vec![0.0],
        };
        let tol = Tolerance::default();
        let result = generate_sweep_surface(&profile, &frames, &ScaleMethod::None, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn empty_profile_returns_error() {
        let profile = AnalyzedProfile {
            spans: vec![],
            is_closed: false,
            centroid: Point3::new(0.0, 0.0, 0.0),
        };
        let frames = z_axis_frames(5, 3.0);
        let tol = Tolerance::default();
        let result = generate_sweep_surface(&profile, &frames, &ScaleMethod::None, &tol);
        assert!(result.is_err());
    }

    #[test]
    fn caps_are_generated() {
        let profile = AnalyzedProfile {
            spans: vec![circle_span()],
            is_closed: true,
            centroid: Point3::new(0.0, 0.0, 0.0),
        };
        let frames = z_axis_frames(5, 3.0);
        let tol = Tolerance::default();
        let result = generate_sweep_surface(&profile, &frames, &ScaleMethod::None, &tol);
        assert!(result.is_ok());
        let result = result.unwrap(); // ok in tests
        assert_eq!(result.start_cap.len(), 1);
        assert_eq!(result.end_cap.len(), 1);

        // Start cap should be at z=0, end cap at z=3.
        let start_p = result.start_cap[0].point_at(0.0);
        assert!(start_p.z.abs() < 0.01, "start cap z≈0, got {}", start_p.z);

        let end_p = result.end_cap[0].point_at(0.0);
        assert!((end_p.z - 3.0).abs() < 0.01, "end cap z≈3, got {}", end_p.z);
    }

    #[test]
    fn clamped_knots_degree_1() {
        let knots = clamped_uniform_knots(3, 1);
        assert_eq!(knots.len(), 5); // n + degree + 1 = 3 + 1 + 1 = 5
        assert_eq!(knots[0], 0.0);
        assert_eq!(knots[1], 0.0);
        assert_eq!(knots[3], 1.0);
        assert_eq!(knots[4], 1.0);
    }

    #[test]
    fn clamped_knots_degree_3() {
        let knots = clamped_uniform_knots(6, 3);
        // Length = 6 + 3 + 1 = 10
        assert_eq!(knots.len(), 10);
        // First 4 = 0.0
        for k in &knots[..4] {
            assert_eq!(*k, 0.0);
        }
        // Last 4 = 1.0
        for k in &knots[6..] {
            assert_eq!(*k, 1.0);
        }
        // 2 interior knots
        assert!(knots[4] > 0.0 && knots[4] < 1.0);
        assert!(knots[5] > 0.0 && knots[5] < 1.0);
    }

    #[test]
    fn interpolate_pairs_basic() {
        let pairs = vec![(0.0, 1.0), (0.5, 2.0), (1.0, 3.0)];
        assert!((interpolate_pairs(&pairs, 0.0) - 1.0).abs() < 1e-10);
        assert!((interpolate_pairs(&pairs, 0.25) - 1.5).abs() < 1e-10);
        assert!((interpolate_pairs(&pairs, 0.5) - 2.0).abs() < 1e-10);
        assert!((interpolate_pairs(&pairs, 1.0) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn interpolate_pairs_empty() {
        assert!((interpolate_pairs(&[], 0.5) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn interpolate_pairs_single() {
        let pairs = vec![(0.5, 7.0)];
        assert!((interpolate_pairs(&pairs, 0.0) - 7.0).abs() < 1e-10);
        assert!((interpolate_pairs(&pairs, 1.0) - 7.0).abs() < 1e-10);
    }

    #[test]
    fn v_degree_selection() {
        assert_eq!(compute_v_degree(2), 1);
        assert_eq!(compute_v_degree(3), 2);
        assert_eq!(compute_v_degree(4), 3);
        assert_eq!(compute_v_degree(5), 3);
        assert_eq!(compute_v_degree(100), 3);
    }
}
