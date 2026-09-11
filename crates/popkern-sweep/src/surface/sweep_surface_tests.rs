//! Tests for sweep-surface skinning (split per file-size rule).

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
        wires: vec![1],
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
        wires: vec![1],
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
        wires: vec![1],
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
        wires: vec![1],
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
        wires: vec![1],
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
        wires: vec![1],
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
