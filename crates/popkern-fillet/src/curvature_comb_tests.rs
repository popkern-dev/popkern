use super::*;
use crate::center_curve::fit_nurbs_through_points;
use popkern_core::geometry::surface::SurfaceGeometry;
use std::f64::consts::PI;

/// A NURBS-fit line has zero curvature everywhere (to within fit
/// noise). Comb max curvature should be effectively zero.
#[test]
fn straight_line_comb_has_near_zero_curvature() {
    let pts: Vec<Point3> = (0..20)
        .map(|i| {
            let t = i as f64 / 19.0;
            Point3::new(t, 2.0 * t, -t)
        })
        .collect();
    let curve = fit_nurbs_through_points(&pts).expect("line fit should succeed");

    let comb = curve_curvature_comb(&curve, 10);
    assert_eq!(comb.samples.len(), 10);
    // Fit noise gives tiny but non-zero curvature near the endpoints.
    assert!(
        comb.max_curvature().unwrap() < 1e-3,
        "straight line should have max κ ≈ 0, got {}",
        comb.max_curvature().unwrap()
    );
}

/// A NURBS-fit circle of radius `r` should have curvature `≈ 1/r`
/// at every interior sample (the fit is interpolating, so endpoints
/// inherit some wiggle).
#[test]
fn circle_comb_has_constant_curvature_one_over_r() {
    let r = 2.0;
    let n = 64;
    let pts: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * PI * (i as f64) / (n as f64 - 1.0);
            Point3::new(r * theta.cos(), r * theta.sin(), 0.0)
        })
        .collect();
    let curve = fit_nurbs_through_points(&pts).expect("circle fit should succeed");

    let comb = curve_curvature_comb(&curve, 20);
    let expected = 1.0 / r;
    // Interior samples only — skip the first and last which inherit
    // endpoint fit wiggle. Tolerance is 15% because a cubic NURBS fit
    // through 64 points approximates a circle, and the 3-point
    // osculating-circle estimator measures the *fit*, not the true
    // circle. The fit is accurate to ~10% on interior knots.
    for s in &comb.samples[3..17] {
        assert!(
            (s.curvature - expected).abs() / expected < 0.15,
            "circle κ should be ≈ {expected}, got {} at t={}",
            s.curvature,
            s.t
        );
    }
}

/// A smooth parabola has curvature that rises toward the vertex and
/// falls on the arms. The comb's sign_changes count for the
/// forward-difference should reveal a single peak — i.e. one sign
/// change in dκ/dt.
#[test]
fn parabola_comb_has_one_curvature_peak() {
    let n = 41;
    let pts: Vec<Point3> = (0..n)
        .map(|i| {
            let x = -2.0 + 4.0 * (i as f64) / (n as f64 - 1.0);
            Point3::new(x, x * x, 0.0)
        })
        .collect();
    let curve = fit_nurbs_through_points(&pts).expect("parabola fit should succeed");

    let comb = curve_curvature_comb(&curve, 30);
    // Expect κ > 0 everywhere and max_curvature at the centre.
    assert!(comb.max_curvature().unwrap() > 0.1);
    // Peak curvature should be near the vertex (centre of the comb).
    // Find the sample with max κ and verify it lies near the middle.
    let (peak_idx, _) = comb
        .samples
        .iter()
        .enumerate()
        .max_by(|a, b| {
            a.1.curvature
                .partial_cmp(&b.1.curvature)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap();
    let mid = comb.samples.len() / 2;
    assert!(
        peak_idx.abs_diff(mid) <= 5,
        "parabola peak should be near the middle ({mid}), found at {peak_idx}"
    );
}

/// `curvature_ratio` returns max/min. On a circle fit (≈ constant κ)
/// this should be close to 1.
#[test]
fn circle_comb_curvature_ratio_is_near_one() {
    let n = 64;
    let pts: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * PI * (i as f64) / (n as f64 - 1.0);
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    let curve = fit_nurbs_through_points(&pts).expect("circle fit should succeed");

    let comb = curve_curvature_comb(&curve, 30);
    // Trim endpoints — interior ratio only.
    let interior: Vec<f64> = comb.samples[5..25].iter().map(|s| s.curvature).collect();
    let min_i = interior.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_i = interior.iter().cloned().fold(0.0_f64, f64::max);
    assert!(min_i > 0.0);
    assert!(
        max_i / min_i < 1.1,
        "interior circle κ ratio should be ≈ 1, got {}",
        max_i / min_i
    );
}

/// A flat cubic curve — `jumps()` should report small deltas
/// everywhere, and `max_jump` should be small.
#[test]
fn straight_line_jumps_are_tiny() {
    let pts: Vec<Point3> = (0..20)
        .map(|i| {
            let t = i as f64 / 19.0;
            Point3::new(t, 0.0, 0.0)
        })
        .collect();
    let curve = fit_nurbs_through_points(&pts).expect("line fit should succeed");

    let comb = curve_curvature_comb(&curve, 10);
    let max = comb.max_jump();
    assert!(
        max < 1e-3,
        "straight line max jump should be tiny, got {max}"
    );
}

/// `CurvatureSample::hair_length(scale)` multiplies curvature by
/// the supplied scale factor — the identity when scale = 1.
#[test]
fn hair_length_scales_with_factor() {
    let s = CurvatureSample {
        t: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
        hair_direction: Vec3::new(0.0, 1.0, 0.0),
        curvature: 2.5,
    };
    assert!((s.hair_length(1.0) - 2.5).abs() < 1e-12);
    assert!((s.hair_length(4.0) - 10.0).abs() < 1e-12);
    assert!((s.hair_length(0.0)).abs() < 1e-12);
}

/// Empty comb accessors should return sensible defaults without
/// panicking.
#[test]
fn empty_comb_accessors_are_safe() {
    let empty = CurvatureComb { samples: vec![] };
    assert!(empty.min_curvature().is_none());
    assert!(empty.max_curvature().is_none());
    assert!(empty.curvature_ratio().is_none());
    assert!(empty.jumps().is_empty());
    assert_eq!(empty.sign_changes(), 0);
    assert_eq!(empty.max_jump(), 0.0);
}

/// Surface isoline comb on a cylindrical surface: along `v` (axial),
/// the iso-line is a straight ruling, so curvature is ~0. Along `u`
/// (circumferential), the iso-line is a circle of radius `r`, so
/// curvature is `≈ 1/r`.
#[test]
fn cylinder_isoline_combs_match_theory() {
    let r = 3.0;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    };

    // Along v: axial ruling — straight.
    let axial = surface_isoline_curvature_comb(&cyl, CombAlong::V, 0.5, (0.0, 2.0), 10);
    assert!(
        axial.max_curvature().unwrap() < 1e-3,
        "axial cylinder ruling should be straight, got κ_max = {}",
        axial.max_curvature().unwrap()
    );

    // Along u: circumferential — a circle of radius r, so κ = 1/r.
    let circ = surface_isoline_curvature_comb(&cyl, CombAlong::U, 0.5, (0.0, 2.0 * PI), 30);
    let expected = 1.0 / r;
    // Interior samples only.
    for s in &circ.samples[5..25] {
        assert!(
            (s.curvature - expected).abs() / expected < 0.15,
            "cylinder circumference κ should be ≈ {expected}, got {} at s={}",
            s.curvature,
            s.t
        );
    }
}

/// Small `n_samples` is clamped to 2.
#[test]
fn tiny_n_samples_clamps_to_two() {
    let pts: Vec<Point3> = (0..10)
        .map(|i| {
            let t = i as f64 / 9.0;
            Point3::new(t, t * t, 0.0)
        })
        .collect();
    let curve = fit_nurbs_through_points(&pts).expect("parabola fit should succeed");

    let comb = curve_curvature_comb(&curve, 0);
    assert_eq!(comb.samples.len(), 2);
    let comb_one = curve_curvature_comb(&curve, 1);
    assert_eq!(comb_one.samples.len(), 2);
}
