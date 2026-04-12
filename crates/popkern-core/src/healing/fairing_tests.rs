use super::*;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::{Point3, Vec3};

/// Build a simple degree-2 NURBS representing a straight line from (0,0,0) to (1,0,0).
fn line_nurbs() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// Build a 5-point NURBS where interior points lie on a line.
fn straight_5pt_nurbs() -> NurbsCurve {
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.25, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.75, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
    }
}

#[test]
fn fair_straight_line_no_change() {
    let mut curve = straight_5pt_nurbs();
    let original = curve.control_points.clone();
    let result = fair_nurbs_curve(&mut curve, 10, 0.5, 1.0);

    // A perfectly straight control polygon should not be changed by
    // Laplacian smoothing (the average of neighbors equals the point).
    assert!(
        result.max_displacement < 1e-12,
        "straight line should have ~zero displacement, got {}",
        result.max_displacement
    );
    assert_eq!(
        result.points_adjusted, 0,
        "no points should be adjusted on a straight line"
    );

    // Verify control points unchanged.
    for (i, (orig, faired)) in original.iter().zip(curve.control_points.iter()).enumerate() {
        assert!(
            orig.distance_to(faired) < 1e-12,
            "control point {i} should be unchanged"
        );
    }
}

#[test]
fn fair_too_few_control_points_is_noop() {
    let mut curve = line_nurbs();
    let result = fair_nurbs_curve(&mut curve, 10, 0.5, 1.0);
    assert!(
        result.max_displacement < f64::EPSILON,
        "2-point curve has no interior points to smooth"
    );
    assert_eq!(result.points_adjusted, 0);
}

#[test]
fn fair_zero_iterations_is_noop() {
    let mut curve = straight_5pt_nurbs();
    let result = fair_nurbs_curve(&mut curve, 0, 0.5, 1.0);
    assert!(result.max_displacement < f64::EPSILON);
    assert_eq!(result.points_adjusted, 0);
}

#[test]
fn fair_respects_max_displacement() {
    // Create a curve with a deliberate wrinkle in the middle.
    let mut curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.25, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0), // wrinkle
            Point3::new(0.75, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
    };

    let max_disp = 0.1;
    let result = fair_nurbs_curve(&mut curve, 100, 0.5, max_disp);

    assert!(
        result.max_displacement <= max_disp + 1e-12,
        "displacement {} should not exceed max {}",
        result.max_displacement,
        max_disp
    );
}

#[test]
fn detect_wrinkles_on_smooth_circle_finds_none() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let wrinkles = detect_wrinkles(&circle, (0.0, std::f64::consts::TAU), 50, 3.0);

    // A circle has constant curvature, so no point should exceed
    // 3x the median (which equals the constant curvature).
    assert!(
        wrinkles.is_empty(),
        "smooth circle should have no wrinkles, got {} flagged points",
        wrinkles.len()
    );
}

#[test]
fn detect_wrinkles_on_line_finds_none() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let wrinkles = detect_wrinkles(&line, (0.0, 1.0), 20, 2.0);

    // A line has zero curvature everywhere; no wrinkles possible.
    assert!(wrinkles.is_empty(), "straight line should have no wrinkles");
}

#[test]
fn detect_wrinkles_too_few_samples() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let wrinkles = detect_wrinkles(&circle, (0.0, 1.0), 1, 2.0);
    assert!(
        wrinkles.is_empty(),
        "fewer than 2 samples should return empty"
    );
}

#[test]
fn energy_fair_straight_line_no_change() {
    let mut curve = straight_5pt_nurbs();
    let result = fair_nurbs_curve_energy(&mut curve, 10, 0.1, 1.0);
    // A perfectly straight control polygon: second differences are zero,
    // so energy gradient is zero everywhere.
    assert!(
        result.max_displacement < 1e-12,
        "straight line should have ~zero displacement, got {}",
        result.max_displacement
    );
    assert_eq!(result.points_adjusted, 0);
}

#[test]
fn energy_fair_too_few_control_points_is_noop() {
    let mut curve = line_nurbs();
    let result = fair_nurbs_curve_energy(&mut curve, 10, 0.1, 1.0);
    assert!(result.max_displacement < f64::EPSILON);
    assert_eq!(result.points_adjusted, 0);
}

#[test]
fn energy_fair_zero_iterations_is_noop() {
    let mut curve = straight_5pt_nurbs();
    let result = fair_nurbs_curve_energy(&mut curve, 0, 0.1, 1.0);
    assert!(result.max_displacement < f64::EPSILON);
    assert_eq!(result.points_adjusted, 0);
}

#[test]
fn energy_fair_respects_max_displacement() {
    let mut curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.25, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0), // wrinkle
            Point3::new(0.75, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
    };

    let max_disp = 0.1;
    let result = fair_nurbs_curve_energy(&mut curve, 100, 0.1, max_disp);

    assert!(
        result.max_displacement <= max_disp + 1e-12,
        "displacement {} should not exceed max {}",
        result.max_displacement,
        max_disp
    );
    assert!(
        result.points_adjusted > 0,
        "wrinkled curve should have adjusted points"
    );
}

#[test]
fn energy_fair_reduces_curvature() {
    // A curve with a deliberate kink should get smoother.
    let original_points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.25, 0.0, 0.0),
        Point3::new(0.5, 0.5, 0.0), // kink
        Point3::new(0.75, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ];
    let mut curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: original_points.clone(),
        weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
    };

    let result = fair_nurbs_curve_energy(&mut curve, 50, 0.1, 1.0);
    assert!(result.points_adjusted > 0);
    // The middle point should have moved closer to the line.
    assert!(
        curve.control_points[2].y < original_points[2].y,
        "energy minimization should reduce the kink"
    );
}

#[test]
fn fair_strength_clamped() {
    // Strength > 1.0 should be clamped to 1.0 and not panic.
    let mut curve = straight_5pt_nurbs();
    let result = fair_nurbs_curve(&mut curve, 5, 2.0, 1.0);
    assert!(result.max_displacement < 1e-12);

    // Strength < 0.0 should be clamped to 0.0 and produce no change.
    let mut curve2 = straight_5pt_nurbs();
    let result2 = fair_nurbs_curve(&mut curve2, 5, -1.0, 1.0);
    assert!(result2.max_displacement < 1e-12);
    assert_eq!(result2.points_adjusted, 0);
}
