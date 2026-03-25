use super::*;

#[test]
fn interpolate_constant() {
    let angle = interpolate_angle(0.5, 0.0, 1.0, &[], 0.1);
    assert!((angle - 0.1).abs() < 1e-12);
}

#[test]
fn interpolate_single_pair() {
    let angle = interpolate_angle(0.5, 0.0, 1.0, &[(0.5, 0.2)], 0.1);
    assert!((angle - 0.2).abs() < 1e-12);
}

#[test]
fn interpolate_two_pairs_midpoint() {
    let angle = interpolate_angle(0.5, 0.0, 1.0, &[(0.0, 0.1), (1.0, 0.3)], 0.0);
    assert!((angle - 0.2).abs() < 1e-10, "expected 0.2, got {}", angle);
}

#[test]
fn interpolate_two_pairs_endpoints() {
    let a0 = interpolate_angle(0.0, 0.0, 1.0, &[(0.0, 0.1), (1.0, 0.3)], 0.0);
    let a1 = interpolate_angle(1.0, 0.0, 1.0, &[(0.0, 0.1), (1.0, 0.3)], 0.0);
    assert!((a0 - 0.1).abs() < 1e-10);
    assert!((a1 - 0.3).abs() < 1e-10);
}

#[test]
fn variable_angle_verify_different_at_endpoints() {
    let pairs = vec![(0.0, 0.05), (1.0, 0.15)];
    let a_start = interpolate_angle(0.0, 0.0, 1.0, &pairs, 0.0);
    let a_end = interpolate_angle(1.0, 0.0, 1.0, &pairs, 0.0);
    assert!(
        (a_start - a_end).abs() > 0.05,
        "angles at endpoints should differ"
    );
}

#[test]
fn rodrigues_identity() {
    let v = Vec3::new(0.0, 1.0, 0.0);
    let k = Vec3::new(0.0, 0.0, 1.0);
    let rotated = rodrigues_rotate(&v, &k, 0.0);
    assert!((rotated.x - v.x).abs() < 1e-12);
    assert!((rotated.y - v.y).abs() < 1e-12);
}

#[test]
fn sample_cross_section_produces_correct_count() {
    let center = Point3::new(0.0, 0.0, 0.0);
    let cross_dir = Vec3::new(1.0, 0.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: center,
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let points = sample_cross_section(&plane, &center, &cross_dir, 5);
    assert_eq!(points.len(), 5);
}

#[test]
fn sample_cross_section_single_point_at_center() {
    let center = Point3::new(1.0, 2.0, 3.0);
    let cross_dir = Vec3::new(0.0, 1.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: center,
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let points = sample_cross_section(&plane, &center, &cross_dir, 1);
    assert_eq!(points.len(), 1);
    // Single sample lands at center (t=0.0).
    assert!((points[0].x - center.x).abs() < 1e-10);
    assert!((points[0].y - center.y).abs() < 1e-10);
}

#[test]
fn rodrigues_rotate_ninety_degrees() {
    // Rotate Y axis 90 degrees about Z: should give -X.
    let v = Vec3::new(0.0, 1.0, 0.0);
    let k = Vec3::new(0.0, 0.0, 1.0);
    let rotated = rodrigues_rotate(&v, &k, std::f64::consts::FRAC_PI_2);
    assert!((rotated.x + 1.0).abs() < 1e-10, "x={}", rotated.x);
    assert!(rotated.y.abs() < 1e-10, "y={}", rotated.y);
    assert!(rotated.z.abs() < 1e-10, "z={}", rotated.z);
}

#[test]
fn rodrigues_rotate_preserves_length() {
    let v = Vec3::new(1.0, 2.0, 3.0);
    let k = Vec3::new(0.0, 1.0, 0.0).normalized().unwrap();
    let angle = 0.7;
    let rotated = rodrigues_rotate(&v, &k, angle);
    let original_len = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
    let rotated_len =
        (rotated.x * rotated.x + rotated.y * rotated.y + rotated.z * rotated.z).sqrt();
    assert!((original_len - rotated_len).abs() < 1e-10);
}

#[test]
fn compute_cross_direction_plane_surface() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let center = Point3::new(0.0, 0.0, 0.0);
    let tangent = Vec3::new(1.0, 0.0, 0.0);
    let cross = compute_cross_direction(&plane, &center, &tangent);
    // Cross direction should be perpendicular to tangent.
    assert!(
        cross.dot(&tangent).abs() < 1e-10,
        "dot={}",
        cross.dot(&tangent)
    );
    // Should be a unit vector.
    let len = (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z).sqrt();
    assert!((len - 1.0).abs() < 1e-10);
}

#[test]
fn compute_cross_direction_nurbs_surface() {
    use crate::surface::skinning;
    // Build a tiny NURBS surface and compute cross direction.
    let rows = vec![
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
    ];
    let nurbs = skinning::skin_points(&rows).unwrap();
    let surf = SurfaceGeometry::Nurbs(nurbs);
    let center = Point3::new(0.5, 0.5, 0.0);
    let tangent = Vec3::new(1.0, 0.0, 0.0);
    let cross = compute_cross_direction(&surf, &center, &tangent);
    // Must be a unit vector.
    let len = (cross.x * cross.x + cross.y * cross.y + cross.z * cross.z).sqrt();
    assert!(
        (len - 1.0).abs() < 1e-10,
        "cross direction not unit: len={}",
        len
    );
}

#[test]
fn interpolate_angle_below_first_pair() {
    // u < pairs[0].0: should clamp to first pair's angle.
    let angle = interpolate_angle(-1.0, 0.0, 1.0, &[(0.5, 0.3), (1.0, 0.6)], 0.0);
    assert!((angle - 0.3).abs() < 1e-12, "got {}", angle);
}

#[test]
fn interpolate_angle_above_last_pair() {
    // u > pairs.last().0: should clamp to last pair's angle.
    let angle = interpolate_angle(2.0, 0.0, 1.0, &[(0.0, 0.1), (0.5, 0.2)], 0.0);
    assert!((angle - 0.2).abs() < 1e-12, "got {}", angle);
}

#[test]
fn interpolate_angle_degenerate_range() {
    // t0 == t1: range is zero, u defaults to 0.5.
    let angle = interpolate_angle(5.0, 5.0, 5.0, &[(0.0, 0.1), (1.0, 0.9)], 0.0);
    // u=0.5 → should interpolate to 0.5 between 0.1 and 0.9.
    assert!((angle - 0.5).abs() < 1e-10, "got {}", angle);
}

#[test]
fn draft_general_constant_angle_produces_nurbs() {
    use crate::neutral::NeutralLine;
    use crate::surface::skinning;
    use crate::types::{AngleVariation, DraftOptions};
    use popkern_core::geometry::NurbsCurve;

    let rows = vec![
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
    ];
    let nurbs_surf = skinning::skin_points(&rows).unwrap();
    let surface = SurfaceGeometry::Nurbs(nurbs_surf);

    let neutral = NeutralLine {
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        },
        direction: Vec3::new(1.0, 0.0, 0.0),
        origin: Point3::new(0.5, 0.5, 0.0),
    };

    let pull = Vec3::new(0.0, 0.0, 1.0);
    let opts = DraftOptions {
        nurbs_samples: 3,
        ..DraftOptions::default()
    };
    let tol = popkern_core::tolerance::Tolerance::default();

    let result = draft_general(
        &surface,
        &neutral,
        0.1,
        &AngleVariation::Constant,
        &pull,
        &opts,
        &tol,
    );
    assert!(result.is_ok(), "draft_general failed: {:?}", result.err());
    assert!(matches!(result.unwrap(), SurfaceGeometry::Nurbs(_)));
}

#[test]
fn draft_general_zero_angle_still_produces_nurbs() {
    use crate::neutral::NeutralLine;
    use crate::surface::skinning;
    use crate::types::{AngleVariation, DraftOptions};
    use popkern_core::geometry::NurbsCurve;

    let rows = vec![
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
    ];
    let nurbs_surf = skinning::skin_points(&rows).unwrap();
    let surface = SurfaceGeometry::Nurbs(nurbs_surf);

    let neutral = NeutralLine {
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        },
        direction: Vec3::new(1.0, 0.0, 0.0),
        origin: Point3::new(0.5, 0.5, 0.0),
    };

    let pull = Vec3::new(0.0, 0.0, 1.0);
    let opts = DraftOptions {
        nurbs_samples: 2,
        ..DraftOptions::default()
    };
    let tol = popkern_core::tolerance::Tolerance::default();

    let result = draft_general(
        &surface,
        &neutral,
        0.0,
        &AngleVariation::Constant,
        &pull,
        &opts,
        &tol,
    );
    assert!(
        result.is_ok(),
        "zero angle draft failed: {:?}",
        result.err()
    );
}

#[test]
fn draft_general_variable_angle() {
    use crate::neutral::NeutralLine;
    use crate::surface::skinning;
    use crate::types::{AngleVariation, DraftOptions};
    use popkern_core::geometry::NurbsCurve;

    let rows = vec![
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
    ];
    let nurbs_surf = skinning::skin_points(&rows).unwrap();
    let surface = SurfaceGeometry::Nurbs(nurbs_surf);

    let neutral = NeutralLine {
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        },
        direction: Vec3::new(1.0, 0.0, 0.0),
        origin: Point3::new(0.5, 0.5, 0.0),
    };

    let pull = Vec3::new(0.0, 0.0, 1.0);
    let variation = AngleVariation::Variable(vec![(0.0, 0.05), (1.0, 0.15)]);
    let opts = DraftOptions {
        nurbs_samples: 4,
        ..DraftOptions::default()
    };
    let tol = popkern_core::tolerance::Tolerance::default();

    let result = draft_general(&surface, &neutral, 0.1, &variation, &pull, &opts, &tol);
    assert!(
        result.is_ok(),
        "variable angle draft failed: {:?}",
        result.err()
    );
}

#[test]
fn draft_general_nurbs_samples_minimum_two() {
    use crate::neutral::NeutralLine;
    use crate::surface::skinning;
    use crate::types::{AngleVariation, DraftOptions};
    use popkern_core::geometry::NurbsCurve;

    // nurbs_samples=1 should be clamped to 2 internally.
    let rows = vec![
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
    ];
    let nurbs_surf = skinning::skin_points(&rows).unwrap();
    let surface = SurfaceGeometry::Nurbs(nurbs_surf);

    let neutral = NeutralLine {
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
            weights: vec![1.0, 1.0],
        },
        direction: Vec3::new(1.0, 0.0, 0.0),
        origin: Point3::new(0.5, 0.5, 0.0),
    };

    let pull = Vec3::new(0.0, 0.0, 1.0);
    let opts = DraftOptions {
        nurbs_samples: 1,
        ..DraftOptions::default()
    };
    let tol = popkern_core::tolerance::Tolerance::default();

    let result = draft_general(
        &surface,
        &neutral,
        0.05,
        &AngleVariation::Constant,
        &pull,
        &opts,
        &tol,
    );
    assert!(
        result.is_ok(),
        "n_samples=1 clamp failed: {:?}",
        result.err()
    );
}
