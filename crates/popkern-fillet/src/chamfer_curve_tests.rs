//! Tests for chamfer curve and surface computation.
//! Included only in test mode via `#[cfg(test)] #[path = "chamfer_curve_tests.rs"] mod tests`.

use super::*;
use crate::analysis::analyze_edges;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn make_cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).unwrap()
}

fn first_edge_analysis(brep: &BRep) -> EdgeAnalysis {
    let edge = brep.edges.keys().next().unwrap();
    let mut analyses = analyze_edges(brep, &[edge], false).unwrap();
    analyses.remove(0)
}

#[test]
fn chamfer_equal_cube_edge() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let result = compute_chamfer(&analysis, &cube, &ChamferKind::Equal(0.1), &tolerance);
    assert!(result.is_ok(), "chamfer failed: {:?}", result.err());

    let chamfer = result.unwrap();

    // A straight-edge plane-plane chamfer bevel is FLAT: the surface
    // must come out as an analytic Plane (canonical recognition),
    // with the true 45° diagonal normal.
    let SurfaceGeometry::Plane { normal, .. } = &chamfer.surface else {
        panic!("expected analytic Plane bevel, got {:?}", chamfer.surface);
    };
    let n = normal.normalized().expect("bevel normal");
    assert!(
        (n.x.abs() - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-9
            || (n.y.abs() - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-9
            || (n.z.abs() - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-9,
        "bevel normal must be a 45° diagonal, got ({}, {}, {})",
        n.x,
        n.y,
        n.z
    );

    // Trim curves should have control points.
    assert!(
        !chamfer.trim_curve_1.control_points.is_empty(),
        "trim_curve_1 should have points"
    );
    assert!(
        !chamfer.trim_curve_2.control_points.is_empty(),
        "trim_curve_2 should have points"
    );

    // Pcurves should use (u, v, 0) convention.
    for pt in &chamfer.pcurve_1.control_points {
        assert!(pt.z.abs() < 1e-10, "pcurve z should be 0, got {}", pt.z);
    }
}

#[test]
fn chamfer_two_distance() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let kind = ChamferKind::TwoDistance {
        face_1: analysis.face_1,
        dist_1: 0.1,
        dist_2: 0.2,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "two-distance chamfer failed: {:?}",
        result.err()
    );

    let chamfer = result.unwrap();

    // Both trim curves should exist and have points.
    assert!(!chamfer.trim_curve_1.control_points.is_empty());
    assert!(!chamfer.trim_curve_2.control_points.is_empty());
}

#[test]
fn chamfer_distance_angle() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let kind = ChamferKind::DistanceAngle {
        face: analysis.face_1,
        distance: 0.1,
        angle: std::f64::consts::FRAC_PI_4, // 45 degrees -> tan = 1.0
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "distance-angle chamfer failed: {:?}",
        result.err()
    );
}

#[test]
fn chamfer_invalid_reference_face() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    // A default handle will never match any real face in the analysis.
    let bogus_face = popkern_core::topology::FaceHandle::default();

    let kind = ChamferKind::TwoDistance {
        face_1: bogus_face,
        dist_1: 0.1,
        dist_2: 0.2,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(result.is_err(), "should fail with invalid reference face");
}

#[test]
fn build_ruled_surface_basic() {
    let pts1 = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ];
    let pts2 = vec![
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
    ];

    let surface = build_ruled_surface(&pts1, &pts2);
    assert!(surface.is_ok());

    // Coplanar rows are recognized as an analytic Plane; only
    // genuinely curved bevels stay NURBS.
    match surface {
        Ok(SurfaceGeometry::Plane { .. }) => {}
        Ok(other) => panic!("coplanar ruled rows must recognize as Plane, got {other:?}"),
        Err(e) => panic!("ruled surface failed: {e}"),
    }
}

#[test]
fn build_ruled_surface_curved_stays_nurbs() {
    // Non-coplanar rows (a genuinely twisted/curved bevel) must NOT
    // be flattened by the planar recognition.
    let pts1 = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.2),
        Point3::new(2.0, 0.0, 0.0),
    ];
    let pts2 = vec![
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
    ];
    let surface = build_ruled_surface(&pts1, &pts2).expect("ruled");
    if let SurfaceGeometry::Nurbs(nurbs) = surface {
        assert_eq!(nurbs.num_u, 3);
        assert_eq!(nurbs.num_v, 2);
        assert_eq!(nurbs.degree_v, 1);
    } else {
        panic!("curved rows must stay NURBS");
    }
}

#[test]
fn build_ruled_surface_mismatched_lengths() {
    let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    let pts2 = vec![Point3::new(0.0, 1.0, 0.0)];
    assert!(build_ruled_surface(&pts1, &pts2).is_err());
}

// -----------------------------------------------------------------------
// Gap 3: Chamfer TwoDistance and DistanceAngle edge cases
// -----------------------------------------------------------------------

#[test]
fn two_distance_equal_distances_is_symmetric() {
    // TwoDistance with equal d1==d2 should behave like Equal.
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let equal_result =
        compute_chamfer(&analysis, &cube, &ChamferKind::Equal(0.15), &tolerance).unwrap();

    let two_dist_kind = ChamferKind::TwoDistance {
        face_1: analysis.face_1,
        dist_1: 0.15,
        dist_2: 0.15,
    };
    let two_dist_result = compute_chamfer(&analysis, &cube, &two_dist_kind, &tolerance).unwrap();

    // Both flat bevels come out as analytic Planes.
    assert!(matches!(
        equal_result.surface,
        SurfaceGeometry::Plane { .. }
    ));
    assert!(matches!(
        two_dist_result.surface,
        SurfaceGeometry::Plane { .. }
    ));

    // With equal distances, trim curves should be at similar offsets.
    // Compare the midpoints of trim_curve_1 for both.
    let (t0_e, t1_e) = equal_result.trim_curve_1.parameter_range();
    let mid_equal = equal_result.trim_curve_1.point_at((t0_e + t1_e) * 0.5);

    let (t0_t, t1_t) = two_dist_result.trim_curve_1.parameter_range();
    let mid_two = two_dist_result.trim_curve_1.point_at((t0_t + t1_t) * 0.5);

    let dist = mid_equal.distance_to(&mid_two);
    assert!(
        dist < 0.1,
        "symmetric TwoDistance should produce similar trim to Equal, got dist={dist}",
    );
}

#[test]
fn two_distance_very_different_distances() {
    // TwoDistance with large asymmetry: 0.05 vs 0.4.
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let kind = ChamferKind::TwoDistance {
        face_1: analysis.face_1,
        dist_1: 0.05,
        dist_2: 0.4,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "asymmetric two-distance chamfer failed: {:?}",
        result.err(),
    );

    let chamfer = result.unwrap();
    assert!(!chamfer.trim_curve_1.control_points.is_empty());
    assert!(!chamfer.trim_curve_2.control_points.is_empty());
}

#[test]
fn two_distance_reversed_face_ref() {
    // Pass face_2 as the reference face to verify the distance swap logic.
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let kind = ChamferKind::TwoDistance {
        face_1: analysis.face_2, // reversed
        dist_1: 0.1,
        dist_2: 0.2,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "reversed-face two-distance chamfer failed: {:?}",
        result.err(),
    );
}

#[test]
fn distance_angle_small_angle() {
    // Small angle (10 degrees): other_dist = distance * tan(10deg) ~ 0.0176.
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let small_angle = 10.0_f64.to_radians();
    let kind = ChamferKind::DistanceAngle {
        face: analysis.face_1,
        distance: 0.1,
        angle: small_angle,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "small-angle chamfer failed: {:?}",
        result.err(),
    );
}

#[test]
fn distance_angle_large_angle() {
    // Large angle (80 degrees): other_dist = distance * tan(80deg) ~ 0.567.
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let large_angle = 80.0_f64.to_radians();
    let kind = ChamferKind::DistanceAngle {
        face: analysis.face_1,
        distance: 0.1,
        angle: large_angle,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "large-angle chamfer failed: {:?}",
        result.err(),
    );
}

#[test]
fn distance_angle_reversed_face_ref() {
    // Using face_2 as reference should swap which face gets the direct distance.
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let kind = ChamferKind::DistanceAngle {
        face: analysis.face_2,
        distance: 0.15,
        angle: std::f64::consts::FRAC_PI_4,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_ok(),
        "reversed-face distance-angle chamfer failed: {:?}",
        result.err(),
    );
}

#[test]
fn distance_angle_invalid_face_errors() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let tolerance = Tolerance::default();

    let bogus_face = popkern_core::topology::FaceHandle::default();
    let kind = ChamferKind::DistanceAngle {
        face: bogus_face,
        distance: 0.1,
        angle: std::f64::consts::FRAC_PI_4,
    };

    let result = compute_chamfer(&analysis, &cube, &kind, &tolerance);
    assert!(
        result.is_err(),
        "invalid face for DistanceAngle should error",
    );
}

#[test]
fn compute_offset_distances_equal() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let (d1, d2) = compute_offset_distances(&analysis, &ChamferKind::Equal(0.25)).unwrap();
    assert!((d1 - 0.25).abs() < f64::EPSILON);
    assert!((d2 - 0.25).abs() < f64::EPSILON);
}

#[test]
fn compute_offset_distances_two_distance_face1() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let kind = ChamferKind::TwoDistance {
        face_1: analysis.face_1,
        dist_1: 0.1,
        dist_2: 0.3,
    };
    let (d1, d2) = compute_offset_distances(&analysis, &kind).unwrap();
    assert!((d1 - 0.1).abs() < f64::EPSILON);
    assert!((d2 - 0.3).abs() < f64::EPSILON);
}

#[test]
fn compute_offset_distances_two_distance_face2() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let kind = ChamferKind::TwoDistance {
        face_1: analysis.face_2, // reversed
        dist_1: 0.1,
        dist_2: 0.3,
    };
    // When ref face is face_2, distances are swapped.
    let (d1, d2) = compute_offset_distances(&analysis, &kind).unwrap();
    assert!((d1 - 0.3).abs() < f64::EPSILON);
    assert!((d2 - 0.1).abs() < f64::EPSILON);
}

#[test]
fn compute_offset_distances_distance_angle() {
    let cube = make_cube();
    let analysis = first_edge_analysis(&cube);
    let kind = ChamferKind::DistanceAngle {
        face: analysis.face_1,
        distance: 0.2,
        angle: std::f64::consts::FRAC_PI_4, // tan(45deg) = 1.0
    };
    let (d1, d2) = compute_offset_distances(&analysis, &kind).unwrap();
    assert!((d1 - 0.2).abs() < f64::EPSILON);
    assert!((d2 - 0.2).abs() < 1e-10, "tan(PI/4) * 0.2 = 0.2, got {d2}");
}
