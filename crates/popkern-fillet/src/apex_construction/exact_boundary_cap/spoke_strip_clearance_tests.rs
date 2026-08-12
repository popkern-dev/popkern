//! Unit tests for the construct-time spoke/strip footprint clearance
//! (PK-261).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::super::corner::GeneratedInnerSpoke;
use super::super::sector::ExactBoundaryOuterSide;
use super::footprint::{arc_windows, unwrap_periodic};
use super::{StripFootprint, clear_one_spoke, cubic_bezier_form, worst_violation};
use crate::types::ApexConstructionDiagnostic;

const BAND: f64 = 4.0e-6;

fn unit_z_cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    }
}

fn test_strip() -> StripFootprint {
    StripFootprint {
        face: popkern_core::topology::FaceHandle::default(),
        surface: unit_z_cylinder(),
        u_periodic: true,
        u_range: (0.0, std::f64::consts::PI),
        v_range: (0.0, 2.0),
        // Corner end arc at v = 0 (axis-perpendicular), interior above.
        arc_profile: vec![(0.0, 0.0), (std::f64::consts::PI, 0.0)],
        far_side_positive: true,
    }
}

fn cubic(points: [Point3; 4]) -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: points.to_vec(),
        weights: vec![1.0; 4],
    }
}

fn diving_spoke() -> GeneratedInnerSpoke {
    // Anchor and pole outside the cylinder; the interior handles pull
    // the curve radially inside within the footprint window.
    let anchor = Point3::new(1.3, 0.05, 1.0);
    let pole = Point3::new(1.2 * (1.0f64).cos(), 1.2 * (1.0f64).sin(), 1.0);
    GeneratedInnerSpoke {
        name: "test_spoke".to_string(),
        source_side: ExactBoundaryOuterSide::FilletEnd,
        anchor,
        pole,
        curve: cubic([
            anchor,
            Point3::new(0.4, 0.2, 1.0),
            Point3::new(0.3, 0.5, 1.0),
            pole,
        ]),
    }
}

#[test]
fn violating_spoke_clears_the_strip_footprint() {
    let strips = [test_strip()];
    let sides = [Some(1.0)];
    let mut spoke = diving_spoke();
    let anchor = spoke.curve.control_points[0];
    let pole = spoke.curve.control_points[3];
    let start_dir = (spoke.curve.control_points[1] - anchor)
        .normalized()
        .unwrap();
    let before = worst_violation(&spoke.curve, &strips, &sides, BAND)
        .expect("the synthetic spoke must violate the footprint");
    assert!(
        before.depth > 0.3,
        "synthetic dip should be deep; got {}",
        before.depth
    );

    let diagnostic = clear_one_spoke(&mut spoke, &strips, &sides, BAND, "sector", "FilletEnd")
        .expect("a violating spoke must produce a diagnostic");
    let ApexConstructionDiagnostic::ExactBoundaryCap(record) = diagnostic else {
        panic!("clearance emits exact-boundary diagnostics");
    };
    assert_eq!(record.code, "spoke_footprint_cleared");

    assert!(
        worst_violation(&spoke.curve, &strips, &sides, BAND).is_none(),
        "the corrected spoke must measure clear"
    );
    // Locked quantities: endpoints exactly, anchor tangent direction.
    assert_eq!(spoke.curve.control_points[0], anchor);
    assert_eq!(spoke.curve.control_points[3], pole);
    let corrected_dir = (spoke.curve.control_points[1] - anchor)
        .normalized()
        .unwrap();
    assert!(
        corrected_dir.cross(&start_dir).length() < 1.0e-9,
        "the anchor-side tangent direction is locked by seam-end G1"
    );
}

#[test]
fn clear_spoke_is_left_untouched() {
    let strips = [test_strip()];
    let sides = [Some(1.0)];
    let anchor = Point3::new(1.5, 0.0, 1.0);
    let pole = Point3::new(1.5 * (1.0f64).cos(), 1.5 * (1.0f64).sin(), 1.0);
    let mut spoke = GeneratedInnerSpoke {
        name: "clear_spoke".to_string(),
        source_side: ExactBoundaryOuterSide::FilletEnd,
        anchor,
        pole,
        curve: cubic([
            anchor,
            Point3::new(1.5, 0.4, 1.0),
            Point3::new(1.3, 0.8, 1.0),
            pole,
        ]),
    };
    let before = spoke.curve.clone();
    assert!(
        clear_one_spoke(&mut spoke, &strips, &sides, BAND, "sector", "FilletEnd").is_none(),
        "a clear spoke produces no diagnostic"
    );
    assert_eq!(spoke.curve, before, "a clear spoke is not modified");
}

#[test]
fn twin_spoke_records_stay_bit_identical() {
    let strips = [test_strip()];
    let sides = [Some(1.0)];
    let mut first = diving_spoke();
    let mut second = diving_spoke();
    clear_one_spoke(&mut first, &strips, &sides, BAND, "sector_a", "FaceEdge");
    clear_one_spoke(&mut second, &strips, &sides, BAND, "sector_b", "FilletEnd");
    assert_eq!(
        first.curve, second.curve,
        "shared-spoke twins must correct identically from identical inputs"
    );
}

#[test]
fn negative_side_defends_the_inside_of_the_primitive() {
    // With the rim on the INSIDE (a concave configuration), a spoke
    // that stays inside is legal and the same spoke violates when the
    // rim side flips to positive.
    let strips = [test_strip()];
    let inside_curve = cubic([
        Point3::new(0.5, 0.05, 1.0),
        Point3::new(0.4, 0.2, 1.0),
        Point3::new(0.35, 0.35, 1.0),
        Point3::new(0.3, 0.5, 1.0),
    ]);
    assert!(
        worst_violation(&inside_curve, &strips, &[Some(-1.0)], BAND).is_none(),
        "an inside spoke is legal when the rim side is negative"
    );
    assert!(
        worst_violation(&inside_curve, &strips, &[Some(1.0)], BAND).is_some(),
        "the same spoke violates when the rim side is positive"
    );
}

#[test]
fn straight_spokes_elevate_exactly_to_cubic() {
    let a = Point3::new(1.0, 2.0, 3.0);
    let b = Point3::new(4.0, 5.0, 6.0);
    let line = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0; 2],
    };
    let cubic = cubic_bezier_form(&line).expect("lines elevate");
    assert_eq!(cubic.degree, 3);
    for k in 0..=8 {
        let t = k as f64 / 8.0;
        assert!(
            cubic.point_at(t).distance_to(&line.point_at(t)) < 1.0e-12,
            "elevation must be geometrically exact"
        );
    }
}

#[test]
fn arc_window_handles_the_periodic_seam() {
    // A quarter arc straddling u=0 on the unit cylinder: from -45° to
    // +45°, sampled as 3D points. The window must span ~[-π/4, π/4]
    // (unwrapped), not fold to [0, 2π).
    let surface = unit_z_cylinder();
    let quarter = std::f64::consts::FRAC_PI_4;
    let arc_points: Vec<Point3> = (0..=8)
        .map(|k| {
            let a = -quarter + 2.0 * quarter * k as f64 / 8.0;
            Point3::new(a.cos(), a.sin(), 1.0)
        })
        .collect();
    // Interpolating polyline as a degree-1 NURBS is enough for window
    // sampling.
    let n = arc_points.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    let arc = NurbsCurve {
        degree: 1,
        knots,
        control_points: arc_points,
        weights: vec![1.0; n],
    };
    let ((lo, hi), _) = arc_windows(&arc, &surface, true).expect("arc projects");
    assert!(
        hi - lo < 2.0 * quarter * 1.2,
        "the window must stay a quarter turn, not fold around the seam: [{lo},{hi}]"
    );
}

#[test]
fn unwrap_periodic_selects_the_nearest_branch() {
    let tau = std::f64::consts::TAU;
    assert!((unwrap_periodic(0.1 + tau, 0.0) - 0.1).abs() < 1.0e-12);
    assert!((unwrap_periodic(-0.1, tau) - (tau - 0.1)).abs() < 1.0e-12);
    assert!((unwrap_periodic(1.0, 1.5) - 1.0).abs() < 1.0e-12);
}
