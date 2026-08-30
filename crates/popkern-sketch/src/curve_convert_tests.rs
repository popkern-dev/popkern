//! Tests for 2D-to-3D curve conversion.

use super::*;
use crate::plane::SketchPlane;
use popkern_constraint::Sketch;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use std::f64::consts::FRAC_PI_2;

/// Helper: create a directed edge for a line between two 2D points.
fn line_directed_edge(
    sketch: &Sketch,
    start_pt: popkern_constraint::PointHandle,
    end_pt: popkern_constraint::PointHandle,
    line_handle: popkern_constraint::LineHandle,
    reversed: bool,
) -> DirectedEdge {
    let start = sketch.point_position(start_pt).unwrap();
    let end = sketch.point_position(end_pt).unwrap();
    if reversed {
        DirectedEdge {
            element: ElementHandle::Line(line_handle),
            start: end,
            end: start,
            reversed: true,
        }
    } else {
        DirectedEdge {
            element: ElementHandle::Line(line_handle),
            start,
            end,
            reversed: false,
        }
    }
}

#[test]
fn line_on_xy_plane_correct_3d() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(3.0, 4.0);
    let lh = sk.add_line(p0, p1);
    let plane = SketchPlane::xy();

    let de = line_directed_edge(&sk, p0, p1, lh, false);
    let curve = convert_to_3d(&de, &sk, &plane).unwrap();

    let start = curve.geometry.point_at(curve.t_start);
    let end = curve.geometry.point_at(curve.t_end);

    assert!((start.x - 0.0).abs() < 1e-10);
    assert!((start.y - 0.0).abs() < 1e-10);
    assert!((start.z - 0.0).abs() < 1e-10);

    assert!((end.x - 3.0).abs() < 1e-10);
    assert!((end.y - 4.0).abs() < 1e-10);
    assert!((end.z - 0.0).abs() < 1e-10);
}

#[test]
fn line_on_tilted_plane() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(1.0, 0.0);
    let p1 = sk.add_point(0.0, 1.0);
    let lh = sk.add_line(p0, p1);

    // Plane with u=X, v=Z, origin at (0,5,0)
    let plane = SketchPlane::new(
        Point3::new(0.0, 5.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    )
    .unwrap();

    let de = line_directed_edge(&sk, p0, p1, lh, false);
    let curve = convert_to_3d(&de, &sk, &plane).unwrap();

    let start = curve.geometry.point_at(curve.t_start);
    let end = curve.geometry.point_at(curve.t_end);

    // start: origin + 1*u + 0*v = (1, 5, 0)
    assert!((start.x - 1.0).abs() < 1e-10);
    assert!((start.y - 5.0).abs() < 1e-10);
    assert!((start.z - 0.0).abs() < 1e-10);

    // end: origin + 0*u + 1*v = (0, 5, 1)
    assert!((end.x - 0.0).abs() < 1e-10);
    assert!((end.y - 5.0).abs() < 1e-10);
    assert!((end.z - 1.0).abs() < 1e-10);
}

#[test]
fn arc_3d_correct_center_and_normal() {
    let mut sk = Sketch::new();
    let _arc = sk.add_arc(1.0, 2.0, 5.0, 0.0, FRAC_PI_2);
    let plane = SketchPlane::xy();

    let de = DirectedEdge {
        element: ElementHandle::Arc(sk.arc_handles().next().unwrap()),
        start: (6.0, 2.0), // center + (r, 0)
        end: (1.0, 7.0),   // center + (0, r)
        reversed: false,
    };

    let curve = convert_to_3d(&de, &sk, &plane).unwrap();
    match &curve.geometry {
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            assert!((center.x - 1.0).abs() < 1e-10);
            assert!((center.y - 2.0).abs() < 1e-10);
            assert!((center.z - 0.0).abs() < 1e-10);
            assert!((normal.z - 1.0).abs() < 1e-10, "normal should be +Z");
            assert!((*radius - 5.0).abs() < 1e-10);
        }
        _ => panic!("expected Arc geometry"),
    }
}

#[test]
fn reversed_line_swaps_t() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let lh = sk.add_line(p0, p1);
    let plane = SketchPlane::xy();

    let de = line_directed_edge(&sk, p0, p1, lh, true);
    let curve = convert_to_3d(&de, &sk, &plane).unwrap();

    assert!(
        curve.t_start > curve.t_end,
        "reversed line should have t_start > t_end"
    );
}

#[test]
fn reversed_arc_runs_from_its_directed_start_with_a_clockwise_sweep() {
    // The kernel's Arc evaluator measures angles in its own frame
    // (`circle_basis`), not the sketch plane's, so the contract is on
    // EVALUATED endpoints and the signed sweep, never on raw angles.
    let mut sk = Sketch::new();
    let _arc = sk.add_arc(0.0, 0.0, 1.0, 0.0, FRAC_PI_2);
    let plane = SketchPlane::xy();

    let de = DirectedEdge {
        element: ElementHandle::Arc(sk.arc_handles().next().unwrap()),
        start: (0.0, 1.0),
        end: (1.0, 0.0),
        reversed: true,
    };

    let curve = convert_to_3d(&de, &sk, &plane).unwrap();
    let CurveGeometry::Arc {
        start_angle,
        end_angle,
        ..
    } = &curve.geometry
    else {
        panic!("expected Arc geometry");
    };
    assert!(
        (end_angle - start_angle + FRAC_PI_2).abs() < 1e-10,
        "reversed quarter arc sweeps −π/2, got {}",
        end_angle - start_angle
    );
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    assert!(
        p0.distance_to(&Point3::new(0.0, 1.0, 0.0)) < 1e-10,
        "start {p0:?}"
    );
    assert!(
        p1.distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-10,
        "end {p1:?}"
    );
    // Midpoint on the quarter between them (never the ¾ complement).
    let m = curve.geometry.point_at((curve.t_start + curve.t_end) / 2.0);
    assert!(
        m.x > 0.0 && m.y > 0.0,
        "midpoint {m:?} must stay in the first quadrant"
    );
}

#[test]
fn line_pcurve_has_correct_2d_points() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(2.0, 3.0);
    let p1 = sk.add_point(5.0, 7.0);
    let lh = sk.add_line(p0, p1);

    let de = line_directed_edge(&sk, p0, p1, lh, false);
    let pcurve = convert_to_pcurve(&de, &sk, &SketchPlane::xy()).unwrap();
    let pcurve = pcurve.as_nurbs().expect("line pcurve is NURBS");

    assert_eq!(pcurve.degree, 1);
    assert_eq!(pcurve.control_points.len(), 2);
    assert!((pcurve.control_points[0].x - 2.0).abs() < 1e-10);
    assert!((pcurve.control_points[0].y - 3.0).abs() < 1e-10);
    assert!((pcurve.control_points[0].z - 0.0).abs() < 1e-10);
    assert!((pcurve.control_points[1].x - 5.0).abs() < 1e-10);
    assert!((pcurve.control_points[1].y - 7.0).abs() < 1e-10);
}

/// A periodic sketch spline is spelled for the kernel as the
/// equivalent unclamped open NURBS (extended knots, wrapped control
/// points); both evaluators agree along the whole period, and each
/// quarter-piece edge carries its own parameter sub-range — the lamina
/// validates under the full gate (before: every piece carried the
/// whole spline on `[0, 1]`, and evaluating it indexed past the
/// control points).
#[test]
fn periodic_spline_pieces_carry_their_sub_range_and_match_the_sketch_evaluator() {
    let mut sop = crate::SketchOnPlane::new(SketchPlane::xy());
    let through = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let (sh, _) = sop
        .sketch
        .add_closed_spline_through_points(&through)
        .unwrap();
    let result = sop.to_brep(crate::Tolerance::default(), None).unwrap();
    assert!(
        result.brep.validate_all().is_empty(),
        "{:?}",
        result.brep.validate_all()
    );
    assert_eq!(result.brep.edges.len(), 4);
    let mut covered = 0.0;
    for e in result.brep.edges.values() {
        let c = e.curve.as_ref().unwrap();
        let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        covered += hi - lo;
        for i in 0..=20 {
            let t = lo + (hi - lo) * i as f64 / 20.0;
            let p = c.geometry.point_at(t);
            let (x, y) = sop.sketch.spline_point_at(sh, t).unwrap();
            assert!(
                (p.x - x).abs() < 1e-9 && (p.y - y).abs() < 1e-9,
                "t={t}: kernel ({}, {}) vs sketch ({x}, {y})",
                p.x,
                p.y
            );
        }
    }
    // Four quarters tile the period exactly once.
    assert!(
        (covered - 4.0).abs() < 1e-6,
        "pieces cover {covered} of period 4"
    );
}
