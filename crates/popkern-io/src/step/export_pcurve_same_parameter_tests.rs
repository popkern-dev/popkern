//! PK-557: exported pcurves run over their 3D curve's written
//! parameter range.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::{Point3, Vec3};

use super::*;

/// The STEP CIRCLE the exporter writes for `center`/`normal`/`radius`,
/// evaluated at parameter `s` (angle from the placement's x axis).
fn written_circle(center: Point3, normal: Vec3, radius: f64, s: f64) -> Point3 {
    let x = super::super::export_nurbs::compute_perpendicular(&normal);
    let y = normal
        .normalized()
        .unwrap_or_else(|| panic!("unit normal"))
        .cross(&x);
    center + x * (radius * s.cos()) + y * (radius * s.sin())
}

#[test]
fn a_circle_edge_runs_over_the_written_circle_angle() {
    let (center, radius) = (Point3::new(1.0, -2.0, 0.5), 1.5);
    for normal in [Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.3, -0.8, 0.52)] {
        let normal = normal.normalized().unwrap_or_else(|| panic!("unit"));
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle: 0.0,
                end_angle: 6.0,
            },
            t_start: 4.1,
            t_end: 5.9,
        };
        let (s0, s1) = written_parameter_range(&curve).unwrap_or_else(|| panic!("mapped"));
        assert!((s1 - s0 - 1.8).abs() < 1e-12, "same sweep");
        for k in 0..=10 {
            let f = k as f64 / 10.0;
            let t = curve.t_start + f * (curve.t_end - curve.t_start);
            let s = s0 + f * (s1 - s0);
            let gap = curve
                .geometry
                .point_at(t)
                .distance_to(&written_circle(center, normal, radius, s));
            assert!(gap < 1e-12, "normal {normal:?} f {f}: {gap:e}");
        }
    }
}

#[test]
fn lines_and_splines_keep_their_parameter_and_other_curves_are_left_alone() {
    let line = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(2.0, 0.0, 0.0),
        },
        t_start: 0.25,
        t_end: 3.0,
    };
    assert_eq!(written_parameter_range(&line), Some((0.25, 3.0)));
    let reversed = CurveOnEdge {
        t_start: 3.0,
        t_end: 0.25,
        ..line.clone()
    };
    assert_eq!(written_parameter_range(&reversed), None);
    let ellipse = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center: Point3::new(0.0, 0.0, 0.0),
            u_basis: Vec3::new(2.0, 0.0, 0.0),
            v_basis: Vec3::new(0.0, 1.0, 0.0),
            start_angle: 0.0,
            end_angle: 1.0,
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    assert_eq!(
        written_parameter_range(&ellipse),
        None,
        "written as a rational B-spline"
    );
}

#[test]
fn remapping_knots_keeps_the_shape() {
    let mut nurbs = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(2.0, -1.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 4],
    };
    let before: Vec<Point3> = (0..=8).map(|k| nurbs.point_at(k as f64 / 8.0)).collect();
    remap_knots(&mut nurbs, (4.1, 5.9));
    assert_eq!((nurbs.knots[0], nurbs.knots[6]), (4.1, 5.9));
    for (k, p) in before.iter().enumerate() {
        let q = nurbs.point_at(4.1 + 1.8 * k as f64 / 8.0);
        assert!(p.distance_to(&q) < 1e-12);
    }
}
