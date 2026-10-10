//! Plane × circle closed form (PK-758).

use super::{PlaneFrame, intersect_plane_circle};
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn points(r: &[SurfaceCurveIntersection]) -> Vec<(f64, Point3)> {
    r.iter()
        .filter_map(|i| match i {
            SurfaceCurveIntersection::Point { t, point, .. } => Some((*t, *point)),
            SurfaceCurveIntersection::Overlap { .. } => None,
        })
        .collect()
}

#[test]
fn an_arc_on_a_window_a_turn_below_its_stored_angles_meets_a_plane_twice() {
    // The ring's notch rim: radius 49 about (-817, y, 162.5), axis y,
    // crossing the plane x = -770 twice; the edge window is one turn
    // below the arc's stored angles.
    let arc = CurveGeometry::Arc {
        center: Point3::new(-817.0, -25.0, 162.5),
        normal: Vec3::new(0.0, 1.0, 0.0),
        radius: 49.0,
        start_angle: 0.0,
        end_angle: std::f64::consts::TAU,
    };
    let tol = Tolerance::default();
    // The window covering the +x extreme, a turn below.
    let probe = |t: f64| arc.point_at(t);
    let t_max_x = (0..3600)
        .map(|k| k as f64 * std::f64::consts::TAU / 3600.0)
        .max_by(|a, b| probe(*a).x.total_cmp(&probe(*b).x))
        .expect("samples");
    let window = (
        t_max_x - 0.5 - std::f64::consts::TAU,
        t_max_x + 0.5 - std::f64::consts::TAU,
    );
    let (o, n, u, v) = (
        Point3::new(-770.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    let hits = intersect_plane_circle(
        PlaneFrame {
            origin: &o,
            normal: &n,
            u_axis: &u,
            v_axis: &v,
        },
        &arc,
        window,
        ((-1e3, 1e3), (-1e3, 1e3)),
        &tol,
    )
    .expect("circular");
    let pts = points(&hits);
    assert_eq!(pts.len(), 2, "{pts:?}");
    for (t, p) in pts {
        assert!(t >= window.0 && t <= window.1, "t {t} outside {window:?}");
        assert!((p.x + 770.0).abs() < 1e-9, "{p:?}");
        let dz = (p.z - 162.5).abs();
        assert!(
            (dz - (49.0f64 * 49.0 - 47.0 * 47.0).sqrt()).abs() < 1e-9,
            "{p:?}"
        );
    }
}

#[test]
fn a_tangent_plane_gives_one_root_and_an_in_plane_circle_an_overlap() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let tol = Tolerance::default();
    let (u, v) = (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
    let o = Point3::new(2.0, 0.0, 0.0);
    let n = Vec3::new(1.0, 0.0, 0.0);
    let frame = PlaneFrame {
        origin: &o,
        normal: &n,
        u_axis: &u,
        v_axis: &v,
    };
    let tangent = intersect_plane_circle(
        frame,
        &circle,
        (0.0, std::f64::consts::TAU),
        ((-9.0, 9.0), (-9.0, 9.0)),
        &tol,
    )
    .expect("circular");
    assert_eq!(points(&tangent).len(), 1, "{tangent:?}");
    let (o2, n2, u2, v2) = (
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    );
    let in_plane = intersect_plane_circle(
        PlaneFrame {
            origin: &o2,
            normal: &n2,
            u_axis: &u2,
            v_axis: &v2,
        },
        &circle,
        (0.0, 1.0),
        ((-9.0, 9.0), (-9.0, 9.0)),
        &tol,
    )
    .expect("circular");
    assert!(
        matches!(
            in_plane.as_slice(),
            [SurfaceCurveIntersection::Overlap { .. }]
        ),
        "{in_plane:?}"
    );
}
