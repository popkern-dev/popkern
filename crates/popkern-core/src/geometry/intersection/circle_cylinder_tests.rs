//! Circle × cylinder closed form with a tolerance-band tangency (PK-887).

use super::{CylinderFrame, intersect_circle_cylinder};
use crate::geometry::curve::CurveGeometry;
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// A radius-49 circle about (-817, y, 162.5) in the plane y = -17
/// (axis y) against a radius-33 cylinder along x whose top generator
/// lies `proud` past that plane (the web's ring arc against the boss).
fn hits(proud: f64) -> Vec<Point3> {
    let circle = CurveGeometry::Circle {
        center: Point3::new(-817.0, -17.0, 162.5),
        normal: Vec3::new(0.0, 1.0, 0.0),
        radius: 49.0,
    };
    let origin = Point3::new(0.0, -50.0 + proud, 162.5);
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let found = intersect_circle_cylinder(
        CylinderFrame {
            origin: &origin,
            axis: &axis,
            radius: 33.0,
        },
        &circle,
        (0.0, std::f64::consts::TAU),
        ((0.0, std::f64::consts::TAU), (-1e3, 1e3)),
        &Tolerance::default(),
    )
    .expect("circular");
    found
        .iter()
        .filter_map(|h| match h {
            SurfaceCurveIntersection::Point { point, .. } => Some(*point),
            SurfaceCurveIntersection::Overlap { .. } => None,
        })
        .collect()
}

/// The circle meets the (unbounded) cylinder where it is nearest the
/// axis: its two points at z = 162.5, x = −768 and x = −866.
#[test]
fn a_circle_past_a_cylinder_within_tolerance_touches_it_once_at_each_extreme() {
    for proud in [0.0, 1e-7, 5e-7, 9e-7] {
        let mut pts = hits(proud);
        pts.sort_by(|p, q| q.x.total_cmp(&p.x));
        assert_eq!(pts.len(), 2, "proud {proud:e}: {pts:?}");
        for (p, x) in pts.iter().zip([-768.0, -866.0]) {
            let extreme = Point3::new(x, -17.0, 162.5);
            assert!((*p - extreme).length() < 1e-9, "proud {proud:e}: {pts:?}");
        }
    }
}

#[test]
fn a_circle_a_whole_tolerance_past_a_cylinder_crosses_it_twice_at_each_extreme() {
    let tol = Tolerance::default().linear();
    for proud in [tol, 1e-4, 1.0] {
        let pts = hits(proud);
        assert_eq!(pts.len(), 4, "proud {proud:e}: {pts:?}");
        for p in &pts {
            let r = ((p.y - (-50.0 + proud)).powi(2) + (p.z - 162.5).powi(2)).sqrt();
            assert!(
                (r - 33.0).abs() < 1e-9,
                "proud {proud:e}: {p:?} off the wall"
            );
        }
    }
}

#[test]
fn a_circle_short_of_a_cylinder_misses_it() {
    assert!(hits(-1e-7).is_empty());
    assert!(hits(-1.0).is_empty());
}

#[test]
fn a_circle_on_the_wall_defers_to_the_general_solver() {
    let rim = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 3.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let origin = Point3::new(0.0, 0.0, 0.0);
    let axis = Vec3::new(0.0, 0.0, 1.0);
    let found = intersect_circle_cylinder(
        CylinderFrame {
            origin: &origin,
            axis: &axis,
            radius: 2.0,
        },
        &rim,
        (0.0, std::f64::consts::TAU),
        ((0.0, std::f64::consts::TAU), (-1e3, 1e3)),
        &Tolerance::default(),
    );
    assert!(found.is_none());
}

/// An in-band dip at the window's own seam is one touch, not a root at
/// each end of the window: the full-turn window starts at the
/// circle's extreme nearest the boss.
#[test]
fn a_dip_at_the_window_seam_touches_once() {
    use crate::geometry::curve::CurveEval as _;
    let circle = CurveGeometry::Circle {
        center: Point3::new(-817.0, -17.0, 162.5),
        normal: Vec3::new(0.0, 1.0, 0.0),
        radius: 49.0,
    };
    let t_extreme = (0..36000)
        .map(|k| k as f64 * std::f64::consts::TAU / 36000.0)
        .max_by(|a, b| circle.point_at(*a).x.total_cmp(&circle.point_at(*b).x))
        .expect("samples");
    let origin = Point3::new(0.0, -50.0 + 5e-7, 162.5);
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let found = intersect_circle_cylinder(
        CylinderFrame {
            origin: &origin,
            axis: &axis,
            radius: 33.0,
        },
        &circle,
        (t_extreme, t_extreme + std::f64::consts::TAU),
        ((0.0, std::f64::consts::TAU), (-1e3, 1e3)),
        &Tolerance::default(),
    )
    .expect("circular");
    assert_eq!(found.len(), 2, "{found:?}");
}
