//! Tangential contact consolidation (PK-732): a curve touching a surface
//! is one root at the touch; two genuine crossings close together stay
//! two; a touch at a closed curve's own seam is one root, not one at
//! each end of the parameter window.

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::surface_curve::{
    SurfaceCurveIntersection, intersect_surface_curve_dense,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// The equal-radius tee's main cylinder: r = 40 about z.
fn main_cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 40.0,
    }
}

/// A circle of radius 40 about x through (x, 0, 60): at x = 0 the
/// branch end circle, touching the main cylinder at (0, ±40, 60).
fn end_circle(x: f64) -> CurveGeometry {
    CurveGeometry::Circle {
        center: Point3::new(x, 0.0, 60.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        radius: 40.0,
    }
}

fn roots(curve: &CurveGeometry, t_range: (f64, f64)) -> Vec<(f64, Point3)> {
    let tau = std::f64::consts::TAU;
    let found = intersect_surface_curve_dense(
        &main_cylinder(),
        (0.0, tau),
        (-50.0, 170.0),
        curve,
        t_range,
        &Tolerance::default(),
    )
    .expect("intersection");
    found
        .intersections
        .iter()
        .filter_map(|r| match *r {
            SurfaceCurveIntersection::Point { t, point, .. } => Some((t, point)),
            SurfaceCurveIntersection::Overlap { .. } => None,
        })
        .collect()
}

#[test]
fn a_touching_circle_meets_the_cylinder_once_at_each_touch() {
    let tau = std::f64::consts::TAU;
    let found = roots(&end_circle(0.0), (0.0, tau));
    assert_eq!(found.len(), 2, "one root per touch: {found:?}");
    for want in [Point3::new(0.0, 40.0, 60.0), Point3::new(0.0, -40.0, 60.0)] {
        let d = found
            .iter()
            .map(|(_, p)| p.distance_to(&want))
            .fold(f64::INFINITY, f64::min);
        assert!(d < 1e-6, "the touch at {want:?} is {d} off");
    }
}

#[test]
fn close_genuine_crossings_stay_apart() {
    // Offset 0.05: the circle crosses the cylinder at z = 60 ± 0.05
    // beside each former touch, bulging 3e-5 outside it between —
    // beyond the tolerance, so four crossings.
    let tau = std::f64::consts::TAU;
    let found = roots(&end_circle(0.05), (0.0, tau));
    assert_eq!(found.len(), 4, "{found:?}");
    for (_, p) in &found {
        assert!(((p.z - 60.0).abs() - 0.05).abs() < 1e-6, "{p:?}");
    }
}

#[test]
fn a_touch_at_the_closed_curves_seam_is_one_root() {
    // Start the window at the touch (0, 40, 60), as the tee's end
    // circle edge does: Newton converged just after the start and just
    // before the end — one contact.
    let circle = end_circle(0.0);
    let tau = std::f64::consts::TAU;
    let start = (0..3600)
        .map(|k| -std::f64::consts::PI + k as f64 * tau / 3600.0)
        .min_by(|a, b| {
            let d = |t: f64| {
                circle
                    .point_at(t)
                    .distance_to(&Point3::new(0.0, 40.0, 60.0))
            };
            d(*a).total_cmp(&d(*b))
        })
        .expect("start");
    let found = roots(&circle, (start, start + tau));
    assert_eq!(found.len(), 2, "start {start}: {found:?}");
    let at_seam = found
        .iter()
        .filter(|(_, p)| p.distance_to(&Point3::new(0.0, 40.0, 60.0)) < 1e-6)
        .count();
    assert_eq!(at_seam, 1, "{found:?}");
}

/// A ball's seam meridian touching an ice-cream cone at its rim, a
/// thousand units from the origin: one root, at the touch to 1e-9, not
/// the 1.7e-6 off it a minimum located from distance values gave
/// (PK-1133).
#[test]
fn a_meridian_touching_a_cone_far_from_the_origin_lands_on_the_touch() {
    use crate::geometry::Transform3;
    let tol = Tolerance::default();
    let o = Vec3::new(1234.5, -678.25, 910.125);
    let ball = crate::topology::BRep::sphere(15.0, tol)
        .and_then(|b| b.transform(&Transform3::translation(o.x, o.y, o.z)))
        .expect("ball");
    let seam = ball
        .edges
        .values()
        .find_map(|e| e.curve.clone())
        .expect("seam");
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(o.x, o.y, o.z - 25.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.75f64.atan(),
    };
    let found = intersect_surface_curve_dense(
        &cone,
        (0.0, std::f64::consts::TAU),
        (-8.0, 24.0),
        &seam.geometry,
        (seam.t_start, seam.t_end),
        &tol,
    )
    .expect("intersects");
    let touch = Point3::new(o.x + 12.0, o.y, o.z - 9.0);
    let points: Vec<Point3> = found
        .intersections
        .iter()
        .filter_map(|i| match i {
            SurfaceCurveIntersection::Point { point, .. } => Some(*point),
            SurfaceCurveIntersection::Overlap { .. } => None,
        })
        .collect();
    assert_eq!(points.len(), 1, "{points:?}");
    assert!(
        points[0].distance_to(&touch) < 1e-9,
        "{:e}",
        points[0].distance_to(&touch)
    );
}

/// A NURBS circle (the general solver's, not the analytic arm's) of
/// radius 4 poking past the plane x = 4 − δ: half a tolerance past it,
/// one touch at the extreme, as the analytic circle × plane arm decides
/// it; two tolerances past, two crossings (PK-1128: a ball's seam
/// meridian half a tolerance past a plate's side face read as two
/// crossings 4e-3 apart and the boolean declined non-manifold).
#[test]
fn a_sub_tolerance_bulge_of_a_nurbs_circle_is_one_touch() {
    let tol = Tolerance::default();
    let circle = CurveGeometry::Nurbs(super::super::circle_of_revolution::nurbs_circle(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        4.0,
    ));
    let range = (0.0, 1.0);
    for (k, want) in [(0.5, 1usize), (2.0, 2)] {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(4.0 - k * tol.linear(), 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let found = intersect_surface_curve_dense(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &circle,
            range,
            &tol,
        )
        .expect("intersects");
        let points: Vec<Point3> = found
            .intersections
            .iter()
            .filter_map(|i| match i {
                SurfaceCurveIntersection::Point { point, .. } => Some(*point),
                SurfaceCurveIntersection::Overlap { .. } => None,
            })
            .collect();
        assert_eq!(points.len(), want, "{k} tol past: {points:?}");
        if want == 1 {
            assert!(
                points[0].distance_to(&Point3::new(4.0, 0.0, 0.0)) < 1e-9,
                "{points:?}"
            );
        }
    }
}

/// Viviani's cylinder: r = 5 about the line x = 5, y = 0.
fn viviani_cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(5.0, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    }
}

/// A ball's seam meridian (radius 10, a rational half circle) from its
/// pole at (10, 0, 0), where Viviani's cylinder touches the ball: the
/// touch is the curve's own start, found at t = 0 exactly (PK-1139: no
/// sign change of C′ · N brackets a root at a window's end, and the
/// golden-section seed sat 2.6e-6 along, 7.5e-5 off the pole).
#[test]
fn a_touch_at_a_curves_own_end_lands_on_the_end() {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let meridian = CurveGeometry::Nurbs(crate::geometry::nurbs_curve::NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, -10.0, 0.0),
            Point3::new(0.0, -10.0, 0.0),
            Point3::new(-10.0, -10.0, 0.0),
            Point3::new(-10.0, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0],
    });
    let tol = Tolerance::default();
    for range in [(0.0, 1.0), (0.0, 0.05)] {
        let found = intersect_surface_curve_dense(
            &viviani_cylinder(),
            (0.0, std::f64::consts::TAU),
            (0.0, 24.0),
            &meridian,
            range,
            &tol,
        )
        .expect("intersects");
        let points: Vec<(f64, Point3)> = found
            .intersections
            .iter()
            .filter_map(|i| match i {
                SurfaceCurveIntersection::Point { t, point, .. } => Some((*t, *point)),
                _ => None,
            })
            .collect();
        assert_eq!(points.len(), 1, "{range:?}: {points:?}");
        assert!(
            points[0].0.abs() < 1e-15,
            "{range:?}: the touch is the curve's start, t {}",
            points[0].0
        );
        assert!(points[0].1.distance_to(&Point3::new(10.0, 0.0, 0.0)) < 1e-12);
    }
}

/// A transversal crossing a ten-millionth of a unit inside a line's
/// start, at 30 degrees: the band reaches the window's end but the
/// crossing keeps its place (C′ · N is the crossing's sine at both of
/// the band's edges, never vanishing at the end).
#[test]
fn a_crossing_just_inside_a_window_end_keeps_its_place() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let d = Vec3::new(0.75_f64.sqrt(), 0.0, 0.5);
    let start = Point3::new(0.0, 0.0, 0.0) - d * 1e-7;
    let line = CurveGeometry::Line {
        origin: start,
        direction: d,
    };
    let found = intersect_surface_curve_dense(
        &plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &line,
        (0.0, 5.0),
        &Tolerance::default(),
    )
    .expect("intersects");
    let t: Vec<f64> = found
        .intersections
        .iter()
        .filter_map(|i| match i {
            SurfaceCurveIntersection::Point { t, .. } => Some(*t),
            _ => None,
        })
        .collect();
    assert_eq!(t.len(), 1, "{t:?}");
    assert!((t[0] - 1e-7).abs() < 1e-12, "the crossing at {:.3e}", t[0]);
}
