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
