//! The edge clearance proof (`incidence::clear_of`): a point is cleared
//! of an edge only when no point of the edge lies within reach.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::{Point3, Vec3};

use super::incidence::clear_of;

const REACH: f64 = 1e-5;

fn rims() -> [CurveGeometry; 2] {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 30.0,
    };
    let nurbs = CurveGeometry::Nurbs(circle.to_nurbs_in_range((0.0, std::f64::consts::TAU)));
    [circle, nurbs]
}

/// A point inside a disc's rim, a few millimetres in, is clear of it,
/// both as a circle and as the rational NURBS a primitive builds.
#[test]
fn a_point_well_inside_a_rim_is_clear_of_it() {
    for rim in rims() {
        let window = (0.0, std::f64::consts::TAU);
        assert!(
            clear_of(&rim, window, &Point3::new(24.5, 0.0, 0.0), REACH),
            "{rim:?}"
        );
        assert!(
            clear_of(&rim, window, &Point3::new(0.0, 0.0, 0.0), REACH),
            "{rim:?}"
        );
    }
}

/// Negative controls: a point on the rim, or within reach of it, is
/// never cleared, whatever the curve's form.
#[test]
fn a_point_on_or_within_reach_of_a_rim_is_not_clear_of_it() {
    for rim in rims() {
        let window = (0.0, std::f64::consts::TAU);
        let a = 1.234_f64;
        for r in [30.0, 30.0 + 0.5 * REACH, 30.0 - 0.5 * REACH] {
            let p = Point3::new(r * a.cos(), r * a.sin(), 0.0);
            assert!(!clear_of(&rim, window, &p, REACH), "{rim:?} r={r}");
        }
    }
}

/// Only the edge's window counts: a point of the other half of the
/// circle is clear of a quarter-rim edge, and not of the half holding
/// it.
#[test]
fn a_point_beside_the_rest_of_the_circle_is_clear_of_the_edge() {
    use popkern_core::geometry::curve::CurveEval;
    use std::f64::consts::PI;
    let [_, nurbs] = rims();
    let p = nurbs.point_at(PI);
    assert!(clear_of(&nurbs, (0.0, 0.5 * PI), &p, REACH));
    assert!(clear_of(&nurbs, (1.5 * PI, 2.0 * PI), &p, REACH));
    assert!(!clear_of(&nurbs, (0.5 * PI, 1.5 * PI), &p, REACH));
}
