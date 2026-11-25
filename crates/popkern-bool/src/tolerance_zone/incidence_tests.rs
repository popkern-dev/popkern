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

/// A point of an edge carrying its own tolerance is placed only to it
/// along the edge: its zone takes the edge's tolerance (a thousandfold
/// NURBS edge of 9.1e-4, 1440). Control: on the same edge without one,
/// the zone stays the radius.
#[test]
fn a_point_of_an_imprecise_edge_has_the_edges_tolerance_as_its_zone() {
    use super::candidates::{CandidateVertex, Incidence, SourceKind};
    use super::incidence::add_crossing_zones;
    use super::tests::{TOL, cube};
    let mut body = cube();
    let (edge, at) = body
        .edges
        .iter()
        .find_map(|(h, e)| {
            let c = e.curve.as_ref()?;
            use popkern_core::geometry::curve::CurveEval;
            Some((h, c.geometry.point_at(0.5 * (c.t_start + c.t_end))))
        })
        .expect("an edge");
    let cand = || CandidateVertex {
        position: at,
        radius: TOL,
        zone: TOL,
        source: SourceKind::RecordEndpoint,
        incidence: vec![Incidence::Edge {
            operand: 0,
            edge,
            t: 0.5,
        }],
        key: (SourceKind::RecordEndpoint, 0, 0, 0, 0),
    };
    let mut precise = [cand()];
    add_crossing_zones(&mut precise, [&body, &body], 1e-3);
    assert_eq!(precise[0].zone, TOL);
    body.edges[edge].tolerance = Some(9.1e-4);
    let mut imprecise = [cand()];
    add_crossing_zones(&mut imprecise, [&body, &body], 1e-3);
    assert_eq!(imprecise[0].zone, 9.1e-4);
}
