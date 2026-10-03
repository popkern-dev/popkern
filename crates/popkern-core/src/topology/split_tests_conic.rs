//! PK-646: splitting a circle or ellipse edge keeps each child on its
//! own span of the parent's exact geometry. The rational form is one
//! turn on [0, 2π] whatever the trim, with parameter ≠ angle inside
//! each quarter, so routing the halves through it left an edge
//! trimmed past one turn holding the whole turn (measured: a closed
//! ellipse rim on jack_in_the_box trimmed [2.98, 6.34] split at 2π, a
//! 3.44 shell gap) and cut off-quarter splits off their vertex.

use super::*;
use crate::geometry::Vec3;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};

/// A one-edge body carrying `geometry` trimmed to `(t0, t1)`.
fn conic_brep(
    geometry: CurveGeometry,
    (t0, t1): (f64, f64),
) -> (BRep, crate::topology::EdgeHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(geometry.point_at(t0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, geometry.point_at(t1))
        .unwrap();
    brep.edges[mev.edge].curve = Some(CurveOnEdge {
        geometry,
        t_start: t0,
        t_end: t1,
    });
    (brep, mev.edge)
}

/// Split at `fraction` and require both children's curves to run
/// between their own vertices along the parent's geometry.
fn check_split(geometry: CurveGeometry, range: (f64, f64), fraction: f64) {
    let (mut brep, edge) = conic_brep(geometry.clone(), range);
    let split = brep.split_edge(edge, fraction).expect("split");
    let expected = geometry.point_at(range.0 + fraction * (range.1 - range.0));
    assert!(
        brep.vertices[split.vertex].position.distance_to(&expected) < 1e-12,
        "{range:?} at {fraction}: split vertex off the curve"
    );
    for child in [split.edge_a, split.edge_b] {
        let curve = brep.edges[child].curve.as_ref().expect("child curve");
        assert!(
            matches!(
                curve.geometry,
                CurveGeometry::Circle { .. } | CurveGeometry::Ellipse { .. }
            ),
            "{range:?} at {fraction}: the child keeps the exact conic"
        );
        let he = &brep.half_edges[brep.edges[child].half_edge];
        let ends = [
            brep.vertices[he.origin].position,
            brep.vertices[brep.half_edges[he.twin].origin].position,
        ];
        let (s, e) = (
            curve.geometry.point_at(curve.t_start),
            curve.geometry.point_at(curve.t_end),
        );
        let gap = (s.distance_to(&ends[0]).max(e.distance_to(&ends[1])))
            .min(s.distance_to(&ends[1]).max(e.distance_to(&ends[0])));
        assert!(
            gap < 1e-12,
            "{range:?} at {fraction}: child {child:?} gap {gap}"
        );
        assert!(
            (curve.t_end - curve.t_start).abs() < (range.1 - range.0).abs(),
            "{range:?} at {fraction}: child {child:?} spans its share, not the turn"
        );
    }
}

#[test]
fn an_ellipse_trimmed_past_one_turn_splits_at_the_turn() {
    let ellipse = CurveGeometry::Ellipse {
        center: Point3::new(1.0, 2.0, 3.0),
        u_basis: Vec3::new(3.0, 0.0, 0.0),
        v_basis: Vec3::new(0.0, 1.5, 0.0),
        start_angle: 0.0,
        end_angle: std::f64::consts::TAU,
    };
    let range = (2.9826190236611296, 6.343564382764563);
    let at_turn = (std::f64::consts::TAU - range.0) / (range.1 - range.0);
    check_split(ellipse.clone(), range, at_turn);
    check_split(ellipse, range, 0.5);
}

#[test]
fn a_circle_splits_off_its_quarter_knots_on_its_vertex() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    check_split(circle.clone(), (0.3, 2.9), 0.37);
    check_split(circle, (4.71238898038469, 10.995574287564276), 0.25);
}
