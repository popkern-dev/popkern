use crate::geometry::Point3;
#[cfg(test)]
#[allow(unused_imports)]
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[path = "split_tests_face.rs"]
mod split_tests_face;

// ==========================================================================
// split_edge tests
// ==========================================================================

#[test]
fn split_edge_simple() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    let result = brep.split_edge(mev1.edge, 0.5).unwrap();

    assert_eq!(brep.vertices.len(), 3);
    assert_eq!(brep.edges.len(), 2);
    assert_eq!(brep.half_edges.len(), 4);
    let v = &brep.vertices[result.vertex];
    assert!((v.position.x - 0.5).abs() < 1e-10);
}

#[test]
fn split_edge_on_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let edge = brep.edges.keys().next().unwrap();

    brep.split_edge(edge, 0.5).unwrap();

    assert_eq!(brep.vertices.len(), v_before + 1);
    assert_eq!(brep.edges.len(), e_before + 1);
}

#[test]
fn split_edge_invalid_handle() {
    let mut brep = BRep::new(Tolerance::default());
    assert!(brep.split_edge(Default::default(), 0.5).is_err());
}

#[test]
fn split_edge_half_edge_chain_valid() {
    // After split, walking the loop should still form a cycle
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let _mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();

    brep.split_edge(mev1.edge, 0.5).unwrap();

    // Walk the loop -- should return to start
    let start = brep.loops[mvfs.loop_].half_edge;
    let mut he = start;
    let mut count = 0;
    loop {
        he = brep.half_edges[he].next;
        count += 1;
        if he == start || count > 100 {
            break;
        }
    }
    assert_eq!(he, start, "Loop should be closed after split_edge");
}

#[test]
fn split_edge_t_out_of_range() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    assert!(brep.split_edge(mev1.edge, 0.0).is_err());
    assert!(brep.split_edge(mev1.edge, 1.0).is_err());
    assert!(brep.split_edge(mev1.edge, -0.1).is_err());
    assert!(brep.split_edge(mev1.edge, 1.1).is_err());
}

// ==========================================================================
// Curve-propagation through split_edge (β-1-A per
// 2026-06-21-split-edge-surface-derived-vertex-placement.md)
// ==========================================================================

#[test]
fn split_edge_propagates_parent_curve_to_both_children() {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};

    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    // Author a line curve on the edge from (0,0,0) to (1,0,0).
    brep.edges[mev.edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    });

    let result = brep.split_edge(mev.edge, 0.4).unwrap();

    let edge_a_curve = brep.edges[result.edge_a]
        .curve
        .as_ref()
        .expect("β-1-A: original edge should retain its A→V curve portion");
    let edge_b_curve = brep.edges[result.edge_b]
        .curve
        .as_ref()
        .expect("β-1-A: new edge should inherit the V→B curve portion");

    // Edge A covers [0.0, 0.4] of the parent curve; edge B covers [0.4, 1.0].
    assert!((edge_a_curve.t_start - 0.0).abs() < 1.0e-12);
    assert!((edge_a_curve.t_end - 0.4).abs() < 1.0e-12);
    assert!((edge_b_curve.t_start - 0.4).abs() < 1.0e-12);
    assert!((edge_b_curve.t_end - 1.0).abs() < 1.0e-12);
}

#[test]
fn split_edge_preserves_none_curve_for_curveless_edge() {
    // Edges without curves stay None on both children; the
    // surface-derived path or chord fallback is unaffected for
    // curveless inputs.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    assert!(brep.edges[mev.edge].curve.is_none());

    let result = brep.split_edge(mev.edge, 0.5).unwrap();
    assert!(brep.edges[result.edge_a].curve.is_none());
    assert!(brep.edges[result.edge_b].curve.is_none());
}

// ==========================================================================
// split_edge_at_point tests (per ADR §14.2 / Stage 0b)
// ==========================================================================

/// Helper: build a one-edge BRep with a circular-arc curve on the edge,
/// returning the BRep and the edge handle.
#[cfg(test)]
fn arc_brep(
    center: Point3,
    normal: crate::geometry::Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> (BRep, crate::topology::EdgeHandle) {
    use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};

    let arc = CurveGeometry::Arc {
        center,
        normal,
        radius,
        start_angle,
        end_angle,
    };
    let p_start = arc.point_at(start_angle);
    let p_end = arc.point_at(end_angle);

    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(p_start).unwrap();
    let mev = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, p_end).unwrap();

    // Attach the arc curve to the edge so split_edge_at_point can invert.
    brep.edges[mev.edge].curve = Some(CurveOnEdge {
        geometry: arc,
        t_start: start_angle,
        t_end: end_angle,
    });

    (brep, mev.edge)
}

#[test]
fn split_edge_at_point_on_arc_midpoint() {
    use crate::geometry::Vec3;
    use crate::geometry::curve::{CurveEval, CurveGeometry};
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    let arc_geom = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let mid = arc_geom.point_at(FRAC_PI_4);

    let (mut brep, edge) = arc_brep(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );

    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();

    let result = brep
        .split_edge_at_point(edge, mid, &Tolerance::default())
        .expect("split should succeed at arc midpoint");

    assert_eq!(brep.vertices.len(), v_before + 1);
    assert_eq!(brep.edges.len(), e_before + 1);
    let v_pos = brep.vertices[result.vertex].position;
    assert!(
        v_pos.distance_to(&mid) < 1e-9,
        "new vertex should sit at the requested point: got {v_pos:?}, want {mid:?}"
    );
}

#[test]
fn split_edge_at_point_on_nurbs_arc() {
    use crate::geometry::Vec3;
    use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    // Build a NURBS curve representing a quarter arc, attach it to an edge.
    let arc_geom = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let nurbs = arc_geom.to_nurbs();
    let p_start = arc_geom.point_at(0.0);
    let p_end = arc_geom.point_at(FRAC_PI_2);
    let mid = arc_geom.point_at(FRAC_PI_4);

    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(p_start).unwrap();
    let mev = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, p_end).unwrap();
    brep.edges[mev.edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: FRAC_PI_2,
    });

    let result = brep
        .split_edge_at_point(mev.edge, mid, &Tolerance::default())
        .expect("split should succeed via NURBS-detected-as-arc tier 1b");

    let v_pos = brep.vertices[result.vertex].position;
    // NURBS detection involves a fit — accept slightly looser tolerance
    // (1e-3 m for the quarter-arc fit).
    assert!(
        v_pos.distance_to(&mid) < 1e-3,
        "new vertex on NURBS arc should sit ≈at target: got {v_pos:?}, want {mid:?}"
    );
}

#[test]
fn split_edge_at_point_off_curve_returns_error() {
    use crate::geometry::Vec3;
    use std::f64::consts::FRAC_PI_2;

    let (mut brep, edge) = arc_brep(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );

    // A point well off the unit-radius arc.
    let off_point = Point3::new(0.5, 0.5, 1.5);
    let result = brep.split_edge_at_point(edge, off_point, &Tolerance::default());
    assert!(
        matches!(result, Err(crate::error::CoreError::PreconditionFailed(_))),
        "off-curve point should fail: {result:?}"
    );
}

#[test]
fn split_edge_at_point_at_endpoint_returns_error() {
    use crate::geometry::Vec3;
    use crate::geometry::curve::CurveEval;
    use crate::geometry::curve::CurveGeometry;
    use std::f64::consts::FRAC_PI_2;

    let (mut brep, edge) = arc_brep(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );

    // Attempt to split at the start endpoint — should be rejected with
    // a "too close to existing endpoint" precondition error.
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let endpoint = arc.point_at(0.0);
    let result = brep.split_edge_at_point(edge, endpoint, &Tolerance::default());
    assert!(
        matches!(result, Err(crate::error::CoreError::PreconditionFailed(_))),
        "splitting at endpoint should fail: {result:?}"
    );
}

#[test]
fn split_edge_at_point_no_curve_returns_error() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    // No curve attached — split_edge_at_point requires one.
    let result =
        brep.split_edge_at_point(mev.edge, Point3::new(0.5, 0.0, 0.0), &Tolerance::default());
    assert!(
        matches!(result, Err(crate::error::CoreError::PreconditionFailed(_))),
        "edge without curve should fail: {result:?}"
    );
}
