//! Tests for check 9 (curve-surface consistency), including the
//! PK-235 collect-all contract: every failing edge is reported with
//! its own per-edge worst deviation.

use crate::error::ValidationError;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

fn make_cube() -> BRep {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    brep
}

/// Replace `edge`'s curve with a quadratic whose mid control point is
/// lifted `bulge` off the chord along (1,1,1)/√3 — off every cube
/// face plane by construction.
fn bulge_edge(brep: &mut BRep, edge: EdgeHandle, bulge: f64) {
    let half_edge = brep.edges[edge].half_edge;
    let start = brep.vertices[brep.half_edges[half_edge].origin].position;
    let twin = brep.half_edges[half_edge].twin;
    let end = brep.vertices[brep.half_edges[twin].origin].position;
    let s = bulge / 3.0_f64.sqrt();
    let lifted = Point3::new(
        0.5 * (start.x + end.x) + s,
        0.5 * (start.y + end.y) + s,
        0.5 * (start.z + end.z) + s,
    );
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(crate::geometry::NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![start, lifted, end],
            weights: vec![1.0, 1.0, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    });
}

#[test]
fn topology_only_cube_passes_curve_surface() {
    let brep = make_cube();
    assert!(brep.validate_curve_surface().is_ok());
    assert!(brep.validate_curve_surface_all().is_empty());
}

#[test]
fn geometric_cube_passes_curve_surface() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    assert!(brep.validate_curve_surface().is_ok());
    assert!(brep.validate_curve_surface_all().is_empty());
}

#[test]
fn off_surface_edge_curve_is_rejected() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    bulge_edge(&mut brep, edge, 0.3 * 3.0_f64.sqrt());
    assert!(
        brep.validate_curve_surface().is_err(),
        "an edge curve bulging off its face planes must fail check 9"
    );
}

/// PK-235 — the collect-all form reports EVERY failing edge, each
/// with its own per-edge worst deviation. The early-exit form once
/// fed `validate_all`, which capped the report at one error carrying
/// the first-found (not worst) deviation.
#[test]
fn collect_all_reports_every_failing_edge_with_its_own_worst() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edges: Vec<EdgeHandle> = brep.edges.keys().take(2).collect();
    let bulge_small = 0.05;
    let bulge_large = 0.2;
    bulge_edge(&mut brep, edges[0], bulge_small);
    bulge_edge(&mut brep, edges[1], bulge_large);

    let errors = brep.validate_curve_surface_all();
    assert_eq!(
        errors.len(),
        2,
        "both broken edges must be reported; got {errors:?}"
    );
    let mut by_edge = std::collections::BTreeMap::new();
    for err in &errors {
        let ValidationError::CurveSurfaceDeviation {
            edge,
            max_deviation,
        } = err
        else {
            panic!("unexpected error kind from check 9: {err:?}");
        };
        by_edge.insert(*edge, *max_deviation);
    }
    let d_small = by_edge[&edges[0]];
    let d_large = by_edge[&edges[1]];
    // Per-edge worst: the larger bulge must report the larger
    // deviation, and each must be within its own bulge magnitude
    // (the deviation is at most the mid-CP lift).
    assert!(
        d_large > d_small,
        "per-edge worst deviations must be independent: small={d_small:.3e} large={d_large:.3e}"
    );
    assert!(
        d_small > 0.5 * bulge_small * 0.5 && d_small <= bulge_small,
        "small-bulge deviation out of range: {d_small:.3e}"
    );
    assert!(
        d_large > 0.5 * bulge_large * 0.5 && d_large <= bulge_large,
        "large-bulge deviation out of range: {d_large:.3e}"
    );

    // The early-exit form still reports exactly one of them.
    assert!(brep.validate_curve_surface().is_err());

    // validate_all carries the full set (both edges present).
    let all = brep.validate_all();
    let csd_count = all
        .iter()
        .filter(|e| matches!(e, ValidationError::CurveSurfaceDeviation { .. }))
        .count();
    assert_eq!(
        csd_count, 2,
        "validate_all must surface every failing edge; got {all:?}"
    );
}
