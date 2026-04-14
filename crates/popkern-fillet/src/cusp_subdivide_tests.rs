//! Tests for the §13.1(1) cusp-subdivide repair.

use super::*;
use crate::types::FilletRadius;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// find_interior_cusp: negative path (clean input)
// ---------------------------------------------------------------------------

#[test]
fn find_interior_cusp_returns_none_on_straight_cube_edge() {
    // Every cube edge is a straight line with constant non-zero
    // derivative. The detector must not fire anywhere.
    let brep = BRep::unit_cube(tol()).unwrap();
    for edge in brep.edges.keys() {
        let result = find_interior_cusp(&brep, edge, &tol());
        assert!(
            result.is_none(),
            "cube edge {edge:?} false-positived interior cusp: t = {result:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// find_interior_cusp: positive path (synthetic cubic Bezier cusp)
// ---------------------------------------------------------------------------

/// A cubic Bezier with control points `(P0, P1, P1, P0)` — where
/// the middle two control points coincide and the endpoints match
/// — has first derivative `B'(t) = 3·(1 − 2t)·(P1 − P0)` which
/// vanishes exactly at `t = 0.5`. This is a genuine cusp in the
/// "vanishing derivative" sense: the curve doubles back on itself
/// at the midpoint. Suitable as a synthetic test input for the
/// interior-cusp detector.
fn cusp_cubic_bezier() -> NurbsCurve {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    // Degree-3 clamped-uniform Bezier: knots [0, 0, 0, 0, 1, 1, 1, 1].
    NurbsCurve {
        degree: 3,
        control_points: vec![p0, p1, p1, p0],
        weights: vec![1.0, 1.0, 1.0, 1.0],
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
    }
}

#[test]
fn cusp_cubic_bezier_has_zero_derivative_at_midpoint() {
    // Sanity check: our synthetic test curve really does have a
    // vanishing first derivative at t = 0.5. If this breaks, the
    // positive-path detector test below would silently become
    // uninformative.
    let curve = cusp_cubic_bezier();
    let d = curve.derivative_at(0.5, 1);
    assert!(
        d.length() < 1e-9,
        "cusp Bezier derivative at t=0.5 should vanish, got magnitude {}",
        d.length()
    );
}

#[test]
fn find_interior_cusp_detects_cubic_bezier_midpoint_cusp() {
    // Synthetic strategy: build a cube, replace one edge's curve
    // with the cuspidal cubic Bezier, then run the detector on
    // that edge. The cusp at t = 0.5 should be reported as the
    // returned normalized parameter.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(cusp_cubic_bezier()),
        t_start: 0.0,
        t_end: 1.0,
    });

    let result = find_interior_cusp(&brep, edge, &tol());
    let t = result.expect("interior cusp detector should fire on cuspidal curve");
    // The 20-sample grid lands exactly on t = 10/20 = 0.5, which
    // is the zero-derivative point. Allow a slight tolerance in
    // case the sweep samples slightly off from the exact midpoint.
    assert!(
        (t - 0.5).abs() < 1e-6,
        "detector should return t ≈ 0.5, got {t}"
    );
}

// ---------------------------------------------------------------------------
// subdivide_cusp_edges: no-op invariant (negative)
// ---------------------------------------------------------------------------

#[test]
fn subdivide_cusp_edges_is_noop_on_healthy_cube() {
    // The cube has no cusps anywhere, so the pre-pass should
    // return:
    //   1. A BRep with the same edge count as the input.
    //   2. A spec list identical to the input (same length, same
    //      edge handles, same radii).
    let brep = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = brep
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();

    let original_edge_count = brep.edges.len();
    let (healed_brep, new_specs) = subdivide_cusp_edges(&brep, &specs, &tol());
    assert_eq!(
        healed_brep.edges.len(),
        original_edge_count,
        "cusp pre-pass should not add edges on a clean cube"
    );
    assert_eq!(
        new_specs.len(),
        specs.len(),
        "cusp pre-pass should not change spec count on a clean cube"
    );
    for (orig, new) in specs.iter().zip(new_specs.iter()) {
        assert_eq!(orig.edge, new.edge);
    }
}

// ---------------------------------------------------------------------------
// subdivide_cusp_edges: positive path
// ---------------------------------------------------------------------------

#[test]
fn subdivide_cusp_edges_splits_cuspidal_bezier_edge_into_two() {
    // Build a cube, swap one edge's curve for the cuspidal cubic,
    // and feed that edge's spec to the pre-pass. The pre-pass
    // should clone the BRep, split the edge at t = 0.5, and
    // return a spec list with TWO entries instead of one — one
    // for each sub-edge created by the split.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(cusp_cubic_bezier()),
        t_start: 0.0,
        t_end: 1.0,
    });
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let original_edge_count = brep.edges.len();
    let (healed_brep, new_specs) = subdivide_cusp_edges(&brep, &specs, &tol());

    // Edge count should increase by 1 (original edge becomes two
    // sub-edges).
    assert_eq!(
        healed_brep.edges.len(),
        original_edge_count + 1,
        "cusp subdivide should add exactly one edge on a single-cusp input, \
         original {original_edge_count} → new {}",
        healed_brep.edges.len()
    );
    // Spec list should have grown from 1 to 2.
    assert_eq!(
        new_specs.len(),
        2,
        "single cuspidal spec should expand to two sub-specs"
    );
    // Both sub-specs should share the original radius.
    for sub in &new_specs {
        match sub.radius {
            FilletRadius::Constant(r) => assert!((r - 0.1).abs() < 1e-12),
            _ => panic!("sub-spec should inherit constant radius"),
        }
    }
    // Both sub-spec edges should exist in the healed BRep.
    for sub in &new_specs {
        assert!(
            healed_brep.edges.contains_key(sub.edge),
            "sub-spec edge {:?} should exist in healed BRep",
            sub.edge
        );
    }
}
