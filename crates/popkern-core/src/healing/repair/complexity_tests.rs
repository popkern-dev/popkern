use super::*;
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

// ---------------------------------------------------------------------------
// Helper: build a simple degree-2 NURBS line (3 control points, no interior
// knots) — this is the minimal representation.
// ---------------------------------------------------------------------------

fn minimal_nurbs_line() -> NurbsCurve {
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    }
}

// ---------------------------------------------------------------------------
// Curve tests
// ---------------------------------------------------------------------------

#[test]
fn minimal_curve_returns_none() {
    let curve = minimal_nurbs_line();
    let result = reduce_curve_complexity(&curve, 1e-6);
    assert!(result.is_none(), "already-minimal curve should return None");
}

#[test]
fn inserted_knots_are_removed() {
    // Start with a simple line, insert several knots (which doesn't change
    // the shape), then verify that reduce_curve_complexity removes them.
    let base = minimal_nurbs_line();
    let inflated = base
        .insert_knot(0.25)
        .and_then(|c| c.insert_knot(0.5))
        .and_then(|c| c.insert_knot(0.75))
        .unwrap_or_else(|_| base.clone());

    assert!(
        inflated.control_points.len() > base.control_points.len(),
        "insertion should add control points"
    );

    let result = reduce_curve_complexity(&inflated, 1e-6);
    assert!(result.is_some(), "should remove inserted knots");

    let r = result.unwrap_or(CurveReductionResult {
        points_before: 0,
        points_after: 0,
        knots_removed: 0,
        max_deviation: 0.0,
    });
    assert!(r.knots_removed > 0, "should have removed at least one knot");
    assert!(
        r.points_after < r.points_before,
        "should have fewer control points after"
    );
    assert!(
        r.max_deviation < 1e-6,
        "deviation {} should be near zero for shape-preserving knots",
        r.max_deviation
    );
}

#[test]
fn curve_reduction_result_fields() {
    let base = minimal_nurbs_line();
    let inflated = base.insert_knot(0.5).unwrap_or_else(|_| base.clone());

    if let Some(r) = reduce_curve_complexity(&inflated, 1e-6) {
        assert_eq!(r.points_before, inflated.control_points.len());
        assert!(r.points_after <= r.points_before);
        assert!(r.max_deviation >= 0.0);
    }
}

// ---------------------------------------------------------------------------
// Surface tests
// ---------------------------------------------------------------------------

/// Build a flat bilinear NURBS patch (degree 1x1, no interior knots).
fn minimal_nurbs_patch() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    }
}

#[test]
fn minimal_surface_returns_none() {
    let surface = minimal_nurbs_patch();
    let result = reduce_surface_complexity(&surface, 1e-6);
    assert!(
        result.is_none(),
        "already-minimal surface should return None"
    );
}

#[test]
fn surface_inserted_knots_are_removed() {
    let base = minimal_nurbs_patch();
    // Insert knots in both directions.
    let inflated = base
        .insert_knot_u(0.5)
        .and_then(|s| s.insert_knot_v(0.5))
        .unwrap_or_else(|_| base.clone());

    assert!(
        inflated.num_u > base.num_u || inflated.num_v > base.num_v,
        "insertion should increase grid dimensions"
    );

    let result = reduce_surface_complexity(&inflated, 1e-6);
    assert!(result.is_some(), "should remove inserted surface knots");

    let r = result.unwrap_or(SurfaceReductionResult {
        points_before: (0, 0),
        points_after: (0, 0),
        u_knots_removed: 0,
        v_knots_removed: 0,
        max_deviation: 0.0,
    });
    let total = r.u_knots_removed + r.v_knots_removed;
    assert!(total > 0, "should have removed at least one knot");
    assert!(
        r.max_deviation < 1e-6,
        "deviation {} should be near zero for shape-preserving knots",
        r.max_deviation
    );
}

// ---------------------------------------------------------------------------
// BRep-wide tests
// ---------------------------------------------------------------------------

#[test]
fn reduce_all_curves_on_cube() {
    // Cube edges are lines, not NURBS — should produce no reductions.
    let mut brep =
        BRep::unit_cube(Tolerance::default()).unwrap_or_else(|_| BRep::new(Tolerance::default()));
    let result = reduce_all_curves(&mut brep, 1e-6);
    assert_eq!(result.curves_reduced, 0);
    assert_eq!(result.total_knots_removed, 0);
}

#[test]
fn reduce_all_curves_on_empty_brep() {
    let mut brep = BRep::new(Tolerance::default());
    let result = reduce_all_curves(&mut brep, 1e-6);
    assert_eq!(result.curves_reduced, 0);
    assert_eq!(result.surfaces_reduced, 0);
    assert_eq!(result.total_knots_removed, 0);
}

// ---------------------------------------------------------------------------
// Interior knot index helper tests
// ---------------------------------------------------------------------------

#[test]
fn interior_indices_empty_for_minimal_knots() {
    // degree 2, knots = [0,0,0,1,1,1] → no interior
    let indices = interior_knot_indices(2, &[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    assert!(indices.is_empty());
}

#[test]
fn interior_indices_finds_one() {
    // degree 2, knots = [0,0,0, 0.5, 1,1,1] → index 3 is interior
    let indices = interior_knot_indices(2, &[0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]);
    assert_eq!(indices, vec![3]);
}

#[test]
fn interior_indices_deduplicates() {
    // degree 2, knots with a repeated interior value
    let indices = interior_knot_indices(2, &[0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]);
    // Should only report the first occurrence of 0.5
    assert_eq!(indices, vec![3]);
}
