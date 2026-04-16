//! NURBS-support tests for the 2D corner helpers. Exercises
//! [`tangent_at_corner_outgoing`], [`edge_arc_length`],
//! [`point_at_arc_length_from_corner`], [`setback_from_corner`],
//! and [`rebuild_edge_curve_clipped`] on edges whose curve is a
//! [`CurveGeometry::Nurbs`].
//!
//! Test vehicles rewrite edges of a `sheet_quad` to carry
//! degree-1 NURBS curves (which are line segments represented
//! in NURBS form) and higher-degree planar NURBS, then run
//! `chamfer_corner_2d` and `chamfer_corner_2d_by_angle` against
//! the result.

use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Find the edge of a `sheet_quad` whose two endpoints are
/// `p_from` and `p_to` (or the reverse), and return its handle.
fn find_edge(brep: &BRep, p_from: Point3, p_to: Point3) -> EdgeHandle {
    for (eh, edge) in &brep.edges {
        let he = edge.half_edge;
        let a = brep.vertices[brep.half_edges[he].origin].position;
        let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
        let ok = (a.distance_to(&p_from) < 1e-9 && b.distance_to(&p_to) < 1e-9)
            || (a.distance_to(&p_to) < 1e-9 && b.distance_to(&p_from) < 1e-9);
        if ok {
            return eh;
        }
    }
    panic!("no edge between {p_from:?} and {p_to:?}");
}

/// Find the vertex at a specific position in a sheet body.
fn find_vertex(brep: &BRep, pos: Point3) -> VertexHandle {
    brep.vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&pos) < 1e-9)
        .map(|(h, _)| h)
        .unwrap()
}

/// Rewrite an edge's curve so its canonical half-edge now
/// points at a degree-1 two-control-point NURBS (a line segment
/// represented in NURBS form). The resulting curve is
/// geometrically identical to the original Line, so chamfer /
/// fillet results on the modified sheet should match the plain
/// Line case.
fn rewrite_edge_as_linear_nurbs(brep: &mut BRep, edge: EdgeHandle) {
    let he = brep.edges[edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let a = brep.vertices[brep.half_edges[he].origin].position;
    let b = brep.vertices[brep.half_edges[twin].origin].position;
    let nurbs = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    });
}

// ---------------------------------------------------------------------------
// Degree-1 NURBS (line-in-NURBS-form) — should match Line case
// ---------------------------------------------------------------------------

/// Chamfer a corner where BOTH adjacent edges are degree-1
/// NURBS. The curves are geometrically identical to the
/// original straight lines, so the chamfer output must validate
/// and have the same vertex / edge counts as the plain-Line
/// version.
#[test]
fn chamfer_corner_2d_with_both_edges_as_linear_nurbs() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    // Rewrite the two edges incident at (0, 0, 0) as degree-1 NURBS.
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    let left = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    );
    rewrite_edge_as_linear_nurbs(&mut sheet, bottom);
    rewrite_edge_as_linear_nurbs(&mut sheet, left);
    sheet.validate().unwrap();

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = chamfer_corner_2d(&sheet, corner, 0.1, 0.1)
        .expect("NURBS-NURBS chamfer should succeed on degree-1 supports");

    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();

    // Both tangent points land where the plain-Line chamfer
    // would place them: (0.1, 0, 0) and (0, 0.1, 0).
    let has_t1 = result.vertices.values().any(|v| {
        (v.position.x - 0.1).abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z.abs() < 1e-9
    });
    let has_t2 = result.vertices.values().any(|v| {
        v.position.x.abs() < 1e-9 && (v.position.y - 0.1).abs() < 1e-9 && v.position.z.abs() < 1e-9
    });
    assert!(has_t1, "missing tangent point at (0.1, 0, 0)");
    assert!(has_t2, "missing tangent point at (0, 0.1, 0)");
}

/// Asymmetric chamfer on NURBS supports (dist_a = 0.2,
/// dist_b = 0.1). Exercises `point_at_arc_length_from_corner`
/// on NURBS with a non-unit walk.
#[test]
fn chamfer_corner_2d_asymmetric_with_linear_nurbs() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    let left = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    );
    rewrite_edge_as_linear_nurbs(&mut sheet, bottom);
    rewrite_edge_as_linear_nurbs(&mut sheet, left);

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = chamfer_corner_2d(&sheet, corner, 0.2, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();
}

/// Distance+angle chamfer on a Line-NURBS corner. The curved-
/// support path handles the NURBS edge via
/// `tangent_at_forward_distance` → `invert_nurbs_arc_length`
/// (when NURBS is edge_a) and via `ray_nurbs_edge_intersection`
/// (when NURBS is edge_b). Since the NURBS is degree-1 (a line),
/// both paths must produce results that match the plain-Line
/// by-angle chamfer.
#[test]
fn chamfer_corner_2d_by_angle_with_linear_nurbs_support() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    rewrite_edge_as_linear_nurbs(&mut sheet, bottom);

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = chamfer_corner_2d_by_angle(&sheet, corner, 0.1, 3.0 * std::f64::consts::FRAC_PI_4)
        .expect("by-angle chamfer should succeed on a Line-NURBS corner");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();
}

/// Distance+angle chamfer on a NURBS-NURBS corner (both adjacent
/// edges rewritten as degree-1 NURBS). The by-angle path uses
/// `tangent_at_forward_distance` on edge_a (NURBS) AND
/// `ray_nurbs_edge_intersection` on edge_b (NURBS). Since the
/// curves are geometrically straight the answer matches the
/// plain-Line case (tangent points at (0.1, 0, 0) and (0, 0.1, 0)).
#[test]
fn chamfer_corner_2d_by_angle_with_both_edges_as_linear_nurbs() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    let left = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    );
    rewrite_edge_as_linear_nurbs(&mut sheet, bottom);
    rewrite_edge_as_linear_nurbs(&mut sheet, left);

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = chamfer_corner_2d_by_angle(&sheet, corner, 0.1, 3.0 * std::f64::consts::FRAC_PI_4)
        .expect("by-angle chamfer should succeed on a NURBS-NURBS corner");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();
}

// ---------------------------------------------------------------------------
// Degree-2 planar NURBS — exercises inverse arc-length + trim
// ---------------------------------------------------------------------------

/// Rewrite an edge as a degree-2 planar NURBS arc approximation.
/// The control polygon is chosen so that the curve starts at
/// `a`, ends at `b`, and bulges toward `midpoint` (a
/// quasi-quarter-circle via a Bezier weight). The resulting
/// curve is genuinely curved, so the NURBS arc-length
/// integration and inverse arc-length paths both get exercised.
fn rewrite_edge_as_quadratic_nurbs(
    brep: &mut BRep,
    edge: EdgeHandle,
    a: Point3,
    b: Point3,
    bulge: Point3,
) {
    let he = brep.edges[edge].half_edge;
    let twin_origin = brep.half_edges[brep.half_edges[he].twin].origin;
    let canonical_from = brep.vertices[brep.half_edges[he].origin].position;
    // Honour the canonical half-edge direction so `validate()`
    // sees matching watertight endpoints.
    let (cp0, cp2) = if canonical_from.distance_to(&a) < 1e-9 {
        (a, b)
    } else {
        (b, a)
    };
    let _ = twin_origin; // suppress unused-variable lint, topology already valid
    let nurbs = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![cp0, bulge, cp2],
        weights: vec![1.0, 1.0, 1.0],
    };
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// Rewrite an edge as a degree-3 S-curve NURBS between its
/// endpoints. The two inner control points bulge in opposite
/// directions, producing a real non-convex shape that causes
/// spurious far-lobe projections for naive global searches and
/// multi-bracket sign changes for ray-curve intersections.
fn rewrite_edge_as_s_curve_nurbs(
    brep: &mut BRep,
    edge: EdgeHandle,
    a: Point3,
    b: Point3,
    bulge1: Point3,
    bulge2: Point3,
) {
    let he = brep.edges[edge].half_edge;
    let canonical_from = brep.vertices[brep.half_edges[he].origin].position;
    let (cp0, cp3) = if canonical_from.distance_to(&a) < 1e-9 {
        (a, b)
    } else {
        (b, a)
    };
    let nurbs = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![cp0, bulge1, bulge2, cp3],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// Build a `sheet_quad`, rewrite the bottom edge as a degree-2
/// NURBS that bulges outward (away from the sheet interior),
/// then chamfer the `(0, 0, 0)` corner. Exercises `edge_arc_length`,
/// `point_at_arc_length_from_corner`, and
/// `rebuild_edge_curve_clipped` on a real NURBS curve.
#[test]
fn chamfer_corner_2d_with_quadratic_nurbs_support() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    // Bulge the bottom edge downward (away from the interior).
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    rewrite_edge_as_quadratic_nurbs(
        &mut sheet,
        bottom,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, -0.2, 0.0),
    );
    // Don't revalidate the sheet — the planar pcurve on the
    // rewritten edge is stale relative to the new curved
    // geometry, which the strict watertight check would catch.
    // The chamfer path re-runs `attach_planar_pcurves` on its
    // output, so its own validate() will exercise the final
    // state.

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result =
        chamfer_corner_2d(&sheet, corner, 0.1, 0.1).expect("NURBS-backed chamfer should succeed");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();

    // After the chamfer, the bottom edge remains a NURBS curve
    // (clipped to a subrange). Count curve kinds to verify.
    let nurbs_count = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Nurbs(_))
            )
        })
        .count();
    assert_eq!(
        nurbs_count, 1,
        "expected exactly 1 NURBS edge (the clipped bottom)"
    );
}

// ---------------------------------------------------------------------------
// Multi-bracket ray-NURBS (follow-up #1)
// ---------------------------------------------------------------------------

/// Chamfer by distance+angle when `edge_b` is an S-curve NURBS.
/// Multiple sign changes of `(C(t) − t1) × ray_dir` can appear
/// along the sampled parameter range, so the multi-bracket
/// ray-NURBS intersector has to pick the smallest positive `r`
/// across all brackets rather than the first one encountered.
/// This test doesn't need to exhibit multiple brackets to pass
/// — it just verifies that the multi-bracket code path
/// produces a valid chamfer on a non-trivial NURBS support.
#[test]
fn chamfer_corner_2d_by_angle_with_s_curve_nurbs_edge_b() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    rewrite_edge_as_s_curve_nurbs(
        &mut sheet,
        bottom,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.333, -0.1, 0.0),
        Point3::new(0.666, 0.1, 0.0),
    );

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result =
        chamfer_corner_2d_by_angle(&sheet, corner, 0.05, 3.0 * std::f64::consts::FRAC_PI_4)
            .expect("by-angle chamfer should succeed on S-curve NURBS support");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();
}
