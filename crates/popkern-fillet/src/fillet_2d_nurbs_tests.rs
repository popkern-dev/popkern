//! NURBS-support tests for `fillet_corner_2d`. Exercises the
//! generic Newton fillet solver in
//! [`crate::fillet_2d_solvers::solve_generic_fillet`], which
//! kicks in whenever at least one adjacent edge carries a
//! `CurveGeometry::Nurbs`.
//!
//! Test vehicles rewrite edges of a `sheet_quad` to carry
//! degree-1 NURBS curves (line segments in NURBS form) so the
//! generic Newton solver's output can be directly compared
//! against the existing closed-form results.

use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

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

fn find_vertex(brep: &BRep, pos: Point3) -> VertexHandle {
    brep.vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&pos) < 1e-9)
        .map(|(h, _)| h)
        .unwrap()
}

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

fn rewrite_edge_as_quadratic_nurbs(brep: &mut BRep, edge: EdgeHandle, bulge: Point3) {
    let he = brep.edges[edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let cp0 = brep.vertices[brep.half_edges[he].origin].position;
    let cp2 = brep.vertices[brep.half_edges[twin].origin].position;
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

// ---------------------------------------------------------------------------
// Degree-1 NURBS: must match the Line-Line closed-form result
// ---------------------------------------------------------------------------

/// Fillet a corner where BOTH adjacent edges are degree-1 NURBS
/// (geometrically straight lines). The generic Newton solver
/// must produce the same tangent points as the Line-Line closed
/// form: `(0.1, 0, 0)` and `(0, 0.1, 0)` on the unit-quad
/// `(0, 0, 0)` corner.
#[test]
fn fillet_corner_2d_with_both_edges_as_linear_nurbs() {
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
    let result = fillet_corner_2d(&sheet, corner, 0.1)
        .expect("NURBS-NURBS fillet via generic Newton should succeed");

    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();

    // The tangent points on linear NURBS supports must sit at
    // the same positions as the Line-Line closed form.
    let has_t1 = result.vertices.values().any(|v| {
        (v.position.x - 0.1).abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z.abs() < 1e-9
    });
    let has_t2 = result.vertices.values().any(|v| {
        v.position.x.abs() < 1e-9 && (v.position.y - 0.1).abs() < 1e-9 && v.position.z.abs() < 1e-9
    });
    assert!(has_t1, "missing tangent point at (0.1, 0, 0)");
    assert!(has_t2, "missing tangent point at (0, 0.1, 0)");

    // The generic solver builds a new fillet arc; verify its
    // centre matches the Line-Line closed-form `(0.1, 0.1, 0)`.
    let fillet_center = result
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { center, radius, .. }) if (*radius - 0.1).abs() < 1e-9 => {
                Some(*center)
            }
            _ => None,
        })
        .expect("no fillet arc with r=0.1");
    assert!(
        (fillet_center.x - 0.1).abs() < 1e-9,
        "fillet centre x = {}",
        fillet_center.x
    );
    assert!(
        (fillet_center.y - 0.1).abs() < 1e-9,
        "fillet centre y = {}",
        fillet_center.y
    );
}

/// Mixed supports: edge_a is a degree-1 NURBS, edge_b is the
/// plain Line left over from `sheet_quad`. The generic Newton
/// solver takes over (because one edge is NURBS) and must still
/// match the Line-Line closed form.
#[test]
fn fillet_corner_2d_with_one_edge_as_linear_nurbs() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    rewrite_edge_as_linear_nurbs(&mut sheet, bottom);

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = fillet_corner_2d(&sheet, corner, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();
}

// ---------------------------------------------------------------------------
// Degree-2 NURBS: genuinely curved support
// ---------------------------------------------------------------------------

/// Fillet with a degree-2 NURBS on one side. The quadratic
/// bulges downward from the bottom edge of the unit quad; the
/// fillet centre and tangent points land somewhere off the
/// straight-line answer but the chamfer must still produce a
/// topologically valid result with the fillet radius preserved.
#[test]
fn fillet_corner_2d_with_quadratic_nurbs_support() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    rewrite_edge_as_quadratic_nurbs(&mut sheet, bottom, Point3::new(0.5, -0.2, 0.0));

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = fillet_corner_2d(&sheet, corner, 0.08)
        .expect("generic Newton should converge on a quadratic NURBS support");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();

    // The fillet arc's radius must match the requested 0.08 by
    // construction — the generic solver only terminates when the
    // residuals |F − closest_point| − r are below tolerance on
    // both supports.
    let fillet_radius = result
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { radius, .. }) => Some(*radius),
            _ => None,
        })
        .expect("no fillet arc found");
    assert!(
        (fillet_radius - 0.08).abs() < 1e-9,
        "fillet radius = {fillet_radius}, expected 0.08"
    );

    // The clipped NURBS edge is still a NURBS.
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
// Robustness: non-convex NURBS + local projection (follow-up #3)
// ---------------------------------------------------------------------------

/// Build a degree-3 NURBS with **two** internal control points
/// bulging alternately up and down, producing an S-shaped curve
/// between `a` and `b`. The curve is non-convex near one end,
/// so the global closest-point search for a fillet centre could
/// find a spurious match on the far lobe if the solver used
/// `CurveEval::closest_point` directly. The local-projection
/// path seeded from the corner parameter avoids this pitfall.
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

/// Fillet a corner where one support is an S-shaped degree-3
/// NURBS. Without the local-projection seed the generic solver
/// could find a "spurious" closest point on the far lobe of the
/// S; with local projection it tracks the near-corner branch
/// and converges to the correct fillet centre.
#[test]
fn fillet_corner_2d_with_s_curve_nurbs_support() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let bottom = find_edge(
        &sheet,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    );
    // S-curve: bulge down near x=1/3, bulge up near x=2/3.
    rewrite_edge_as_s_curve_nurbs(
        &mut sheet,
        bottom,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.333, -0.15, 0.0),
        Point3::new(0.666, 0.15, 0.0),
    );

    let corner = find_vertex(&sheet, Point3::new(0.0, 0.0, 0.0));
    let result = fillet_corner_2d(&sheet, corner, 0.05)
        .expect("fillet with S-curve NURBS support should succeed via local projection");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();

    // Fillet radius preserved to tolerance.
    let fillet_radius = result
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { radius, .. }) => Some(*radius),
            _ => None,
        })
        .expect("no fillet arc");
    assert!(
        (fillet_radius - 0.05).abs() < 1e-9,
        "fillet radius = {fillet_radius}, expected 0.05"
    );
}
