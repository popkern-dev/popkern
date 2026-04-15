use super::*;
use crate::geometry::surface::SurfaceGeometry;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// sheet_quad
// ---------------------------------------------------------------------------

#[test]
fn sheet_quad_has_expected_topology_counts() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    // 2 faces (front + backside), 4 edges, 4 vertices, 1 shell, 1 solid.
    assert_eq!(sheet.faces.len(), 2);
    assert_eq!(sheet.edges.len(), 4);
    assert_eq!(sheet.vertices.len(), 4);
    assert_eq!(sheet.shells.len(), 1);
    assert_eq!(sheet.solids.len(), 1);
    // 8 half-edges total (4 edges × 2 halves).
    assert_eq!(sheet.half_edges.len(), 8);
}

#[test]
fn sheet_quad_validates() {
    let sheet = BRep::sheet_quad(0.0, 2.0, 0.0, 3.0, 0.5, tol()).unwrap();
    // Full validate passes — implies manifold edges, loop closure,
    // half-edge pairing, vertex valence, and Euler identity all hold.
    sheet.validate().unwrap();
}

#[test]
fn sheet_quad_euler_identity() {
    // For a topological sphere (which a 2-face planar paper-thin
    // shell is, topologically: two hemispheres joined at all 4 edges),
    // V - E + F = 2. With 4 vertices, 4 edges, 2 faces: 4 − 4 + 2 = 2. ✓
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let euler = sheet.vertices.len() as i64 - sheet.edges.len() as i64 + sheet.faces.len() as i64;
    assert_eq!(euler, 2);
}

#[test]
fn sheet_quad_front_face_has_plane_surface() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    // At least one face has a Plane surface attached.
    let with_plane: usize = sheet
        .faces
        .values()
        .filter(|f| {
            f.surface
                .as_ref()
                .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
        })
        .count();
    assert!(
        with_plane >= 1,
        "at least one face should have a Plane surface"
    );
}

#[test]
fn sheet_quad_rejects_degenerate_extents() {
    // x_min == x_max → zero-width quad.
    assert!(BRep::sheet_quad(0.0, 0.0, 0.0, 1.0, 0.0, tol()).is_err());
    // y_min > y_max → inverted.
    assert!(BRep::sheet_quad(0.0, 1.0, 1.0, 0.5, 0.0, tol()).is_err());
}

#[test]
fn sheet_quad_vertices_are_at_requested_corners() {
    let sheet = BRep::sheet_quad(1.0, 3.0, -1.0, 2.0, 5.0, tol()).unwrap();
    let mut positions: Vec<(f64, f64, f64)> = sheet
        .vertices
        .values()
        .map(|v| (v.position.x, v.position.y, v.position.z))
        .collect();
    positions.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut expected = vec![
        (1.0, -1.0, 5.0),
        (3.0, -1.0, 5.0),
        (3.0, 2.0, 5.0),
        (1.0, 2.0, 5.0),
    ];
    expected.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(positions, expected);
}

// ---------------------------------------------------------------------------
// sheet_triangle
// ---------------------------------------------------------------------------

#[test]
fn sheet_triangle_validates() {
    let sheet = BRep::sheet_triangle(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
        tol(),
    )
    .unwrap();
    assert_eq!(sheet.vertices.len(), 3);
    assert_eq!(sheet.edges.len(), 3);
    assert_eq!(sheet.faces.len(), 2);
    sheet.validate().unwrap();
}

#[test]
fn sheet_triangle_rejects_colinear_points() {
    // Three colinear points.
    let err = BRep::sheet_triangle(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        tol(),
    );
    assert!(err.is_err(), "colinear points should reject");
}

#[test]
fn sheet_triangle_euler_identity() {
    let sheet = BRep::sheet_triangle(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        tol(),
    )
    .unwrap();
    // V=3, E=3, F=2 → 3 − 3 + 2 = 2 ✓
    let euler = sheet.vertices.len() as i64 - sheet.edges.len() as i64 + sheet.faces.len() as i64;
    assert_eq!(euler, 2);
}

// ---------------------------------------------------------------------------
// sheet_polygon
// ---------------------------------------------------------------------------

#[test]
fn sheet_polygon_pentagon_validates() {
    // Regular pentagon in the xy plane.
    use std::f64::consts::PI;
    let n = 5;
    let pts: Vec<Point3> = (0..n)
        .map(|i| {
            let a = 2.0 * PI * (i as f64) / (n as f64);
            Point3::new(a.cos(), a.sin(), 0.0)
        })
        .collect();
    let sheet = BRep::sheet_polygon(&pts, tol()).unwrap();
    assert_eq!(sheet.vertices.len(), 5);
    assert_eq!(sheet.edges.len(), 5);
    assert_eq!(sheet.faces.len(), 2);
    sheet.validate().unwrap();
}

#[test]
fn sheet_polygon_rejects_too_few_vertices() {
    let two: [Point3; 2] = [Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    assert!(BRep::sheet_polygon(&two, tol()).is_err());
    let empty: [Point3; 0] = [];
    assert!(BRep::sheet_polygon(&empty, tol()).is_err());
}

// ---------------------------------------------------------------------------
// sheet_arched_triangle
// ---------------------------------------------------------------------------

#[test]
fn sheet_arched_triangle_has_expected_topology() {
    let sheet = BRep::sheet_arched_triangle(1.0, 0.0, tol()).unwrap();
    // Triangle: 3 vertices, 3 edges, 2 faces (front + backside paper-thin shell).
    assert_eq!(sheet.vertices.len(), 3);
    assert_eq!(sheet.edges.len(), 3);
    assert_eq!(sheet.faces.len(), 2);
    // Euler = V − E + F = 3 − 3 + 2 = 2 ✓
    let euler = sheet.vertices.len() as i64 - sheet.edges.len() as i64 + sheet.faces.len() as i64;
    assert_eq!(euler, 2);
}

#[test]
fn sheet_arched_triangle_validates() {
    let sheet = BRep::sheet_arched_triangle(2.0, 0.5, tol()).unwrap();
    sheet.validate().unwrap();
}

#[test]
fn sheet_arched_triangle_vertices_are_at_requested_positions() {
    let r = 1.5_f64;
    let z = 0.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, z, tol()).unwrap();
    let has_v0 = sheet.vertices.values().any(|v| {
        (v.position.x - (-r)).abs() < 1e-12
            && v.position.y.abs() < 1e-12
            && (v.position.z - z).abs() < 1e-12
    });
    let has_v1 = sheet.vertices.values().any(|v| {
        v.position.x.abs() < 1e-12
            && (v.position.y - r).abs() < 1e-12
            && (v.position.z - z).abs() < 1e-12
    });
    let has_v2 = sheet.vertices.values().any(|v| {
        (v.position.x - r).abs() < 1e-12
            && v.position.y.abs() < 1e-12
            && (v.position.z - z).abs() < 1e-12
    });
    assert!(has_v0 && has_v1 && has_v2);
}

#[test]
fn sheet_arched_triangle_has_two_arc_edges_and_one_line_edge() {
    use crate::geometry::curve::CurveGeometry;
    let sheet = BRep::sheet_arched_triangle(1.0, 0.0, tol()).unwrap();
    let arc_count = sheet
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Arc { .. })
            )
        })
        .count();
    let line_count = sheet
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Line { .. })
            )
        })
        .count();
    assert_eq!(arc_count, 2);
    assert_eq!(line_count, 1);
}

#[test]
fn sheet_arched_triangle_arc_curves_pass_through_midpoint() {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::curve::CurveGeometry;
    let r = 1.0;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();

    // Each arc edge's curve should evaluate to a point on the
    // circle of radius r centred at the origin.
    for (_, edge) in &sheet.edges {
        let curve = edge.curve.as_ref().unwrap();
        if matches!(&curve.geometry, CurveGeometry::Arc { .. }) {
            // Mid-parameter point.
            let t_mid = 0.5 * (curve.t_start + curve.t_end);
            let p = curve.geometry.point_at(t_mid);
            let dist_from_origin = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!(
                (dist_from_origin - r).abs() < 1e-9,
                "arc midpoint not on circle: |{p:?}| = {dist_from_origin}, expected {r}"
            );
        }
    }
}

#[test]
fn sheet_arched_triangle_rejects_bad_radius() {
    assert!(BRep::sheet_arched_triangle(0.0, 0.0, tol()).is_err());
    assert!(BRep::sheet_arched_triangle(-1.0, 0.0, tol()).is_err());
    assert!(BRep::sheet_arched_triangle(f64::NAN, 0.0, tol()).is_err());
    assert!(BRep::sheet_arched_triangle(f64::INFINITY, 0.0, tol()).is_err());
}
