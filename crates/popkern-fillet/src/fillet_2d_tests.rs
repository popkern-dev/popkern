use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Fillet the `(0, 0, 0)` corner of a `sheet_quad` with a 0.1
/// radius. Expect the output to have:
/// - 5 vertices (4 original - 1 corner + 2 tangent)
/// - 5 edges (4 original + 1 arc - 0 deletions)
/// - 2 faces (unchanged sheet topology)
#[test]
fn fillet_corner_2d_on_sheet_quad_expected_counts() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();

    let result = fillet_corner_2d(&sheet, corner, 0.1).expect("2D corner fillet should succeed");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    assert_eq!(result.faces.len(), 2);
}

/// The fillet operation must preserve the Euler identity for the
/// sphere-topology sheet body: `V - E + F = 2`.
#[test]
fn fillet_corner_2d_preserves_euler_identity() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, 0.15).unwrap();
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
}

/// On the `sheet_quad` corner at `(0, 0, 0)` with r = 0.1, the
/// new tangent points should be at `(0.1, 0, 0)` and `(0, 0.1, 0)`
/// (within tolerance), and the arc centre should be at
/// `(0.1, 0.1, 0)`.
#[test]
fn fillet_corner_2d_tangent_points_and_arc_center_match_theory() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, 0.1).unwrap();

    let tol_eps = 1e-9;

    // Expected tangent positions.
    let has_tangent_a = result.vertices.values().any(|v| {
        (v.position.x - 0.1).abs() < tol_eps
            && v.position.y.abs() < tol_eps
            && v.position.z.abs() < tol_eps
    });
    let has_tangent_b = result.vertices.values().any(|v| {
        v.position.x.abs() < tol_eps
            && (v.position.y - 0.1).abs() < tol_eps
            && v.position.z.abs() < tol_eps
    });
    assert!(has_tangent_a, "missing tangent point at (0.1, 0, 0)");
    assert!(has_tangent_b, "missing tangent point at (0, 0.1, 0)");

    // The arc edge has a Circle/Arc curve with centre (0.1, 0.1, 0)
    // and radius 0.1.
    let arc_edges: Vec<_> = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Arc { .. })
            )
        })
        .collect();
    assert_eq!(arc_edges.len(), 1, "expected exactly 1 arc edge");
    if let CurveGeometry::Arc { center, radius, .. } =
        &arc_edges[0].curve.as_ref().unwrap().geometry
    {
        assert!((center.x - 0.1).abs() < tol_eps);
        assert!((center.y - 0.1).abs() < tol_eps);
        assert!(center.z.abs() < tol_eps);
        assert!((radius - 0.1).abs() < tol_eps);
    }
}

/// Fillet a corner of a triangular sheet body.
#[test]
fn fillet_corner_2d_on_sheet_triangle() {
    let sheet = BRep::sheet_triangle(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 2.0, 0.0),
        tol(),
    )
    .unwrap();
    let corner = sheet.vertices.keys().next().unwrap(); // pick any
    let result = fillet_corner_2d(&sheet, corner, 0.2).expect("triangle corner fillet should work");
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
}

/// Negative, zero, NaN, infinite radii are rejected at the API
/// boundary.
#[test]
fn fillet_corner_2d_rejects_invalid_radii() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();

    for bad in [-0.1_f64, 0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let err = fillet_corner_2d(&sheet, corner, bad).unwrap_err();
        assert!(
            matches!(err, BlendError::InvalidParameter { ref parameter, .. } if parameter == "radius"),
            "radius={bad} should produce InvalidParameter, got {err:?}"
        );
    }
}

/// A radius so large that the tangent setback exceeds the
/// adjacent edge length is rejected with `InvalidTopology`
/// (wrapped in `PartialFailure`).
#[test]
fn fillet_corner_2d_rejects_radius_too_large() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    // radius 2.0 on a unit-width quad has setback 2.0 which is
    // larger than the adjacent edge length of 1.0.
    let err = fillet_corner_2d(&sheet, corner, 2.0).unwrap_err();
    match err {
        BlendError::PartialFailure { failures, .. } => {
            let (_, f) = &failures[0];
            match f {
                BlendFailure::InvalidTopology { reason, .. } => {
                    assert!(reason.contains("too large"), "got reason: {reason}");
                }
                other => panic!("expected InvalidTopology, got {other:?}"),
            }
        }
        other => panic!("expected PartialFailure, got {other:?}"),
    }
}

/// Invalid vertex handle returns `InvalidParameter`.
#[test]
fn fillet_corner_2d_rejects_invalid_vertex() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let err = fillet_corner_2d(&sheet, VertexHandle::default(), 0.1).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

/// The result must still pass `validate` (manifold edges, loop
/// closure, watertight geometry, etc.). This test is redundant
/// with `fillet_corner_2d`'s internal validate call, but pins the
/// behaviour as a regression guard.
#[test]
fn fillet_corner_2d_result_validates() {
    let sheet = BRep::sheet_quad(0.0, 2.0, 0.0, 2.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, 0.25).unwrap();
    result.validate().unwrap();
}

// ---------------------------------------------------------------------------
// Batch (§5.7 + §5.8)
// ---------------------------------------------------------------------------

/// Fillet all 4 corners of a unit `sheet_quad` with radius 0.1
/// in a single batch call. Expect 8 vertices, 8 edges, 2 faces.
#[test]
fn fillet_corners_2d_all_four_sheet_quad_corners() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let specs: Vec<_> = sheet.vertices.keys().map(|v| (v, 0.1_f64)).collect();
    assert_eq!(specs.len(), 4);

    let result = fillet_corners_2d(&sheet, &specs).expect("all-4 fillet should succeed");
    // Each corner: +1 vertex (2 new tangent - 1 deleted = +1), +1 edge (arc).
    // 4 × (+1) = 4 new vertices + 4 new edges.
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);
    assert_eq!(result.faces.len(), 2);
    // Euler identity.
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// Fillet all 3 corners of a triangle in one batch call.
#[test]
fn fillet_corners_2d_all_three_triangle_corners() {
    let sheet = BRep::sheet_triangle(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
        Point3::new(1.5, 3.0, 0.0),
        tol(),
    )
    .unwrap();
    let specs: Vec<_> = sheet.vertices.keys().map(|v| (v, 0.2)).collect();
    let result = fillet_corners_2d(&sheet, &specs).unwrap();
    // 3 corners × (+1 vertex +1 edge) = 6 vertices, 6 edges.
    assert_eq!(result.vertices.len(), 6);
    assert_eq!(result.edges.len(), 6);
    assert_eq!(result.faces.len(), 2);
    result.validate().unwrap();
}

/// Empty spec list returns the input unchanged.
#[test]
fn fillet_corners_2d_empty_specs_is_identity() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let result = fillet_corners_2d(&sheet, &[]).unwrap();
    assert_eq!(result.vertices.len(), sheet.vertices.len());
    assert_eq!(result.edges.len(), sheet.edges.len());
    assert_eq!(result.faces.len(), sheet.faces.len());
}

/// Collision detection: two adjacent corner fillets whose
/// combined setbacks exceed the shared edge length should fail
/// on the second call with a "setback too large" error.
///
/// Sheet quad with edge length 1.0. Filleting the first corner
/// with radius 0.6 leaves the shared edge at length 0.4.
/// Filleting the adjacent corner with radius 0.6 would need a
/// 0.6 setback on an edge of length 0.4 — rejected.
#[test]
fn fillet_corners_2d_collision_rejects_second_corner() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    // Pick two vertices that share an edge. On a unit_quad both
    // (0,0) and (1,0) share the edge y=0.
    let v0 = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let v1 = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 1.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();

    // Collision: 0.6 + 0.6 > 1.0 shared edge length.
    let err = fillet_corners_2d(&sheet, &[(v0, 0.6), (v1, 0.6)]).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

/// Each batch operation must be independently valid. Passing an
/// invalid radius in the middle of a batch surfaces the error
/// without a partial result leaking out.
#[test]
fn fillet_corners_2d_aborts_on_invalid_spec() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    // First spec valid, second spec has a NaN radius.
    let err = fillet_corners_2d(&sheet, &[(corner, 0.1), (corner, f64::NAN)]).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

/// The result's front face still has a planar surface (we didn't
/// accidentally lose it during the surgery).
#[test]
fn fillet_corner_2d_front_face_plane_is_preserved() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, 0.1).unwrap();
    assert!(result.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
    }));
}
