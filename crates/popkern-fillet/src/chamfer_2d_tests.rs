use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Symmetric 45° chamfer on a unit `sheet_quad` corner.
#[test]
fn chamfer_corner_2d_symmetric_on_sheet_quad() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.1, 0.1).expect("chamfer should succeed");

    // 5 vertices (3 original + 2 tangent), 5 edges (4 original + 1 new), 2 faces.
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    assert_eq!(result.faces.len(), 2);
    // Euler = 2 for the sphere-topology sheet.
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
}

/// Asymmetric chamfer: 0.1 along one edge, 0.2 along the other.
/// The new edge should be a straight line between the two
/// distance-determined tangent points.
#[test]
fn chamfer_corner_2d_asymmetric_tangent_points_match_distances() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.1, 0.2).unwrap();

    let tol_eps = 1e-9;

    // One tangent at distance 0.1 along the +x edge and another
    // at distance 0.2 along the +y edge. The ordering of
    // `dist_a` and `dist_b` is slot-map dependent, so both
    // permutations are acceptable.
    let has_01_along_x = result.vertices.values().any(|v| {
        (v.position.x - 0.1).abs() < tol_eps
            && v.position.y.abs() < tol_eps
            && v.position.z.abs() < tol_eps
    });
    let has_02_along_y = result.vertices.values().any(|v| {
        v.position.x.abs() < tol_eps
            && (v.position.y - 0.2).abs() < tol_eps
            && v.position.z.abs() < tol_eps
    });
    let has_02_along_x = result.vertices.values().any(|v| {
        (v.position.x - 0.2).abs() < tol_eps
            && v.position.y.abs() < tol_eps
            && v.position.z.abs() < tol_eps
    });
    let has_01_along_y = result.vertices.values().any(|v| {
        v.position.x.abs() < tol_eps
            && (v.position.y - 0.1).abs() < tol_eps
            && v.position.z.abs() < tol_eps
    });

    let permutation_1 = has_01_along_x && has_02_along_y;
    let permutation_2 = has_02_along_x && has_01_along_y;
    assert!(
        permutation_1 || permutation_2,
        "expected tangent points at distances 0.1 and 0.2 along the two axes"
    );
}

/// The new chamfer edge carries a `Line` curve between the two
/// tangent points, not an `Arc`.
#[test]
fn chamfer_corner_2d_new_edge_is_a_line_curve() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.1, 0.1).unwrap();

    let line_edges = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Line { .. })
            )
        })
        .count();
    // 5 total edges, all of which are `Line` (3 original unchanged,
    // 2 original rebuilt after endpoint moves, 1 new chamfer).
    assert_eq!(line_edges, 5);
    // Zero arc edges.
    let arc_edges = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Arc { .. })
            )
        })
        .count();
    assert_eq!(arc_edges, 0);
}

/// Chamfer on a triangular sheet body with asymmetric distances.
#[test]
fn chamfer_corner_2d_on_sheet_triangle() {
    let sheet = BRep::sheet_triangle(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
        Point3::new(1.5, 3.0, 0.0),
        tol(),
    )
    .unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.3, 0.4).unwrap();
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    result.validate().unwrap();
}

/// Negative, zero, NaN, infinite distances are rejected at the
/// API boundary.
#[test]
fn chamfer_corner_2d_rejects_invalid_distances() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();

    for bad in [-0.1_f64, 0.0, f64::NAN, f64::INFINITY] {
        let err = chamfer_corner_2d(&sheet, corner, bad, 0.1).unwrap_err();
        assert!(
            matches!(err, BlendError::InvalidParameter { ref parameter, .. } if parameter == "dist_a"),
            "dist_a={bad} should produce InvalidParameter, got {err:?}"
        );
        let err = chamfer_corner_2d(&sheet, corner, 0.1, bad).unwrap_err();
        assert!(
            matches!(err, BlendError::InvalidParameter { ref parameter, .. } if parameter == "dist_b"),
            "dist_b={bad} should produce InvalidParameter, got {err:?}"
        );
    }
}

/// Chamfer distances larger than the adjacent edge length are
/// rejected with `InvalidTopology`.
#[test]
fn chamfer_corner_2d_rejects_oversized_distance() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    // Edge length is 1.0; distance 1.5 exceeds.
    let err = chamfer_corner_2d(&sheet, corner, 1.5, 0.1).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

/// Invalid vertex handle returns `InvalidParameter`.
#[test]
fn chamfer_corner_2d_rejects_invalid_vertex() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let err = chamfer_corner_2d(&sheet, VertexHandle::default(), 0.1, 0.1).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

// ---------------------------------------------------------------------------
// Batch (§5.7 + §5.8)
// ---------------------------------------------------------------------------

/// Chamfer all 4 corners of a `sheet_quad` with equal 0.1×0.1
/// distances in a single batch call.
#[test]
fn chamfer_corners_2d_all_four_sheet_quad_corners() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let specs: Vec<_> = sheet.vertices.keys().map(|v| (v, 0.1, 0.1)).collect();
    assert_eq!(specs.len(), 4);

    let result = chamfer_corners_2d(&sheet, &specs).expect("all-4 chamfer should succeed");
    // Each corner: +1 vertex, +1 edge.
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);
    assert_eq!(result.faces.len(), 2);
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// Empty spec list returns the input unchanged.
#[test]
fn chamfer_corners_2d_empty_specs_is_identity() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let result = chamfer_corners_2d(&sheet, &[]).unwrap();
    assert_eq!(result.vertices.len(), 4);
}

/// Collision: two adjacent chamfers whose combined distances
/// exceed the shared edge length fail on the second call.
#[test]
fn chamfer_corners_2d_collision_rejects_second_corner() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
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
    // 0.7 + 0.7 > 1.0 shared edge length.
    let err = chamfer_corners_2d(&sheet, &[(v0, 0.7, 0.1), (v1, 0.7, 0.1)]).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

/// Mixed batch: chamfer a corner with asymmetric distances
/// followed by another with symmetric distances.
#[test]
fn chamfer_corners_2d_mixed_symmetric_asymmetric() {
    let sheet = BRep::sheet_quad(0.0, 2.0, 0.0, 2.0, 0.0, tol()).unwrap();
    let mut iter = sheet.vertices.iter().map(|(h, v)| (h, v.position));
    let (v0, _) = iter.next().unwrap();
    let (v1, _) = iter.next().unwrap();
    let result =
        chamfer_corners_2d(&sheet, &[(v0, 0.3, 0.5), (v1, 0.2, 0.2)]).expect("mixed batch works");
    // 2 corners × (+1 vertex +1 edge) = 6 vertices, 6 edges.
    assert_eq!(result.vertices.len(), 6);
    assert_eq!(result.edges.len(), 6);
    result.validate().unwrap();
}

/// Chamfer result passes full BRep validation.
#[test]
fn chamfer_corner_2d_result_validates() {
    let sheet = BRep::sheet_quad(0.0, 2.0, 0.0, 2.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.25, 0.5).unwrap();
    result.validate().unwrap();
    // Plane surface preserved.
    assert!(result.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
    }));
}
