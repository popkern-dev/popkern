use super::*;
use crate::chamfer_2d_by_angle::chamfer_corner_2d_by_angle;
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
// Distance + angle mode (§2.3 / §5.7)
// ---------------------------------------------------------------------------

/// 45-45 symmetric chamfer via distance+angle: on a 90° corner
/// with dist_a = 0.1 and angle = 135° (3π/4), the result should
/// be identical to the symmetric distance-based chamfer
/// `chamfer_corner_2d(b, v, 0.1, 0.1)`.
#[test]
fn chamfer_by_angle_45_45_matches_symmetric_distance_chamfer() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();

    let by_distance =
        chamfer_corner_2d(&sheet, corner, 0.1, 0.1).expect("distance chamfer should succeed");
    let by_angle =
        chamfer_corner_2d_by_angle(&sheet, corner, 0.1, 3.0 * std::f64::consts::PI / 4.0)
            .expect("angle chamfer should succeed");

    // Same topology counts.
    assert_eq!(by_angle.vertices.len(), by_distance.vertices.len());
    assert_eq!(by_angle.edges.len(), by_distance.edges.len());
    assert_eq!(by_angle.faces.len(), by_distance.faces.len());

    // Same tangent-point positions — both methods should produce
    // vertices at (0.1, 0, 0) and (0, 0.1, 0).
    for (_, v) in &by_angle.vertices {
        let eq_to_t1 = (v.position.x - 0.1).abs() < 1e-9
            && v.position.y.abs() < 1e-9
            && v.position.z.abs() < 1e-9;
        let eq_to_t2 = v.position.x.abs() < 1e-9
            && (v.position.y - 0.1).abs() < 1e-9
            && v.position.z.abs() < 1e-9;
        let eq_to_orig_v1 = (v.position.x - 1.0).abs() < 1e-9
            && v.position.y.abs() < 1e-9
            && v.position.z.abs() < 1e-9;
        let eq_to_orig_v2 = (v.position.x - 1.0).abs() < 1e-9
            && (v.position.y - 1.0).abs() < 1e-9
            && v.position.z.abs() < 1e-9;
        let eq_to_orig_v3 = v.position.x.abs() < 1e-9
            && (v.position.y - 1.0).abs() < 1e-9
            && v.position.z.abs() < 1e-9;
        assert!(
            eq_to_t1 || eq_to_t2 || eq_to_orig_v1 || eq_to_orig_v2 || eq_to_orig_v3,
            "unexpected vertex position {:?}",
            v.position
        );
    }
}

/// Asymmetric distance+angle: 0.1 along the first edge at 120°
/// yields a non-equal second distance via law of sines.
#[test]
fn chamfer_by_angle_asymmetric_uses_law_of_sines() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();

    let dist_a = 0.2;
    let angle_rad = 2.0 * std::f64::consts::PI / 3.0; // 120°
    // Corner angle on a square is π/2 = 90°. With angle_rad = 120°
    // the expected dist_b via law of sines is:
    // dist_b = 0.2 · sin(120°) / sin(120° − 90°)
    //        = 0.2 · (√3/2) / (1/2)
    //        = 0.2 · √3
    //        ≈ 0.3464
    let expected_dist_b = 0.2 * 3.0_f64.sqrt();

    let result = chamfer_corner_2d_by_angle(&sheet, corner, dist_a, angle_rad).unwrap();
    // Verify a vertex at dist_b along the y-axis or x-axis (slot
    // map ordering may pick either as "dist_a side").
    let expected_found = result.vertices.values().any(|v| {
        // Either (dist_a, 0, 0) + (0, expected_dist_b, 0) matches
        // the (+x, +y) pair, or the permuted version.
        ((v.position.x - expected_dist_b).abs() < 1e-6 && v.position.y.abs() < 1e-6)
            || (v.position.x.abs() < 1e-6 && (v.position.y - expected_dist_b).abs() < 1e-6)
    });
    assert!(
        expected_found,
        "expected a vertex at distance {expected_dist_b} along one axis; got positions: {:?}",
        result
            .vertices
            .values()
            .map(|v| v.position)
            .collect::<Vec<_>>()
    );
    result.validate().unwrap();
}

/// `angle_rad ≤ corner_angle` is rejected — the chamfer line
/// would be parallel to (or away from) the second edge.
#[test]
fn chamfer_by_angle_rejects_angle_at_or_below_corner_angle() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    // Corner angle is π/2; angle_rad = π/2 should be rejected.
    let err =
        chamfer_corner_2d_by_angle(&sheet, corner, 0.1, std::f64::consts::FRAC_PI_2).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

/// `angle_rad ≥ π` is rejected.
#[test]
fn chamfer_by_angle_rejects_angle_at_or_above_pi() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    let err = chamfer_corner_2d_by_angle(&sheet, corner, 0.1, std::f64::consts::PI).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

/// Negative / zero / NaN / infinite distances and angles are
/// rejected at the API boundary.
#[test]
fn chamfer_by_angle_rejects_invalid_inputs() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();

    for bad in [-0.1_f64, 0.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            chamfer_corner_2d_by_angle(&sheet, corner, bad, 2.0).unwrap_err(),
            BlendError::InvalidParameter { .. }
        ));
        assert!(matches!(
            chamfer_corner_2d_by_angle(&sheet, corner, 0.1, bad).unwrap_err(),
            BlendError::InvalidParameter { .. }
        ));
    }
}

/// Derived `dist_b` exceeding the adjacent edge length is
/// rejected with `InvalidTopology`.
#[test]
fn chamfer_by_angle_rejects_derived_dist_b_too_large() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    // dist_a=0.9, angle_rad=100° → dist_b = 0.9 · sin(100°)/sin(10°) ≈ 5.10,
    // well over the 1.0 edge length.
    let err = chamfer_corner_2d_by_angle(&sheet, corner, 0.9, 100.0_f64.to_radians()).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

// ---------------------------------------------------------------------------
// Curved boundary edges — see `chamfer_2d_curved_tests.rs`
// ---------------------------------------------------------------------------

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
