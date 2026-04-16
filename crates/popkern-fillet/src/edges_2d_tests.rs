//! Tests for the edge-sequence 2D corner ops
//! ([`fillet_edges_2d`], [`chamfer_edges_2d`], and their
//! variable-profile siblings).

use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Walk the outer loop of the planar front face of a sheet body
/// and return the first `count` edges in loop order. Shared by
/// several tests below that need a connected-chain input.
fn front_face_first_edges(sheet: &BRep, count: usize) -> Vec<EdgeHandle> {
    let front_face = sheet
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
        })
        .map(|(h, _)| h)
        .unwrap();
    let loop_h = sheet.faces[front_face].outer_loop;
    let start_he = sheet.loops[loop_h].half_edge;
    let mut out = Vec::with_capacity(count);
    let mut current = start_he;
    for _ in 0..count {
        out.push(sheet.half_edges[current].edge);
        current = sheet.half_edges[current].next;
    }
    out
}

// ---------------------------------------------------------------------------
// Constant-radius edge-sequence fillet / chamfer
// ---------------------------------------------------------------------------

#[test]
fn fillet_edges_2d_all_edges_of_quad_fillets_every_corner() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let all_edges: Vec<_> = sheet.edges.keys().collect();
    assert_eq!(all_edges.len(), 4);

    let result = fillet_edges_2d(&sheet, &all_edges, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);
    assert_eq!(result.faces.len(), 2);
    result.validate().unwrap();
}

#[test]
fn fillet_edges_2d_three_edge_chain_fillets_two_interior_corners() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let chain = front_face_first_edges(&sheet, 3);
    let result = fillet_edges_2d(&sheet, &chain, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 6);
    assert_eq!(result.edges.len(), 6);
    result.validate().unwrap();
}

#[test]
fn fillet_edges_2d_single_edge_is_no_op() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let edge: Vec<_> = sheet.edges.keys().take(1).collect();
    let result = fillet_edges_2d(&sheet, &edge, 0.1).unwrap();
    assert_eq!(result.vertices.len(), sheet.vertices.len());
    assert_eq!(result.edges.len(), sheet.edges.len());
}

#[test]
fn fillet_edges_2d_empty_is_no_op() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let result = fillet_edges_2d(&sheet, &[], 0.1).unwrap();
    assert_eq!(result.vertices.len(), sheet.vertices.len());
    assert_eq!(result.edges.len(), sheet.edges.len());
}

#[test]
fn fillet_edges_2d_disconnected_pair_is_no_op() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let mut bottom = None;
    let mut top = None;
    for (eh, edge) in &sheet.edges {
        let he = edge.half_edge;
        let a = sheet.vertices[sheet.half_edges[he].origin].position;
        let b = sheet.vertices[sheet.half_edges[sheet.half_edges[he].twin].origin].position;
        if a.y == 0.0 && b.y == 0.0 {
            bottom = Some(eh);
        } else if a.y == 1.0 && b.y == 1.0 {
            top = Some(eh);
        }
    }
    let edges = [bottom.unwrap(), top.unwrap()];
    let result = fillet_edges_2d(&sheet, &edges, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
}

#[test]
fn chamfer_edges_2d_all_edges_of_quad_chamfers_every_corner() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let all_edges: Vec<_> = sheet.edges.keys().collect();
    let result = chamfer_edges_2d(&sheet, &all_edges, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);
    result.validate().unwrap();
}

#[test]
fn chamfer_edges_2d_single_edge_is_no_op() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let edge: Vec<_> = sheet.edges.keys().take(1).collect();
    let result = chamfer_edges_2d(&sheet, &edge, 0.1).unwrap();
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
}

// ---------------------------------------------------------------------------
// Variable-parameter edge-sequence fillet / chamfer
// ---------------------------------------------------------------------------

/// Linear profile on a 3-edge chain of a quad gives exactly two
/// interior corners with distinct radii — 0.05 at the first,
/// 0.15 at the second.
#[test]
fn fillet_edges_2d_variable_linear_profile_on_quad_chain() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let chain = front_face_first_edges(&sheet, 3);
    let result = fillet_edges_2d_variable(&sheet, &chain, 0.05, 0.15).unwrap();
    assert_eq!(result.vertices.len(), 6);
    assert_eq!(result.edges.len(), 6);
    result.validate().unwrap();

    // Each fillet arc edge picks up its per-corner radius. Two
    // interior corners ⇒ two `Arc` curves in the result with
    // distinct radii (0.05 and 0.15).
    let mut radii: Vec<f64> = result
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(popkern_core::geometry::curve::CurveGeometry::Arc { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(radii.len(), 2);
    assert!((radii[0] - 0.05).abs() < 1e-9, "first radius {}", radii[0]);
    assert!((radii[1] - 0.15).abs() < 1e-9, "last radius {}", radii[1]);
}

/// Constant profile (`r_start == r_end`) must give the same
/// result as the constant-radius path.
#[test]
fn fillet_edges_2d_variable_constant_equals_constant_api() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let all_edges: Vec<_> = sheet.edges.keys().collect();
    let r = 0.1_f64;
    let result_constant = fillet_edges_2d(&sheet, &all_edges, r).unwrap();
    let result_variable = fillet_edges_2d_variable(&sheet, &all_edges, r, r).unwrap();
    assert_eq!(
        result_constant.vertices.len(),
        result_variable.vertices.len()
    );
    assert_eq!(result_constant.edges.len(), result_variable.edges.len());
    assert_eq!(result_constant.faces.len(), result_variable.faces.len());
    result_variable.validate().unwrap();
}

/// Closed loop of 4 edges walks from the smallest-key edge and
/// interpolates across 4 corners (the walk visits every corner).
#[test]
fn fillet_edges_2d_variable_closed_quad_loop() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let all_edges: Vec<_> = sheet.edges.keys().collect();
    let result = fillet_edges_2d_variable(&sheet, &all_edges, 0.05, 0.2).unwrap();
    // 4 corners filleted regardless of interpolation.
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);
    result.validate().unwrap();

    // 4 distinct radii expected from the linear profile:
    // t ∈ {0, 1/3, 2/3, 1} ⇒ r ∈ {0.05, 0.10, 0.15, 0.20}.
    let mut radii: Vec<f64> = result
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(popkern_core::geometry::curve::CurveGeometry::Arc { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(radii.len(), 4);
    let expected = [0.05, 0.10, 0.15, 0.20];
    for (got, want) in radii.iter().zip(expected.iter()) {
        assert!(
            (got - want).abs() < 1e-9,
            "radius mismatch: got {got}, want {want}"
        );
    }
}

/// Non-finite start radius is rejected at the API boundary.
#[test]
fn fillet_edges_2d_variable_rejects_invalid_radii() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let all_edges: Vec<_> = sheet.edges.keys().collect();
    for bad in [-0.1, 0.0, f64::NAN, f64::INFINITY] {
        let err = fillet_edges_2d_variable(&sheet, &all_edges, bad, 0.1).unwrap_err();
        assert!(
            matches!(err, BlendError::InvalidParameter { parameter, .. } if parameter == "r_start")
        );
    }
}

/// Disconnected edge set is rejected with a chain-walker error.
#[test]
fn fillet_edges_2d_variable_rejects_disconnected_set() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let mut bottom = None;
    let mut top = None;
    for (eh, edge) in &sheet.edges {
        let he = edge.half_edge;
        let a = sheet.vertices[sheet.half_edges[he].origin].position;
        let b = sheet.vertices[sheet.half_edges[sheet.half_edges[he].twin].origin].position;
        if a.y == 0.0 && b.y == 0.0 {
            bottom = Some(eh);
        } else if a.y == 1.0 && b.y == 1.0 {
            top = Some(eh);
        }
    }
    let edges = [bottom.unwrap(), top.unwrap()];
    let err = fillet_edges_2d_variable(&sheet, &edges, 0.05, 0.1).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

/// The chamfer variable path walks the chain and interpolates
/// distance between `d_start` and `d_end`; verified end-to-end
/// by running it on the same quad chain and checking counts +
/// validate.
#[test]
fn chamfer_edges_2d_variable_linear_profile_on_quad_chain() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let chain = front_face_first_edges(&sheet, 3);
    let result = chamfer_edges_2d_variable(&sheet, &chain, 0.05, 0.15).unwrap();
    assert_eq!(result.vertices.len(), 6);
    assert_eq!(result.edges.len(), 6);
    result.validate().unwrap();
}

// ---------------------------------------------------------------------------
// Face-based wrappers (closed-loop walker)
// ---------------------------------------------------------------------------

fn front_planar_face(sheet: &BRep) -> popkern_core::topology::FaceHandle {
    sheet
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
        })
        .map(|(h, _)| h)
        .unwrap()
}

/// `fillet_face_corners_2d` on a sheet_quad front face should
/// match the topology counts of `fillet_edges_2d` with all 4
/// edges — it's the same operation, just with auto-derived
/// edges.
#[test]
fn fillet_face_corners_2d_equivalent_to_fillet_edges_2d() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let face = front_planar_face(&sheet);
    let result_face = fillet_face_corners_2d(&sheet, face, 0.1).unwrap();
    let all_edges: Vec<_> = sheet.edges.keys().collect();
    let result_edges = fillet_edges_2d(&sheet, &all_edges, 0.1).unwrap();
    assert_eq!(result_face.vertices.len(), result_edges.vertices.len());
    assert_eq!(result_face.edges.len(), result_edges.edges.len());
    assert_eq!(result_face.faces.len(), result_edges.faces.len());
    result_face.validate().unwrap();
}

/// Front and back faces of a sheet body must give identical
/// fillet counts — the chain walker normalises the walk start,
/// so the output doesn't depend on which face was passed.
#[test]
fn fillet_face_corners_2d_front_and_back_faces_give_same_result() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let faces: Vec<_> = sheet.faces.keys().collect();
    assert_eq!(faces.len(), 2);
    let r_f0 = fillet_face_corners_2d(&sheet, faces[0], 0.1).unwrap();
    let r_f1 = fillet_face_corners_2d(&sheet, faces[1], 0.1).unwrap();
    assert_eq!(r_f0.vertices.len(), r_f1.vertices.len());
    assert_eq!(r_f0.edges.len(), r_f1.edges.len());
}

/// Variable-radius face-loop fillet produces the 4 distinct
/// radii expected from the linear profile (same as the
/// explicit-edge-list version).
#[test]
fn fillet_face_corners_2d_variable_linear_profile_on_quad() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let face = front_planar_face(&sheet);
    let result = fillet_face_corners_2d_variable(&sheet, face, 0.05, 0.2).unwrap();
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);

    let mut radii: Vec<f64> = result
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(popkern_core::geometry::curve::CurveGeometry::Arc { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let expected = [0.05, 0.10, 0.15, 0.20];
    for (got, want) in radii.iter().zip(expected.iter()) {
        assert!(
            (got - want).abs() < 1e-9,
            "radius mismatch: got {got}, want {want}"
        );
    }
}

/// `chamfer_face_corners_2d` on the front face of a triangle
/// chamfers all three corners.
#[test]
fn chamfer_face_corners_2d_on_sheet_triangle() {
    let sheet = BRep::sheet_triangle(
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Point3::new(3.0, 0.0, 0.0),
        popkern_core::geometry::Point3::new(1.5, 3.0, 0.0),
        tol(),
    )
    .unwrap();
    let face = front_planar_face(&sheet);
    let result = chamfer_face_corners_2d(&sheet, face, 0.2).unwrap();
    // 3 corners chamfered: +3 vertices, +3 line edges.
    assert_eq!(result.vertices.len(), 6);
    assert_eq!(result.edges.len(), 6);
    result.validate().unwrap();
}

/// Variable-distance chamfer on a face loop. Smoke-test end-to-end.
#[test]
fn chamfer_face_corners_2d_variable_on_quad() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let face = front_planar_face(&sheet);
    let result = chamfer_face_corners_2d_variable(&sheet, face, 0.05, 0.15).unwrap();
    assert_eq!(result.vertices.len(), 8);
    assert_eq!(result.edges.len(), 8);
    result.validate().unwrap();
}

/// Invalid face handle is rejected as `InvalidParameter`.
#[test]
fn fillet_face_corners_2d_rejects_invalid_handle() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let err = fillet_face_corners_2d(&sheet, popkern_core::topology::FaceHandle::default(), 0.1)
        .unwrap_err();
    assert!(
        matches!(err, BlendError::InvalidParameter { ref parameter, .. } if parameter == "face")
    );
}

/// Error propagation: oversized radius on a unit-quad face
/// surfaces as `PartialFailure`/`InvalidTopology` from the
/// underlying `fillet_corner_2d`.
#[test]
fn fillet_face_corners_2d_propagates_oversize_error() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let face = front_planar_face(&sheet);
    let err = fillet_face_corners_2d(&sheet, face, 2.0).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}
