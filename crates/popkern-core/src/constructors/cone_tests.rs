//! Tests for the cone constructor.

use super::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

#[test]
fn cone_topology() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    assert_eq!(cone.vertices.len(), 2);
    assert_eq!(cone.edges.len(), 2);
    assert_eq!(cone.faces.len(), 2);
}

#[test]
fn cone_passes_validation() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    assert!(cone.validate().is_ok(), "Cone should pass validation");
}

#[test]
fn cone_has_edge_geometry() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    for edge in cone.edges.values() {
        assert!(
            edge.curve.is_some(),
            "Every edge should have curve geometry"
        );
    }
}

#[test]
fn cone_vertex_positions() {
    let cone = BRep::cone(1.0, 3.0, Tolerance::default()).unwrap();
    let positions: Vec<_> = cone.vertices.values().map(|v| v.position).collect();
    let has_base = positions
        .iter()
        .any(|p| (p.x - 1.0).abs() < 1e-10 && p.y.abs() < 1e-10 && p.z.abs() < 1e-10);
    let has_apex = positions
        .iter()
        .any(|p| (p.z - 3.0).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
    assert!(has_base, "Should have base vertex at (radius, 0, 0)");
    assert!(has_apex, "Should have apex vertex at (0, 0, height)");
}

#[test]
fn cone_invalid_params() {
    assert!(BRep::cone(0.0, 1.0, Tolerance::default()).is_err());
    assert!(BRep::cone(-1.0, 1.0, Tolerance::default()).is_err());
    assert!(BRep::cone(1.0, 0.0, Tolerance::default()).is_err());
    assert!(BRep::cone(1.0, -1.0, Tolerance::default()).is_err());
}

#[test]
fn cone_all_faces_have_surface_geometry() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    for face in cone.faces.values() {
        assert!(
            face.surface.is_some(),
            "every cone face should have surface geometry"
        );
    }
}

#[test]
fn cone_has_surface_types() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    let expected_half_angle = (1.0_f64 / 2.0_f64).atan();
    let mut cone_faces = 0;
    let mut plane_faces = 0;
    for face in cone.faces.values() {
        let surface = face
            .surface
            .as_ref()
            .expect("cone face should have surface geometry");
        match &surface.geometry {
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } => {
                assert!(
                    (apex.z - 2.0).abs() < 1e-10,
                    "cone apex z should be at height"
                );
                assert!(
                    apex.x.abs() < 1e-10 && apex.y.abs() < 1e-10,
                    "cone apex should be on z-axis"
                );
                assert!(
                    (axis.z + 1.0).abs() < 1e-10 && axis.x.abs() < 1e-10 && axis.y.abs() < 1e-10,
                    "cone axis should point in -Z direction"
                );
                assert!(
                    (*half_angle - expected_half_angle).abs() < 1e-10,
                    "cone half-angle should match atan(radius/height)"
                );
                cone_faces += 1;
            }
            SurfaceGeometry::Plane { .. } => {
                plane_faces += 1;
            }
            other => panic!(
                "expected Cone or Plane, got {:?}",
                std::mem::discriminant(other)
            ),
        }
    }
    assert_eq!(cone_faces, 1, "Should have 1 conical face");
    assert_eq!(plane_faces, 1, "Should have 1 planar cap face");
}

#[test]
fn cone_surface_geometry_matches_params() {
    let radius = 3.0_f64;
    let height = 4.0_f64;
    let cone = BRep::cone(radius, height, Tolerance::default()).unwrap();
    let expected_half_angle = (radius / height).atan();
    for face in cone.faces.values() {
        if let Some(ref sof) = face.surface
            && let SurfaceGeometry::Cone { half_angle, .. } = sof.geometry
        {
            assert!(
                (half_angle - expected_half_angle).abs() < 1e-10,
                "half_angle mismatch: got {half_angle}, expected {expected_half_angle}"
            );
        }
    }
}

#[test]
fn cone_seam_edge_both_half_edges_on_lateral_face() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    let mut found_seam = false;
    for edge in cone.edges.values() {
        let he = &cone.half_edges[edge.half_edge];
        let twin = &cone.half_edges[he.twin];
        let face_he = cone.loops[he.loop_].face;
        let face_twin = cone.loops[twin.loop_].face;
        if face_he == face_twin {
            found_seam = true;
        }
    }
    assert!(
        found_seam,
        "Should have a seam edge with both halves on the lateral face"
    );
}

#[test]
fn cone_cap_loop_is_self_loop() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    let mut single_he_loops = 0;
    for loop_ in cone.loops.values() {
        let he = &cone.half_edges[loop_.half_edge];
        if he.next == loop_.half_edge {
            single_he_loops += 1;
        }
    }
    assert_eq!(
        single_he_loops, 1,
        "Should have 1 single-half-edge cap loop"
    );
}
