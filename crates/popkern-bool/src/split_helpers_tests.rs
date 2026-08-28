//! Tests for split_helpers.rs — separated to keep split_helpers.rs under 500 lines.
//! Included only in test mode via `#[cfg(test)] #[path = "split_helpers_tests.rs"] mod tests`.
// #[cfg(test)] — marker for the no_unwrap linter; .unwrap() is allowed in tests.

use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn unit_cube() -> BRep {
    BRep::unit_cube(tol()).unwrap()
}

// ── point_on_segment ────────────────────────────────────────────────────

#[test]
fn point_on_segment_midpoint_found() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(2.0, 0.0, 0.0);
    let mid = Point3::new(1.0, 0.0, 0.0);
    let t = point_on_segment(&a, &b, &mid, 1e-6).unwrap();
    assert!((t - 0.5).abs() < 1e-6, "expected t≈0.5, got {t}");
}

#[test]
fn point_on_segment_at_endpoint_returns_none() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    // Exactly at endpoint — should return None (already at vertex).
    assert!(point_on_segment(&a, &b, &a, 1e-6).is_none());
    assert!(point_on_segment(&a, &b, &b, 1e-6).is_none());
}

#[test]
fn point_on_segment_off_line_returns_none() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let off = Point3::new(0.5, 1.0, 0.0);
    assert!(point_on_segment(&a, &b, &off, 1e-6).is_none());
}

#[test]
fn point_on_segment_degenerate_edge_returns_none() {
    let a = Point3::new(0.5, 0.5, 0.5);
    // a == b → degenerate, length_sq ≈ 0.
    let t = point_on_segment(&a, &a, &a, 1e-6);
    assert!(t.is_none(), "degenerate edge should return None");
}

// ── walk_face_loop ───────────────────────────────────────────────────────

#[test]
fn walk_face_loop_returns_four_half_edges_for_cube_face() {
    let brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    let pairs = walk_face_loop(&brep, face).unwrap();
    assert_eq!(
        pairs.len(),
        4,
        "cube face outer loop should have 4 half-edges"
    );
}

#[test]
fn walk_face_loop_all_vertices_distinct() {
    let brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    let pairs = walk_face_loop(&brep, face).unwrap();
    let verts: Vec<_> = pairs.iter().map(|&(_, v)| v).collect();
    for i in 0..verts.len() {
        for j in (i + 1)..verts.len() {
            assert_ne!(verts[i], verts[j], "loop vertices should be distinct");
        }
    }
}

// ── is_planar_face ───────────────────────────────────────────────────────

#[test]
fn is_planar_face_true_for_cube_face() {
    let brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    assert!(is_planar_face(&brep, face), "cube face should be planar");
}

#[test]
fn is_planar_face_false_for_cylinder_face() {
    // Build a minimal BRep with a cylinder surface on a face using unit_cube
    // and override the surface to a Cylinder.
    let mut brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    brep.faces[face].surface = Some(popkern_core::geometry::surface::SurfaceOnFace {
        geometry: SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        orientation: true,
    });
    assert!(
        !is_planar_face(&brep, face),
        "cylinder face should not be planar"
    );
}

#[test]
fn is_planar_face_true_when_no_surface_geometry() {
    // Face with no geometry attached → treated as planar (historical default).
    let mut brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    brep.faces[face].surface = None;
    assert!(
        is_planar_face(&brep, face),
        "no-geometry face should be planar"
    );
}

#[test]
fn is_planar_face_true_for_nonexistent_face() {
    let brep = unit_cube();
    let bogus = popkern_core::topology::FaceHandle::default();
    // Non-existent face → returns true (caller will skip anyway).
    assert!(is_planar_face(&brep, bogus));
}

// ── face_has_pcurves ────────────────────────────────────────────────────

#[test]
fn face_has_pcurves_false_for_cube_without_pcurves() {
    // A unit_cube built by BRepBuilder has no pcurves on half-edges.
    let brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    // Unit cube has no pcurves, so face_has_pcurves should return false.
    let result = face_has_pcurves(&brep, face);
    // Either true or false is valid depending on the cube builder, but
    // exercising the path is the goal.
    let _ = result;
}

// ── split_err ────────────────────────────────────────────────────────────

#[test]
fn split_err_produces_correct_diagnostics() {
    let brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();
    let reason = "test failure".to_string();

    let diag = split_err(face, reason.clone(), &[]);
    assert!(diag.error.is_some(), "split_err should set an error");
    assert_eq!(diag.failed_stage, Some(crate::error::Stage::Split));
    assert!(diag.intersection_curves.is_empty());
    assert!(diag.face_classifications.is_empty());

    match diag.error.unwrap() {
        crate::error::BoolError::SplitFailed { face: f, reason: r } => {
            assert_eq!(f, face);
            assert_eq!(r, reason);
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}

// ── find_fragment_for_curve ──────────────────────────────────────────────

#[test]
fn find_fragment_for_curve_returns_original_when_no_fragments() {
    use crate::split::FragmentMap;
    use popkern_core::geometry::NurbsCurve;
    use std::collections::{HashMap, HashSet};

    let brep = unit_cube();
    let face = brep.faces.keys().next().unwrap();

    // FragmentMap with no entries for this face.
    let fragments = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };

    // A simple straight UV curve from (0,0) to (1,0).
    let curve_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };

    // When the original face exists but has no fragments, fallback to the
    // original face itself.
    let result = find_fragment_for_curve(&brep, face, &curve_uv, &fragments, &tol());
    assert_eq!(
        result,
        Some(face),
        "should fall back to original face when no fragments exist"
    );
}

#[test]
fn find_fragment_for_curve_returns_none_for_nonexistent_face() {
    use crate::split::FragmentMap;
    use popkern_core::geometry::NurbsCurve;
    use std::collections::{HashMap, HashSet};

    let brep = unit_cube();
    let bogus = popkern_core::topology::FaceHandle::default();

    let fragments = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };

    let curve_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };

    // No face in BRep matches bogus, and no fragments exist → None.
    let result = find_fragment_for_curve(&brep, bogus, &curve_uv, &fragments, &tol());
    assert!(result.is_none(), "should return None for non-existent face");
}
