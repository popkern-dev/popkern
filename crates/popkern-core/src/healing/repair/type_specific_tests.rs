//! Tests for geometry-type-specific healing.

use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn detect_surface_type_on_cube_returns_general() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    for (fh, _) in &brep.faces {
        let result = detect_surface_type(&brep, fh, &tolerance);
        assert_eq!(
            result.surface_type,
            SurfaceType::General,
            "cube face should be General, got {:?}",
            result.surface_type
        );
        assert_eq!(result.face, fh);
    }
}

#[test]
fn detect_all_surface_types_on_cube_all_general() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    let results = detect_all_surface_types(&brep, &tolerance);
    assert_eq!(results.len(), brep.faces.len());
    for result in &results {
        assert_eq!(result.surface_type, SurfaceType::General);
    }
}

#[test]
fn repair_fillet_on_general_face_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    let face = brep.faces.keys().next();
    if let Some(fh) = face {
        let result = repair_fillet(&mut brep, fh, &tolerance, 1.0);
        assert!(result.invariant_preserved);
        assert!((result.max_displacement - 0.0).abs() < 1e-15);
        assert_eq!(result.surface_type, SurfaceType::Fillet);
    }
}

#[test]
fn repair_offset_on_general_face_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    let face = brep.faces.keys().next();
    if let Some(fh) = face {
        let result = repair_offset(&mut brep, fh, &tolerance, 1.0);
        assert!(result.invariant_preserved);
        assert!((result.max_displacement - 0.0).abs() < 1e-15);
        assert_eq!(result.surface_type, SurfaceType::Offset);
    }
}

#[test]
fn repair_swept_on_general_face_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    let face = brep.faces.keys().next();
    if let Some(fh) = face {
        let result = repair_swept(&mut brep, fh, &tolerance, 1.0);
        assert!(result.invariant_preserved);
        assert!((result.max_displacement - 0.0).abs() < 1e-15);
        assert_eq!(result.surface_type, SurfaceType::Swept);
    }
}

#[test]
fn repair_all_type_specific_on_valid_cube_returns_empty() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    let results = repair_all_type_specific(&mut brep, &tolerance, 1.0);
    assert!(
        results.is_empty(),
        "expected no type-specific repairs on a valid cube, got {}",
        results.len()
    );
}

#[test]
fn detect_surface_type_on_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let tolerance = Tolerance::default();

    let results = detect_all_surface_types(&brep, &tolerance);
    assert!(results.is_empty());
}

#[test]
fn repair_all_type_specific_on_empty_brep_returns_empty() {
    let mut brep = BRep::new(Tolerance::default());
    let tolerance = Tolerance::default();

    let results = repair_all_type_specific(&mut brep, &tolerance, 1.0);
    assert!(results.is_empty());
}

#[test]
fn surface_type_enum_equality() {
    assert_eq!(SurfaceType::Fillet, SurfaceType::Fillet);
    assert_eq!(SurfaceType::Offset, SurfaceType::Offset);
    assert_eq!(SurfaceType::Swept, SurfaceType::Swept);
    assert_eq!(SurfaceType::General, SurfaceType::General);
    assert_ne!(SurfaceType::Fillet, SurfaceType::Offset);
    assert_ne!(SurfaceType::Offset, SurfaceType::Swept);
    assert_ne!(SurfaceType::Swept, SurfaceType::General);
}

#[test]
fn adjacent_faces_on_cube() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    for (fh, _) in &brep.faces {
        let neighbors = adjacent_faces(&brep, fh);
        // A cube face has 4 edges, each shared with a different face.
        assert!(!neighbors.is_empty(), "cube face should have neighbors");
        // No self-references.
        assert!(
            !neighbors.contains(&fh),
            "adjacent faces should not include self"
        );
    }
}

#[test]
fn measure_cross_section_uniformity_on_flat_surface() {
    // A bilinear NURBS where u is translation-only: all cross-sections
    // are identical lines at different u positions.
    use crate::geometry::Point3;
    use crate::geometry::nurbs_surface::NurbsSurface;

    let nurbs = NurbsSurface {
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
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };

    let (var_ratio, _max_shape_dev) = measure_cross_section_uniformity(&nurbs, 5);
    // Length variance should be near zero: all cross-sections have the
    // same length (they are parallel lines of equal length).
    assert!(
        var_ratio < 0.01,
        "flat surface should have near-zero length variance ratio, got {var_ratio}"
    );
    // Shape deviation can be non-trivial because cross-sections at different
    // u values are translated from each other, but the length ratio is the
    // primary swept-surface indicator.
}

#[test]
fn type_detection_result_fields() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    if let Some((fh, _)) = brep.faces.iter().next() {
        let result = detect_surface_type(&brep, fh, &tolerance);
        assert_eq!(result.face, fh);
        assert!(result.confidence >= 0.0 && result.confidence <= 1.0);
        assert!(!result.description.is_empty());
    }
}

#[test]
fn type_specific_repair_result_fields() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tolerance = Tolerance::default();

    if let Some(fh) = brep.faces.keys().next() {
        let result = repair_fillet(&mut brep, fh, &tolerance, 1.0);
        assert_eq!(result.face, fh);
        assert!(result.max_displacement >= 0.0);
        assert!(!result.description.is_empty());
    }
}

// -----------------------------------------------------------------------
// Tests with NURBS surfaces that have special type properties
// -----------------------------------------------------------------------

use crate::constructors::builder::BRepBuilder;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};

/// Build a BRep with a single NURBS face that has uniform cross-sections
/// (swept-like). Uses a 4x4 grid where each u-row is identical.
fn make_swept_like_brep() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);

    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();

    let mut cps = Vec::with_capacity(16);
    let mut weights = Vec::with_capacity(16);
    for i in 0..4 {
        for j in 0..4 {
            cps.push(Point3::new(i as f64 / 3.0, j as f64 / 3.0, 0.0));
            weights.push(1.0);
        }
    }

    let nurbs = NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps,
        num_u: 4,
        num_v: 4,
        weights,
    };

    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(nurbs),
            orientation: true,
        },
    );

    b.build_unchecked()
}

#[test]
fn detect_type_on_flat_nurbs_isolated_face() {
    // A flat NURBS face with no neighbors: swept detection checks
    // cross-section uniformity. A flat surface has uniform lengths
    // but high shape deviation (different u-positions), so it may
    // be detected as swept or general depending on normalization.
    let brep = make_swept_like_brep();
    let tolerance = Tolerance::default();
    let face = brep.faces.keys().next().unwrap();
    let result = detect_surface_type(&brep, face, &tolerance);
    // Isolated NURBS face: no neighbors for fillet/offset, so either
    // swept (if cross-sections pass) or general.
    assert!(
        result.surface_type == SurfaceType::Swept || result.surface_type == SurfaceType::General,
    );
}

#[test]
fn measure_cross_section_uniformity_on_swept_like() {
    // The make_swept_like_brep creates a flat 4x4 surface.
    // Verify the measurement functions work correctly.
    let brep = make_swept_like_brep();
    let face_data = brep.faces.values().next().unwrap();
    let sof = face_data.surface.as_ref().unwrap();
    if let SurfaceGeometry::Nurbs(ref nurbs) = sof.geometry {
        let (var_ratio, max_shape_dev) = measure_cross_section_uniformity(nurbs, 5);
        // Length variance should be near zero (all cross-sections are
        // the same length).
        assert!(var_ratio < 0.01, "var_ratio={var_ratio}");
        // Shape deviation measures absolute position displacement
        // which is non-zero for translated sections.
        assert!(max_shape_dev >= 0.0);
    }
}

#[test]
fn repair_swept_on_non_swept_surface_is_noop() {
    // Calling repair_swept on a face not detected as swept should
    // return "not a swept surface" immediately.
    let mut brep = make_swept_like_brep();
    let tolerance = Tolerance::default();
    let face = brep.faces.keys().next().unwrap();
    let det = detect_surface_type(&brep, face, &tolerance);
    if det.surface_type != SurfaceType::Swept {
        let result = repair_swept(&mut brep, face, &tolerance, 1.0);
        assert!(result.invariant_preserved);
        assert!(result.description.contains("not a swept surface"));
    }
}

#[test]
fn measure_offset_consistency_identical_surfaces() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let (median, max_dev, mean_dev) = measure_offset_consistency(&plane, &plane, 8);
    assert!(median < 1e-10);
    assert!(max_dev < 1e-10);
    assert!(mean_dev < 1e-10);
}

#[test]
fn measure_cross_section_uniformity_nonuniform() {
    let mut cps = Vec::with_capacity(9);
    let mut weights = Vec::with_capacity(9);
    for i in 0..3 {
        for j in 0..3 {
            let y = j as f64 / 2.0 * (1.0 + i as f64);
            cps.push(Point3::new(i as f64 / 2.0, y, 0.0));
            weights.push(1.0);
        }
    }
    let nurbs = NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: cps,
        num_u: 3,
        num_v: 3,
        weights,
    };
    let (var_ratio, _) = measure_cross_section_uniformity(&nurbs, 5);
    assert!(var_ratio > 0.01);
}

#[test]
fn edges_shared_with_neighbors_on_cube() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    for (fh, _) in &brep.faces {
        let edges = edges_shared_with_neighbors(&brep, fh);
        assert_eq!(edges.len(), 4);
    }
}

#[test]
fn detect_type_on_nurbs_face_runs_all_checks() {
    // Build a BRep with a NURBS face and run detection.
    // Even if nothing is detected, all detection branches execute.
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    // Convert face surface to NURBS.
    if let Some(ref surf) = brep.faces[fh].surface {
        let nurbs = surf.geometry.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0));
        brep.faces[fh].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(nurbs),
            orientation: true,
        });
    }
    let tol = Tolerance::default();
    let result = detect_surface_type(&brep, fh, &tol);
    assert_eq!(result.surface_type, SurfaceType::General);
}

#[test]
fn repair_fillet_on_nurbs_face() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Convert all face surfaces to NURBS so fillet detection runs.
    let faces: Vec<_> = brep.faces.keys().collect();
    for &fh in &faces {
        if let Some(ref surf) = brep.faces[fh].surface {
            let nurbs = surf.geometry.to_nurbs_in_range((-1.0, 2.0), (-1.0, 2.0));
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(nurbs),
                orientation: true,
            });
        }
    }
    let tol = Tolerance::default();
    let result = repair_fillet(&mut brep, faces[0], &tol, 1.0);
    // Even if not detected as fillet, the function should run without error.
    assert!(result.max_displacement >= 0.0);
}

#[test]
fn repair_offset_on_nurbs_face() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces: Vec<_> = brep.faces.keys().collect();
    for &fh in &faces {
        if let Some(ref surf) = brep.faces[fh].surface {
            let nurbs = surf.geometry.to_nurbs_in_range((-1.0, 2.0), (-1.0, 2.0));
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(nurbs),
                orientation: true,
            });
        }
    }
    let tol = Tolerance::default();
    let result = repair_offset(&mut brep, faces[0], &tol, 1.0);
    assert!(result.max_displacement >= 0.0);
}

#[test]
fn repair_swept_on_nurbs_face() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces: Vec<_> = brep.faces.keys().collect();
    for &fh in &faces {
        if let Some(ref surf) = brep.faces[fh].surface {
            let nurbs = surf.geometry.to_nurbs_in_range((-1.0, 2.0), (-1.0, 2.0));
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(nurbs),
                orientation: true,
            });
        }
    }
    let tol = Tolerance::default();
    let result = repair_swept(&mut brep, faces[0], &tol, 1.0);
    assert!(result.max_displacement >= 0.0);
}

#[test]
fn repair_all_type_specific_on_nurbs_cube() {
    // Convert all surfaces to NURBS, run repair_all.
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces: Vec<_> = brep.faces.keys().collect();
    for &fh in &faces {
        if let Some(ref surf) = brep.faces[fh].surface {
            let nurbs = surf.geometry.to_nurbs_in_range((-1.0, 2.0), (-1.0, 2.0));
            brep.faces[fh].surface = Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(nurbs),
                orientation: true,
            });
        }
    }
    let tol = Tolerance::default();
    let results = repair_all_type_specific(&mut brep, &tol, 1.0);
    // Plane-converted NURBS typically not detected as special types,
    // so results should be empty or contain only type-specific repairs.
    for r in &results {
        assert!(r.max_displacement >= 0.0);
    }
}
