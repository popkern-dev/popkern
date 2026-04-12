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
