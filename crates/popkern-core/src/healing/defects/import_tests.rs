//! Tests for import-specific defect detection.

use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_no_import_defects() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    // The unit cube has planar faces (non-periodic), so no lost
    // periodic topology, and valid orientation.
    assert!(
        report.is_empty(),
        "Valid cube should have no import defects, got {}",
        report.len()
    );
}

#[test]
fn orientation_inconsistency_detected() {
    use crate::geometry::Point3;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    // Attach a degenerate surface to trigger orientation failure.
    brep.faces[face_handle].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.0,
        },
        orientation: true,
    });

    let mut report = HealingReport::new();
    check_orientation_inconsistencies(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Degenerate surface should trigger orientation defect"
    );
    assert_eq!(report.defects()[0].category, DefectCategory::Import);
}

#[test]
fn empty_brep_no_import_defects() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    assert!(report.is_empty());
}

#[test]
fn non_periodic_face_skipped() {
    // A cube has only planar faces; none should trigger periodic checks.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_lost_periodic_topology(&brep, &mut report);
    assert!(report.is_empty());
}

#[test]
fn lost_singularity_topology_on_cube_produces_no_defect() {
    // A cube has planar faces; no sphere/cone surfaces.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_lost_singularity_topology(&brep, &mut report);
    assert!(
        report.is_empty(),
        "Cube with planar faces should have no singularity topology defects"
    );
}

#[test]
fn coordinate_artifacts_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_coordinate_artifacts(&brep, &tol, &mut report);
    // Unit cube vertices are at 0.0 and 1.0 — none in the suspicious range.
    assert!(
        report.is_empty(),
        "Valid cube should have no coordinate artifact defects"
    );
}

#[test]
fn empty_brep_no_singularity_or_artifact_defects() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_lost_singularity_topology(&brep, &mut report);
    check_coordinate_artifacts(&brep, &tol, &mut report);
    assert!(report.is_empty());
}

/// Knot multiplicities on a knot vector with no interior knots.
#[test]
fn knot_multiplicities_no_interior() {
    // Degree 2, knot vector [0,0,0,1,1,1] has no interior knots.
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let result = knot_multiplicities(&knots, 2);
    assert!(
        result.is_empty(),
        "Should have no interior knots, got {:?}",
        result
    );
}

/// Knot multiplicities detect a C0 break (multiplicity == degree).
#[test]
fn knot_multiplicities_c0_break() {
    // Degree 3, knot vector [0,0,0,0, 0.5,0.5,0.5, 1,1,1,1]
    // Interior knot 0.5 with multiplicity 3 == degree 3 => C0.
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0];
    let result = knot_multiplicities(&knots, 3);
    assert_eq!(
        result.len(),
        1,
        "Should find one distinct interior knot group"
    );
    assert_eq!(result[0].1, 3, "Multiplicity should be 3");
}

/// Knot multiplicities on a simple cubic with one interior knot.
#[test]
fn knot_multiplicities_single_interior() {
    // Degree 3, knot vector [0,0,0,0, 0.5, 1,1,1,1]
    // Interior knot 0.5 with multiplicity 1 => C2 (no warning).
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
    let result = knot_multiplicities(&knots, 3);
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].1, 1, "Multiplicity should be 1");
}

/// Cube with no NURBS geometry produces no representation discontinuities.
#[test]
fn cube_no_representation_discontinuities() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_representation_discontinuities(&brep, &tol, &mut report);
    assert!(
        report.is_empty(),
        "Cube with analytic geometry should have no representation discontinuities"
    );
}
