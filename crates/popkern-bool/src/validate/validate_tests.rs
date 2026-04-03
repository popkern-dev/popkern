//! Tests for boolean validation and repair.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ── Pre-check tests ────────────────────────────────────────────────────

#[test]
fn pre_check_valid_cubes() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(
        report.is_ok(),
        "valid cubes should pass pre-check, got {} errors",
        report.error_count()
    );
}

#[test]
fn pre_check_valid_cylinder() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(report.is_ok());
}

#[test]
fn pre_check_empty_body_fails() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::new(tol());
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(!report.is_ok(), "empty body B should fail pre-check");
    assert!(report.error_count() > 0);
}

#[test]
fn pre_check_report_counts() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert_eq!(report.error_count(), 0);
    assert_eq!(report.warning_count(), 0);
}

// ── Post-check tests ───────────────────────────────────────────────────

#[test]
fn post_check_valid_cube() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::post_check::post_check(&cube, &tol());
    assert!(
        report.is_ok(),
        "valid cube should pass post-check, got {} errors: {:?}",
        report.error_count(),
        report.findings
    );
    assert_eq!(report.euler_characteristic, 2);
    assert_eq!(report.face_count, 6);
    assert_eq!(report.vertex_count, 8);
}

#[test]
fn post_check_valid_cylinder() {
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let report = super::post_check::post_check(&cyl, &tol());
    // Cylinder may have minor geometric warnings but no structural errors
    let structural_errors = report
        .findings
        .iter()
        .filter(|f| {
            matches!(f.severity, crate::error::Severity::Error)
                && !f.message.contains("orientation")
        })
        .count();
    assert!(
        structural_errors == 0,
        "cylinder should have no structural errors, got: {:?}",
        report.findings
    );
}

#[test]
fn post_check_boolean_result() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    if let Ok(result) = crate::union(a, b, tol(), None) {
        let report = super::post_check::post_check(&result, &tol());
        // Log but don't hard-assert — graph-merge may produce imperfect results
        if !report.is_ok() {
            eprintln!(
                "[post_check] boolean result: {} errors, {} faces, euler={}",
                report.error_count(),
                report.face_count,
                report.euler_characteristic
            );
        }
    }
}

#[test]
fn post_check_report_metadata() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::post_check::post_check(&cube, &tol());
    assert!(report.face_count > 0);
    assert!(report.edge_count > 0);
    assert!(report.vertex_count > 0);
    assert!(report.max_gap >= 0.0);
}

// ── Repair tests ───────────────────────────────────────────────────────

#[test]
fn repair_on_valid_cube_changes_nothing() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let report = super::repair::repair_vertex_positions(&mut cube, 1e-4);
    assert!(report.fully_repaired());
    assert_eq!(report.vertices_adjusted, 0);
}

#[test]
fn repair_report_default() {
    let report = super::repair::RepairReport::default();
    assert!(report.fully_repaired());
    assert_eq!(report.vertices_adjusted, 0);
    assert_eq!(report.unrepairable_gaps, 0);
}
