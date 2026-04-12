//! Tests for the top-level healing API (healing/mod.rs).
//!
//! Covers: analyze, analyze_with, heal, heal_with, preview_repairs,
//! analyze_with_telemetry, heal_with_telemetry on both valid and empty BReps.

use crate::healing::options::{AnalysisOptions, CheckSet, HealOptions, HealStrategy};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn analyze_valid_cube_returns_clean_report() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let report = super::analyze(&brep, &Tolerance::default());
    let summary = report.summary();
    assert_eq!(summary.critical, 0);
    assert_eq!(summary.error, 0);
}

#[test]
fn analyze_with_custom_options() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let opts = AnalysisOptions {
        checks: CheckSet::FAST,
        ..AnalysisOptions::new(Tolerance::default())
    };
    let report = super::analyze_with(&brep, &opts);
    assert_eq!(report.summary().critical, 0);
}

#[test]
fn heal_valid_cube_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let face_count = brep.faces.len();
    let result = super::heal(&mut brep, &Tolerance::default());
    assert!(result.is_ok());
    let result = result.unwrap_or_else(|_| panic!("heal should succeed on a valid cube"));
    assert_eq!(brep.faces.len(), face_count);
    assert!(result.defects_after <= result.defects_before);
}

#[test]
fn heal_with_conservative_strategy() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let opts = HealOptions {
        strategy: HealStrategy::Conservative,
        ..HealOptions::new(Tolerance::default())
    };
    let result = super::heal_with(&mut brep, &opts);
    assert!(result.is_ok());
    let result = result.unwrap_or_else(|_| panic!("heal_with should succeed on a valid cube"));
    assert!(result.defects_after <= result.defects_before);
}

#[test]
fn preview_repairs_does_not_modify_brep() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let vertex_count = brep.vertices.len();
    let opts = HealOptions::new(Tolerance::default());
    let preview = super::preview_repairs(&brep, &opts);
    assert_eq!(brep.vertices.len(), vertex_count);
    // Preview may have proposals (or not for a valid cube).
    let _ = preview.proposals.len();
}

#[test]
fn analyze_with_telemetry_returns_timing() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let opts = AnalysisOptions::new(Tolerance::default());
    let (report, telemetry) = super::analyze_with_telemetry(&brep, &opts);
    assert_eq!(report.summary().critical, 0);
    assert!(telemetry.operation_count() > 0);
    assert!(telemetry.total_duration().as_nanos() > 0);
}

#[test]
fn heal_with_telemetry_returns_timing() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let opts = HealOptions::new(Tolerance::default());
    let result = super::heal_with_telemetry(&mut brep, &opts);
    assert!(result.is_ok());
    let (result, telemetry) =
        result.unwrap_or_else(|_| panic!("heal_with_telemetry should succeed"));
    assert!(result.defects_after <= result.defects_before);
    assert!(telemetry.operation_count() > 0);
}

#[test]
fn analyze_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let report = super::analyze(&brep, &Tolerance::default());
    assert!(report.is_empty());
}

#[test]
fn heal_empty_brep() {
    let mut brep = BRep::new(Tolerance::default());
    let result = super::heal(&mut brep, &Tolerance::default());
    assert!(result.is_ok());
    let result = result.unwrap_or_else(|_| panic!("heal should succeed on empty brep"));
    assert_eq!(result.defects_before, 0);
    assert_eq!(result.defects_after, 0);
}

#[test]
fn analyze_delegates_to_analyze_with() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let tol = Tolerance::default();
    let direct = super::analyze(&brep, &tol);
    let via_opts = super::analyze_with(&brep, &AnalysisOptions::new(tol));
    assert_eq!(direct.len(), via_opts.len());
}

#[test]
fn heal_delegates_to_heal_with() {
    let mut brep1 = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let mut brep2 = brep1.clone();
    let tol = Tolerance::default();

    let r1 = super::heal(&mut brep1, &tol);
    let r2 = super::heal_with(&mut brep2, &HealOptions::new(tol));

    assert!(r1.is_ok());
    assert!(r2.is_ok());
    let r1 = r1.unwrap_or_else(|_| panic!("heal should succeed"));
    let r2 = r2.unwrap_or_else(|_| panic!("heal_with should succeed"));
    assert_eq!(r1.defects_before, r2.defects_before);
    assert_eq!(r1.defects_after, r2.defects_after);
}

#[test]
fn preview_repairs_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let opts = HealOptions::new(Tolerance::default());
    let preview = super::preview_repairs(&brep, &opts);
    assert!(preview.is_empty());
    assert_eq!(preview.estimated_defects_after, 0);
}

#[test]
fn analyze_with_geometric_only_checks() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let opts = AnalysisOptions {
        checks: CheckSet::GEOMETRIC_ONLY,
        ..AnalysisOptions::new(Tolerance::default())
    };
    let report = super::analyze_with(&brep, &opts);
    assert_eq!(report.summary().critical, 0);
}

#[test]
fn heal_with_aggressive_strategy() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let opts = HealOptions {
        strategy: HealStrategy::Aggressive,
        ..HealOptions::new(Tolerance::default())
    };
    let result = super::heal_with(&mut brep, &opts);
    assert!(result.is_ok());
    let result = result.unwrap_or_else(|_| panic!("aggressive heal should succeed on valid cube"));
    assert!(result.defects_after <= result.defects_before);
}
