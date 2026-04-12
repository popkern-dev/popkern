//! Tests for the validation test oracle.

use super::*;
use crate::tolerance::Tolerance;

#[test]
fn clean_cube_has_no_planted_defects() {
    let model = clean_cube();
    assert!(
        model.planted_defects.is_empty(),
        "clean cube should have no planted defects"
    );
}

#[test]
fn gap_model_has_one_planted_defect() {
    let model = cube_with_gap();
    assert_eq!(model.planted_defects.len(), 1);
    assert_eq!(
        model.planted_defects[0].category,
        DefectCategory::Topological
    );
}

#[test]
fn tolerance_violation_model_has_one_planted_defect() {
    let model = cube_with_tolerance_violation();
    assert_eq!(model.planted_defects.len(), 1);
    assert_eq!(model.planted_defects[0].category, DefectCategory::Tolerance);
}

#[test]
fn missing_geometry_model_has_one_planted_defect() {
    let model = cube_with_missing_geometry();
    assert_eq!(model.planted_defects.len(), 1);
    assert_eq!(
        model.planted_defects[0].category,
        DefectCategory::Degenerate
    );
}

#[test]
fn small_edge_model_has_one_planted_defect() {
    let model = cube_with_small_edge();
    assert_eq!(model.planted_defects.len(), 1);
    assert_eq!(
        model.planted_defects[0].category,
        DefectCategory::Degenerate
    );
}

#[test]
fn all_test_models_returns_five() {
    let models = all_test_models();
    assert_eq!(models.len(), 5);
}

#[test]
fn oracle_clean_cube_passes() {
    let tol = Tolerance::default();
    let results = run_oracle(&tol);
    let clean = results
        .iter()
        .find(|r| r.model_name == "clean_cube")
        .expect("clean_cube should be in results");

    assert_eq!(clean.planted, 0);
    assert!(
        clean.passed,
        "clean cube should pass: found {} defects, {} false positives",
        clean.found, clean.false_positives
    );
}

#[test]
fn oracle_detects_gap() {
    let tol = Tolerance::default();
    let results = run_oracle(&tol);
    let result = results
        .iter()
        .find(|r| r.model_name == "cube_with_gap")
        .expect("cube_with_gap should be in results");

    assert!(
        result.passed,
        "gap model should pass oracle: planted={}, found={}, missed={:?}",
        result.planted, result.found, result.missed
    );
}

#[test]
fn oracle_detects_tolerance_violation() {
    let tol = Tolerance::default();
    let results = run_oracle(&tol);
    let result = results
        .iter()
        .find(|r| r.model_name == "cube_with_tolerance_violation")
        .expect("cube_with_tolerance_violation should be in results");

    assert!(
        result.passed,
        "tolerance violation model should pass oracle: planted={}, found={}, missed={:?}",
        result.planted, result.found, result.missed
    );
}

#[test]
fn oracle_detects_missing_geometry() {
    let tol = Tolerance::default();
    let results = run_oracle(&tol);
    let result = results
        .iter()
        .find(|r| r.model_name == "cube_with_missing_geometry")
        .expect("cube_with_missing_geometry should be in results");

    assert!(
        result.passed,
        "missing geometry model should pass oracle: planted={}, found={}, missed={:?}",
        result.planted, result.found, result.missed
    );
}

#[test]
fn oracle_detects_small_edge() {
    let tol = Tolerance::default();
    let results = run_oracle(&tol);
    let result = results
        .iter()
        .find(|r| r.model_name == "cube_with_small_edge")
        .expect("cube_with_small_edge should be in results");

    assert!(
        result.passed,
        "small edge model should pass oracle: planted={}, found={}, missed={:?}",
        result.planted, result.found, result.missed
    );
}

#[test]
fn oracle_all_models_pass() {
    let tol = Tolerance::default();
    let results = run_oracle(&tol);
    assert_eq!(results.len(), 5);
    for result in &results {
        assert!(
            result.passed,
            "model '{}' failed: planted={}, found={}, missed={}, false_positives={}",
            result.model_name,
            result.planted,
            result.found,
            result.missed.len(),
            result.false_positives
        );
    }
}
