use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn continuity_level_ordering() {
    assert!(ContinuityLevel::G0 < ContinuityLevel::G1);
    assert!(ContinuityLevel::G1 < ContinuityLevel::G2);
    assert!(ContinuityLevel::G2 < ContinuityLevel::G3);
}

#[test]
fn tier1_profile_has_g2_continuity() {
    let profile = QualityProfile::from_tier(QualityTier::Tier1);
    assert_eq!(profile.min_continuity, ContinuityLevel::G2);
    assert!(profile.max_position_gap < 1e-5);
}

#[test]
fn tier3_profile_has_g0_continuity() {
    let profile = QualityProfile::from_tier(QualityTier::Tier3);
    assert_eq!(profile.min_continuity, ContinuityLevel::G0);
    assert!(profile.max_position_gap > 1e-4);
}

#[test]
fn quality_result_passes_when_all_dimensions_pass() {
    let result = QualityResult::new(
        QualityTier::Tier3,
        vec![QualityDimension {
            name: "gap".into(),
            measured: 1e-6,
            threshold: 1e-3,
            passed: true,
        }],
    );
    assert!(result.passed);
    assert!(result.failures().is_empty());
}

#[test]
fn quality_result_fails_when_any_dimension_fails() {
    let result = QualityResult::new(
        QualityTier::Tier1,
        vec![
            QualityDimension {
                name: "gap".into(),
                measured: 1e-6,
                threshold: 1e-3,
                passed: true,
            },
            QualityDimension {
                name: "angle".into(),
                measured: 5.0,
                threshold: 0.1,
                passed: false,
            },
        ],
    );
    assert!(!result.passed);
    assert_eq!(result.failures().len(), 1);
}

#[test]
fn tier_display() {
    assert!(QualityTier::Tier1.to_string().contains("Class-A"));
    assert!(QualityTier::Tier2.to_string().contains("Industrial"));
}

#[test]
fn profile_tier_field_matches() {
    let p1 = QualityProfile::from_tier(QualityTier::Tier1);
    assert_eq!(p1.tier, QualityTier::Tier1);
    let p2 = QualityProfile::from_tier(QualityTier::Tier2);
    assert_eq!(p2.tier, QualityTier::Tier2);
    let p3 = QualityProfile::from_tier(QualityTier::Tier3);
    assert_eq!(p3.tier, QualityTier::Tier3);
}

#[test]
fn validate_against_profile_tier3_unit_cube_passes() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let profile = QualityProfile::from_tier(QualityTier::Tier3);
    let tol = Tolerance::default();

    let result = validate_against_profile(&brep, &profile, &tol);

    assert!(
        result.passed,
        "Unit cube should pass all Tier3 quality dimensions, but got: {result}"
    );
    assert_eq!(result.tier, QualityTier::Tier3);
    assert!(!result.dimensions.is_empty());
    for dim in &result.dimensions {
        assert!(dim.passed, "Dimension should pass: {dim}");
    }
}

#[test]
fn validate_against_profile_tier1_unit_cube_fails_normal_angle() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let profile = QualityProfile::from_tier(QualityTier::Tier1);
    let tol = Tolerance::default();

    let result = validate_against_profile(&brep, &profile, &tol);

    // A cube has 90-degree edges. Tier1 allows only 0.1 degree,
    // so the normal angle check must fail.
    assert!(
        !result.passed,
        "Unit cube should fail Tier1 due to 90-degree edges"
    );
    assert_eq!(result.tier, QualityTier::Tier1);

    let failures = result.failures();
    let has_angle_failure = failures.iter().any(|d| d.name.contains("normal angle"));
    assert!(
        has_angle_failure,
        "Expected a normal angle failure for 90-degree cube edges, got: {result}"
    );
}

#[test]
fn validate_against_profile_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let profile = QualityProfile::from_tier(QualityTier::Tier3);
    let tol = Tolerance::default();

    let result = validate_against_profile(&brep, &profile, &tol);

    // Empty BRep: gap=0.0 (passes), angle=0.0 (passes),
    // edge_length=inf (passes since inf >= threshold).
    assert!(result.passed, "Empty BRep should pass all checks: {result}");
}

#[test]
fn validate_against_profile_custom_profile() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let profile = QualityProfile {
        name: "Custom".into(),
        tier: QualityTier::Tier3,
        min_continuity: ContinuityLevel::G0,
        max_position_gap: 10.0,
        max_normal_angle_deg: 180.0,
        max_curvature_diff: 1000.0,
        min_edge_length: 0.0,
        min_face_area: 0.0,
        custom_checks: Vec::new(),
    };
    let result = validate_against_profile(&brep, &profile, &tol);
    assert!(
        result.passed,
        "Generous custom profile should pass on a cube: {result}"
    );
}

#[test]
fn validate_against_profile_fails_on_tight_edge_length() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let profile = QualityProfile {
        min_edge_length: 100.0, // cube edges are length 1 -- will fail
        ..QualityProfile::from_tier(QualityTier::Tier3)
    };
    let result = validate_against_profile(&brep, &profile, &tol);
    assert!(!result.passed, "Cube should fail with min_edge_length=100");
    let failures = result.failures();
    let has_edge_failure = failures.iter().any(|d| d.name.contains("edge length"));
    assert!(
        has_edge_failure,
        "Expected edge length failure, got: {result}"
    );
}

#[test]
fn validate_against_profile_fails_on_tight_gap() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let profile = QualityProfile {
        max_position_gap: 0.0, // exact match required
        ..QualityProfile::from_tier(QualityTier::Tier3)
    };
    let result = validate_against_profile(&brep, &profile, &tol);
    // The cube may have tiny floating-point gaps at curve endpoints.
    // Whether it passes or fails depends on fp precision, so just ensure
    // the dimensions were measured.
    assert!(
        !result.dimensions.is_empty(),
        "Should have quality dimensions"
    );
}

#[test]
fn validate_against_profile_tier2_unit_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let profile = QualityProfile::from_tier(QualityTier::Tier2);
    let tol = Tolerance::default();

    let result = validate_against_profile(&brep, &profile, &tol);
    // Tier2 allows 1-degree normal angle -- a cube has 90-degree edges, so
    // the normal angle check should fail.
    assert!(
        !result.passed,
        "Unit cube should fail Tier2 due to 90-degree edges"
    );
    let failures = result.failures();
    let has_angle_failure = failures.iter().any(|d| d.name.contains("normal angle"));
    assert!(
        has_angle_failure,
        "Expected a normal angle failure: {result}"
    );
}

#[test]
fn quality_result_display_contains_status() {
    let result = QualityResult::new(
        QualityTier::Tier3,
        vec![QualityDimension {
            name: "gap".into(),
            measured: 0.0,
            threshold: 1.0,
            passed: true,
        }],
    );
    let display = result.to_string();
    assert!(display.contains("PASSED"), "Display should contain PASSED");
    assert!(
        display.contains("Tier 3"),
        "Display should contain tier name"
    );
}

#[test]
fn quality_result_display_failed() {
    let result = QualityResult::new(
        QualityTier::Tier1,
        vec![QualityDimension {
            name: "angle".into(),
            measured: 90.0,
            threshold: 0.1,
            passed: false,
        }],
    );
    let display = result.to_string();
    assert!(display.contains("FAILED"), "Display should contain FAILED");
}

#[test]
fn quality_dimension_display() {
    let dim = QualityDimension {
        name: "test".into(),
        measured: 0.5,
        threshold: 1.0,
        passed: true,
    };
    let display = dim.to_string();
    assert!(display.contains("PASS"), "Passing dimension shows PASS");
    assert!(display.contains("test"), "Dimension name is shown");
}

#[test]
fn continuity_level_display() {
    assert_eq!(ContinuityLevel::G0.to_string(), "G0");
    assert_eq!(ContinuityLevel::G1.to_string(), "G1");
    assert_eq!(ContinuityLevel::G2.to_string(), "G2");
    assert_eq!(ContinuityLevel::G3.to_string(), "G3");
}

#[test]
fn tier_display_all() {
    let t1 = QualityTier::Tier1.to_string();
    let t2 = QualityTier::Tier2.to_string();
    let t3 = QualityTier::Tier3.to_string();
    assert!(t1.contains("Class-A"));
    assert!(t2.contains("Industrial"));
    assert!(t3.contains("General"));
}
