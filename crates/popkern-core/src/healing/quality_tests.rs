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
