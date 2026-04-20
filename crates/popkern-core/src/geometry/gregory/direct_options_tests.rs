//! Tests for [`super::DirectConstructionOptions`].

use super::*;

#[test]
fn defaults_are_standard_validation_g1_target() {
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    assert_eq!(opts.validation_level, ValidationLevel::Standard);
    assert_eq!(opts.target_continuity, ContinuityLevel::G1);
}

#[test]
fn characteristic_length_scales_absolute_thresholds() {
    let a = DirectConstructionOptions::with_characteristic_length(1.0);
    let b = DirectConstructionOptions::with_characteristic_length(1000.0);
    // `min_absolute_radius` should scale linearly.
    assert!((b.min_absolute_radius() / a.min_absolute_radius() - 1000.0).abs() < 1e-12);
    // So should the max.
    assert!((b.max_absolute_radius() / a.max_absolute_radius() - 1000.0).abs() < 1e-12);
}

#[test]
fn with_continuity_overrides_target() {
    let opts = DirectConstructionOptions::with_characteristic_length(1.0)
        .with_continuity(ContinuityLevel::G2);
    assert_eq!(opts.target_continuity, ContinuityLevel::G2);
}

#[test]
fn with_validation_overrides_level() {
    let opts = DirectConstructionOptions::with_characteristic_length(1.0)
        .with_validation(ValidationLevel::Full);
    assert_eq!(opts.validation_level, ValidationLevel::Full);
    // Other fields stay at defaults.
    assert_eq!(opts.target_continuity, ContinuityLevel::G1);
}

#[test]
fn defaults_satisfy_radius_ratio_invariants() {
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    // min < max, and both positive.
    assert!(opts.min_radius_ratio > 0.0);
    assert!(opts.max_radius_ratio > opts.min_radius_ratio);
    // Asymmetry ceiling is > 1 (degenerate if <= 1).
    assert!(opts.max_radius_asymmetry > 1.0);
}

#[test]
fn tolerance_knobs_are_positive() {
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    assert!(opts.tangent_corner_rad > 0.0);
    assert!(opts.boundary_tangent_eps > 0.0);
    assert!(opts.corner_singularity_eps > 0.0);
    assert!(opts.self_check_g1_tol > 0.0);
    assert!(opts.self_check_g2_tol > 0.0);
}
