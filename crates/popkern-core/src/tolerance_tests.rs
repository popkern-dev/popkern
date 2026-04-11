use super::*;

// -----------------------------------------------------------------------
// Original tests (migrated to getter API)
// -----------------------------------------------------------------------

#[test]
fn default_tolerance_values() {
    let tol = Tolerance::default();
    assert_eq!(tol.linear(), 1e-6);
    assert_eq!(tol.angular(), 1e-4);
    assert_eq!(tol.parametric_base(), 1e-8);
}

#[test]
fn builder_pattern() {
    let tol = Tolerance::default().with_linear(1e-3).with_angular(1e-2);
    assert_eq!(tol.linear(), 1e-3);
    assert_eq!(tol.angular(), 1e-2);
    assert_eq!(tol.parametric_base(), 1e-8);
}

#[test]
fn points_coincident_within_tolerance() {
    let tol = Tolerance::default();
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1e-7, 0.0, 0.0);
    assert!(tol.points_coincident(&a, &b));
}

#[test]
fn points_not_coincident_outside_tolerance() {
    let tol = Tolerance::default();
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1e-5, 0.0, 0.0);
    assert!(!tol.points_coincident(&a, &b));
}

#[test]
fn angles_equal_within_tolerance() {
    let tol = Tolerance::default();
    assert!(tol.angles_equal(1.0, 1.0 + 1e-5));
}

#[test]
fn angles_not_equal_outside_tolerance() {
    let tol = Tolerance::default();
    assert!(!tol.angles_equal(1.0, 1.0 + 1e-3));
}

#[test]
#[allow(deprecated)]
fn params_equal_within_tolerance() {
    let tol = Tolerance::default();
    assert!(tol.params_equal(0.5, 0.5 + 1e-9));
}

#[test]
fn distances_equal_within_tolerance() {
    let tol = Tolerance::default();
    assert!(tol.distances_equal(5.0, 5.0 + 1e-7));
}

#[test]
fn with_parametric_builder() {
    let tol = Tolerance::default().with_parametric(1e-3);
    assert_eq!(tol.parametric_base(), 1e-3);
    assert_eq!(tol.linear(), 1e-6);
    assert_eq!(tol.angular(), 1e-4);
}

#[test]
fn with_angular_builder() {
    let tol = Tolerance::default().with_angular(1e-1);
    assert_eq!(tol.angular(), 1e-1);
    assert_eq!(tol.linear(), 1e-6);
    assert_eq!(tol.parametric_base(), 1e-8);
}

#[test]
fn chained_builders() {
    let tol = Tolerance::default()
        .with_linear(1e-3)
        .with_angular(1e-2)
        .with_parametric(1e-1);
    assert_eq!(tol.linear(), 1e-3);
    assert_eq!(tol.angular(), 1e-2);
    assert_eq!(tol.parametric_base(), 1e-1);
}

// -----------------------------------------------------------------------
// New field tests
// -----------------------------------------------------------------------

#[test]
fn default_new_fields() {
    let tol = Tolerance::default();
    assert_eq!(tol.numerical_floor(), 1e-14);
    assert_eq!(tol.solver(), 1e-10);
    assert_eq!(tol.approximation(), 1e-5);
}

#[test]
fn with_solver_builder() {
    let tol = Tolerance::default().with_solver(1e-8);
    assert_eq!(tol.solver(), 1e-8);
    assert_eq!(tol.linear(), 1e-6);
}

#[test]
fn with_numerical_floor_builder() {
    let tol = Tolerance::default().with_numerical_floor(1e-12);
    assert_eq!(tol.numerical_floor(), 1e-12);
}

#[test]
fn with_approximation_builder() {
    let tol = Tolerance::default().with_approximation(1e-3);
    assert_eq!(tol.approximation(), 1e-3);
}

// -----------------------------------------------------------------------
// parametric_at tests
// -----------------------------------------------------------------------

#[test]
fn parametric_at_normal_derivative() {
    let tol = Tolerance::default();
    // derivative_mag = 1.0 -> linear / 1.0 = 1e-6, max(1e-8) = 1e-6
    let p = tol.parametric_at(1.0);
    assert!((p - 1e-6).abs() < 1e-15);
}

#[test]
fn parametric_at_large_derivative() {
    let tol = Tolerance::default();
    // derivative_mag = 1000 -> linear / 1000 = 1e-9, max(1e-8) = 1e-8
    let p = tol.parametric_at(1000.0);
    assert_eq!(p, tol.parametric_base());
}

#[test]
fn parametric_at_degenerate() {
    let tol = Tolerance::default();
    // derivative_mag below numerical_floor -> returns parametric_base
    let p = tol.parametric_at(1e-16);
    assert_eq!(p, tol.parametric_base());
}

// -----------------------------------------------------------------------
// from_model_scale tests
// -----------------------------------------------------------------------

#[test]
fn from_model_scale_one_meter() {
    let tol = Tolerance::from_model_scale(1.0);
    assert!((tol.linear() - 1e-6).abs() < 1e-15);
    assert_eq!(tol.angular(), 1e-4);
}

#[test]
fn from_model_scale_tiny() {
    // diagonal = 1e-8 -> linear = max(1e-14, 1e-12) = 1e-12
    let tol = Tolerance::from_model_scale(1e-8);
    assert_eq!(tol.linear(), 1e-12);
}

#[test]
fn from_model_scale_large() {
    // diagonal = 1000 -> linear = 1e-3
    let tol = Tolerance::from_model_scale(1000.0);
    assert!((tol.linear() - 1e-3).abs() < 1e-15);
}

// -----------------------------------------------------------------------
// scale tests
// -----------------------------------------------------------------------

#[test]
fn scale_preserves_invariants() {
    let tol = Tolerance::default();
    let scaled = tol.scale(1000.0);
    // Invariant fields unchanged
    assert_eq!(scaled.angular(), 1e-4);
    assert_eq!(scaled.parametric_base(), 1e-8);
    assert_eq!(scaled.numerical_floor(), 1e-14);
    // Scaled fields
    assert!((scaled.linear() - 1e-3).abs() < 1e-15);
    assert!((scaled.solver() - 1e-7).abs() < 1e-15);
    assert!((scaled.approximation() - 1e-2).abs() < 1e-15);
}

// -----------------------------------------------------------------------
// Preset tests
// -----------------------------------------------------------------------

#[test]
fn preset_standard_matches_default() {
    let tol = Tolerance::preset(TolerancePreset::Standard, None);
    assert_eq!(tol.linear(), 1e-6);
    assert_eq!(tol.solver(), 1e-10);
    assert_eq!(tol.approximation(), 1e-5);
}

#[test]
fn preset_fine() {
    let tol = Tolerance::preset(TolerancePreset::Fine, None);
    assert_eq!(tol.linear(), 1e-8);
    assert_eq!(tol.solver(), 1e-12);
    assert_eq!(tol.approximation(), 1e-7);
}

#[test]
fn preset_draft() {
    let tol = Tolerance::preset(TolerancePreset::Draft, None);
    assert_eq!(tol.linear(), 1e-4);
    assert_eq!(tol.solver(), 1e-8);
    assert_eq!(tol.approximation(), 1e-3);
}

// -----------------------------------------------------------------------
// Comparison method tests
// -----------------------------------------------------------------------

#[test]
fn nearly_zero() {
    let tol = Tolerance::default();
    assert!(tol.nearly_zero(1e-15));
    assert!(!tol.nearly_zero(1e-13));
}

#[test]
fn solver_converged() {
    let tol = Tolerance::default();
    assert!(tol.solver_converged(1e-11));
    assert!(!tol.solver_converged(1e-9));
}

#[test]
fn params_equal_at_method() {
    let tol = Tolerance::default();
    // derivative_mag = 1.0 -> threshold = 1e-6
    assert!(tol.params_equal_at(0.5, 0.5 + 1e-7, 1.0));
    assert!(!tol.params_equal_at(0.5, 0.5 + 1e-5, 1.0));
}

#[test]
fn approx_adequate() {
    let tol = Tolerance::default();
    assert!(tol.approx_adequate(1e-6));
    assert!(!tol.approx_adequate(1e-4));
}

// -----------------------------------------------------------------------
// Serde round-trip
// -----------------------------------------------------------------------

#[test]
fn serde_round_trip() {
    let tol = Tolerance::default().with_linear(1e-3);
    let json = serde_json::to_string(&tol).expect("serialize");
    let back: Tolerance = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(tol, back);
}

#[test]
fn serde_old_format_compat() {
    // Old format used "parametric" instead of "parametric_base"
    let json = r#"{"linear":1e-6,"angular":1e-4,"parametric":1e-8}"#;
    let tol: Tolerance = serde_json::from_str(json).expect("deserialize old format");
    assert_eq!(tol.parametric_base(), 1e-8);
    // New fields should get defaults
    assert_eq!(tol.numerical_floor(), 1e-14);
    assert_eq!(tol.solver(), 1e-10);
    assert_eq!(tol.approximation(), 1e-5);
}
