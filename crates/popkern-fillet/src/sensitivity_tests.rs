use super::*;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Sanity: running sensitivity on a plain cube edge produces a
/// report with both base and step quality filled in, and the
/// `radius_base` matches the input spec.
#[test]
fn sensitivity_on_cube_edge_produces_report() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let report = quality_sensitivity(&cube, &specs, &FilletOptions::default(), tol(), 0.01)
        .expect("sensitivity should compute for a valid spec");

    assert!((report.radius_base - 0.1).abs() < 1e-12);
    assert!((report.epsilon - 0.01).abs() < 1e-12);
    // Both qualities were computed — `composite_score` in [0, 1].
    assert!(report.quality_base.composite_score >= 0.0);
    assert!(report.quality_base.composite_score <= 1.0);
    assert!(report.quality_step.composite_score >= 0.0);
    assert!(report.quality_step.composite_score <= 1.0);
}

/// On a unit cube, the blend surface area should grow with the
/// radius (larger torus strip ⇒ more area). Expect
/// `d_area_d_r > 0`.
#[test]
fn d_area_d_r_is_positive_on_cube_fillet() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let report =
        quality_sensitivity(&cube, &specs, &FilletOptions::default(), tol(), 0.01).unwrap();
    assert!(
        report.d_area_d_r > 0.0,
        "expected d(area)/dr > 0 on a simple cube fillet, got {}",
        report.d_area_d_r
    );
}

/// `d_r_min_d_r` should be close to 1 on a cube fillet where the
/// surface is an exact torus and the radius flows through linearly.
#[test]
fn d_r_min_d_r_is_near_unity_on_cube_fillet() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let report =
        quality_sensitivity(&cube, &specs, &FilletOptions::default(), tol(), 0.01).unwrap();
    // Torus minor radius IS the fillet radius — the numerical
    // derivative should be 1.0 exactly modulo FD noise.
    assert!(
        (report.d_r_min_d_r - 1.0).abs() < 0.1,
        "expected d(r_min)/dr ≈ 1 on a torus-producing fillet, got {}",
        report.d_r_min_d_r
    );
}

/// Empty spec list returns `InvalidParameter`.
#[test]
fn empty_specs_returns_invalid_parameter() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let err = quality_sensitivity(&cube, &[], &FilletOptions::default(), tol(), 0.01).unwrap_err();
    match err {
        BlendError::InvalidParameter { parameter, .. } => {
            assert_eq!(parameter, "specs");
        }
        other => panic!("expected InvalidParameter, got {other:?}"),
    }
}

/// Non-positive epsilon returns `InvalidParameter`.
#[test]
fn non_positive_epsilon_returns_invalid_parameter() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    for bad in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        let err =
            quality_sensitivity(&cube, &specs, &FilletOptions::default(), tol(), bad).unwrap_err();
        match err {
            BlendError::InvalidParameter { parameter, .. } => assert_eq!(parameter, "epsilon"),
            other => panic!("expected InvalidParameter for epsilon={bad}, got {other:?}"),
        }
    }
}

/// `radius_for_perfect_score` returns `Some(radius)` when the
/// derivative is positive, and `None` when non-positive. Exercise
/// both branches via a synthetic report.
#[test]
fn radius_for_perfect_score_handles_both_signs() {
    use crate::precheck::BlendQuality;

    let mut report = SensitivityReport {
        radius_base: 0.1,
        epsilon: 0.01,
        quality_base: BlendQuality {
            r_min: 0.1,
            patch_count: 1,
            surface_area: 1.0,
            composite_score: 0.8,
            achieved_continuity: crate::ContinuityLevel::G1,
        },
        quality_step: BlendQuality {
            r_min: 0.11,
            patch_count: 1,
            surface_area: 1.1,
            composite_score: 0.9,
            achieved_continuity: crate::ContinuityLevel::G1,
        },
        d_composite_d_r: 10.0,
        d_r_min_d_r: 1.0,
        d_area_d_r: 10.0,
    };
    // Positive derivative: linear extrapolation `(1 - 0.8) / 10 = 0.02` → `0.12`.
    let perfect = report.radius_for_perfect_score().unwrap();
    assert!((perfect - 0.12).abs() < 1e-12);

    // Zero derivative: returns None.
    report.d_composite_d_r = 0.0;
    assert!(report.radius_for_perfect_score().is_none());

    // Negative derivative: returns None.
    report.d_composite_d_r = -1.0;
    assert!(report.radius_for_perfect_score().is_none());

    // Derivative positive but already at perfect score: return the
    // base radius unchanged.
    report.d_composite_d_r = 10.0;
    report.quality_base.composite_score = 1.0;
    assert_eq!(report.radius_for_perfect_score().unwrap(), 0.1);
}

/// `grows_with_radius` is a simple wrapper around the sign of the
/// composite derivative.
#[test]
fn grows_with_radius_matches_derivative_sign() {
    let base = SensitivityReport {
        radius_base: 0.1,
        epsilon: 0.01,
        quality_base: crate::precheck::BlendQuality {
            r_min: 0.1,
            patch_count: 1,
            surface_area: 1.0,
            composite_score: 0.8,
            achieved_continuity: crate::ContinuityLevel::G1,
        },
        quality_step: crate::precheck::BlendQuality {
            r_min: 0.11,
            patch_count: 1,
            surface_area: 1.1,
            composite_score: 0.9,
            achieved_continuity: crate::ContinuityLevel::G1,
        },
        d_composite_d_r: 10.0,
        d_r_min_d_r: 1.0,
        d_area_d_r: 10.0,
    };
    assert!(base.grows_with_radius());

    let mut flat = base;
    flat.d_composite_d_r = 0.0;
    assert!(!flat.grows_with_radius());

    let mut down = base;
    down.d_composite_d_r = -0.5;
    assert!(!down.grows_with_radius());
}
