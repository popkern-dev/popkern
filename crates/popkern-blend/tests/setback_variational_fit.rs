//! PK-574: the variational corner fit serves beneath the KP20 setback
//! cap. Row 877 (the class KP20 folds on), row 876 (every KP20-served
//! corner with the KP20 cap forced to decline), row 878 (structured
//! declines: `setback_variational_fit_controls.rs`).

#[path = "support/setback_fit_fixtures.rs"]
mod setback_fit_fixtures;

use popkern_blend::dispatch_overrides::with_kp20_cap_forced_to_decline;
use popkern_blend::{
    BlendResult, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, SetbackFitMode,
    SetbackFitOutcome, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use setback_fit_fixtures::*;

#[test]
fn the_fit_serves_the_cube_corner_kp20_folds_on() {
    // Row 877: the C1 census class (KP20 folds on the wide mixed-radius
    // cube corner) is served by the fit. At 0.05 / 0.2 / 0.3 the
    // curvature-weighted fit certifies G2 but its skinny, sheared sector
    // beside the thin strip bends to 30 / r_min (chart sine 0.12), so
    // the fairer tangent fit (2.2) serves at G1 (PK-614).
    let result = cube_corner(&[0.05, 0.2, 0.3]);
    let record = only_record("cube 0.05/0.2/0.3", &result);
    assert!(
        record.kp20_decline.contains("fold"),
        "{}",
        record.kp20_decline
    );
    assert_fit_serves(
        "cube 0.05/0.2/0.3",
        &result,
        SetbackFitMode::Tangent,
        ContinuityLevel::G1,
    );
    let first = &record.attempts[0];
    assert_eq!(first.mode, SetbackFitMode::CurvatureWeighted);
    assert!(
        first.interior.is_some_and(|k| k > 4.0)
            && matches!(&first.outcome, SetbackFitOutcome::Rejected { certified, reason }
                if certified.level == ContinuityLevel::G2 && reason.contains("fairer")),
        "{record:#?}"
    );
    assert!(
        record.attempts[1].interior.is_some_and(|k| k < 4.0),
        "{record:#?}"
    );
    // At 0.08 / 0.2 / 0.3 the fan bends to 468 / r_min and the G2 fit
    // serves fair (1.75).
    let result = cube_corner(&[0.08, 0.2, 0.3]);
    assert_fit_serves(
        "cube 0.08/0.2/0.3",
        &result,
        SetbackFitMode::CurvatureWeighted,
        ContinuityLevel::G2,
    );
}

#[test]
fn row_876_triangular_and_square_apexes_serve_with_kp20_forced_off() {
    let mut fixtures: Vec<Fixture> = Vec::new();
    for (name, n, h) in [
        ("pyramid3 h0.8", 3, 0.8),
        ("pyramid3 h1.2", 3, 1.2),
        ("pyramid3 h3", 3, 3.0),
        ("pyramid4 h0.8", 4, 0.8),
        ("pyramid4 h1.2", 4, 1.2),
        ("pyramid4 h3", 4, 3.0),
    ] {
        // PK-595: the corner-cell fit certifies G2 on every one (the
        // tangent fit served them at G1 before).
        fixtures.push((
            name,
            Box::new(move || pyramid_apex(n, h, &[0.16], ContinuityLevel::G2)),
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        ));
    }
    census(fixtures);
}

#[test]
fn row_876_wide_apexes_serve_with_kp20_forced_off() {
    // PK-595: the corner-cell fit certifies G2 on the whole rim family
    // (rim g2 bounds 1.9e-3..5.0e-2, spokes 7.9e-3..4.9e-2).
    let curvature_weighted = |name, build: Box<dyn Fn() -> BlendResult>| {
        (
            name,
            build,
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        )
    };
    census(vec![
        curvature_weighted(
            "pentagon flagship",
            Box::new(|| pyramid_apex(5, 1.2, &[0.16], ContinuityLevel::G2)),
        ),
        curvature_weighted(
            "pentagon5 blunt",
            Box::new(|| pyramid_apex(5, 3.0, &[0.05], ContinuityLevel::G2)),
        ),
        curvature_weighted(
            "hexagon6 blunt",
            Box::new(|| pyramid_apex(6, 3.0, &[0.05], ContinuityLevel::G2)),
        ),
        curvature_weighted(
            "hexagon6 mixed",
            Box::new(|| {
                pyramid_apex(
                    6,
                    0.8,
                    &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
                    ContinuityLevel::G2,
                )
            }),
        ),
        curvature_weighted(
            "pentagon5 variable",
            Box::new(|| pyramid_apex(5, 0.8, &[0.02, 0.04, 0.06, 0.04, 0.02], ContinuityLevel::G2)),
        ),
        curvature_weighted(
            "octagon8 mixed",
            Box::new(|| {
                pyramid_apex(
                    8,
                    0.8,
                    &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
                    ContinuityLevel::G2,
                )
            }),
        ),
    ]);
}

#[test]
fn row_876_cube_corners_and_g1_requests_serve_with_kp20_forced_off() {
    census(vec![
        (
            "cube corner 0.2",
            Box::new(|| cube_corner(&[0.2, 0.2, 0.2])),
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        ),
        (
            "cube corner 0.1/0.15/0.2",
            Box::new(|| cube_corner(&[0.1, 0.15, 0.2])),
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        ),
        (
            "pentagon G1 request",
            Box::new(|| pyramid_apex(5, 1.2, &[0.16], ContinuityLevel::G1)),
            SetbackFitMode::Tangent,
            ContinuityLevel::G1,
        ),
    ]);
}

/// Row 876 (E5): the saddle carve's bridge sides (degree-12 curves on
/// the reflex support plane) serve through the same fit — at certified
/// G2 where requested, the saddle's strips being circular.
#[test]
fn row_876_saddles_and_variable_radii_serve_with_kp20_forced_off() {
    let (body, saddle) = notched_saddle();
    let saddle_at = |radii: [f64; 3], level: ContinuityLevel| -> Box<dyn Fn() -> BlendResult> {
        let body = body.clone();
        Box::new(move || fillet(&body, saddle, &radii, level))
    };
    census(vec![
        (
            "saddle G1",
            saddle_at([0.06, 0.08, 0.10], ContinuityLevel::G1),
            SetbackFitMode::Tangent,
            ContinuityLevel::G1,
        ),
        (
            "saddle G2",
            saddle_at([0.06, 0.08, 0.10], ContinuityLevel::G2),
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        ),
        (
            "saddle uniform G2",
            saddle_at([0.08, 0.08, 0.08], ContinuityLevel::G2),
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        ),
        (
            "hexagon6 variable",
            Box::new(|| {
                pyramid_apex(
                    6,
                    0.8,
                    &[0.02, 0.04, 0.06, 0.04, 0.02, 0.03],
                    ContinuityLevel::G2,
                )
            }),
            SetbackFitMode::CurvatureWeighted,
            ContinuityLevel::G2,
        ),
    ]);
}

/// Row 876: pyramid bodies with every edge filleted — the apex and each
/// base corner is a setback corner, and every one is served.
#[test]
fn row_876_all_edge_pyramids_serve_every_corner_with_kp20_forced_off() {
    let tol = Tolerance::default();
    for (n, h) in [(3usize, 0.8), (4, 1.2), (5, 0.8)] {
        let name = format!("pyramid{n} h{h} all edges");
        let body = BRep::n_gonal_pyramid(n, 1.0, h, tol).expect("pyramid");
        let specs: Vec<FilletSpec> = body
            .edges
            .keys()
            .map(|edge| FilletSpec {
                edge,
                radius: FilletRadius::Constant(0.1),
            })
            .collect();
        let options = FilletOptions {
            continuity: ContinuityLevel::G2,
            ..Default::default()
        };
        let result =
            with_kp20_cap_forced_to_decline(|| fillet_edges(&body, &specs, &options, tol, None))
                .expect("fillet");
        let records = fit_records(&result);
        assert_eq!(records.len(), n + 1, "{name}: one record per corner");
        for record in records {
            // PK-595: every corner certifies G2 on the corner-cell fit.
            assert_record_serves(
                &name,
                record,
                SetbackFitMode::CurvatureWeighted,
                ContinuityLevel::G2,
            );
        }
        result
            .brep
            .validate()
            .unwrap_or_else(|e| panic!("{name}: {e}"));
    }
}

#[test]
fn the_fit_serves_the_extreme_cube_corner_ahead_of_the_exact_ring() {
    // Row 877, second C1 fixture: the 1:10 cube corner. KP20's cap folds;
    // before PK-574 the corner fell through to the exact-boundary ring
    // (G1 by measurement). The tangent fit now serves it, certified G1
    // (the curvature-weighted attempt is declined on its own rim
    // witness, a curvature jump ~80 against the 5e-2 gate).
    let result = cube_corner(&[0.02, 0.10, 0.20]);
    let record = only_record("cube 0.02/0.1/0.2", &result);
    assert!(
        record.kp20_decline.contains("fold"),
        "{}",
        record.kp20_decline
    );
    assert!(
        matches!(&record.attempts[0].outcome, SetbackFitOutcome::Declined(r) if r.contains("G2 cannot certify")),
        "{record:#?}"
    );
    assert_fit_serves(
        "cube 0.02/0.1/0.2",
        &result,
        SetbackFitMode::Tangent,
        ContinuityLevel::G1,
    );
}
