//! PK-603: the fillet corner-cap preference. A stated construction
//! reorders the setback rung (fitted cap first, or the setback cap
//! with direct apex held off), `classic` skips the rung, and a strict
//! preference serves its construction or fails typed — never a silent
//! fallback. Design: docs/architecture/fillet-corner-cap-preference.md.

#[path = "support/setback_fit_fixtures.rs"]
mod setback_fit_fixtures;

use popkern_blend::{
    BlendError, BlendFailure, BlendResult, ContinuityLevel, CornerCapPreference,
    CornerConstruction, CornerMethod, CornerSubPatchFamily, FilletOptions, FilletRadius,
    FilletSpec, fillet_edges,
};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use setback_fit_fixtures::{edges_at, fit_records};

fn cube_corner_with(
    radii: &[f64],
    continuity: ContinuityLevel,
    construction: CornerConstruction,
    strict: bool,
) -> Result<BlendResult, BlendError> {
    cube_corner_with_retry(radii, continuity, construction, strict, true)
}

/// [`cube_corner_with`] with the serving retry on or off (off observes
/// the construction, before any exclusion re-run).
fn cube_corner_with_retry(
    radii: &[f64],
    continuity: ContinuityLevel,
    construction: CornerConstruction,
    strict: bool,
    serving_self_intersection_retry: bool,
) -> Result<BlendResult, BlendError> {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9)
        .map(|(h, _)| h)
        .expect("corner");
    let specs: Vec<_> = edges_at(&cube, corner)
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity,
        corner_cap: CornerCapPreference {
            construction,
            strict,
        },
        serving_self_intersection_retry,
        ..Default::default()
    };
    fillet_edges(&cube, &specs, &options, Tolerance::default(), None)
}

fn fitted_sectors(result: &BlendResult) -> usize {
    result
        .corner_patches
        .iter()
        .filter(|p| {
            p.sub_patch_provenance
                .as_ref()
                .is_some_and(|s| s.family == CornerSubPatchFamily::VariationalFit)
        })
        .count()
}

fn setback_kp20_patches(result: &BlendResult) -> usize {
    result
        .corner_patches
        .iter()
        .filter(|p| p.setback_cap && p.method != CornerMethod::PlateMin)
        .count()
}

fn strict_failure(error: BlendError) -> (CornerConstruction, String) {
    let BlendError::PartialFailure { failures, .. } = error else {
        panic!("a strict decline is a transactional failure, got {error:?}");
    };
    failures
        .iter()
        .find_map(|(_, f)| match f {
            BlendFailure::CornerCapUnavailable {
                construction,
                reason,
                location,
                ..
            } => {
                let at = location.expect("the corner is located");
                assert!(at.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9, "{at:?}");
                Some((*construction, reason.clone()))
            }
            _ => None,
        })
        .unwrap_or_else(|| panic!("no CornerCapUnavailable in {failures:?}"))
}

#[test]
fn fitted_serves_the_fitted_cap_where_auto_serves_the_setback_cap() {
    let auto = cube_corner_with(&[0.2], ContinuityLevel::G2, CornerConstruction::Auto, false)
        .expect("auto");
    assert!(
        setback_kp20_patches(&auto) >= 3,
        "{:?}",
        auto.corner_patches
    );
    assert_eq!(fitted_sectors(&auto), 0);

    let fitted = cube_corner_with(
        &[0.2],
        ContinuityLevel::G2,
        CornerConstruction::Fitted,
        false,
    )
    .expect("fitted");
    assert!(fitted_sectors(&fitted) >= 3, "{:?}", fitted.corner_patches);
    assert_eq!(setback_kp20_patches(&fitted), 0);
    let record = fit_records(&fitted);
    assert_eq!(record.len(), 1, "{record:#?}");
    assert!(
        record[0].kp20_decline.contains("not attempted"),
        "{record:#?}"
    );
    assert!(fitted.brep.validate().is_ok());
}

#[test]
fn a_strict_setback_cap_declines_typed_where_kp20_folds() {
    let error = cube_corner_with(
        &[0.05, 0.2, 0.3],
        ContinuityLevel::G2,
        CornerConstruction::Setback,
        true,
    )
    .expect_err("KP20 folds on the wide mixed corner; strict must not fall back");
    let (construction, reason) = strict_failure(error);
    assert_eq!(construction, CornerConstruction::Setback);
    assert!(reason.contains("fold"), "{reason}");

    // Non-strict: the same preference falls back to the fitted cap.
    let fallback = cube_corner_with(
        &[0.05, 0.2, 0.3],
        ContinuityLevel::G2,
        CornerConstruction::Setback,
        false,
    )
    .expect("setback falls back");
    assert!(
        fitted_sectors(&fallback) >= 3,
        "{:?}",
        fallback.corner_patches
    );
}

#[test]
fn a_strict_fitted_cap_serves_or_declines_but_never_takes_kp20() {
    let served = cube_corner_with(
        &[0.05, 0.2, 0.3],
        ContinuityLevel::G2,
        CornerConstruction::Fitted,
        true,
    )
    .expect("the fit serves the wide mixed corner");
    assert!(fitted_sectors(&served) >= 3);
    let record = fit_records(&served);
    assert!(record[0].kp20_decline.contains("strict"), "{record:#?}");
}

#[test]
fn classic_skips_the_setback_rung() {
    let classic = cube_corner_with(
        &[0.2],
        ContinuityLevel::G2,
        CornerConstruction::Classic,
        false,
    )
    .expect("classic");
    assert!(
        classic.corner_patches.iter().all(|p| !p.setback_cap),
        "{:?}",
        classic.corner_patches
    );
    assert!(fit_records(&classic).is_empty());
    assert!(
        !classic.corner_patches.is_empty(),
        "a classic cap still closes the corner"
    );
}

#[test]
fn a_stated_preference_holds_direct_apex_off_rung_corners() {
    // G1 mixed-radius cylindrical corners are a promoted direct-apex
    // class: under `auto` direct apex claims the corner before the
    // rung (observed at construction; its sector seams crease, PK-609,
    // so the serve then takes the exclusion re-run onto the rung). A
    // stated `fitted` preference reaches the rung without that re-run.
    let auto = cube_corner_with_retry(
        &[0.05, 0.2, 0.3],
        ContinuityLevel::G1,
        CornerConstruction::Auto,
        false,
        false,
    )
    .expect("auto");
    assert_eq!(
        fitted_sectors(&auto),
        0,
        "control: auto's construction is not the fit: {:?}",
        auto.corner_patches
    );
    assert!(
        auto.repair_report
            .apex_routes
            .iter()
            .any(|r| r.route == popkern_blend::ApexConstructionRoute::DirectApex),
        "control: direct apex claims the corner under auto: {:?}",
        auto.repair_report.apex_routes
    );
    let fitted = cube_corner_with(
        &[0.05, 0.2, 0.3],
        ContinuityLevel::G1,
        CornerConstruction::Fitted,
        false,
    )
    .expect("fitted");
    assert!(fitted_sectors(&fitted) >= 3, "{:?}", fitted.corner_patches);
    assert!(
        fitted
            .repair_report
            .serving_self_intersection_retries
            .is_empty(),
        "the preference reaches the rung directly: {:?}",
        fitted.repair_report.serving_self_intersection_retries
    );
    assert!(fitted.brep.validate().is_ok());
}
