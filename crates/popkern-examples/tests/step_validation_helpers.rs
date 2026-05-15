use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::step_validation::{
    StepRoundTripValidationOptions, export_and_validate_breps,
};

#[test]
fn shared_step_validation_accepts_non_pns_cube_body() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let (_bytes, report) = export_and_validate_breps(
        "unit_cube",
        &[cube],
        &StepRoundTripValidationOptions::default(),
    )
    .unwrap();

    assert_eq!(report.imported_bodies, 1);
    assert_eq!(report.body_reports[0].faces, 6);
    assert_eq!(report.body_reports[0].edge_curve_blockers, 0);
    assert_eq!(report.body_reports[0].pcurve_blockers, 0);
    assert!(report.body_reports[0].max_pcurve_deviation < 1e-6);
}
