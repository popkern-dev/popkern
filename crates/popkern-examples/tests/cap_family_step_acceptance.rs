use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::step_validation::{
    StepRoundTripValidationOptions, StepStrictAcceptanceStatus, validate_strict_step_acceptance,
};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
fn direct_apex_gregory_cap_runs_shared_strict_step_acceptance() {
    let blend = mixed_radius_gregory_cube();
    let route = blend
        .repair_report
        .apex_routes
        .first()
        .expect("fixture should produce one apex route");
    assert_eq!(route.exact_boundary_report.schema_version, 1);

    let report = validate_strict_step_acceptance(
        "direct_apex_gregory_cap",
        &[blend.brep],
        &StepRoundTripValidationOptions::default(),
    );

    assert_eq!(report.label, "direct_apex_gregory_cap");
    assert_eq!(report.pre_export_body_reports.len(), 1);
    let pre_export = &report.pre_export_body_reports[0];
    match report.status {
        StepStrictAcceptanceStatus::Accepted => {
            let round_trip = report
                .round_trip_report
                .as_ref()
                .expect("accepted strict STEP report should include round-trip evidence");
            assert_eq!(round_trip.imported_bodies, 1);
            assert_eq!(round_trip.body_reports[0].pcurve_blockers, 0);
            assert_eq!(round_trip.body_reports[0].edge_surface_blockers, 0);
        }
        StepStrictAcceptanceStatus::Declined => {
            let reason = report
                .failure_reason
                .as_ref()
                .expect("declined strict STEP report should explain the gate");
            assert!(
                reason.contains("export-quality blocked") || reason.contains("STEP import"),
                "strict STEP decline should identify the failed shared gate: {reason}"
            );
            assert!(
                pre_export.pcurve_blockers > 0 || pre_export.edge_surface_blockers > 0,
                "declined cap-family strict STEP report should expose fidelity blockers: {pre_export:#?}"
            );
            assert!(
                !pre_export.pcurve_evidence.is_empty()
                    || !pre_export.edge_surface_evidence.is_empty(),
                "declined cap-family strict STEP report should include per-edge evidence: {pre_export:#?}"
            );
        }
    }
}

fn mixed_radius_gregory_cube() -> popkern_fillet::BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|half_edge| half_edge.origin == vertex)
        .map(|half_edge| half_edge.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<_> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let options = FilletOptions {
        use_rational_gregory: true,
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &options, Tolerance::default(), None)
        .expect("mixed-radius Gregory cube must build")
}
