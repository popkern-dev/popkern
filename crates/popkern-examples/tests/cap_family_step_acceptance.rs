use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::step_validation::{
    StepImportedValidationFailureKind, StepRoundTripValidationOptions, StepStrictAcceptanceStatus,
    validate_strict_step_acceptance,
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
        std::slice::from_ref(&blend.brep),
        &StepRoundTripValidationOptions::default(),
    );

    assert_eq!(report.label, "direct_apex_gregory_cap");
    assert_eq!(report.pre_export_body_reports.len(), 1);
    let pre_export = &report.pre_export_body_reports[0];
    assert_eq!(
        pre_export.topology_blockers, 0,
        "raw exact-boundary cap-family export must have no topology blockers: {pre_export:#?}"
    );
    assert_eq!(
        pre_export.pcurve_blockers,
        0,
        "raw exact-boundary cap-family export must have no pcurve blockers: route={:#?} blocker_roles={:#?} pre_export={pre_export:#?}",
        route.exact_boundary_report,
        blocker_roles(&blend.brep, pre_export)
    );
    assert_eq!(
        pre_export.edge_surface_blockers, 0,
        "raw exact-boundary cap-family export must have no edge/surface blockers: {pre_export:#?}"
    );
    if pre_export.edge_curve_blockers > 0 {
        assert!(
            !pre_export.edge_curve_evidence.is_empty(),
            "raw exact-boundary cap-family missing-curve blockers must expose per-edge evidence: {pre_export:#?}"
        );
    }
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
            if reason.contains("export-quality blocked") {
                let failure = report
                    .failure
                    .as_ref()
                    .expect("STEP export-quality decline should carry structured failure evidence");
                let imported_quality = failure
                    .export_quality
                    .as_ref()
                    .expect("STEP export-quality decline should include body evidence");
                assert!(
                    imported_quality.trim_blockers > 0
                        || imported_quality.edge_curve_blockers > 0
                        || imported_quality.pcurve_blockers > 0
                        || imported_quality.edge_surface_blockers > 0,
                    "imported export-quality decline should expose blockers: {imported_quality:#?}"
                );
                assert!(
                    !imported_quality.trim_evidence.is_empty()
                        || !imported_quality.edge_curve_evidence.is_empty()
                        || !imported_quality.pcurve_evidence.is_empty()
                        || !imported_quality.edge_surface_evidence.is_empty(),
                    "imported export-quality decline should include per-entity evidence: {imported_quality:#?}"
                );
                assert!(
                    failure.imported_validation.is_none(),
                    "line endpoint repair should move the old shell-gap failure out of import validation"
                );
            } else {
                let failure = report
                    .failure
                    .as_ref()
                    .expect("STEP import decline should carry structured failure evidence");
                assert_eq!(pre_export.pcurve_blockers, 0);
                assert_eq!(pre_export.edge_surface_blockers, 0);
                let imported_validation = failure
                    .imported_validation
                    .as_ref()
                    .expect("STEP import decline should identify the imported validation gate");
                match &imported_validation.kind {
                    StepImportedValidationFailureKind::ShellGap(shell_gap) => {
                        assert!(shell_gap.gap > 0.0);
                        assert!(
                            shell_gap.natural_gap > Tolerance::default().linear()
                                || shell_gap.swapped_gap > Tolerance::default().linear(),
                            "shell-gap evidence should include endpoint disagreement: {shell_gap:#?}"
                        );
                    }
                    other => panic!(
                        "remaining STEP import decline should expose shell-gap evidence: {other:#?}"
                    ),
                }
            }
        }
    }
}

fn blocker_roles(
    brep: &BRep,
    report: &popkern_examples::step_validation::StepRoundTripBodyReport,
) -> Vec<String> {
    report
        .pcurve_evidence
        .iter()
        .flat_map(|edge| {
            edge.half_edges.iter().map(move |half_edge| {
                let edge_role = brep
                    .edges
                    .get(edge.edge)
                    .map(|edge| format!("{:?}", edge.role))
                    .unwrap_or_else(|| "missing_edge".to_string());
                let face_role = brep
                    .faces
                    .get(half_edge.face)
                    .map(|face| format!("{:?}", face.role))
                    .unwrap_or_else(|| "missing_face".to_string());
                let face_surface = brep
                    .faces
                    .get(half_edge.face)
                    .and_then(|face| face.surface.as_ref())
                    .map(|surface| surface_kind(&surface.geometry))
                    .unwrap_or("missing_surface");
                let pcurve_kind = brep
                    .half_edges
                    .get(half_edge.half_edge)
                    .and_then(|half_edge| half_edge.pcurve.as_ref())
                    .map(pcurve_kind)
                    .unwrap_or("missing_pcurve");
                let curve_kind = brep
                    .edges
                    .get(edge.edge)
                    .and_then(|edge| edge.curve.as_ref())
                    .map(|curve| edge_curve_kind(&curve.geometry))
                    .unwrap_or("missing_curve");
                format!(
                    "edge={:?} edge_role={} curve={} half_edge={:?} face={:?} face_role={} surface={} pcurve={} max={:.6e} same={:.6e} trace={:.6e}",
                    edge.edge,
                    edge_role,
                    curve_kind,
                    half_edge.half_edge,
                    half_edge.face,
                    face_role,
                    face_surface,
                    pcurve_kind,
                    half_edge.max_deviation,
                    half_edge.same_fraction_deviation,
                    half_edge.trace_deviation
                )
            })
        })
        .collect()
}

fn edge_curve_kind(curve: &popkern_core::geometry::curve::CurveGeometry) -> &'static str {
    match curve {
        popkern_core::geometry::curve::CurveGeometry::Line { .. } => "line",
        popkern_core::geometry::curve::CurveGeometry::Circle { .. } => "circle",
        popkern_core::geometry::curve::CurveGeometry::Arc { .. } => "arc",
        popkern_core::geometry::curve::CurveGeometry::Ellipse { .. } => "ellipse",
        popkern_core::geometry::curve::CurveGeometry::Parabola { .. } => "parabola",
        popkern_core::geometry::curve::CurveGeometry::Hyperbola { .. } => "hyperbola",
        popkern_core::geometry::curve::CurveGeometry::Nurbs(_) => "nurbs",
    }
}

fn pcurve_kind(pcurve: &popkern_core::geometry::pcurve::Pcurve) -> &'static str {
    match &pcurve.geometry {
        popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => "linear",
        popkern_core::geometry::pcurve::PcurveGeometry::Elliptical { .. } => "elliptical",
        popkern_core::geometry::pcurve::PcurveGeometry::Sinusoidal { .. } => "sinusoidal",
        popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(_) => "nurbs",
        popkern_core::geometry::pcurve::PcurveGeometry::CylinderArcViaRationalQuad { .. } => {
            "cylinder_arc_via_rational_quad"
        }
    }
}

fn surface_kind(surface: &popkern_core::geometry::SurfaceGeometry) -> &'static str {
    match surface {
        popkern_core::geometry::SurfaceGeometry::Plane { .. } => "plane",
        popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "cylinder",
        popkern_core::geometry::SurfaceGeometry::Sphere { .. } => "sphere",
        popkern_core::geometry::SurfaceGeometry::Cone { .. } => "cone",
        popkern_core::geometry::SurfaceGeometry::Torus { .. } => "torus",
        popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "nurbs",
        popkern_core::geometry::SurfaceGeometry::Gregory(_) => "gregory",
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
