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
        "raw exact-boundary cap-family export must have no topology blockers: topology_evidence={:#?} boundary_evidence={:#?} report={pre_export:#?}",
        pre_export.topology_evidence, pre_export.boundary_edge_evidence
    );
    // Census provenance: 21 pcurve-blocked EDGES when authored
    // (21b95c60); c3eb03ea (PK-095 tessellation hardening —
    // pcurve_recompute + setback topology) repaired the six
    // setback-cluster edges (40v1–42v1, 48v1–50v1). Verified as a
    // strict subset drop: every current blocker existed in the old
    // set, none were added (PK-116 ledger holds the edge-set diff).
    // Dropped 15 -> 14 with PK-124's exact-boundary seam fidelity
    // package (canonical shared-trim curve prefers the measured
    // on-support cap curve; the support pcurve is reprojected from it
    // with periodic unwrapping), which reconciled one more
    // cap/cylinder seam pcurve.
    // Dropped 14 -> 9 with PK-125 (ribbon cross-derivative
    // orientation continuity + corner-compatible spokes): five cap
    // seam pcurve residuals resolved outright.
    // Note the two counters use different granularity by design:
    // pcurve_blockers counts EDGES, edge_surface_blockers counts
    // incompatible HALF-EDGES (unchanged by the pcurve repair).
    assert_eq!(
        pre_export.pcurve_blockers,
        9,
        "raw exact-boundary cap-family export should expose the current structured pcurve blocker class: route={:#?} blocker_roles={:#?} pre_export={pre_export:#?}",
        route.exact_boundary_report,
        blocker_roles(&blend.brep, pre_export)
    );
    // Census provenance: 21 through the PK-095 era; dropped to 19 when
    // the pcurve recompute's degree-1 shortcut started verifying
    // parameterization at the quarter points (PK-074 closure) — two
    // previously fraction-skewed linear pcurves now take the dense
    // interpolated path and reconcile with their edge curves. Dropped
    // 19 -> 11 with PK-124: measured-gated transactional relock (no
    // more destructive full-row rewrites) plus the canonical
    // shared-trim curve + reprojected support pcurves reconciled eight
    // more cap-seam half-edges with their edge curves. 11 -> 9 with
    // PK-125 (ribbon orientation continuity + corner-compatible
    // spokes).
    assert_eq!(
        pre_export.edge_surface_blockers, 9,
        "raw exact-boundary cap-family export should expose the current structured edge/surface blocker class: {pre_export:#?}"
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
                let imported_validation = failure
                    .imported_validation
                    .as_ref()
                    .expect("STEP import decline should identify the imported validation gate");
                match &imported_validation.kind {
                    StepImportedValidationFailureKind::EulerPoincare(euler) => {
                        // Census provenance: authored at actual=23 with
                        // boundary_edges=0 (a closed accounting failure).
                        // The import class drifted in the same PK-095-era
                        // window as the blocker counts above: the imported
                        // model carried 38 unpaired half-edges (an
                        // OPEN-shell decline, actual=36, V=80 E=63 F=19).
                        // PK-124's seam-fidelity package (canonical
                        // shared-trim curves + reprojected support pcurves
                        // + measured-gated relock) paired more imported
                        // half-edges: actual=28, boundary_edges=32. Still a
                        // structured, deterministic decline of the
                        // known-blocked raw export (PK-116).
                        assert_eq!(euler.expected, 2);
                        assert_eq!(euler.actual, 28);
                        assert_eq!(
                            euler.boundary_edges, 32,
                            "Euler-Poincare evidence should expose the imported open-shell census: {euler:#?}"
                        );
                        assert_eq!(
                            euler.solid_shell_faces, euler.faces,
                            "the failing solid shell should own the imported faces: {euler:#?}"
                        );
                    }
                    StepImportedValidationFailureKind::ShellGap(shell_gap) => {
                        assert!(shell_gap.gap > 0.0);
                        assert!(
                            shell_gap.natural_gap > Tolerance::default().linear()
                                || shell_gap.swapped_gap > Tolerance::default().linear(),
                            "shell-gap evidence should include endpoint disagreement: {shell_gap:#?}"
                        );
                    }
                    StepImportedValidationFailureKind::OpenLoop(open_loop) => {
                        assert!(open_loop.gap > Tolerance::default().linear());
                        assert!(
                            open_loop.face.is_some(),
                            "open-loop evidence should identify the owning face: {open_loop:#?}"
                        );
                        assert!(
                            open_loop.half_edges > 0,
                            "open-loop evidence should include loop membership: {open_loop:#?}"
                        );
                    }
                    other => panic!(
                        "remaining STEP import decline should expose structured topology evidence: {other:#?}"
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
