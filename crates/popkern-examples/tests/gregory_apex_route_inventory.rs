use popkern_blend::{
    ApexConstructionRoute, ExactBoundaryCandidateOutcome, ExactBoundaryReportTransactionStatus,
    FilletOptions, FilletRadius, FilletSpec, build_blend_repair_action_artifact, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
    analyze_brep_edge_surface_compatibility, analyze_brep_export_quality,
};

#[test]
fn mixed_radius_gregory_records_default_direct_apex_route() {
    let result = mixed_radius_gregory_cube(FilletOptions::default());
    let routes = &result.repair_report.apex_routes;

    assert_eq!(routes.len(), 1, "fixture should have one smooth apex route");
    assert_eq!(routes[0].route, ApexConstructionRoute::DirectApex);
    assert_eq!(routes[0].incident_edges.len(), 3);
    assert_eq!(routes[0].built_faces.len(), 3);
    assert_eq!(routes[0].exact_boundary_sector_count, 3);
    assert_eq!(routes[0].exact_boundary_prepared_curve_count, 6);
    assert_eq!(routes[0].exact_boundary_inner_spoke_count, 6);
    assert_eq!(routes[0].exact_boundary_ribbon_count, 6);
    assert_eq!(routes[0].exact_boundary_solved_surface_count, 3);
    assert_eq!(routes[0].exact_boundary_pcurve_pair_count, 6);
    assert_exact_boundary_performance_report(&routes[0]);
    assert_exact_boundary_ux_report(&routes[0]);
    assert_exact_boundary_accepted(&result);
    assert_exact_boundary_artifact_has_no_decline_actions(&result);
    assert!(routes[0].fallback_reason.is_none());
}

#[test]
fn mixed_radius_gregory_records_direct_apex_requested_route() {
    let result = mixed_radius_gregory_cube(FilletOptions {
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    });
    let routes = &result.repair_report.apex_routes;

    assert_eq!(routes.len(), 1, "fixture should have one smooth apex route");
    assert_eq!(routes[0].route, ApexConstructionRoute::DirectApex);
    assert_eq!(routes[0].incident_edges.len(), 3);
    assert_eq!(routes[0].built_faces.len(), 3);
    assert_eq!(routes[0].exact_boundary_sector_count, 3);
    assert_eq!(routes[0].exact_boundary_prepared_curve_count, 6);
    assert_eq!(routes[0].exact_boundary_inner_spoke_count, 6);
    assert_eq!(routes[0].exact_boundary_ribbon_count, 6);
    assert_eq!(routes[0].exact_boundary_solved_surface_count, 3);
    assert_eq!(routes[0].exact_boundary_pcurve_pair_count, 6);
    assert_exact_boundary_performance_report(&routes[0]);
    assert_exact_boundary_ux_report(&routes[0]);
    assert_exact_boundary_accepted(&result);
    assert_exact_boundary_artifact_has_no_decline_actions(&result);
    assert!(routes[0].fallback_reason.is_none());
}

fn assert_exact_boundary_performance_report(route: &popkern_blend::ApexConstructionRouteDecision) {
    let performance = &route.exact_boundary_performance;
    assert_eq!(performance.sector_count, route.exact_boundary_sector_count);
    assert_eq!(
        performance.prepared_curve_count,
        route.exact_boundary_prepared_curve_count
    );
    assert_eq!(
        performance.inner_spoke_count,
        route.exact_boundary_inner_spoke_count
    );
    assert_eq!(
        performance.support_ribbon_count,
        route.exact_boundary_ribbon_count
    );
    assert_eq!(
        performance.solved_surface_count,
        route.exact_boundary_solved_surface_count
    );
    assert_eq!(
        performance.pcurve_pair_count,
        route.exact_boundary_pcurve_pair_count
    );
    assert!(
        performance.elapsed_micros > 0,
        "exact-boundary performance report should carry elapsed budget evidence"
    );
    assert!(performance.deterministic_skip_code.is_none());
}

fn assert_exact_boundary_ux_report(route: &popkern_blend::ApexConstructionRouteDecision) {
    let report = &route.exact_boundary_report;
    assert_eq!(report.schema_version, 1);
    assert_eq!(report.route_label, "direct_apex");
    assert_eq!(
        report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted
    );
    assert_eq!(
        report.provenance.transaction_status,
        ExactBoundaryReportTransactionStatus::Committed
    );
    assert!(!report.provenance.topology_rollback_performed);
    assert_eq!(
        report.committed_faces.len(),
        route.exact_boundary_sector_count
    );
    assert_eq!(report.sectors.len(), route.exact_boundary_sector_count);
    assert_eq!(report.performance, route.exact_boundary_performance);
    assert!(report.preview_available);
    assert_eq!(
        report.tolerance_evidence.deterministic_skip_code,
        route.exact_boundary_performance.deterministic_skip_code
    );
    assert!(report.tolerance_evidence.model_linear_tolerance.is_some());
    // Residual evidence is recorded per EXCEEDANCE, not per sector:
    // through the PK-095 era every sector exceeded a residual gate so
    // the count floored at sector_count; the PK-074 closure left a
    // single support-scope exceedance (2.1e-5 vs 1e-6 on
    // cap_sector.1). Re-pinned 2026-07-17 (PK-125): the ribbon
    // cross-derivative orientation-continuity fix + corner-compatible
    // spokes resolved the last exceedance — an accepted mixed-radius
    // route now records NO residual diagnostics.
    assert!(
        report
            .diagnostics
            .iter()
            .all(|d| popkern_blend::types::diagnostic_code_is_informational(&d.code)),
        "the mixed-radius fixture regressed a residual exceedance: {:#?}",
        report.diagnostics
    );
    assert!(report.diagnostics.iter().all(|diagnostic| {
        diagnostic.schema_version == 1
            && (popkern_blend::types::diagnostic_code_is_informational(&diagnostic.code)
                || (diagnostic.residual.is_some_and(|residual| residual > 0.0)
                    && ((diagnostic.stage == "exact_boundary_cap.support_scope"
                        && diagnostic.code
                            == "exact_boundary_cap.support_scope.residual_exceeded"
                        && diagnostic.tolerance == Some(1.0e-6))
                        || (diagnostic.stage == "exact_boundary_cap.cap_outer"
                            && diagnostic.code
                                == "exact_boundary_cap.cap_outer.pair_residual_exceeded"
                            && diagnostic.tolerance == Some(5.0e-3)))))
    }));
    assert!(
        report
            .diagnostics
            .iter()
            .all(|diagnostic| diagnostic.schema_version == 1)
    );
    assert!(report.source_fact_ids.is_empty());
    assert_eq!(report.sectors[0].sector_name, "direct_apex.cap_sector.0");
    assert_eq!(report.sectors[0].cap_half_edges.len(), 2);

    let json = serde_json::to_string(report).expect("exact-boundary report serializes to JSON");
    assert!(json.contains("\"schema_version\":1"));
    assert!(json.contains("\"candidate_status\":\"Accepted\""));
    assert!(json.contains("\"transaction_status\":\"Committed\""));
    assert!(json.contains("\"source_fact_ids\""));
}

fn assert_exact_boundary_accepted(result: &popkern_blend::BlendResult) {
    assert!(
        result.repair_report.apex_routes.iter().all(|route| route
            .exact_boundary_report
            .diagnostics
            .iter()
            .all(|diagnostic| {
                popkern_blend::types::diagnostic_code_is_informational(&diagnostic.code)
                    || (diagnostic.stage == "exact_boundary_cap.support_scope"
                        && diagnostic.code == "exact_boundary_cap.support_scope.residual_exceeded")
                    || (diagnostic.stage == "exact_boundary_cap.cap_outer"
                        && diagnostic.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded")
            })),
        "accepted exact-boundary mixed-radius Gregory should only promote structured cap/support diagnostics"
    );
}

fn assert_exact_boundary_artifact_has_no_decline_actions(result: &popkern_blend::BlendResult) {
    let artifact = build_blend_repair_action_artifact(&result.repair_report);
    assert!(
        artifact
            .actions
            .iter()
            .all(|action| !action.stage.starts_with("exact_boundary_cap.")
                || action.stage == "exact_boundary_cap.route_decision"
                || action.stage == "exact_boundary_cap.support_scope"
                || action.stage == "exact_boundary_cap.cap_outer"),
        "accepted exact-boundary candidates should not emit declined construction actions: {artifact:#?}"
    );
}

#[test]
fn mixed_radius_gregory_default_matches_requested_direct_apex_quality() {
    let default = mixed_radius_gregory_cube(FilletOptions::default());
    let requested_direct_apex = mixed_radius_gregory_cube(FilletOptions {
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    });
    let default_quality = quality(&default.brep);
    let requested_direct_apex_quality = quality(&requested_direct_apex.brep);
    let _evidence_rows =
        default_quality.evidence_rows() + requested_direct_apex_quality.evidence_rows();

    assert_eq!(
        default.repair_report.apex_routes[0].route,
        ApexConstructionRoute::DirectApex
    );
    assert_eq!(
        requested_direct_apex.repair_report.apex_routes[0].route,
        ApexConstructionRoute::DirectApex
    );
    assert!(
        default_quality.pcurve_blockers == requested_direct_apex_quality.pcurve_blockers
            && default_quality.edge_surface_blockers
                == requested_direct_apex_quality.edge_surface_blockers
            && default_quality.max_pcurve_deviation
                == requested_direct_apex_quality.max_pcurve_deviation
            && default_quality.max_edge_surface_deviation
                == requested_direct_apex_quality.max_edge_surface_deviation,
        "default mixed-radius Gregory should use the same STEP-fidelity path as requested DirectApex: default={default_quality:?}, requested_direct_apex={requested_direct_apex_quality:?}"
    );
}

impl RouteQuality {
    fn evidence_rows(&self) -> usize {
        let pcurve_rows = self
            .worst_pcurve_edges
            .iter()
            .map(|edge| {
                edge.edge.len()
                    + (edge.max_deviation > 0.0) as usize
                    + edge
                        .half_edges
                        .iter()
                        .map(|half_edge| {
                            half_edge.half_edge.len()
                                + half_edge.face.len()
                                + (half_edge.max_deviation > 0.0) as usize
                        })
                        .sum::<usize>()
            })
            .sum::<usize>();
        let edge_surface_rows = self
            .worst_edge_surface_edges
            .iter()
            .map(|edge| {
                edge.edge.len()
                    + edge.half_edge.len()
                    + edge.face.len()
                    + (edge.max_deviation > 0.0) as usize
            })
            .sum::<usize>();
        pcurve_rows + edge_surface_rows
    }
}

fn mixed_radius_gregory_cube(options: FilletOptions) -> popkern_blend::BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges = cube
        .half_edges
        .values()
        .filter(|half_edge| half_edge.origin == vertex)
        .map(|half_edge| half_edge.edge)
        .collect::<Vec<_>>();
    let radii = [0.08, 0.10, 0.12];
    let specs = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    fillet_edges(&cube, &specs, &options, Tolerance::default(), None)
        .expect("mixed-radius Gregory cube must build")
}

#[derive(Debug)]
struct RouteQuality {
    pcurve_blockers: usize,
    max_pcurve_deviation: f64,
    edge_surface_blockers: usize,
    max_edge_surface_deviation: f64,
    worst_pcurve_edges: Vec<EdgePcurveEvidence>,
    worst_edge_surface_edges: Vec<EdgeSurfaceEvidence>,
}

#[derive(Debug)]
struct EdgePcurveEvidence {
    edge: String,
    max_deviation: f64,
    half_edges: Vec<HalfEdgePcurveEvidence>,
}

#[derive(Debug)]
struct HalfEdgePcurveEvidence {
    half_edge: String,
    face: String,
    max_deviation: f64,
}

#[derive(Debug)]
struct EdgeSurfaceEvidence {
    edge: String,
    half_edge: String,
    face: String,
    max_deviation: f64,
}

fn quality(brep: &BRep) -> RouteQuality {
    let export = analyze_brep_export_quality(
        brep,
        popkern_core::topology::BrepExportQualityOptions {
            topology_output: BrepTopologyOutputPolicy::OpenSheet,
            trim: BrepTrimRepairOptions {
                endpoint_tolerance: 1e-6,
                sampled_tolerance: 1e-6,
                samples: 9,
                validate_output: true,
            },
            require_edge_curves: true,
            allow_collapsed_trim_edge_curves: true,
            max_pcurve_deviation: Some(1e-6),
            max_edge_surface_deviation: Some(1e-6),
            validate_brep: true,
        },
    );
    let edge_surface = analyze_brep_edge_surface_compatibility(
        brep,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.0,
            samples: 9,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    RouteQuality {
        pcurve_blockers: export.pcurve_blockers,
        max_pcurve_deviation: export.max_pcurve_deviation,
        edge_surface_blockers: edge_surface.incompatible_half_edges,
        max_edge_surface_deviation: edge_surface.max_deviation,
        worst_pcurve_edges: worst_pcurve_edges(brep, 30),
        worst_edge_surface_edges: edge_surface
            .findings
            .iter()
            .take(30)
            .map(|finding| EdgeSurfaceEvidence {
                edge: format!("{:?}", finding.edge),
                half_edge: format!("{:?}", finding.half_edge),
                face: format!("{:?}", finding.face),
                max_deviation: finding.max_deviation,
            })
            .collect(),
    }
}

fn worst_pcurve_edges(brep: &BRep, limit: usize) -> Vec<EdgePcurveEvidence> {
    let mut reports = brep
        .edges
        .keys()
        .filter_map(|edge| brep.measure_edge_pcurve_deviation(edge))
        .filter(|report| report.max_deviation > 1e-12)
        .collect::<Vec<_>>();
    reports.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    reports
        .into_iter()
        .take(limit)
        .map(|report| EdgePcurveEvidence {
            edge: format!("{:?}", report.edge),
            max_deviation: report.max_deviation,
            half_edges: report
                .half_edges
                .iter()
                .map(|half_edge| HalfEdgePcurveEvidence {
                    half_edge: format!("{:?}", half_edge.half_edge),
                    face: format!("{:?}", half_edge.face),
                    max_deviation: half_edge.max_deviation,
                })
                .collect(),
        })
        .collect()
}
