//! Exact-boundary export-readiness repair without broad support-face demotion.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::{pcurve_measured_repair, surface_driven_edge_repair};

pub(super) fn repair_exact_boundary_export_blockers(
    result: &mut BRep,
    repair_report: &mut crate::types::BlendRepairReport,
) {
    let exact_boundary_committed = repair_report.apex_routes.iter().any(|route| {
        route.exact_boundary_report.candidate_status
            == crate::types::ExactBoundaryCandidateOutcome::Accepted
    });
    if !exact_boundary_committed {
        return;
    }
    if let Ok(repaired) = popkern_core::topology::repair_brep_boundary_edges_transactionally(
        result,
        &[],
        exact_boundary_export_repair_options(),
    ) {
        *result = repaired.brep;
        result.repair_line_edge_curve_endpoint_coherence(1.0e-6);
        if let Ok(edge_curves_completed) =
            popkern_core::topology::complete_missing_edge_curves_transactionally(
                result,
                exact_boundary_edge_curve_completion_options(),
            )
        {
            *result = edge_curves_completed.brep;
        }
        if let Ok(trim_repaired) =
            popkern_core::topology::reconstruct_curve_backed_trims_transactionally(
                result,
                exact_boundary_post_surface_trim_options(),
            )
        {
            *result = trim_repaired.brep;
        }
    }
    result.fix_stale_pcurves();
    pcurve_measured_repair::refit_measured_pcurve_blockers(result, 1.0e-6);
    pcurve_measured_repair::repair_measured_plane_cylinder_pcurves(result, 1.0e-6);
    surface_driven_edge_repair::repair_measured_edges_from_single_surface_pcurves(
        result, 1.0e-6, 0.06,
    );
    pcurve_measured_repair::refit_measured_pcurve_blockers(result, 1.0e-6);
    record_exact_boundary_support_scope_declines(result, repair_report, 1.0e-6);
}

fn record_exact_boundary_support_scope_declines(
    result: &BRep,
    repair_report: &mut crate::types::BlendRepairReport,
    tolerance: f64,
) {
    let mut apex_diagnostics = Vec::new();
    for route in &mut repair_report.apex_routes {
        if route.exact_boundary_report.candidate_status
            != crate::types::ExactBoundaryCandidateOutcome::Accepted
        {
            continue;
        }
        for sector in route.exact_boundary_report.sectors.clone() {
            for (edge, support_half_edge, support_face) in [
                (
                    sector.fillet_end_edge,
                    sector.support_half_edges[0],
                    sector.fillet_support_face,
                ),
                (
                    sector.face_edge,
                    sector.support_half_edges[1],
                    sector.bridge_support_face,
                ),
            ] {
                let (code, residual) = exact_boundary_support_scope_decline(
                    result,
                    support_half_edge,
                    support_face,
                    tolerance,
                );
                let Some(code) = code else {
                    continue;
                };
                let diagnostic = crate::types::ExactBoundaryCapDiagnostic {
                    stage: "exact_boundary_cap.support_scope".to_string(),
                    code,
                    sector_name: sector.sector_name.clone(),
                    edge: Some(edge),
                    support_half_edge: Some(support_half_edge),
                    support_face: Some(support_face),
                    residual,
                    shared_fact_ids: Vec::new(),
                };
                route.exact_boundary_report.diagnostics.push(
                    crate::types::ExactBoundaryReportDiagnostic {
                        schema_version: 1,
                        stage: diagnostic.stage.clone(),
                        code: diagnostic.code.clone(),
                        sector_name: diagnostic.sector_name.clone(),
                        edge: diagnostic.edge,
                        support_half_edge: diagnostic.support_half_edge,
                        support_face: diagnostic.support_face,
                        residual: diagnostic.residual,
                        tolerance: Some(tolerance),
                        shared_fact_ids: Vec::new(),
                    },
                );
                apex_diagnostics.push(crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                    diagnostic,
                ));
            }
        }
    }
    repair_report.apex_diagnostics.extend(apex_diagnostics);
}

fn exact_boundary_support_scope_decline(
    result: &BRep,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    tolerance: f64,
) -> (Option<String>, Option<f64>) {
    let Some(half_edge) = result.half_edges.get(support_half_edge) else {
        return (
            Some("exact_boundary_cap.support_scope.topology_missing".to_string()),
            None,
        );
    };
    let Some(loop_) = result.loops.get(half_edge.loop_) else {
        return (
            Some("exact_boundary_cap.support_scope.topology_missing".to_string()),
            None,
        );
    };
    if loop_.face != support_face {
        return (
            Some("exact_boundary_cap.support_scope.ownership_mismatch".to_string()),
            None,
        );
    }
    if half_edge.pcurve.is_none() {
        return (
            Some("exact_boundary_cap.support_scope.pcurve_missing".to_string()),
            None,
        );
    }
    let Some(report) = result.measure_half_edge_pcurve_deviation(support_half_edge) else {
        return (
            Some("exact_boundary_cap.support_scope.measurement_unavailable".to_string()),
            None,
        );
    };
    if report.max_deviation > tolerance {
        return (
            Some("exact_boundary_cap.support_scope.residual_exceeded".to_string()),
            Some(report.max_deviation),
        );
    }
    (None, None)
}

fn exact_boundary_edge_curve_completion_options()
-> popkern_core::topology::BrepEdgeCurveCompletionOptions {
    popkern_core::topology::BrepEdgeCurveCompletionOptions {
        tolerance: 1e-6,
        samples: 9,
        allow_vertex_segments_without_surface: true,
        require_all_completed: false,
        validate_output: false,
    }
}

fn exact_boundary_export_repair_options() -> popkern_core::topology::BrepRepairPipelineOptions {
    let clean_sewing = popkern_core::topology::BrepSewingOptions {
        tolerance: 1e-6,
        samples: 5,
        require_surface_agreement: false,
        validate_output: false,
        allow_nonmanifold_edges: false,
    };
    popkern_core::topology::BrepRepairPipelineOptions {
        classification: popkern_core::topology::BrepSewingClassificationOptions {
            exact_tolerance: 1e-9,
            snap_tolerance: 1e-6,
            samples: 5,
            require_surface_agreement: false,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: popkern_core::topology::BrepSnapSewingOptions {
            max_vertex_motion: 1e-6,
            clean_sewing,
        },
        trim_repair: popkern_core::topology::BrepTrimRepairOptions {
            endpoint_tolerance: 1e-6,
            sampled_tolerance: 1e-6,
            samples: 9,
            validate_output: false,
        },
        edge_curve_completion: popkern_core::topology::BrepEdgeCurveCompletionOptions {
            tolerance: 1e-6,
            samples: 9,
            allow_vertex_segments_without_surface: true,
            require_all_completed: false,
            validate_output: false,
        },
        curve_trim_reconstruction: popkern_core::topology::BrepCurveTrimReconstructionOptions {
            tolerance: 1e-6,
            samples: 9,
            replace_existing: false,
            replace_stale_existing: true,
            expand_bounded_surface_domains: true,
            require_all_reconstructed: false,
            validate_output: false,
        },
        edge_surface_compatibility: popkern_core::topology::BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.06,
            samples: 9,
            require_all_compatible: false,
            validate_output: false,
        },
        face_use_trim_authority: popkern_core::topology::BrepFaceUseTrimAuthorityOptions {
            tolerance: 1e-6,
            samples: 9,
        },
        degenerate_edges: popkern_core::topology::BrepDegenerateEdgeOptions {
            tolerance: 1e-6,
            samples: 9,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: false,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: true,
        repair_face_use_trim_authority: true,
        repair_trims: false,
        diagnose_loop_shell: false,
        topology_output: popkern_core::topology::BrepTopologyOutputPolicy::OpenSheet,
    }
}

fn exact_boundary_post_surface_trim_options()
-> popkern_core::topology::BrepCurveTrimReconstructionOptions {
    popkern_core::topology::BrepCurveTrimReconstructionOptions {
        tolerance: 1e-6,
        samples: 9,
        replace_existing: true,
        replace_stale_existing: true,
        expand_bounded_surface_domains: true,
        require_all_reconstructed: false,
        validate_output: false,
    }
}
