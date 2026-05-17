use popkern_core::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepDegenerateEdgeOptions,
    BrepEdgeCurveCompletionOptions, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairReport, BrepEdgeSurfaceKind,
    BrepEdgeSurfaceVertexCandidateOutcome, BrepEdgeSurfaceVertexCandidateSource,
    BrepRepairPipelineOptions, BrepSewingClassificationOptions, BrepSewingOptions,
    BrepSnapSewingOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
    analyze_brep_edge_surface_compatibility, analyze_brep_export_quality,
    repair_brep_boundary_edges_transactionally,
};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq)]
struct FidelityState {
    label: &'static str,
    edges: usize,
    faces: usize,
    trim_findings: usize,
    pcurve_blockers: usize,
    max_pcurve_deviation: f64,
    edge_surface_blockers: usize,
    max_edge_surface_deviation: f64,
    worst_pcurve_edges: Vec<EdgePcurveEvidence>,
    worst_edge_surface_edges: Vec<EdgeSurfaceEvidence>,
}

#[derive(Debug, Clone, PartialEq)]
struct EdgePcurveEvidence {
    edge: String,
    max_deviation: f64,
    half_edges: Vec<HalfEdgePcurveEvidence>,
}

#[derive(Debug, Clone, PartialEq)]
struct HalfEdgePcurveEvidence {
    half_edge: String,
    face: String,
    max_deviation: f64,
    same_fraction_deviation: f64,
    trace_deviation: f64,
    worst_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct EdgeSurfaceEvidence {
    edge: String,
    half_edge: String,
    face: String,
    max_deviation: f64,
    max_curve_parameter: f64,
}

#[test]
fn gregory_step_fidelity_probe_classifies_blocker_source() {
    let construction = mixed_radius_gregory_cube();
    let repaired = repair_brep_boundary_edges_transactionally(
        &construction,
        &[],
        gregory_step_repair_options(false),
    )
    .expect("diagnostic repair pipeline should commit export candidate");
    if std::env::var_os("POPKERN_GREGORY_FIDELITY_TRACE").is_some() {
        eprintln!(
            "{:#?}",
            edge_surface_repair_summary(&repaired.report.edge_surface_compatibility)
        );
        eprintln!(
            "{:#?}",
            vertex_star_repair_summary(&repaired.report.edge_surface_compatibility)
        );
    }
    let export_input = repaired.brep;
    let bytes = export_step(&step_document(&export_input), &StepExportOptions::default())
        .expect("Gregory export candidate should serialize");
    let imported = import_step(
        &bytes,
        &StepImportOptions {
            healing: popkern_io::HealingLevel::None,
            ..Default::default()
        },
    )
    .expect("Gregory export candidate should import");
    assert!(
        imported.diagnostics.errors().next().is_none(),
        "STEP import diagnostics contain errors: {:?}",
        imported.diagnostics.entries()
    );
    let imported_brep = &imported.document.products[0].shapes[0].solids[0];

    let states = [
        fidelity_state("construction", &construction),
        fidelity_state("export_input", &export_input),
        fidelity_state("imported", imported_brep),
    ];
    if std::env::var_os("POPKERN_GREGORY_FIDELITY_TRACE").is_some() {
        for state in &states {
            eprintln!("{state:#?}");
        }
    }

    let construction_state = &states[0];
    let export_state = &states[1];
    let imported_state = &states[2];
    assert!(
        export_state.pcurve_blockers <= construction_state.pcurve_blockers,
        "repair should not increase construction pcurve blockers: {states:#?}"
    );
    assert!(
        export_state.edge_surface_blockers < construction_state.edge_surface_blockers,
        "repair should reduce construction edge/surface blockers: {states:#?}"
    );
    assert!(
        imported_state.pcurve_blockers > 0 || imported_state.edge_surface_blockers > 0,
        "probe should keep the remaining STEP fidelity blocker observable: {states:#?}"
    );

    let classification =
        classify_fidelity_blocker(construction_state, export_state, imported_state);
    assert_eq!(
        classification, "construction_geometry_issue",
        "unexpected Gregory STEP fidelity blocker classification; states: {states:#?}"
    );
}

fn vertex_star_repair_summary(
    report: &Option<BrepEdgeSurfaceCompatibilityRepairReport>,
) -> Vec<String> {
    let Some(report) = report else {
        return vec!["vertex_star_repair: not_run".to_string()];
    };
    report
        .vertex_stars
        .decisions
        .iter()
        .filter(|star| {
            star.candidate_diagnostics.iter().any(|diagnostic| {
                matches!(
                    diagnostic.source,
                    BrepEdgeSurfaceVertexCandidateSource::SurfaceIntersection
                        | BrepEdgeSurfaceVertexCandidateSource::SideAgreement
                )
            })
        })
        .map(|star| {
            let accepted = star
                .candidate_diagnostics
                .iter()
                .filter(|diagnostic| diagnostic.point.is_some())
                .count();
            let failures = star
                .candidate_diagnostics
                .iter()
                .filter(|diagnostic| {
                    !matches!(
                        diagnostic.outcome,
                        BrepEdgeSurfaceVertexCandidateOutcome::Accepted
                            | BrepEdgeSurfaceVertexCandidateOutcome::ReconciledEndpoint
                    )
                })
                .map(|diagnostic| {
                    format!(
                        "{:?}:{:?}:{:?}/motion={:.6}/spread={:.6}/point={:?}",
                        diagnostic.edge,
                        diagnostic.source,
                        diagnostic.outcome,
                        diagnostic.anchor_motion,
                        diagnostic.side_disagreement,
                        diagnostic.point
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            format!(
                "{:?} {:?} incompatible={:?} incident={:?} accepted_candidates={} motion={:.6} spread={:.6} preserved_gap={:.6} failures=[{}]",
                star.vertex,
                star.outcome,
                star.incompatible_edges,
                star.incident_edges,
                accepted,
                star.max_anchor_motion,
                star.max_candidate_spread,
                star.max_preserved_incident_gap,
                failures
            )
        })
        .collect()
}

fn edge_surface_repair_summary(
    report: &Option<BrepEdgeSurfaceCompatibilityRepairReport>,
) -> Vec<String> {
    let Some(report) = report else {
        return vec!["edge_surface_repair: not_run".to_string()];
    };
    report
        .decisions
        .iter()
        .filter(|decision| {
            matches!(
                (decision.evidence.surface_a, decision.evidence.surface_b),
                (Some(BrepEdgeSurfaceKind::Gregory), Some(BrepEdgeSurfaceKind::Cylinder))
                    | (Some(BrepEdgeSurfaceKind::Cylinder), Some(BrepEdgeSurfaceKind::Gregory))
            )
        })
        .map(|decision| {
            let declines = decision
                .evidence
                .boundary_declines
                .iter()
                .map(|decline| {
                    format!(
                        "{:?}:{:?}/best={:?}/dev={:.6}",
                        decline.side,
                        decline.reason,
                        decline.best_boundary,
                        decline.best_boundary_deviation
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            format!(
                "{:?} {:?} endpoint_gap={:.6} residual={:.6} boundary={} intersection={} anchor={} boundary_failure={:?} intersection_failure={:?} declines=[{}]",
                decision.edge,
                decision.outcome,
                decision.endpoint_gap,
                decision.residual_deviation,
                decision.evidence.attempted_boundary_reconstruction,
                decision.evidence.attempted_surface_intersection,
                decision.evidence.attempted_anchor_repair,
                decision.evidence.boundary_failure,
                decision.evidence.intersection_failure,
                declines
            )
        })
        .collect()
}

fn classify_fidelity_blocker(
    construction: &FidelityState,
    export_input: &FidelityState,
    imported: &FidelityState,
) -> &'static str {
    if construction.edge_surface_blockers > 0
        && export_input.edge_surface_blockers > 0
        && imported.edge_surface_blockers > 0
    {
        "construction_geometry_issue"
    } else if export_input.pcurve_blockers == 0 && imported.pcurve_blockers > 0 {
        "export_pcurve_issue"
    } else if export_input.max_pcurve_deviation <= 1e-6 && imported.max_pcurve_deviation > 1e-6 {
        "import_association_issue"
    } else if construction.pcurve_blockers > 0 && export_input.pcurve_blockers == 0 {
        "repair_pipeline_resolved_construction_issue"
    } else {
        "validator_gap"
    }
}

fn fidelity_state(label: &'static str, brep: &BRep) -> FidelityState {
    let quality = analyze_brep_export_quality(brep, quality_options());
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
    FidelityState {
        label,
        edges: brep.edges.len(),
        faces: brep.faces.len(),
        trim_findings: quality.trims.findings.len(),
        pcurve_blockers: quality.pcurve_blockers,
        max_pcurve_deviation: quality.max_pcurve_deviation,
        edge_surface_blockers: quality.edge_surface_blockers,
        max_edge_surface_deviation: quality.max_edge_surface_deviation,
        worst_pcurve_edges: worst_pcurve_edges(brep, 5),
        worst_edge_surface_edges: edge_surface
            .findings
            .iter()
            .take(8)
            .map(|finding| EdgeSurfaceEvidence {
                edge: format!("{:?}", finding.edge),
                half_edge: format!("{:?}", finding.half_edge),
                face: format!("{:?}", finding.face),
                max_deviation: finding.max_deviation,
                max_curve_parameter: finding.max_curve_parameter,
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
                    same_fraction_deviation: half_edge.max_same_fraction_deviation,
                    trace_deviation: half_edge.max_trace_deviation,
                    worst_fraction: half_edge.worst_fraction,
                })
                .collect(),
        })
        .collect()
}

fn mixed_radius_gregory_cube() -> BRep {
    let cube = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|half_edge| half_edge.origin == vertex)
        .map(|half_edge| half_edge.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let specs = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    fillet_edges(
        &cube,
        &specs,
        &FilletOptions {
            use_rational_gregory: true,
            ..FilletOptions::default()
        },
        popkern_core::tolerance::Tolerance::default(),
        None,
    )
    .expect("mixed-radius Gregory cube must build")
    .brep
}

fn quality_options() -> popkern_core::topology::BrepExportQualityOptions {
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
    }
}

fn gregory_step_repair_options(require_all_compatible: bool) -> BrepRepairPipelineOptions {
    let tolerance = 1e-6;
    let samples = 9;
    BrepRepairPipelineOptions {
        classification: BrepSewingClassificationOptions {
            exact_tolerance: tolerance,
            snap_tolerance: tolerance,
            samples,
            require_surface_agreement: false,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: BrepSnapSewingOptions {
            max_vertex_motion: 0.0,
            clean_sewing: BrepSewingOptions {
                tolerance,
                samples,
                require_surface_agreement: false,
                validate_output: true,
                allow_nonmanifold_edges: false,
            },
        },
        trim_repair: BrepTrimRepairOptions {
            endpoint_tolerance: tolerance,
            sampled_tolerance: tolerance,
            samples,
            validate_output: true,
        },
        edge_curve_completion: BrepEdgeCurveCompletionOptions {
            tolerance,
            samples,
            allow_vertex_segments_without_surface: true,
            require_all_completed: false,
            validate_output: true,
        },
        curve_trim_reconstruction: BrepCurveTrimReconstructionOptions {
            tolerance,
            samples,
            replace_existing: false,
            replace_stale_existing: true,
            expand_bounded_surface_domains: true,
            require_all_reconstructed: false,
            validate_output: true,
        },
        edge_surface_compatibility: BrepEdgeSurfaceCompatibilityOptions {
            tolerance,
            max_anchor_motion: 0.06,
            samples,
            require_all_compatible,
            validate_output: true,
        },
        degenerate_edges: BrepDegenerateEdgeOptions { tolerance, samples },
        require_all_pairs_repaired: true,
        complete_edge_curves: true,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: true,
        repair_trims: false,
        diagnose_loop_shell: true,
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
    }
}

fn step_document(brep: &BRep) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: "GregoryFidelityProbe".into(),
            id: "GREGORY-FIDELITY-PROBE".into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}
