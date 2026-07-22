use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepDegenerateEdgeOptions,
    BrepEdgeCurveCompletionOptions, BrepEdgeCurveCompletionOutcome,
    BrepEdgeSurfaceCompatibilityOptions, BrepEdgeSurfaceKind, BrepFaceUseTrimAuthorityOptions,
    BrepRepairPipelineOptions, BrepSewingClassificationOptions, BrepSewingOptions,
    BrepSnapSewingOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
    analyze_brep_edge_surface_compatibility, analyze_brep_export_quality,
    classify_brep_fidelity_blocker_ownership, repair_brep_boundary_edges_transactionally,
};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
    step::diagnostics::DiagnosticEntry,
};
use std::sync::Arc;
#[path = "gregory_step_fidelity/boundary_knowledge.rs"]
mod gregory_step_fidelity_boundary_knowledge;
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
    surface_kind: &'static str,
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
    surface_kind: &'static str,
    max_deviation: f64,
    max_curve_parameter: f64,
}

#[test]
fn gregory_step_fidelity_probe_classifies_blocker_source() {
    let blend = mixed_radius_gregory_cube_result();
    let construction = blend.brep;
    let repaired = repair_brep_boundary_edges_transactionally(
        &construction,
        &[],
        gregory_step_repair_options(false),
    )
    .expect("diagnostic repair pipeline should commit export candidate");
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
    // Re-pinned 2026-07-17 (PK-125): the cap-seam pcurve residual
    // blockers this diagnostic used to classify are fully resolved —
    // a clean import must not report any.
    assert!(
        imported.diagnostics.entries().iter().all(|entry| !entry
            .message
            .contains("STEP import pcurve fidelity blocker")),
        "Gregory STEP roundtrip regressed pcurve residual blockers: {:?}",
        imported.diagnostics.entries()
    );
    let imported_brep = &imported.document.products[0].shapes[0].solids[0];

    let construction_state = fidelity_state("construction", &construction);
    let export_state = fidelity_state("export_input", &export_input);
    let imported_state =
        imported_fidelity_state("imported", imported_brep, imported.diagnostics.entries());
    let states = [
        construction_state.clone(),
        export_state.clone(),
        imported_state.clone(),
    ];
    let construction_state = &states[0];
    let export_state = &states[1];
    let imported_state = &states[2];
    // Re-pinned 2026-07-21 (PK-153): the V-free DefaultOn n>=3
    // consolidated serving retires the cube-edge stub class this
    // probe used to classify — the fidelity blocker set is empty at
    // construction and must stay empty through export repair.
    let construction_ownership = classify_brep_fidelity_blocker_ownership(&construction, 1e-6);
    assert_eq!(
        construction_ownership.blockers, 0,
        "consolidated V-free construction should expose no fidelity blockers: {construction_ownership:#?}"
    );
    let export_ownership = classify_brep_fidelity_blocker_ownership(&export_input, 1e-6);
    assert_eq!(
        export_ownership.blockers, 0,
        "export repair must not invent fidelity blockers: {export_ownership:#?}"
    );

    assert!(
        export_state.pcurve_blockers <= construction_state.pcurve_blockers,
        "repair should not increase construction pcurve blockers: {states:#?}"
    );
    assert!(
        export_state.edge_surface_blockers <= construction_state.edge_surface_blockers,
        "repair should not increase construction edge/surface blockers: {states:#?}"
    );
    assert_eq!(
        export_state.pcurve_blockers, 0,
        "consolidated V-free export input should carry no pcurve blockers: {states:#?}"
    );
    let edge_surface_report = repaired
        .report
        .edge_surface_compatibility
        .as_ref()
        .expect("edge/surface repair report should be present");
    let edge_curve_report = repaired
        .report
        .edge_curve_completion
        .as_ref()
        .expect("edge curve completion report should be present");
    let edge_curve_decline_reasons = edge_curve_report
        .decisions
        .iter()
        .filter_map(|decision| match decision.outcome {
            BrepEdgeCurveCompletionOutcome::Declined(reason) => Some(reason),
            _ => None,
        })
        .collect::<Vec<_>>();
    // Re-pinned 2026-07-21 (PK-153): the consolidated V-free serving
    // authors every edge curve, so the completion stage has nothing
    // to classify — a decline reappearing here is a regression.
    assert!(
        edge_curve_decline_reasons.is_empty(),
        "consolidated V-free serving should leave no edge-curve completion declines: {edge_curve_decline_reasons:#?}"
    );
    assert!(
        edge_surface_report.decisions.iter().all(|decision| {
            !matches!(
                (decision.evidence.surface_a, decision.evidence.surface_b),
                (
                    Some(BrepEdgeSurfaceKind::Gregory),
                    Some(BrepEdgeSurfaceKind::Cylinder)
                ) | (
                    Some(BrepEdgeSurfaceKind::Cylinder),
                    Some(BrepEdgeSurfaceKind::Gregory)
                )
            ) || decision
                .evidence
                .repair_knowledge_diagnostics
                .iter()
                .any(|code| code.starts_with("brep.surface.boundary"))
        }),
        "Gregory/cylinder blocker decisions must expose promoted surface-boundary knowledge diagnostics when such decisions remain"
    );
    // Re-pinned 2026-07-17 (PK-125): the imported model is now fully
    // clean (0 pcurve / 0 edge-surface blockers) — the residual
    // pcurve class no longer survives the roundtrip. The remaining
    // 9-blocker census lives on the construction-side analyzers only
    // (the validator_gap classification below).
    assert_eq!(
        imported_state.pcurve_blockers, 0,
        "STEP import regressed pcurve blockers: {states:#?}"
    );

    let classification =
        classify_fidelity_blocker(construction_state, export_state, imported_state);
    // Census provenance: "construction_geometry_issue" through PK-116
    // (edge/surface blockers survived construction, export, AND import).
    // PK-124's seam fidelity package reconciled every imported
    // edge/surface pairing (imported edge_surface_blockers 0; only the
    // residual imported pcurve class remained — "validator_gap").
    // PK-153's consolidated V-free serving empties every class at
    // construction, export, and import.
    assert_eq!(
        classification, "no_fidelity_blocker",
        "unexpected Gregory STEP fidelity classification; states: {states:#?}"
    );
}

fn classify_fidelity_blocker(
    construction: &FidelityState,
    export_input: &FidelityState,
    imported: &FidelityState,
) -> &'static str {
    if construction.pcurve_blockers == 0
        && construction.edge_surface_blockers == 0
        && export_input.pcurve_blockers == 0
        && export_input.edge_surface_blockers == 0
        && imported.pcurve_blockers == 0
        && imported.edge_surface_blockers == 0
    {
        "no_fidelity_blocker"
    } else if construction.edge_surface_blockers > 0
        && export_input.edge_surface_blockers > 0
        && imported.edge_surface_blockers > 0
    {
        "construction_geometry_issue"
    } else if export_input.pcurve_blockers == 0
        && export_input.edge_surface_blockers == 0
        && (imported.pcurve_blockers > 0 || imported.edge_surface_blockers > 0)
    {
        "step_import_mapping_issue"
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
                surface_kind: face_surface_kind(brep, finding.face),
                max_deviation: finding.max_deviation,
                max_curve_parameter: finding.max_curve_parameter,
            })
            .collect(),
    }
}

fn imported_fidelity_state(
    label: &'static str,
    brep: &BRep,
    diagnostics: &[DiagnosticEntry],
) -> FidelityState {
    FidelityState {
        label,
        edges: brep.edges.len(),
        faces: brep.faces.len(),
        trim_findings: 0,
        pcurve_blockers: diagnostics
            .iter()
            .filter(|entry| {
                entry
                    .message
                    .contains("STEP import pcurve fidelity blocker")
            })
            .count(),
        max_pcurve_deviation: 0.0,
        edge_surface_blockers: diagnostics
            .iter()
            .filter(|entry| {
                entry
                    .message
                    .contains("STEP import edge/surface fidelity blocker")
            })
            .count(),
        max_edge_surface_deviation: 0.0,
        worst_pcurve_edges: Vec::new(),
        worst_edge_surface_edges: Vec::new(),
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
                    surface_kind: face_surface_kind(brep, half_edge.face),
                    max_deviation: half_edge.max_deviation,
                    same_fraction_deviation: half_edge.max_same_fraction_deviation,
                    trace_deviation: half_edge.max_trace_deviation,
                    worst_fraction: half_edge.worst_fraction,
                })
                .collect(),
        })
        .collect()
}

fn face_surface_kind(brep: &BRep, face: popkern_core::topology::FaceHandle) -> &'static str {
    match brep
        .faces
        .get(face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    {
        Some(SurfaceGeometry::Plane { .. }) => "plane",
        Some(SurfaceGeometry::Cylinder { .. }) => "cylinder",
        Some(SurfaceGeometry::Cone { .. }) => "cone",
        Some(SurfaceGeometry::Sphere { .. }) => "sphere",
        Some(SurfaceGeometry::Torus { .. }) => "torus",
        Some(SurfaceGeometry::Nurbs(_)) => "nurbs",
        Some(SurfaceGeometry::Gregory(_)) => "gregory",
        None => "missing",
    }
}

pub(crate) fn mixed_radius_gregory_cube() -> BRep {
    mixed_radius_gregory_cube_result().brep
}

fn mixed_radius_gregory_cube_result() -> popkern_fillet::BlendResult {
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
        face_use_trim_authority: BrepFaceUseTrimAuthorityOptions { tolerance, samples },
        degenerate_edges: BrepDegenerateEdgeOptions { tolerance, samples },
        require_all_pairs_repaired: true,
        complete_edge_curves: true,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: true,
        repair_face_use_trim_authority: false,
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
