use popkern_core::geometry::{
    Point3,
    surface::{SurfaceEval, SurfaceGeometry},
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepDegenerateEdgeOptions,
    BrepEdgeCurveCompletionOptions, BrepEdgeSurfaceCompatibilityOptions, BrepRepairPipelineOptions,
    BrepSewingClassificationOptions, BrepSewingOptions, BrepSnapSewingOptions,
    BrepTopologyOutputPolicy, BrepTrimRepairOptions, analyze_brep_edge_surface_compatibility,
    repair_brep_boundary_edges_transactionally,
};
use popkern_fillet::{
    ApexConstructionDiagnostic, BlendResult, ExactBoundaryCandidateOutcome, FilletOptions,
    FilletRadius, FilletSpec, fillet_edges,
};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};
use std::sync::Arc;

fn mixed_radius_gregory_cube_result() -> BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
        .expect("mixed-radius Gregory cube must build")
}

fn mixed_radius_gregory_cube() -> BRep {
    mixed_radius_gregory_cube_result().brep
}

struct CornerPatchFaceSamples {
    samples: Vec<(Point3, f64, f64)>,
}

fn sample_native_cap_faces(result: &BlendResult) -> Vec<CornerPatchFaceSamples> {
    let exact_faces: Vec<_> = result
        .repair_report
        .apex_routes
        .iter()
        .flat_map(|route| route.exact_boundary_report.committed_faces.iter().copied())
        .collect();
    sample_faces_by_handle(&result.brep, &exact_faces)
}

fn sample_faces_by_handle(
    brep: &BRep,
    faces: &[popkern_core::topology::FaceHandle],
) -> Vec<CornerPatchFaceSamples> {
    let mut out = Vec::new();
    for face in faces {
        let Some(sof) = brep.faces.get(*face).and_then(|face| face.surface.as_ref()) else {
            continue;
        };
        let (SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)) = &sof.geometry else {
            continue;
        };
        let mut samples = Vec::new();
        for i in 1..=3 {
            for j in 1..=3 {
                let u = i as f64 / 4.0;
                let v = j as f64 / 4.0;
                samples.push((sof.geometry.point_at(u, v), u, v));
            }
        }
        out.push(CornerPatchFaceSamples { samples });
    }
    out
}

/// STEP round-trip on a mixed-radius Gregory fillet. Post-round-trip
/// the BRep must:
/// (1) re-import without errors,
/// (2) carry ≥1 solid with ≥ the original face count,
/// (3) validate without structural errors.
#[test]
fn gregory_fillet_step_roundtrip_preserves_topology() {
    let brep = mixed_radius_gregory_cube();
    let orig_face_count = brep.faces.len();

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "GregoryFilletedCube".into(),
            id: "GFC001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default())
        .expect("export_step must succeed on Gregory-carrying BRep");
    assert!(
        bytes.len() > 1000,
        "Gregory-faced BRep STEP encoding should produce non-trivial output; got {} bytes",
        bytes.len()
    );

    let import_result =
        import_step(&bytes, &StepImportOptions::default()).expect("import_step must succeed");
    let reimported = &import_result.document.products[0].shapes[0].solids[0];
    reimported
        .validate_loop_closure()
        .expect("re-imported Gregory BRep loops must close");

    assert!(
        reimported.faces.len() >= orig_face_count,
        "re-imported BRep must have ≥ {orig_face_count} faces; got {}",
        reimported.faces.len()
    );
    assert_eq!(
        reimported.solids.len(),
        1,
        "re-imported BRep must have 1 solid"
    );

    // Structural validation: re-imported BRep must not have broken
    // cycles, broken twins, or non-manifold edges.
    use popkern_core::error::ValidationError;
    let errors = reimported.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                ValidationError::BrokenCycle { .. }
                    | ValidationError::BrokenTwin { .. }
                    | ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(
        structural.is_empty(),
        "re-imported BRep has structural errors: {structural:?}"
    );
}

/// Native cap faces in the original BRep are exported as STEP NURBS geometry.
/// Geometric fidelity must be preserved: sample the original cap faces at
/// interior UV points, locate the corresponding face in the re-imported BRep,
/// and assert the position deviation is within tolerance.
///
/// Older Gregory wrappers use `to_nurbs_exact`; newer macro-cap rungs
/// may already be NURBS before export. In both cases STEP should carry
/// the same-UV surface samples to < 1e-6 after parser re-encoding.
#[test]
fn gregory_face_geometry_survives_step_roundtrip() {
    let result = mixed_radius_gregory_cube_result();
    let orig_gregory_faces = sample_native_cap_faces(&result);
    let brep = result.brep;
    assert!(
        !orig_gregory_faces.is_empty(),
        "test precondition: original BRep must contain native cap patch faces"
    );

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "GregoryFilletedCube".into(),
            id: "GFC001".into(),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default()).unwrap();
    // Use healing=None so the re-imported NURBS is byte-for-byte
    // the exported form; `Basic` healing can simplify / re-fit
    // geometry which would introduce additional deviation beyond
    // the rational→NURBS export path that item #2 targets.
    let import_opts = popkern_io::StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let reimport = import_step(&bytes, &import_opts).unwrap();
    assert!(
        reimport.diagnostics.errors().next().is_none(),
        "STEP reimport must not produce errors: {:?}",
        reimport.diagnostics.entries()
    );
    let reimport_brep = &reimport.document.products[0].shapes[0].solids[0];
    reimport_brep
        .validate_loop_closure()
        .expect("re-imported Gregory BRep loops must close");
    let imported_nurbs: Vec<_> = reimport_brep
        .faces
        .values()
        .filter_map(|face| face.surface.as_ref())
        .filter(|sof| matches!(sof.geometry, SurfaceGeometry::Nurbs(_)))
        .collect();
    assert!(
        imported_nurbs.len() >= orig_gregory_faces.len(),
        "reimported BRep must contain at least {} NURBS faces converted from Gregory faces; got {}",
        orig_gregory_faces.len(),
        imported_nurbs.len()
    );

    // STEP exports Gregory-method cap faces as NURBS with the same
    // full-range parameterisation. Match each original cap face to the
    // imported NURBS face that minimises same-UV sample error.
    let mut max_err: f64 = 0.0;
    for original_face in &orig_gregory_faces {
        let mut best_face_err = f64::INFINITY;
        for imported in &imported_nurbs {
            let mut face_err: f64 = 0.0;
            for (orig_pt, u, v) in &original_face.samples {
                let p = imported.geometry.point_at(*u, *v);
                face_err = face_err.max(p.distance_to(orig_pt));
            }
            best_face_err = best_face_err.min(face_err);
        }
        max_err = max_err.max(best_face_err);
    }

    assert!(
        max_err < 1e-6,
        "§5.3 item #2: Gregory-method STEP round-trip same-UV position deviation {max_err} exceeds 1e-6 tolerance (corner patch faces checked: {})",
        orig_gregory_faces.len()
    );
}

#[test]
fn gregory_mixed_radius_edge_surface_pipeline_diagnostics_are_structured() {
    let blend = mixed_radius_gregory_cube_result();
    let brep = &blend.brep;
    let support_report = &blend.repair_report.support_seam_relocation;
    assert_eq!(
        support_report.candidates,
        support_report.authored + support_report.declined
    );
    assert_eq!(support_report.decisions.len(), support_report.candidates);
    let artifact = popkern_fillet::build_blend_repair_action_artifact(&blend.repair_report);
    let apex_diagnostic_actions = blend
        .repair_report
        .apex_diagnostics
        .iter()
        .filter(|diagnostic| {
            matches!(
                diagnostic,
                ApexConstructionDiagnostic::ExactBoundaryCap(_)
                    | ApexConstructionDiagnostic::CapSupportContinuity(_)
            )
        })
        .count();
    assert_eq!(
        artifact.summary.attempted_actions,
        support_report.candidates + apex_diagnostic_actions
    );
    assert_gregory_boundary_evidence_or_exact_boundary_acceptance(&blend);
    let before = analyze_brep_edge_surface_compatibility(
        brep,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.0,
            samples: 9,
            require_all_compatible: false,
            validate_output: false,
        },
    );

    let result = repair_brep_boundary_edges_transactionally(
        brep,
        &[],
        gregory_edge_surface_pipeline_options(false),
    )
    .expect("diagnostic repair pipeline should commit partial edge/surface repair");

    let report = result
        .report
        .edge_surface_compatibility
        .as_ref()
        .expect("edge/surface compatibility stage must run");
    assert_eq!(
        report.before.incompatible_half_edges,
        before.incompatible_half_edges
    );
    assert!(
        report.before.incompatible_half_edges > 0,
        "real mixed-radius Gregory fixture must exercise edge/surface diagnostics"
    );
    assert_eq!(
        report.attempted_edges,
        report.repaired_edges + report.declined_edges
    );
    assert!(
        report.after.incompatible_half_edges <= report.before.incompatible_half_edges,
        "pipeline repair must not increase edge/surface incompatibility"
    );
    assert!(
        report.repaired_edges > 0 || report.declined_edges > 0,
        "pipeline must produce per-edge repair or decline decisions"
    );
    assert!(
        report.declined_edges > 0,
        "real mixed-radius Gregory fixture must preserve unresolved edge/surface declines"
    );
    assert!(
        brep.faces.len() == result.brep.faces.len(),
        "edge/surface repair must preserve face ownership"
    );
}

fn assert_gregory_boundary_evidence_or_exact_boundary_acceptance(blend: &BlendResult) {
    let exact_report = &blend.repair_report.apex_routes[0].exact_boundary_report;
    assert!(
        !exact_report.committed_faces.is_empty(),
        "mixed-radius fixture should use accepted exact-boundary cap faces"
    );
    assert_eq!(
        exact_report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted
    );
    assert_eq!(exact_report.performance.pcurve_pair_count, 6);
    assert!(exact_report.diagnostics.is_empty());
}

#[test]
fn gregory_mixed_radius_strict_edge_surface_pipeline_accepts_exact_boundary_route() {
    let brep = mixed_radius_gregory_cube();
    let result = repair_brep_boundary_edges_transactionally(
        &brep,
        &[],
        gregory_edge_surface_pipeline_options(true),
    )
    .expect("strict edge/surface compatibility should accept the exact-boundary cube_mixed route");

    let report = result
        .report
        .edge_surface_compatibility
        .as_ref()
        .expect("strict repair must include edge/surface report");
    assert_eq!(
        report.after.incompatible_half_edges, 0,
        "strict exact-boundary cube_mixed route must close all edge/surface incompatibilities: {report:#?}"
    );
}

fn gregory_edge_surface_pipeline_options(
    require_all_compatible: bool,
) -> BrepRepairPipelineOptions {
    let clean_sewing = BrepSewingOptions {
        tolerance: 1e-6,
        samples: 5,
        require_surface_agreement: false,
        validate_output: false,
        allow_nonmanifold_edges: false,
    };
    BrepRepairPipelineOptions {
        classification: BrepSewingClassificationOptions {
            exact_tolerance: 1e-9,
            snap_tolerance: 1e-6,
            samples: 5,
            require_surface_agreement: false,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: BrepSnapSewingOptions {
            max_vertex_motion: 1e-6,
            clean_sewing,
        },
        trim_repair: BrepTrimRepairOptions {
            endpoint_tolerance: 1e-6,
            sampled_tolerance: 1e-6,
            samples: 9,
            validate_output: false,
        },
        edge_curve_completion: BrepEdgeCurveCompletionOptions {
            tolerance: 1e-6,
            samples: 9,
            allow_vertex_segments_without_surface: true,
            require_all_completed: false,
            validate_output: false,
        },
        curve_trim_reconstruction: BrepCurveTrimReconstructionOptions {
            tolerance: 1e-6,
            samples: 9,
            replace_existing: false,
            replace_stale_existing: true,
            expand_bounded_surface_domains: true,
            require_all_reconstructed: false,
            validate_output: false,
        },
        edge_surface_compatibility: BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.06,
            samples: 9,
            require_all_compatible,
            validate_output: false,
        },
        degenerate_edges: BrepDegenerateEdgeOptions {
            tolerance: 1e-6,
            samples: 9,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: false,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: true,
        repair_trims: false,
        diagnose_loop_shell: false,
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
    }
}
