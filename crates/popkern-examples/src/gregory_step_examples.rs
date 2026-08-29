//! STEP export and manifest example for Gregory-carrying fillets.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use popkern_blend::{
    FilletOptions, FilletRadius, FilletSpec, build_blend_repair_action_artifact, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepDegenerateEdgeOptions,
    BrepEdgeCurveCompletionOptions, BrepEdgeSurfaceCompatibilityOptions, BrepExportQualityOptions,
    BrepFaceUseTrimAuthorityOptions, BrepRepairDiagnosticExport, BrepRepairPipelineOptions,
    BrepRepairPipelineReport, BrepSewingClassificationOptions, BrepSewingOptions,
    BrepSnapSewingOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
    analyze_brep_export_quality, build_brep_export_quality_artifact,
    repair_brep_boundary_edges_transactionally,
};
use popkern_examples::example_manifests::{ExampleRunManifestInput, write_example_run_manifest};
use popkern_examples::repair_artifacts::{
    write_repair_action_artifact, write_whole_brep_repair_knowledge_artifact,
};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};

struct GregoryExampleReport {
    step_path: PathBuf,
    quality_path: PathBuf,
    repair_path: PathBuf,
    generic_repair_path: PathBuf,
    generic_repair_knowledge_path: PathBuf,
    step_bytes: usize,
    validated_bodies: usize,
    faces: usize,
    edges: usize,
    corner_patches: usize,
    repair_actions: usize,
    generic_repair_actions: usize,
}

fn main() {
    let output_dir = output_dir();
    std::fs::create_dir_all(&output_dir).expect("create Gregory STEP output dir");
    let report = export_mixed_radius_gregory_cube(&output_dir)
        .expect("mixed-radius Gregory STEP example must export and validate");
    println!(
        "mixed_radius_gregory_cube -> {} bytes, {} validated body/bodies",
        report.step_bytes, report.validated_bodies
    );
    println!("  step {}", report.step_path.display());
    println!("  quality {}", report.quality_path.display());
    println!("  repair {}", report.repair_path.display());
    println!("  generic repair {}", report.generic_repair_path.display());
    println!(
        "  repair knowledge {}",
        report.generic_repair_knowledge_path.display()
    );
    let manifest = write_manifest(&output_dir, &report).expect("write Gregory manifest");
    println!("manifest {}", manifest.display());
}

fn export_mixed_radius_gregory_cube(output_dir: &Path) -> Result<GregoryExampleReport, String> {
    let blend = mixed_radius_gregory_cube()?;
    let (brep, generic_repair) = repair_for_step_export(&blend.brep)?;
    let step_path = output_dir.join("mixed_radius_gregory_cube.step");
    let bytes = export_step(&step_document(&brep), &StepExportOptions::default())
        .map_err(|err| err.to_string())?;
    if bytes.len() < 1024 {
        return Err(format!(
            "STEP output is unexpectedly small: {} bytes",
            bytes.len()
        ));
    }
    std::fs::write(&step_path, &bytes).map_err(|err| err.to_string())?;
    let validated_bodies = validate_gregory_step_bytes(&bytes)?;
    let quality_path = output_dir.join("mixed_radius_gregory_cube_quality.json");
    write_quality_artifact(&brep, &quality_path)?;
    let repair_path = output_dir.join("mixed_radius_gregory_cube_repair.json");
    let repair_artifact = build_blend_repair_action_artifact(&blend.repair_report);
    write_repair_action_artifact(&repair_path, "mixed_radius_gregory_cube", &repair_artifact)?;
    let generic_repair_path = output_dir.join("mixed_radius_gregory_cube_generic_repair.json");
    write_generic_repair_artifact(&generic_repair_path, &generic_repair)?;
    let generic_repair_knowledge_path =
        output_dir.join("mixed_radius_gregory_cube_repair_knowledge.json");
    write_whole_brep_repair_knowledge_artifact(
        &generic_repair_knowledge_path,
        "mixed_radius_gregory_cube",
        &brep,
        9,
        BrepRepairDiagnosticExport::FullJson,
    )?;

    Ok(GregoryExampleReport {
        step_path,
        quality_path,
        repair_path,
        generic_repair_path,
        generic_repair_knowledge_path,
        step_bytes: bytes.len(),
        validated_bodies,
        faces: brep.faces.len(),
        edges: brep.edges.len(),
        corner_patches: blend.corner_patches.len(),
        repair_actions: repair_artifact.summary.attempted_actions,
        generic_repair_actions: generic_repair.actions.attempted_actions,
    })
}

fn write_generic_repair_artifact(
    path: &Path,
    report: &BrepRepairPipelineReport,
) -> Result<(), String> {
    let edge_surface = report
        .edge_surface_compatibility
        .as_ref()
        .map(|stage| {
            serde_json::json!({
                "attempted_edges": stage.attempted_edges,
                "repaired_edges": stage.repaired_edges,
                "declined_edges": stage.declined_edges,
                "max_endpoint_gap": stage.max_endpoint_gap,
                "shared_fact_ids": stage.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
                "before": {
                    "incompatible_half_edges": stage.before.incompatible_half_edges,
                    "max_deviation": stage.before.max_deviation,
                },
                "after": {
                    "incompatible_half_edges": stage.after.incompatible_half_edges,
                    "max_deviation": stage.after.max_deviation,
                },
                "vertex_stars": {
                    "inspected_vertices": stage.vertex_stars.inspected_vertices,
                    "coherent_vertices": stage.vertex_stars.coherent_vertices,
                    "declined_vertices": stage.vertex_stars.declined_vertices,
                    "shared_fact_ids": stage.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
                    "decisions": stage.vertex_stars.decisions.iter().map(|decision| {
                        serde_json::json!({
                            "vertex": format!("{:?}", decision.vertex),
                            "outcome": format!("{:?}", decision.outcome),
                            "incompatible_edges": decision.incompatible_edges.iter().map(|edge| format!("{edge:?}")).collect::<Vec<_>>(),
                            "incident_edges": decision.incident_edges.iter().map(|edge| format!("{edge:?}")).collect::<Vec<_>>(),
                            "max_anchor_motion": decision.max_anchor_motion,
                            "max_candidate_spread": decision.max_candidate_spread,
                            "max_preserved_incident_gap": decision.max_preserved_incident_gap,
                            "proposed_position": decision.proposed_position.map(|point| serde_json::json!([point.x, point.y, point.z])),
                            "shared_fact_ids": decision.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
                            "candidate_diagnostics": decision.candidate_diagnostics.iter().map(|diagnostic| serde_json::json!({
                                "edge": format!("{:?}", diagnostic.edge),
                                "half_edge": diagnostic.half_edge.map(|half_edge| format!("{half_edge:?}")),
                                "source": format!("{:?}", diagnostic.source),
                                "outcome": format!("{:?}", diagnostic.outcome),
                                "point": diagnostic.point.map(|point| serde_json::json!([point.x, point.y, point.z])),
                                "anchor_motion": diagnostic.anchor_motion,
                                "side_disagreement": diagnostic.side_disagreement,
                                "shared_fact_ids": diagnostic.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
                            })).collect::<Vec<_>>(),
                        })
                    }).collect::<Vec<_>>(),
                },
                "decisions": stage.decisions.iter().map(|decision| {
                    serde_json::json!({
                        "edge": format!("{:?}", decision.edge),
                        "outcome": format!("{:?}", decision.outcome),
                        "endpoint_gap": decision.endpoint_gap,
                        "residual_deviation": decision.residual_deviation,
                        "shared_fact_ids": decision.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
                        "evidence": {
                            "surface_a": decision.evidence.surface_a.map(|kind| format!("{kind:?}")),
                            "surface_b": decision.evidence.surface_b.map(|kind| format!("{kind:?}")),
                            "boundary_hint": decision.evidence.boundary_hint.map(|hint| serde_json::json!({
                                "side": format!("{:?}", hint.side),
                                "boundary": format!("{:?}", hint.boundary),
                            })),
                            "boundary_declines": decision.evidence.boundary_declines.iter().map(|decline| serde_json::json!({
                                "side": format!("{:?}", decline.side),
                                "reason": format!("{:?}", decline.reason),
                                "best_boundary": decline.best_boundary.map(|boundary| format!("{boundary:?}")),
                                "best_boundary_deviation": decline.best_boundary_deviation,
                                "endpoint_gap": decline.endpoint_gap,
                            })).collect::<Vec<_>>(),
                            "attempted_boundary_reconstruction": decision.evidence.attempted_boundary_reconstruction,
                            "attempted_surface_intersection": decision.evidence.attempted_surface_intersection,
                            "attempted_anchor_repair": decision.evidence.attempted_anchor_repair,
                            "anchor_motion": decision.evidence.anchor_motion,
                            "max_allowed_anchor_motion": decision.evidence.max_allowed_anchor_motion,
                            "max_incident_endpoint_gap": decision.evidence.max_incident_endpoint_gap,
                            "boundary_failure": decision.evidence.boundary_failure.as_deref(),
                            "intersection_failure": decision.evidence.intersection_failure.as_deref(),
                        },
                    })
                }).collect::<Vec<_>>(),
            })
        });
    let curve_trim_reconstruction = report.curve_trim_reconstruction.as_ref().map(|stage| {
        serde_json::json!({
            "inspected_half_edges": stage.inspected_half_edges,
            "targeted_half_edges": stage.targeted_half_edges,
            "reconstructed_half_edges": stage.reconstructed_half_edges,
            "declined_half_edges": stage.declined_half_edges,
            "max_sampled_disagreement": stage.max_sampled_disagreement,
            "shared_fact_ids": stage.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
            "decisions": stage.decisions.iter().map(|decision| serde_json::json!({
                "half_edge": format!("{:?}", decision.half_edge),
                "outcome": format!("{:?}", decision.outcome),
                "sampled_disagreement": decision.sampled_disagreement,
                "shared_fact_ids": decision.shared_fact_ids.iter().map(|id| id.0).collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
        })
    });
    let payload = serde_json::json!({
        "schema": "popkern.generic_brep_repair",
        "schema_version": 1,
        "label": "mixed_radius_gregory_cube",
        "curve_trim_reconstruction": curve_trim_reconstruction,
        "edge_surface_compatibility": edge_surface,
    });
    let json = serde_json::to_string_pretty(&payload).map_err(|err| err.to_string())?;
    std::fs::write(path, json).map_err(|err| err.to_string())
}

fn repair_for_step_export(brep: &BRep) -> Result<(BRep, BrepRepairPipelineReport), String> {
    repair_brep_boundary_edges_transactionally(brep, &[], gregory_repair_options())
        .map(|result| (result.brep, result.report))
        .map_err(|failure| {
            format!(
                "generic BRep repair failed at {:?}: {:?}",
                failure.stage, failure.reason
            )
        })
}

fn gregory_repair_options() -> BrepRepairPipelineOptions {
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
            require_all_compatible: false,
            validate_output: true,
        },
        face_use_trim_authority: BrepFaceUseTrimAuthorityOptions { tolerance, samples },
        degenerate_edges: BrepDegenerateEdgeOptions { tolerance, samples },
        require_all_pairs_repaired: true,
        complete_edge_curves: true,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: true,
        repair_face_use_trim_authority: true,
        repair_trims: false,
        diagnose_loop_shell: true,
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
    }
}

fn mixed_radius_gregory_cube() -> Result<popkern_blend::BlendResult, String> {
    let cube = BRep::unit_cube(Tolerance::default()).map_err(|err| err.to_string())?;
    let vertex = cube
        .vertices
        .keys()
        .next()
        .ok_or_else(|| "unit cube has no vertices".to_string())?;
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
    fillet_edges(
        &cube,
        &specs,
        &FilletOptions {
            use_rational_gregory: true,
            ..FilletOptions::default()
        },
        Tolerance::default(),
        None,
    )
    .map_err(|err| err.to_string())
}

fn step_document(brep: &BRep) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: "MixedRadiusGregoryCube".into(),
            id: "GREGORY-MIXED-RADIUS-CUBE".into(),
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

fn write_quality_artifact(brep: &BRep, path: &Path) -> Result<(), String> {
    let quality = analyze_brep_export_quality(
        brep,
        BrepExportQualityOptions {
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
    let artifact = build_brep_export_quality_artifact(brep, &quality);
    let json = serde_json::to_string_pretty(&artifact).map_err(|err| err.to_string())?;
    std::fs::write(path, json).map_err(|err| err.to_string())
}

fn write_manifest(output_dir: &Path, report: &GregoryExampleReport) -> Result<PathBuf, String> {
    let selected = ["mixed_radius_gregory_cube"];
    write_example_run_manifest(
        output_dir,
        "run_manifest.json",
        ExampleRunManifestInput {
            runner: "gregory_step",
            preset: "smoke",
            selected_fixtures: &selected,
            metadata: serde_json::json!({
                "domain": "Gregory fillet STEP examples",
                "step_validation": "full_import_brep_validation_with_export_quality_artifact",
                "repair_artifact": true,
            }),
            fixtures: vec![serde_json::json!({
                "stem": "mixed_radius_gregory_cube",
                "step": {
                    "path": report.step_path.display().to_string(),
                    "bytes": report.step_bytes,
                    "validated_bodies": report.validated_bodies,
                },
                "quality_artifact": report.quality_path.display().to_string(),
                "repair_artifact": report.repair_path.display().to_string(),
                "generic_repair_artifact": report.generic_repair_path.display().to_string(),
                "generic_repair_knowledge_artifact": report.generic_repair_knowledge_path.display().to_string(),
                "topology": {
                    "faces": report.faces,
                    "edges": report.edges,
                },
                "gregory": {
                    "corner_patches": report.corner_patches,
                    "repair_actions": report.repair_actions,
                    "generic_repair_actions": report.generic_repair_actions,
                },
            })],
            failures: Vec::new(),
        },
    )
}

fn validate_gregory_step_bytes(bytes: &[u8]) -> Result<usize, String> {
    let imported =
        import_step(bytes, &StepImportOptions::default()).map_err(|err| err.to_string())?;
    let errors: Vec<_> = imported
        .diagnostics
        .errors()
        .map(|entry| format!("{:?}: {}", entry.phase, entry.message))
        .collect();
    if !errors.is_empty() {
        return Err(format!(
            "STEP import diagnostics contain errors: {}",
            errors.join("; ")
        ));
    }
    let bodies: Vec<_> = imported
        .document
        .products
        .iter()
        .flat_map(|product| &product.shapes)
        .flat_map(|shape| &shape.solids)
        .collect();
    if bodies.len() != 1 {
        return Err(format!(
            "imported body count {} does not match expected 1",
            bodies.len()
        ));
    }
    bodies[0]
        .validate()
        .map_err(|err| format!("imported body validation failed: {err}"))?;
    Ok(bodies.len())
}

fn output_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("examples crate has workspace parent")
        .parent()
        .expect("workspace parent has repo parent")
        .join("examples")
        .join("output")
        .join("gregory_step")
}
