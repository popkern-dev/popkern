//! STEP export showcase for selected ACM TOMS Algorithm 1032 PnS fixtures.
//!
//! Outputs are written under `examples/output/pns_toms/`.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use popkern_core::geometry::polyhedral_spline::{
    ControlNet, PolyhedralCoverageOptions, analyze_surface_complex_sewing,
    build_control_net_surface_complex, weld_clean_surface_complex_boundaries,
};
use popkern_core::topology::{
    BRep, BrepExportQualityOptions, BrepSewingClassificationBucket, BrepTopologyOutputPolicy,
    BrepTrimRepairOptions, analyze_brep_export_quality, brep_sewing_classification_bucket_label,
    build_brep_export_quality_artifact,
};
use popkern_examples::repair_artifacts::write_repair_action_artifact;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step_file,
};

#[path = "pns_toms_step_manifest.rs"]
mod pns_toms_step_manifest;
#[path = "pns_toms_step_report.rs"]
mod pns_toms_step_report;
#[path = "pns_toms_step_runner.rs"]
mod pns_toms_step_runner;
#[path = "pns_toms_step_types.rs"]
mod pns_toms_step_types;
#[path = "pns_toms_step_validation.rs"]
mod pns_toms_step_validation;
use pns_toms_step_manifest::{FixtureFailure, FixtureSuccess, write_run_manifest};
use pns_toms_step_report::write_sidecar_report;
use pns_toms_step_runner::{
    RunnerCommand, fixture_list, parse_runner_command, selected_fixture_stems,
};
use pns_toms_step_types::FixtureReport;
use pns_toms_step_validation::validate_step_file;

struct Fixture {
    name: &'static str,
    file_stem: &'static str,
    obj: &'static str,
    mode: FixtureMode,
}

#[derive(Clone, Copy)]
enum FixtureMode {
    Strict,
    Partial,
}

const FIXTURES: &[Fixture] = &[
    Fixture {
        name: "PnS T0 Transition",
        file_stem: "t0",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS T1 Transition",
        file_stem: "t1",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS T2 Transition",
        file_stem: "t2",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Polar Five-Sector",
        file_stem: "polar5sct",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar5sct.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Extraordinary Five-Star",
        file_stem: "eop5sct",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop5sct.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Pentagon N-gon",
        file_stem: "ngon5",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon5.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Tee",
        file_stem: "tee",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/tee.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS T012 Combined",
        file_stem: "t012",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T012.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Cube",
        file_stem: "cube",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cube.obj"),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Suzanne Configurations",
        file_stem: "suzanne_all_configurations",
        obj: include_str!(
            "../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/suzanne_all_configurations.obj"
        ),
        mode: FixtureMode::Strict,
    },
    Fixture {
        name: "PnS Cylinder Partial",
        file_stem: "cylinder_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cylinder.obj"),
        mode: FixtureMode::Partial,
    },
    Fixture {
        name: "PnS Hand Partial",
        file_stem: "hand_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/hand.obj"),
        mode: FixtureMode::Partial,
    },
    Fixture {
        name: "PnS Airplane Partial",
        file_stem: "airplane_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/airplane.obj"),
        mode: FixtureMode::Partial,
    },
];

pub(crate) fn run() {
    println!("=== PnS TOMS STEP Examples ===");
    let output_dir = output_dir();
    std::fs::create_dir_all(&output_dir).expect("create PnS TOMS output dir");
    let all_stems: Vec<_> = FIXTURES.iter().map(|fixture| fixture.file_stem).collect();
    let command = parse_runner_command(std::env::args()).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(2);
    });
    let selection = match command {
        RunnerCommand::Run(selection) => selection,
        RunnerCommand::ListFixtures => {
            println!("{}", fixture_list(&all_stems));
            return;
        }
    };
    let selected_stems = selected_fixture_stems(&selection, &all_stems).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(2);
    });
    println!(
        "preset: {}, fixtures: {}",
        selection.preset.label(),
        selected_stems.join(", ")
    );

    let mut reports = Vec::new();
    let mut failures = Vec::new();
    for fixture in FIXTURES
        .iter()
        .filter(|fixture| selected_stems.contains(&fixture.file_stem))
    {
        match export_fixture(fixture, &output_dir) {
            Ok(report) => {
                print_report(fixture, &report);
                reports.push((fixture, report));
            }
            Err(err) => {
                eprintln!("  FAILED {}: {err}", fixture.file_stem);
                failures.push(FixtureFailure {
                    stem: fixture.file_stem.to_string(),
                    reason: err,
                });
            }
        }
    }
    let successes: Vec<_> = reports
        .iter()
        .map(|(fixture, report)| FixtureSuccess { fixture, report })
        .collect();
    let manifest_path = write_run_manifest(
        &output_dir,
        selection.preset.label(),
        &selected_stems,
        &successes,
        &failures,
    )
    .expect("write PnS TOMS run manifest");
    println!("manifest {}", manifest_path.display());

    if !failures.is_empty() {
        panic!("{} PnS TOMS fixture export(s) failed", failures.len());
    }
    println!("=== Done ===");
}

fn export_fixture(fixture: &Fixture, output_dir: &Path) -> Result<FixtureReport, String> {
    let net = ControlNet::from_obj_str(fixture.obj).map_err(|err| err.to_string())?;
    let options = PolyhedralCoverageOptions::default();
    let complex =
        build_control_net_surface_complex(&net, options).map_err(|err| err.to_string())?;
    if complex.candidate_splines == 0 {
        return Err("control net did not expose any supported PnS spline family".into());
    }
    let rejections: Vec<String> = complex
        .rejections
        .iter()
        .map(|rejection| format!("spline {} {}", rejection.index, rejection.reason))
        .collect();
    let breps: Vec<BRep> = complex
        .components
        .iter()
        .map(|component| component.brep.clone())
        .collect();
    if breps.is_empty() {
        return Err(format!(
            "all {} candidate splines were rejected",
            complex.candidate_splines
        ));
    }
    if matches!(fixture.mode, FixtureMode::Strict) && !rejections.is_empty() {
        return Err(format!(
            "{} strict rejection(s): {}",
            rejections.len(),
            rejections.join("; ")
        ));
    }
    let sewing = analyze_surface_complex_sewing(&complex, options.c0_tolerance);

    let path = output_dir.join(format!("{}.step", fixture.file_stem));
    let expected_bodies = breps.len();
    export_step_file(
        &step_document(fixture, breps),
        &path,
        &StepExportOptions::default(),
    )
    .map_err(|err| err.to_string())?;
    let bytes = std::fs::metadata(&path)
        .map_err(|err| err.to_string())?
        .len();
    if bytes < 1024 {
        return Err(format!("STEP output is unexpectedly small: {bytes} bytes"));
    }
    let validation = validate_step_file(fixture.file_stem, &path, expected_bodies)?;
    let welded = weld_clean_surface_complex_boundaries(&complex, options)
        .map_err(|err| format!("welded export failed: {err}"))?;
    let (
        welded_path,
        welded_bytes,
        welded_quality_path,
        welded_repair_path,
        welded_validated_bodies,
    ) = if welded.report.welded_pairs > 0 {
        let welded_path = output_dir.join(format!("{}_welded.step", fixture.file_stem));
        let welded_quality_path =
            output_dir.join(format!("{}_welded_quality.json", fixture.file_stem));
        let welded_repair_path =
            output_dir.join(format!("{}_welded_repair.json", fixture.file_stem));
        write_quality_artifact(&welded.brep, options, &welded_quality_path)?;
        write_repair_action_artifact(
            &welded_repair_path,
            &format!("{}_welded", fixture.file_stem),
            &welded.repair_actions,
        )?;
        export_step_file(
            &step_document(fixture, vec![welded.brep]),
            &welded_path,
            &StepExportOptions::default(),
        )
        .map_err(|err| err.to_string())?;
        let welded_bytes = std::fs::metadata(&welded_path)
            .map_err(|err| err.to_string())?
            .len();
        let welded_validation =
            validate_step_file(&format!("{}_welded", fixture.file_stem), &welded_path, 1)?;
        (
            Some(welded_path),
            Some(welded_bytes),
            Some(welded_quality_path),
            Some(welded_repair_path),
            Some(welded_validation.imported_bodies),
        )
    } else {
        (None, None, None, None, None)
    };
    let report = FixtureReport {
        source_vertices: complex.source_vertices,
        source_faces: complex.source_faces,
        total_splines: complex.candidate_splines,
        exported_splines: complex.components.len(),
        rejected_splines: rejections.len(),
        patches: complex
            .components
            .iter()
            .map(|component| component.patches)
            .sum(),
        faces: complex
            .components
            .iter()
            .map(|component| component.brep.faces.len())
            .sum(),
        edges: complex
            .components
            .iter()
            .map(|component| component.brep.edges.len())
            .sum(),
        boundary_edges: complex
            .components
            .iter()
            .map(|component| component.boundary_edges)
            .sum(),
        matched_boundary_pairs: sewing.matched_boundary_pairs,
        unmatched_boundary_edges: sewing.unmatched_boundary_edges,
        max_matched_endpoint_gap: sewing.max_matched_endpoint_gap,
        welded_pairs: welded.report.welded_pairs,
        welded_remaining_boundary_edges: welded.report.remaining_boundary_edges,
        welded_export_ready: welded.report.export_ready,
        welded_expected_open_sheet_boundaries: welded.report.expected_open_sheet_boundaries,
        welded_export_blockers: welded.report.export_blockers,
        welded_edge_curve_repairs: welded.report.edge_curve_repairs,
        welded_edge_curve_declines: welded.report.edge_curve_declines,
        welded_edge_curve_decline_reasons: welded
            .report
            .edge_curve_decline_reasons
            .iter()
            .map(|(reason, count)| format!("{reason:?}:{count}"))
            .collect::<Vec<_>>()
            .join(","),
        welded_classification_buckets: classification_bucket_labels(
            &welded.report.classification_buckets,
        ),
        welded_topology_warnings: welded.report.topology_warnings,
        welded_topology_errors: welded.report.topology_errors,
        welded_missing_edge_curves: welded.report.missing_edge_curves,
        welded_trim_blockers: welded.report.trim_blockers,
        welded_scope_components: welded.report.scope.components.len(),
        welded_scope_edges: welded.report.scope.edges.len(),
        welded_scope_vertices: welded.report.scope.vertices.len(),
        welded_scope_faces: welded.report.scope.faces.len(),
        welded_timing: welded.report.timing,
        welded_bytes,
        welded_path,
        welded_quality_path,
        welded_repair_path,
        welded_validated_bodies,
        bytes,
        path,
        validated_bodies: validation.imported_bodies,
    };
    write_sidecar_report(fixture, &report, &rejections)?;

    Ok(report)
}

fn classification_bucket_labels(
    buckets: &[(BrepSewingClassificationBucket, usize)],
) -> Vec<(String, usize)> {
    buckets
        .iter()
        .map(|(bucket, count)| {
            (
                brep_sewing_classification_bucket_label(*bucket).to_string(),
                *count,
            )
        })
        .collect()
}

fn step_document(fixture: &Fixture, breps: Vec<BRep>) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: fixture.name.into(),
            id: format!("PNS-TOMS-{}", fixture.file_stem.to_uppercase()),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: breps,
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

fn output_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output")
        .join("pns_toms")
}

fn print_report(fixture: &Fixture, report: &FixtureReport) {
    let label = match fixture.mode {
        FixtureMode::Strict => "strict",
        FixtureMode::Partial => "partial",
    };
    println!(
        "  {} ({}) -> {}/{} splines, {} patches, {} faces, {} edges, {} boundary edges, {} bytes, validated {} body/bodies",
        fixture.file_stem,
        label,
        report.exported_splines,
        report.total_splines,
        report.patches,
        report.faces,
        report.edges,
        report.boundary_edges,
        report.bytes,
        report.validated_bodies
    );
    println!("    {}", report.path.display());
    if let Some(path) = &report.welded_path {
        println!(
            "    welded {} (validated={}, ready={}, open_sheet={}, blockers={}, edge_curves +{}/-{} [{}])",
            path.display(),
            report.welded_validated_bodies.unwrap_or(0),
            report.welded_export_ready,
            report.welded_expected_open_sheet_boundaries,
            report.welded_export_blockers,
            report.welded_edge_curve_repairs,
            report.welded_edge_curve_declines,
            report.welded_edge_curve_decline_reasons
        );
    }
    if let Some(path) = &report.welded_quality_path {
        println!("    quality {}", path.display());
    }
    if let Some(path) = &report.welded_repair_path {
        println!("    repair {}", path.display());
    }
}

fn write_quality_artifact(
    brep: &BRep,
    options: PolyhedralCoverageOptions,
    path: &Path,
) -> Result<(), String> {
    let quality = analyze_brep_export_quality(
        brep,
        BrepExportQualityOptions {
            topology_output: BrepTopologyOutputPolicy::OpenSheet,
            trim: BrepTrimRepairOptions {
                endpoint_tolerance: options.c0_tolerance,
                sampled_tolerance: options.c0_tolerance,
                samples: options.samples,
                validate_output: true,
            },
            require_edge_curves: true,
            allow_collapsed_trim_edge_curves: true,
            max_pcurve_deviation: Some(options.c0_tolerance),
            max_edge_surface_deviation: Some(options.c0_tolerance),
            validate_brep: true,
        },
    );
    let artifact = build_brep_export_quality_artifact(brep, &quality);
    let json = serde_json::to_string_pretty(&artifact).map_err(|err| err.to_string())?;
    std::fs::write(path, json).map_err(|err| err.to_string())
}
