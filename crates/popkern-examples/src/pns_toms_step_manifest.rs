use std::path::{Path, PathBuf};

use popkern_examples::example_manifests::{ExampleRunManifestInput, write_example_run_manifest};

use super::{Fixture, FixtureMode, FixtureReport};

const RUNNER_ID: &str = "pns_toms_step";

pub(super) struct FixtureSuccess<'a> {
    pub fixture: &'a Fixture,
    pub report: &'a FixtureReport,
}

pub(super) struct FixtureFailure {
    pub stem: String,
    pub reason: String,
}

pub(super) fn write_run_manifest(
    output_dir: &Path,
    preset: &str,
    selected_stems: &[&str],
    successes: &[FixtureSuccess<'_>],
    failures: &[FixtureFailure],
) -> Result<PathBuf, String> {
    write_example_run_manifest(
        output_dir,
        "run_manifest.json",
        ExampleRunManifestInput {
            runner: RUNNER_ID,
            preset,
            selected_fixtures: selected_stems,
            fixtures: successes.iter().map(success_payload).collect(),
            failures: failures.iter().map(failure_payload).collect(),
            metadata: serde_json::json!({
                "domain": "ACM TOMS Algorithm 1032 PnS STEP examples",
                "welded_outputs": true,
                "step_validation": true,
            }),
        },
    )
}

fn success_payload(success: &FixtureSuccess<'_>) -> serde_json::Value {
    let report = success.report;
    serde_json::json!({
        "name": success.fixture.name,
        "stem": success.fixture.file_stem,
        "mode": fixture_mode_label(success.fixture.mode),
        "source": {
            "vertices": report.source_vertices,
            "faces": report.source_faces,
            "candidate_splines": report.total_splines,
            "exported_splines": report.exported_splines,
            "rejected_splines": report.rejected_splines,
        },
        "surface_complex": {
            "patches": report.patches,
            "faces": report.faces,
            "edges": report.edges,
            "boundary_edges": report.boundary_edges,
            "matched_boundary_pairs": report.matched_boundary_pairs,
            "unmatched_boundary_edges": report.unmatched_boundary_edges,
            "max_matched_endpoint_gap": report.max_matched_endpoint_gap,
        },
        "step": {
            "path": path_label(&report.path),
            "bytes": report.bytes,
            "validated_bodies": report.validated_bodies,
        },
        "welded": {
            "step_path": optional_path_label(&report.welded_path),
            "quality_artifact": optional_path_label(&report.welded_quality_path),
            "repair_artifact": optional_path_label(&report.welded_repair_path),
            "bytes": report.welded_bytes,
            "validated_bodies": report.welded_validated_bodies,
            "pairs": report.welded_pairs,
            "remaining_boundary_edges": report.welded_remaining_boundary_edges,
            "export_ready": report.welded_export_ready,
            "expected_open_sheet_boundaries": report.welded_expected_open_sheet_boundaries,
            "export_blockers": report.welded_export_blockers,
            "edge_curve_repairs": report.welded_edge_curve_repairs,
            "edge_curve_declines": report.welded_edge_curve_declines,
            "edge_curve_decline_reasons": report.welded_edge_curve_decline_reasons,
            "topology_warnings": report.welded_topology_warnings,
            "topology_errors": report.welded_topology_errors,
            "missing_edge_curves": report.welded_missing_edge_curves,
            "trim_blockers": report.welded_trim_blockers,
            "scope": {
                "components": report.welded_scope_components,
                "edges": report.welded_scope_edges,
                "vertices": report.welded_scope_vertices,
                "faces": report.welded_scope_faces,
            },
            "timing_ms": {
                "assemble": report.welded_timing.assemble_ms,
                "candidate_match": report.welded_timing.candidate_match_ms,
                "repair_pipeline": report.welded_timing.repair_pipeline_ms,
                "export_quality": report.welded_timing.export_quality_ms,
                "repair_artifact": report.welded_timing.repair_artifact_ms,
                "total": report.welded_timing.total_ms,
            },
        },
    })
}

fn failure_payload(failure: &FixtureFailure) -> serde_json::Value {
    serde_json::json!({
        "stem": failure.stem,
        "reason": failure.reason,
    })
}

fn fixture_mode_label(mode: FixtureMode) -> &'static str {
    match mode {
        FixtureMode::Strict => "strict",
        FixtureMode::Partial => "partial",
    }
}

fn path_label(path: &Path) -> String {
    path.display().to_string()
}

fn optional_path_label(path: &Option<PathBuf>) -> Option<String> {
    path.as_ref().map(|path| path.display().to_string())
}
