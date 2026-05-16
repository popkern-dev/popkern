use super::{Fixture, FixtureMode, FixtureReport};

pub(super) fn write_sidecar_report(
    fixture: &Fixture,
    report: &FixtureReport,
    rejections: &[String],
) -> Result<(), String> {
    let mut out = String::new();
    out.push_str(&format!("fixture: {}\n", fixture.name));
    out.push_str(&format!("file_stem: {}\n", fixture.file_stem));
    out.push_str(&format!("mode: {}\n", fixture_mode_label(fixture.mode)));
    out.push_str(&format!("source_vertices: {}\n", report.source_vertices));
    out.push_str(&format!("source_faces: {}\n", report.source_faces));
    out.push_str(&format!("total_splines: {}\n", report.total_splines));
    out.push_str(&format!("exported_splines: {}\n", report.exported_splines));
    out.push_str(&format!("rejected_splines: {}\n", report.rejected_splines));
    out.push_str(&format!("patches: {}\n", report.patches));
    out.push_str(&format!("faces: {}\n", report.faces));
    out.push_str(&format!("edges: {}\n", report.edges));
    out.push_str(&format!("boundary_edges: {}\n", report.boundary_edges));
    out.push_str(&format!(
        "matched_boundary_pairs: {}\n",
        report.matched_boundary_pairs
    ));
    out.push_str(&format!(
        "unmatched_boundary_edges: {}\n",
        report.unmatched_boundary_edges
    ));
    out.push_str(&format!(
        "max_matched_endpoint_gap: {}\n",
        report.max_matched_endpoint_gap
    ));
    out.push_str(&format!("welded_pairs: {}\n", report.welded_pairs));
    out.push_str(&format!(
        "welded_remaining_boundary_edges: {}\n",
        report.welded_remaining_boundary_edges
    ));
    out.push_str(&format!(
        "welded_export_ready: {}\n",
        report.welded_export_ready
    ));
    out.push_str(&format!(
        "welded_expected_open_sheet_boundaries: {}\n",
        report.welded_expected_open_sheet_boundaries
    ));
    out.push_str(&format!(
        "welded_export_blockers: {}\n",
        report.welded_export_blockers
    ));
    out.push_str(&format!(
        "welded_edge_curve_repairs: {}\n",
        report.welded_edge_curve_repairs
    ));
    out.push_str(&format!(
        "welded_edge_curve_declines: {}\n",
        report.welded_edge_curve_declines
    ));
    out.push_str(&format!(
        "welded_edge_curve_decline_reasons: {}\n",
        report.welded_edge_curve_decline_reasons
    ));
    out.push_str(&format!(
        "welded_topology_warnings: {}\n",
        report.welded_topology_warnings
    ));
    out.push_str(&format!(
        "welded_topology_errors: {}\n",
        report.welded_topology_errors
    ));
    out.push_str(&format!(
        "welded_missing_edge_curves: {}\n",
        report.welded_missing_edge_curves
    ));
    out.push_str(&format!(
        "welded_trim_blockers: {}\n",
        report.welded_trim_blockers
    ));
    out.push_str(&format!(
        "welded_scope: components={} edges={} vertices={} faces={}\n",
        report.welded_scope_components,
        report.welded_scope_edges,
        report.welded_scope_vertices,
        report.welded_scope_faces
    ));
    out.push_str(&format!(
        "welded_timing_ms: assemble={} candidate_match={} repair_pipeline={} export_quality={} repair_artifact={} total={}\n",
        report.welded_timing.assemble_ms,
        report.welded_timing.candidate_match_ms,
        report.welded_timing.repair_pipeline_ms,
        report.welded_timing.export_quality_ms,
        report.welded_timing.repair_artifact_ms,
        report.welded_timing.total_ms
    ));
    if let Some(bytes) = report.welded_bytes {
        out.push_str(&format!("welded_step_bytes: {bytes}\n"));
    }
    if let Some(path) = &report.welded_quality_path {
        out.push_str(&format!("welded_quality_artifact: {}\n", path.display()));
    }
    if let Some(path) = &report.welded_repair_path {
        out.push_str(&format!("welded_repair_artifact: {}\n", path.display()));
    }
    if let Some(imported) = report.welded_validated_bodies {
        out.push_str(&format!("welded_validated_bodies: {imported}\n"));
    }
    out.push_str(&format!("step_bytes: {}\n", report.bytes));
    out.push_str(&format!("validated_bodies: {}\n", report.validated_bodies));
    out.push_str("rejections:\n");
    if rejections.is_empty() {
        out.push_str("  none\n");
    } else {
        for rejection in rejections {
            out.push_str(&format!("  - {rejection}\n"));
        }
    }
    std::fs::write(report.path.with_extension("report.txt"), out).map_err(|err| err.to_string())
}

fn fixture_mode_label(mode: FixtureMode) -> &'static str {
    match mode {
        FixtureMode::Strict => "strict",
        FixtureMode::Partial => "partial",
    }
}
