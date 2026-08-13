use popkern_core::geometry::polyhedral_spline::{
    ControlNet, PolyhedralCoverageOptions, analyze_surface_complex_sewing,
    build_control_net_surface_complex, weld_clean_surface_complex_boundaries,
};
use popkern_core::topology::BrepSewingClassificationBucket;

struct FixtureExpectation {
    name: &'static str,
    obj: &'static str,
    candidates: usize,
    accepted: usize,
    rejected: usize,
    patches: usize,
}

const SMOKE_FIXTURES: &[FixtureExpectation] = &[
    FixtureExpectation {
        name: "t0",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
        candidates: 1,
        accepted: 1,
        rejected: 0,
        patches: 4,
    },
    FixtureExpectation {
        name: "t1",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
        candidates: 1,
        accepted: 1,
        rejected: 0,
        patches: 8,
    },
    FixtureExpectation {
        name: "t2",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
        candidates: 1,
        accepted: 1,
        rejected: 0,
        patches: 16,
    },
    FixtureExpectation {
        name: "tee",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/tee.obj"),
        candidates: 24,
        accepted: 24,
        rejected: 0,
        // 88 since the PK-198 tiling changes: regular-quad components
        // ship 4 quarter patches (was 1) so mixed complexes weld 1:1.
        patches: 88,
    },
    FixtureExpectation {
        name: "cylinder_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cylinder.obj"),
        candidates: 20,
        accepted: 19,
        rejected: 1,
        patches: 78,
    },
];

const HEAVY_FIXTURES: &[FixtureExpectation] = &[
    FixtureExpectation {
        name: "suzanne",
        obj: include_str!(
            "../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/suzanne_all_configurations.obj"
        ),
        candidates: 471,
        accepted: 471,
        rejected: 0,
        patches: 7648,
    },
    FixtureExpectation {
        name: "hand_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/hand.obj"),
        candidates: 380,
        accepted: 370,
        rejected: 10,
        patches: 5947,
    },
    FixtureExpectation {
        name: "airplane_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/airplane.obj"),
        candidates: 868,
        accepted: 866,
        rejected: 2,
        patches: 3496,
    },
];

#[test]
fn pns_toms_fixture_matrix_matches_expected_coverage() {
    assert_fixture_matrix(SMOKE_FIXTURES);
}

#[test]
#[ignore = "heavy: suzanne/hand/airplane full PnS matrix (1719 candidate splines, 17091 patches after the PK-198 tiling unification) + weld/sew analysis; nightly tier, PK-190 owns the CI wiring — the SMOKE matrix covers the same code paths by default"]
fn pns_toms_heavy_fixture_matrix_matches_expected_coverage() {
    assert_fixture_matrix(HEAVY_FIXTURES);
}

fn assert_fixture_matrix(fixtures: &[FixtureExpectation]) {
    let options = PolyhedralCoverageOptions::default();
    for fixture in fixtures {
        let net = ControlNet::from_obj_str(fixture.obj).unwrap();
        let complex = build_control_net_surface_complex(&net, options).unwrap();
        assert_eq!(
            complex.candidate_splines, fixture.candidates,
            "{} candidates",
            fixture.name
        );
        assert_eq!(
            complex.components.len(),
            fixture.accepted,
            "{} accepted",
            fixture.name
        );
        assert_eq!(
            complex.rejections.len(),
            fixture.rejected,
            "{} rejected",
            fixture.name
        );
        let patches: usize = complex
            .components
            .iter()
            .map(|component| component.patches)
            .sum();
        assert_eq!(patches, fixture.patches, "{} patches", fixture.name);
        let sewing = analyze_surface_complex_sewing(&complex, options.c0_tolerance);
        assert_eq!(
            sewing.components, fixture.accepted,
            "{} sewing components",
            fixture.name
        );
        if fixture.name == "tee" {
            let welded = weld_clean_surface_complex_boundaries(&complex, options).unwrap();
            assert!(welded.report.export_ready, "tee welded export readiness");
            assert!(
                welded.repair_actions.summary.authored_actions >= welded.report.welded_pairs,
                "tee repair actions expose at least welded pair authorship"
            );
            assert_eq!(
                welded.repair_actions.actions.len(),
                welded.repair_actions.summary.attempted_actions,
                "tee repair artifact action count"
            );
            assert_weld_timing_is_recorded("tee", &welded.report);
            assert_weld_scope_covers_candidates("tee", &welded.report);
            assert_repair_artifact_writes_json("tee_welded", &welded.repair_actions);
            // Since the PK-198 tiling split the tee's mixed
            // regular/star interfaces weld 1:1 and the cage closes
            // fully — its historical open boundaries were the
            // dual-cell/quarter-cell mismatch.
            assert_eq!(
                welded.report.remaining_boundary_edges, 0,
                "tee closes fully"
            );
            assert!(
                !welded.report.expected_open_sheet_boundaries,
                "tee no longer expects open-sheet boundaries"
            );
            assert_eq!(welded.report.export_blockers, 0, "tee export blockers");
            assert_eq!(
                welded.report.missing_edge_curves, 0,
                "tee missing edge curves"
            );
            assert_eq!(welded.report.trim_blockers, 0, "tee trim blockers");
        }
        if fixture.name == "suzanne" {
            let welded = weld_clean_surface_complex_boundaries(&complex, options).unwrap();
            assert!(
                welded.report.export_ready,
                "suzanne welded export readiness"
            );
            assert_eq!(
                welded.report.edge_curve_repairs, 16,
                "suzanne completed missing edge curves"
            );
            assert_weld_timing_is_recorded("suzanne", &welded.report);
            assert_weld_scope_covers_candidates("suzanne", &welded.report);
            assert_eq!(
                welded.report.missing_edge_curves, 0,
                "suzanne missing edge curves"
            );
            assert_eq!(welded.report.export_blockers, 0, "suzanne blockers");
        }
        if fixture.name == "cylinder_partial" {
            let welded = weld_clean_surface_complex_boundaries(&complex, options).unwrap();
            assert!(
                welded.report.export_ready,
                "cylinder collapsed-trim export readiness"
            );
            assert_eq!(welded.report.export_blockers, 0, "cylinder blockers");
            assert_weld_timing_is_recorded("cylinder", &welded.report);
            assert_weld_scope_covers_candidates("cylinder", &welded.report);
            assert!(
                welded.report.scope.components.len() < complex.components.len(),
                "cylinder weld scope should exclude unrelated components"
            );
            assert_eq!(
                welded.report.missing_edge_curves, 4,
                "cylinder remaining missing edge curves"
            );
            assert_eq!(
                welded.report.degenerate_edges, 6,
                "cylinder total degenerate edge classification"
            );
            assert!(
                welded.report.collapsed_surface_trims >= welded.report.missing_edge_curves,
                "cylinder collapsed surface trims explain missing edge curves"
            );
            assert_eq!(
                welded.report.ambiguous_degenerate_edges, 0,
                "cylinder ambiguous degenerate edges"
            );
        }
    }
}

fn assert_weld_scope_covers_candidates(
    label: &str,
    report: &popkern_core::geometry::polyhedral_spline::PolyhedralBoundaryWeldReport,
) {
    assert_eq!(
        report.scope.candidate_pairs, report.candidate_pairs,
        "{label} weld scope candidate count must match report"
    );
    assert_eq!(
        report.scope.edges.len(),
        report.candidate_pairs * 2,
        "{label} weld scope should include both edges of every candidate pair"
    );
    assert!(
        !report.scope.components.is_empty(),
        "{label} weld scope should name touched source components"
    );
    assert!(
        !report.scope.faces.is_empty() && !report.scope.loops.is_empty(),
        "{label} weld scope should include adjacent loops and faces"
    );
    assert_eq!(
        report.pre_repair_scoped_quality.edges,
        report.scope.edges.len(),
        "{label} scoped quality should measure every scoped edge"
    );
    assert_eq!(
        report.pre_repair_scoped_quality.unavailable_edges, 0,
        "{label} scoped quality should not reference stale pre-repair edges"
    );
    assert!(
        report.pre_repair_scoped_quality.pcurve_blockers <= report.pre_repair_scoped_quality.edges,
        "{label} scoped pcurve blockers cannot exceed scoped edges"
    );
    assert!(
        report.post_repair_scoped_quality.edges > 0,
        "{label} post-repair scoped quality should retain touched repaired edges"
    );
    assert_classification_bucket(
        report,
        BrepSewingClassificationBucket::ExactClean,
        report.candidate_pairs,
    );
    assert_eq!(
        report.post_repair_scoped_quality.unavailable_edges, 0,
        "{label} post-repair scoped quality should not include deleted sewn edges"
    );
    assert_eq!(
        report.post_repair_scoped_quality.pcurve_blockers, 0,
        "{label} post-repair scoped pcurve blockers"
    );
}

fn assert_classification_bucket(
    report: &popkern_core::geometry::polyhedral_spline::PolyhedralBoundaryWeldReport,
    bucket: BrepSewingClassificationBucket,
    count: usize,
) {
    assert_eq!(
        report
            .classification_buckets
            .iter()
            .find(|(candidate_bucket, _)| *candidate_bucket == bucket)
            .map(|(_, bucket_count)| *bucket_count),
        Some(count),
        "weld classification bucket count"
    );
}

fn assert_weld_timing_is_recorded(
    label: &str,
    report: &popkern_core::geometry::polyhedral_spline::PolyhedralBoundaryWeldReport,
) {
    assert!(
        report.timing.total_ms
            >= report.timing.repair_pipeline_ms + report.timing.export_quality_ms,
        "{label} weld timing must include repair and export-quality stages"
    );
    assert!(
        report.timing.repair_pipeline_ms > 0 || report.timing.export_quality_ms > 0,
        "{label} weld timing should expose nonzero expensive stages"
    );
}

fn assert_repair_artifact_writes_json(
    label: &str,
    artifact: &popkern_core::topology::BrepRepairActionArtifact,
) {
    let path = std::env::temp_dir().join(format!(
        "popkern_{label}_repair_{}.json",
        std::process::id()
    ));
    popkern_examples::repair_artifacts::write_repair_action_artifact(&path, label, artifact)
        .unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let payload: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    let _ = std::fs::remove_file(path);
    assert_eq!(payload["schema"], "popkern.repair_actions");
    assert_eq!(payload["schema_version"], 1);
    assert_eq!(payload["label"], label);
    assert_eq!(
        payload["artifact"]["summary"]["attempted_actions"],
        artifact.summary.attempted_actions
    );
    assert!(
        payload["artifact"]["actions"]
            .as_array()
            .unwrap()
            .iter()
            .any(|action| action["affected_handles"]
                .as_array()
                .is_some_and(|handles| {
                    handles
                        .iter()
                        .filter_map(|handle| handle.as_str())
                        .filter(|handle| handle.starts_with("edge:"))
                        .count()
                        >= 2
                })),
        "repair artifact should expose affected edge handles"
    );
}
