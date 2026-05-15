use popkern_core::geometry::polyhedral_spline::{
    ControlNet, PolyhedralCoverageOptions, analyze_surface_complex_sewing,
    build_control_net_surface_complex, weld_clean_surface_complex_boundaries,
};

struct FixtureExpectation {
    name: &'static str,
    obj: &'static str,
    candidates: usize,
    accepted: usize,
    rejected: usize,
    patches: usize,
}

const FIXTURES: &[FixtureExpectation] = &[
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
        patches: 64,
    },
    FixtureExpectation {
        name: "suzanne",
        obj: include_str!(
            "../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/suzanne_all_configurations.obj"
        ),
        candidates: 471,
        accepted: 471,
        rejected: 0,
        patches: 844,
    },
    FixtureExpectation {
        name: "cylinder_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cylinder.obj"),
        candidates: 20,
        accepted: 19,
        rejected: 1,
        patches: 24,
    },
    FixtureExpectation {
        name: "hand_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/hand.obj"),
        candidates: 380,
        accepted: 370,
        rejected: 10,
        patches: 532,
    },
    FixtureExpectation {
        name: "airplane_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/airplane.obj"),
        candidates: 868,
        accepted: 866,
        rejected: 2,
        patches: 1045,
    },
];

#[test]
fn pns_toms_fixture_matrix_matches_expected_coverage() {
    let options = PolyhedralCoverageOptions::default();
    for fixture in FIXTURES {
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
            assert_repair_artifact_writes_json("tee_welded", &welded.repair_actions);
            assert!(
                welded.report.expected_open_sheet_boundaries,
                "tee expected open-sheet boundaries"
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
}
