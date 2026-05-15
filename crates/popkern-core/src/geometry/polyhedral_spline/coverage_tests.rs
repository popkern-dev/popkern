use super::*;

#[test]
fn coverage_report_accepts_reference_t0_fixture() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"
    ))
    .unwrap();
    let report = analyze_control_net_coverage(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(report.candidate_splines, 1);
    assert_eq!(report.accepted_splines, 1);
    assert_eq!(report.rejected_splines, 0);
    assert_eq!(report.accepted_patches, 4);
    assert_eq!(report.covered_anchor_faces.len(), 1);
}

#[test]
fn coverage_report_records_partial_cylinder_rejection() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cylinder.obj"
    ))
    .unwrap();
    let report = analyze_control_net_coverage(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(report.candidate_splines, 20);
    assert_eq!(report.accepted_splines, 19);
    assert_eq!(report.rejected_splines, 1);
    assert!(
        report
            .candidates
            .iter()
            .any(|candidate| matches!(candidate.status, PolyhedralCoverageStatus::Rejected { .. }))
    );
}
