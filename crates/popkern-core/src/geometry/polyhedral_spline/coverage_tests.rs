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
fn coverage_report_accepts_full_cylinder() {
    // Historically one Polar{6} candidate rejected here with a 0.297 rad
    // seam NORMAL angle at a gap of 1.3e-10 — a pole-degenerate normal
    // artifact of the reference package's ~1e-9 table noise, not real
    // geometry. The regenerated closed-form polar masks (PK-317) are
    // exactly rotation-symmetric, the degenerate normals agree, and the
    // candidate legitimately accepts. Rejection *recording* stays covered
    // by `surface_complex_preserves_partial_hand_rejections`.
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cylinder.obj"
    ))
    .unwrap();
    let report = analyze_control_net_coverage(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(report.candidate_splines, 20);
    assert_eq!(report.accepted_splines, 20);
    assert_eq!(report.rejected_splines, 0);
}
