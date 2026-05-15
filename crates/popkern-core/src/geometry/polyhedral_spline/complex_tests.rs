use super::*;

#[test]
fn surface_complex_accepts_t012_components() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T012.obj"
    ))
    .unwrap();
    let complex =
        build_control_net_surface_complex(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(complex.candidate_splines, 15);
    assert_eq!(complex.components.len(), 15);
    assert!(complex.rejections.is_empty());
    let sewing = analyze_surface_complex_sewing(&complex, 1.0e-6);
    assert_eq!(sewing.components, 15);
    assert!(sewing.boundary_edges > 0);
}

#[test]
fn surface_complex_preserves_partial_hand_rejections() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/hand.obj"
    ))
    .unwrap();
    let complex =
        build_control_net_surface_complex(&net, PolyhedralCoverageOptions::default()).unwrap();
    assert_eq!(complex.candidate_splines, 380);
    assert_eq!(complex.components.len(), 370);
    assert_eq!(complex.rejections.len(), 10);
    assert!(
        complex
            .rejections
            .iter()
            .all(|rejection| !rejection.reason.is_empty())
    );
}

#[test]
fn weld_clean_boundaries_reduces_tee_open_edges_transactionally() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/tee.obj"
    ))
    .unwrap();
    let options = PolyhedralCoverageOptions::default();
    let complex = build_control_net_surface_complex(&net, options).unwrap();
    let before = analyze_surface_complex_sewing(&complex, options.c0_tolerance);
    let welded = weld_clean_surface_complex_boundaries(&complex, options).unwrap();
    assert_eq!(welded.report.welded_pairs, before.matched_boundary_pairs);
    assert_eq!(
        welded.report.remaining_boundary_edges,
        before.unmatched_boundary_edges
    );
    assert!(welded.report.welded_pairs > 0);
    assert_eq!(welded.report.trim_repairs, 0);
    assert_eq!(welded.report.topology_errors, 0);
    assert!(welded.report.topology_warnings > 0);
    assert!(welded.report.export_ready);
    assert!(welded.report.expected_open_sheet_boundaries);
    assert_eq!(welded.report.export_blockers, 0);
    assert_eq!(welded.report.missing_edge_curves, 0);
    assert_eq!(welded.report.trim_blockers, 0);
    assert!(welded.brep.validate().is_ok());
    let original_boundary_edges: usize = complex
        .components
        .iter()
        .map(|component| component.brep.boundary_edges.len())
        .sum();
    assert_eq!(original_boundary_edges, before.boundary_edges);
}
