use super::*;

#[test]
fn debug_artifact_reports_topology_and_obj_edges() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    let artifact = build_polyhedral_spline_debug_artifact(&splines[0], 11, 1e-8).unwrap();
    assert!(artifact.summary.contains("family=T1"));
    assert!(artifact.summary.contains("seams="));
    assert!(artifact.summary.contains("boundary_edges="));
    assert!(artifact.summary.contains("edge_pair_tests="));
    assert!(artifact.obj.contains("o patch_0"));
    assert!(artifact.obj.contains("l "));
}

#[test]
fn work_estimate_is_quadratic_in_patch_pairs_and_bounded_for_reference_matrix() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon8.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    let estimate = estimate_polyhedral_spline_work(&splines[0], 11);
    assert_eq!(estimate.patch_count, 32);
    assert_eq!(estimate.patch_edge_count, 128);
    assert_eq!(estimate.edge_pair_tests, 32 * 31 / 2 * 16);
    assert!(estimate.sampled_edge_evaluations < 100_000);
}

#[test]
fn adversarial_scaled_reference_fixture_still_validates_and_exports() {
    let mut net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop5sct.obj"
    ))
    .unwrap();
    for vertex in &mut net.vertices {
        vertex.point.x *= 1.0e4;
        vertex.point.y = vertex.point.y * -2.5e3 + 17.0;
        vertex.point.z = vertex.point.z * 750.0 - 3.0;
    }
    let splines = construct_from_control_net(&net).unwrap();
    validate_polyhedral_spline(&splines[0], 11, 1e-4, 1e-6, 1e-12).unwrap();
    build_polyhedral_spline_brep(
        &splines[0],
        PolyhedralBrepExportOptions {
            c0_tolerance: 1e-4,
            ..Default::default()
        },
    )
    .unwrap();
}
