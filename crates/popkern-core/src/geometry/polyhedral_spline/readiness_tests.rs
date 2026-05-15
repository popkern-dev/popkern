use super::*;
use crate::geometry::Point3;

fn point(i: usize) -> Point3 {
    Point3::new(i as f64, (i % 5) as f64 * 0.25, (i % 7) as f64 * 0.125)
}

#[test]
fn brep_readiness_assigns_every_reference_patch_edge_one_role() {
    for fixture in reference_fixture_matrix() {
        let net = ControlNet::from_obj_str(fixture).unwrap();
        let splines = construct_from_control_net(&net).unwrap();
        let report = validate_brep_readiness(&splines[0], 11, 1e-8)
            .unwrap_or_else(|err| panic!("{:?}: {err}", splines[0].family));
        assert_eq!(report.patch_count, splines[0].patches.len());
        assert_eq!(report.patch_domains.len(), splines[0].patches.len());
        assert_eq!(report.edge_uses.len(), splines[0].patches.len() * 4);
        assert!(
            report
                .edge_uses
                .iter()
                .any(|edge_use| matches!(edge_use.role, PatchEdgeRole::Boundary)),
            "{:?} should expose trim-boundary candidate edges",
            splines[0].family
        );
        assert!(
            report
                .edge_uses
                .iter()
                .any(|edge_use| matches!(edge_use.role, PatchEdgeRole::Seam { .. })),
            "{:?} should expose internal seam edges",
            splines[0].family
        );
    }
}

#[test]
fn brep_readiness_domains_are_unit_square_and_ordered() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    let report = validate_brep_readiness(&splines[0], 11, 1e-8).unwrap();
    for (index, domain) in report.patch_domains.iter().enumerate() {
        assert_eq!(domain.patch_index, index);
        assert_eq!((domain.u_min, domain.u_max), (0.0, 1.0));
        assert_eq!((domain.v_min, domain.v_max), (0.0, 1.0));
        domain.validate().unwrap();
    }
}

#[test]
fn brep_readiness_rejects_non_manifold_edge_role() {
    let mut spline =
        construct_regular_quad((0..9).map(point).collect::<Vec<_>>().try_into().unwrap()).unwrap();
    let patch = spline.patches[0].clone();
    spline.patches.push(patch.clone());
    spline.patches.push(patch);
    let err = validate_brep_readiness(&spline, 5, 1e-12).unwrap_err();
    assert!(matches!(
        err,
        PolyhedralSplineError::AtlasTopology { reason, .. }
            if reason.contains("multiple roles")
    ));
}

fn reference_fixture_matrix() -> [&'static str; 19] {
    [
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar3sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar4sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar5sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar6sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar7sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar8sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop3sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop5sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop6sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop7sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop8sct.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon3.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon5.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon6.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon7.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon8.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
        include_str!("../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
    ]
}
