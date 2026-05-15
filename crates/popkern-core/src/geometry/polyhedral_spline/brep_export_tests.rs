use super::*;
use crate::diagnostics::Diagnostics;
use crate::naming::OperationId;
use crate::tolerance::Tolerance;

#[test]
fn brep_export_builds_valid_regular_quad_sheet_with_pcurves() {
    let stencil: [crate::geometry::Point3; 9] = (0..9)
        .map(|i| crate::geometry::Point3::new(i as f64, (i / 3) as f64, 0.1 * i as f64))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let spline = construct_regular_quad(stencil).unwrap();
    let export = build_polyhedral_spline_brep(
        &spline,
        PolyhedralBrepExportOptions {
            tolerance: Tolerance::default(),
            operation: OperationId(42),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(export.faces.len(), 1);
    assert_eq!(export.brep.faces.len(), 1);
    assert_eq!(export.boundary_edge_count, 4);
    assert_eq!(export.seam_edge_count, 0);
    assert!(export.brep.validate().is_ok());
    assert!(
        export
            .brep
            .half_edges
            .values()
            .all(|he| he.pcurve.is_some())
    );
    let report = Diagnostics::run(&export.brep);
    assert!(report.defects.pcurve_escapes.is_empty());
}

#[test]
fn brep_export_builds_reference_t0_patch_sheet_with_shared_seams() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    let export =
        build_polyhedral_spline_brep(&splines[0], PolyhedralBrepExportOptions::default()).unwrap();
    assert_eq!(export.faces.len(), splines[0].patches.len());
    assert!(export.seam_edge_count >= 4);
    assert!(export.boundary_edge_count > 0);
    assert!(export.brep.validate().is_ok());
    assert!(
        export
            .brep
            .half_edges
            .values()
            .all(|he| he.pcurve.is_some())
    );
}

#[test]
fn brep_export_builds_polar_sheet_with_collapsed_singularities() {
    let net = ControlNet::from_obj_str(include_str!(
        "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar5sct.obj"
    ))
    .unwrap();
    let splines = construct_from_control_net(&net).unwrap();
    let export =
        build_polyhedral_spline_brep(&splines[0], PolyhedralBrepExportOptions::default()).unwrap();
    assert_eq!(export.faces.len(), splines[0].patches.len());
    assert_eq!(export.collapsed_singularity_count, 5);
    assert!(export.boundary_edge_count >= export.collapsed_singularity_count);
    assert!(export.brep.validate().is_ok());
}

#[test]
fn brep_export_declines_without_mutating_input_spline() {
    let mut spline = construct_regular_quad(
        (0..9)
            .map(|i| crate::geometry::Point3::new(i as f64, 0.0, 0.0))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
    .unwrap();
    let original = spline.clone();
    spline.patches[0].control_points[0].x = f64::NAN;
    assert!(build_polyhedral_spline_brep(&spline, PolyhedralBrepExportOptions::default()).is_err());
    assert_eq!(original.patches.len(), spline.patches.len());
    assert!(original.patches[0].control_points[0].x.is_finite());
}
