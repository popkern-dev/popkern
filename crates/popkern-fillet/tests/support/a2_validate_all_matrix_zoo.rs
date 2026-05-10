use crate::a2_validate_helpers::{a2_fallback_reasons, build_under_a2_with_pyramid_height};
use popkern_core::diagnostics::Diagnostics;

#[test]
#[ignore = "diagnostic gate: exact fillet_zoo pyramid_n04 still leaves retired apex/support topology"]
fn a2_diagnostic_pyramid4_zoo_mixed_export_defects() {
    let result = build_under_a2_with_pyramid_height(true, &[0.07, 0.09, 0.11, 0.13], 1.2);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "diagnostic expects native A2 output before export-defect assertions; got:\n{}",
        reasons.join("\n")
    );

    let diagnostics = Diagnostics::run(&result.brep);
    eprintln!("{diagnostics}");
    assert!(
        diagnostics.defects.degenerate_edges.is_empty(),
        "native pyramid4_zoo_mixed has degenerate edges: {:#?}",
        diagnostics.defects.degenerate_edges
    );
    assert!(
        diagnostics.defects.sliver_faces.is_empty(),
        "native pyramid4_zoo_mixed has sliver faces: {:#?}",
        diagnostics.defects.sliver_faces
    );
    assert!(
        diagnostics.defects.face_overlaps.is_empty(),
        "native pyramid4_zoo_mixed has face overlaps: {:#?}",
        diagnostics.defects.face_overlaps
    );
    assert_eq!(
        diagnostics.defects.tess_non_manifold_count, 0,
        "native pyramid4_zoo_mixed has {} tessellation non-manifold seams, max {:.3e}",
        diagnostics.defects.tess_non_manifold_count, diagnostics.defects.tess_non_manifold_max
    );
}
