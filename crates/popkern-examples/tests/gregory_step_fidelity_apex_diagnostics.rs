use popkern_fillet::{
    ApexConstructionDiagnostic, BlendRepairReport, CylinderTrimMaterializationDeclineReason,
};

pub(crate) fn assert_materialized_cylinder_trim_decline(report: &BlendRepairReport) {
    let diagnostic = report
        .apex_diagnostics
        .iter()
        .find_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::CylinderTrimMaterialization(diagnostic)
                if diagnostic.reason
                    == CylinderTrimMaterializationDeclineReason::NoPlanePlaneCylinderTransition =>
            {
                Some(diagnostic)
            }
            _ => None,
        })
        .expect(
            "Gregory apex repair report should diagnose exact cylinder-trim materialization decline",
        );

    assert!(
        !diagnostic.start_support_faces.is_empty(),
        "decline diagnostic should identify start endpoint support ownership: {diagnostic:#?}"
    );
    assert!(
        !diagnostic.end_support_faces.is_empty(),
        "decline diagnostic should identify end endpoint support ownership: {diagnostic:#?}"
    );
    assert!(
        !diagnostic.exact_authorable,
        "declined cylinder trim span should not be reported exact-authorable: {diagnostic:#?}"
    );
}
