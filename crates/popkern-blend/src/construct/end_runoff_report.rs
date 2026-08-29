//! End-runoff execution + structured diagnostic reporting.
//!
//! Split from `construct::mod`'s apply orchestration (file-size
//! ratchet): runs the fillet end-runoff pass and promotes each
//! outcome into the repair report's apex diagnostics.

use popkern_core::topology::BRep;

/// Run the end-runoff pass and record one
/// `fillet_end_runoff_applied` / `_declined` diagnostic per outcome.
pub(super) fn apply_and_report(
    result: &mut BRep,
    setback_blend_present: bool,
    repair_report: &mut crate::types::BlendRepairReport,
) {
    let runoff_outcomes = super::end_runoff::apply_fillet_end_runoff(result, setback_blend_present);
    for o in &runoff_outcomes {
        repair_report.apex_diagnostics.push(
            crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                crate::types::ExactBoundaryCapDiagnostic {
                    stage: "end_runoff".into(),
                    code: if o.error.is_none() {
                        "fillet_end_runoff_applied".into()
                    } else {
                        "fillet_end_runoff_declined".into()
                    },
                    detail_code: o.error.clone(),
                    sector_name: format!("{:?}->{:?}", o.strip_face, o.end_face),
                    edge: None,
                    support_half_edge: None,
                    support_face: Some(o.end_face),
                    residual: Some(o.overshoot),
                    outer_side: None,
                    sample_parameter: None,
                    tolerance: None,
                    shared_fact_ids: Vec::new(),
                },
            ),
        );
    }
}
