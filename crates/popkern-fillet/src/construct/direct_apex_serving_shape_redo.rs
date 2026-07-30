//! Serving-shape gate rung 3 (PK-219): apex-retaining redo with
//! keep-better semantics for the welded shared-corner class.

use popkern_core::topology::BRep;

use crate::apex_construction::constructive_builder::construct_apex_constructively;

use super::{ServedShapeDefect, run_absorb_consolidate_fixpoint, served_shape_defect};

/// Rung 3: when the served in-place body still carries a defect
/// after the handover completion and the aggregate redo, rebuild
/// through the apex-retaining options and serve whichever body has
/// the strictly better serving signature (measured pyramid4_uniform:
/// in-place welded strict 9 surfaceless/20 dead -> redo 0/8).
pub(super) fn apply_apex_retaining_redo_rung(
    brep: &mut BRep,
    snapshot: &BRep,
    plan: &crate::apex_construction::plan::ApexPlan,
    surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    result: &mut crate::apex_construction::construct::ApexBuildResult,
) {
    let Some(current_defect) = served_shape_defect(brep) else {
        return;
    };
    // An explicit strict V-free request fails loudly rather than
    // silently serving the apex-retaining body (same contract as
    // the PK-124 promoted-default fallback).
    if !super::vfree_fallback::v_free_promoted_by_default() {
        return;
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] serving-shape gate rung 3: {current_defect:?}; attempting \
             apex-retaining redo"
        );
    }
    let incumbent_body = brep.clone();
    let incumbent_result = result.clone();
    *brep = snapshot.clone();
    let retaining_options = crate::apex_construction::constructive_builder_options::
        ConstructiveApexExecutionOptions::from_environment_adapter()
            .with_v_free_repeated_apex_disabled();
    match construct_apex_constructively(brep, plan, surfaces, tolerance, None, retaining_options) {
        Ok(mut redo_result) => {
            run_absorb_consolidate_fixpoint(brep, &mut redo_result);
            if served_shape_defect(brep).is_some() {
                let retired =
                    crate::apex_construction::v_free_inplace_handover::retire_mirror_ghost_complex(
                        brep,
                    );
                if retired > 0 {
                    crate::apex_construction::v_free_shell_consolidation::
                        consolidate_vestigial_support_shell_forcing_winding(brep);
                    run_absorb_consolidate_fixpoint(brep, &mut redo_result);
                    crate::apex_construction::v_free_shell_consolidation::
                        complete_partnerless_chord_rings(brep);
                    run_absorb_consolidate_fixpoint(brep, &mut redo_result);
                    let mut sectors = std::mem::take(&mut redo_result.exact_boundary_sectors);
                    let _ = crate::apex_construction::v_free_inplace_handover::
                        remap_exact_boundary_sector_records(
                            brep,
                            &mut sectors,
                            &mut redo_result.exact_boundary_report.sectors,
                        );
                    redo_result.exact_boundary_sectors = sectors;
                }
            }
            let redo_defect = served_shape_defect(brep);
            let redo_wins = match redo_defect.as_ref() {
                None => brep.validate().is_ok(),
                Some(redo) => redo.strictly_better_than(&current_defect),
            };
            if redo_wins {
                redo_result.diagnostics.push(
                    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                        crate::types::ExactBoundaryCapDiagnostic {
                            stage: "constructive_serving_shape".to_string(),
                            code: "constructive_serving_shape.apex_retaining_redo".to_string(),
                            detail_code: Some(format!(
                                "in-place: {}; apex-retaining redo: {}",
                                current_defect.describe(),
                                redo_defect
                                    .as_ref()
                                    .map(ServedShapeDefect::describe)
                                    .unwrap_or_else(|| "clean".to_string())
                            )),
                            sector_name: "direct_apex".to_string(),
                            edge: None,
                            support_half_edge: None,
                            support_face: None,
                            residual: None,
                            outer_side: None,
                            sample_parameter: None,
                            tolerance: None,
                            shared_fact_ids: Vec::new(),
                        },
                    ),
                );
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] serving-shape gate rung 3: apex-retaining redo \
                         serves ({})",
                        redo_defect
                            .as_ref()
                            .map(ServedShapeDefect::describe)
                            .unwrap_or_else(|| "clean".to_string())
                    );
                }
                *result = redo_result;
            } else {
                *brep = incumbent_body;
                *result = incumbent_result;
                result.diagnostics.push(
                    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                        crate::types::ExactBoundaryCapDiagnostic {
                            stage: "constructive_serving_shape".to_string(),
                            code: "constructive_serving_shape.residual_defect_retained".to_string(),
                            detail_code: Some(format!(
                                "in-place: {}; apex-retaining redo not better: {}",
                                current_defect.describe(),
                                redo_defect
                                    .as_ref()
                                    .map(ServedShapeDefect::describe)
                                    .unwrap_or_else(|| "clean-but-invalid".to_string())
                            )),
                            sector_name: "direct_apex".to_string(),
                            edge: None,
                            support_half_edge: None,
                            support_face: None,
                            residual: None,
                            outer_side: None,
                            sample_parameter: None,
                            tolerance: None,
                            shared_fact_ids: Vec::new(),
                        },
                    ),
                );
            }
        }
        Err(err) => {
            *brep = incumbent_body;
            *result = incumbent_result;
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace] serving-shape gate rung 3: apex-retaining redo \
                     errored ({err}); keeping incumbent"
                );
            }
        }
    }
}
