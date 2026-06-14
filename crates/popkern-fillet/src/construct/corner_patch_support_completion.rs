//! Result-level completion for Gregory cap/support seam ownership.

use std::borrow::Cow;

use popkern_core::topology::BRep;

use super::kp20_polar_support_cell_shared_trim_quality as shared_trim_quality;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo, CornerSubPatchRole};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_completion_apply.rs"]
mod support_cell_completion_apply;
#[path = "kp20_polar_support_cell_completion_rank.rs"]
mod support_cell_completion_rank;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

#[path = "kp20_polar_support_cell_completion_candidates.rs"]
mod support_cell_completion_candidates;
use support_cell_completion_candidates::SupportCellCompletionCandidate;
#[derive(Clone)]
struct SupportCellCompletionSelection {
    brep: BRep,
    candidate: SupportCellCompletionCandidate,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone)]
struct SupportCellCompletionCandidateScore {
    candidate: SupportCellCompletionCandidate,
    status: &'static str,
    operation_status: Cow<'static, str>,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
}

#[derive(Clone)]
pub(super) struct SupportCellCompletionApplication {
    operator_succeeded: bool,
    operation_status: Cow<'static, str>,
}

pub(crate) fn complete_corner_patch_support_pcurves(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    complete_support_cell_pcurves_transactionally(brep, corner_patches, diagnostics);
    let mut completed_cap_faces = Vec::new();
    for patch in corner_patches {
        if patch
            .sub_patch_provenance
            .is_some_and(|provenance| provenance.role == CornerSubPatchRole::SupportCell)
        {
            continue;
        }
        if completed_cap_faces.contains(&patch.face) {
            continue;
        }
        super::three_quad_support_pcurve_completion::complete_missing_cap_support_seams(
            brep, patch.face,
        );
        completed_cap_faces.push(patch.face);
    }
}

fn complete_support_cell_pcurves_transactionally(
    brep: &mut BRep,
    corner_patches: &[CornerPatchInfo],
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) {
    let selected = shared_trim_quality::support_cell_selected_label(diagnostics);
    for group in shared_trim_quality::support_cell_face_groups(corner_patches) {
        let before =
            shared_trim_quality::measure_support_cell_shared_trim_quality(brep, &group.faces, 9);
        let (selection, scores) =
            select_support_cell_completion_candidate(brep, &group.faces, before);
        let (accepted, selected_candidate, after) = if let Some(selection) = selection {
            let candidate = selection.candidate.label();
            let quality = selection.quality;
            *brep = selection.brep;
            (true, candidate, quality)
        } else {
            (false, "none", before)
        };
        for score in scores {
            diagnostics.push(
                shared_trim_quality::support_cell_shared_trim_candidate_fact(
                    group.n_sides,
                    selected,
                    score.candidate.label(),
                    score.status,
                    score.operation_status.as_ref(),
                    before,
                    score.quality,
                ),
            );
        }
        diagnostics.push(shared_trim_quality::support_cell_shared_trim_repair_fact(
            group.n_sides,
            selected,
            selected_candidate,
            accepted,
            before,
            after,
        ));
    }
}

fn select_support_cell_completion_candidate(
    brep: &BRep,
    faces: &[popkern_core::topology::FaceHandle],
    before: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> (
    Option<SupportCellCompletionSelection>,
    Vec<SupportCellCompletionCandidateScore>,
) {
    let mut best: Option<SupportCellCompletionSelection> = None;
    let mut scores = Vec::new();
    for candidate in SupportCellCompletionCandidate::ALL {
        let mut candidate_brep = brep.clone();

        let application =
            apply_support_cell_completion_candidate(&mut candidate_brep, faces, candidate);
        crate::construct::post_pipeline::determinism_trace(
            &format!("post_candidate_{}", candidate.label()),
            &candidate_brep,
        );
        if !application.operator_succeeded {
            scores.push(SupportCellCompletionCandidateScore {
                candidate,
                status: "operator_failed",
                operation_status: application.operation_status,
                quality: before,
            });
            continue;
        }
        let quality = shared_trim_quality::measure_support_cell_shared_trim_quality(
            &candidate_brep,
            faces,
            9,
        );
        let status = if candidate_brep.validate().is_err() {
            "invalid_brep"
        } else if !quality.is_no_worse_than(before) {
            "quality_regressed"
        } else {
            "accepted_not_selected"
        };
        scores.push(SupportCellCompletionCandidateScore {
            candidate,
            status,
            operation_status: application.operation_status,
            quality,
        });
        if status != "accepted_not_selected" {
            continue;
        }
        let selection = SupportCellCompletionSelection {
            brep: candidate_brep,
            candidate,
            quality,
        };
        if support_cell_completion_rank::selection_is_better(
            support_cell_rank(selection.candidate, selection.quality),
            best.as_ref()
                .map(|selection| support_cell_rank(selection.candidate, selection.quality)),
        ) {
            best = Some(selection);
        }
    }
    if let Some(selection) = best.as_ref() {
        for score in &mut scores {
            if score.candidate == selection.candidate {
                score.status = "accepted_selected";
            }
        }
    }
    (best, scores)
}

fn apply_support_cell_completion_candidate(
    brep: &mut BRep,
    faces: &[popkern_core::topology::FaceHandle],
    candidate: SupportCellCompletionCandidate,
) -> SupportCellCompletionApplication {
    match candidate {
        SupportCellCompletionCandidate::PcurveOnly => support_cell_operation_result(
            support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces).is_ok(),
            "pcurve_authored",
            "pcurve_operator_failed",
        ),
        SupportCellCompletionCandidate::BoundaryOnly => support_cell_operation_result(
            support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                brep, faces,
            )
            .is_ok(),
            "boundary_reconciled",
            "boundary_operator_failed",
        ),
        SupportCellCompletionCandidate::PcurveThenBoundary => {
            if support_cell_pcurve_author::author_assigned_support_cell_pcurves(brep, faces)
                .is_err()
            {
                return SupportCellCompletionApplication {
                    operator_succeeded: false,
                    operation_status: Cow::Borrowed("pcurve_operator_failed"),
                };
            }
            support_cell_operation_result(
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                    brep, faces,
                )
                .is_ok(),
                "pcurve_and_boundary_reconciled",
                "boundary_operator_failed",
            )
        }
        SupportCellCompletionCandidate::PcurveThenBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_full_boundary(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenTop3BoundariesIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_top_full_boundaries(
                    brep, faces, 3,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryPrefixIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_ranked_boundary_prefix(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryNeighborhoodIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_boundary_neighborhood(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenWorstEndpointClusterIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_endpoint_cluster_with_status(
                        brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_row_group_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_two_row_groups_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenThreeRowGroupsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_three_row_groups_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_row_group_tangent_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenExactBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_exact_boundary_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenPartialBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_partial_boundary_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenExactAndPartialBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_exact_then_partial_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenExactPartialTangentIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_exact_partial_tangent_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenChainCurveNetworkIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_chain_curve_network_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenChainNetworkCrossTangentRibbonIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_chain_network_cross_tangent_ribbon_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveCellsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_cells_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistCellsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_cells_status(
                    brep, faces,
                )
            })
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistJointLegSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_cells_status(
                        brep, faces,
                    )
                },
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_joint_leg_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryRowGroupIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(brep, faces, |brep, faces| {
                support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_boundary_row_group(
                    brep, faces,
                )
            })
        }
    }
}

fn support_cell_operation_result(
    operator_succeeded: bool,
    success_status: &'static str,
    failure_status: &'static str,
) -> SupportCellCompletionApplication {
    SupportCellCompletionApplication {
        operator_succeeded,
        operation_status: if operator_succeeded {
            Cow::Borrowed(success_status)
        } else {
            Cow::Borrowed(failure_status)
        },
    }
}

fn support_cell_rank(
    candidate: SupportCellCompletionCandidate,
    quality: shared_trim_quality::Kp20SupportCellSharedTrimQuality,
) -> support_cell_completion_rank::SupportCellCompletionRank {
    support_cell_completion_rank::SupportCellCompletionRank {
        boundary_safe: candidate.is_boundary_safe(),
        quality,
    }
}
