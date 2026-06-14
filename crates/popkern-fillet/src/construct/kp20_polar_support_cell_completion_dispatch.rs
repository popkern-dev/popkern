//! Dispatch from measured support-cell completion candidates to operators.

use std::borrow::Cow;

use popkern_core::topology::{BRep, FaceHandle};

use super::support_cell_completion_candidates::SupportCellCompletionCandidate;
use super::{
    SupportCellCompletionApplication, support_cell_boundary_reconcile,
    support_cell_completion_apply, support_cell_pcurve_author,
};

pub(super) fn apply_support_cell_completion_candidate(
    brep: &mut BRep,
    faces: &[FaceHandle],
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
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_full_boundaries(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_full_boundary(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenTop3BoundariesIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_top_full_boundaries(
                        brep, faces, 3,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryPrefixIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_ranked_boundary_prefix(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryNeighborhoodIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_boundary_neighborhood(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenWorstEndpointClusterIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_endpoint_cluster_with_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_row_group_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_two_row_groups_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointClusterThenThreeRowGroupsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_three_row_groups_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_row_group_tangent_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenExactBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_exact_boundary_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenPartialBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_partial_boundary_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenExactAndPartialBoundaryIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_exact_then_partial_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenRankedBoundaryEndpointRowGroupThenExactPartialTangentIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::boundary_composed::reconcile_assigned_ranked_endpoint_exact_partial_tangent_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainCurveNetworkIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_curve_network_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainNetworkCrossTangentRibbonIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_network_cross_tangent_ribbon_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveCellsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_cells_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistCellsIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_cells_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistJointLegSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_joint_leg_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistBetaEnvelopeIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_beta_envelope_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistPairNormalSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_pair_normal_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistActiveNormalSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_active_normal_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistBetaEnvelopeActiveNormalSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_beta_envelope_active_normal_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistWorstG1TangentIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_then_worst_g1_tangent_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistWorstG1NormalSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_then_worst_g1_normal_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenChainConstructiveDualTwistEndpointAdjacentWorstG1NormalSolveIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_with_status_if_g1_safe_against_baseline(
                brep,
                faces,
                dual_twist_cells_status,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_then_endpoint_adjacent_worst_g1_normal_solve_status(
                        brep, faces,
                    )
                },
            )
        }
        SupportCellCompletionCandidate::PcurveThenWorstBoundaryRowGroupIfG0Safe => {
            support_cell_completion_apply::boundary_after_pcurve_if_quality_safe(
                brep,
                faces,
                |brep, faces| {
                    support_cell_boundary_reconcile::reconcile_assigned_support_cell_worst_boundary_row_group(
                        brep, faces,
                    )
                },
            )
        }
    }
}

fn dual_twist_cells_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, popkern_core::error::CoreError> {
    support_cell_boundary_reconcile::reconcile_assigned_chain_constructive_dual_twist_cells_status(
        brep, faces,
    )
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
