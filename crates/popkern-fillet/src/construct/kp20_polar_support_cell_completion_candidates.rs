//! Candidate inventory for KP20 support-cell completion.
//!
//! Each variant names one measured repair strategy of the
//! completion ladder; the labels feed the candidate score facts.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SupportCellCompletionCandidate {
    PcurveOnly,
    BoundaryOnly,
    PcurveThenBoundary,
    PcurveThenBoundaryIfG0Safe,
    PcurveThenWorstBoundaryIfG0Safe,
    PcurveThenTop3BoundariesIfG0Safe,
    PcurveThenRankedBoundaryPrefixIfG0Safe,
    PcurveThenWorstBoundaryNeighborhoodIfG0Safe,
    PcurveThenWorstEndpointClusterIfG0Safe,
    PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe,
    PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe,
    PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe,
    PcurveThenRankedBoundaryEndpointClusterThenThreeRowGroupsIfG0Safe,
    PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe,
    PcurveThenRankedBoundaryEndpointRowGroupThenExactBoundaryIfG0Safe,
    PcurveThenRankedBoundaryEndpointRowGroupThenPartialBoundaryIfG0Safe,
    PcurveThenRankedBoundaryEndpointRowGroupThenExactAndPartialBoundaryIfG0Safe,
    PcurveThenRankedBoundaryEndpointRowGroupThenExactPartialTangentIfG0Safe,
    PcurveThenChainCurveNetworkIfG0Safe,
    PcurveThenChainNetworkCrossTangentRibbonIfG0Safe,
    PcurveThenChainConstructiveCellsIfG0Safe,
    PcurveThenChainConstructiveDualTwistCellsIfG0Safe,
    PcurveThenChainConstructiveDualTwistJointLegSolveIfG0Safe,
    PcurveThenChainConstructiveDualTwistBetaEnvelopeIfG0Safe,
    PcurveThenChainConstructiveDualTwistPairNormalSolveIfG0Safe,
    PcurveThenChainConstructiveDualTwistActiveNormalSolveIfG0Safe,
    PcurveThenChainConstructiveDualTwistWorstG1TangentIfG0Safe,
    PcurveThenChainConstructiveDualTwistWorstG1NormalSolveIfG0Safe,
    PcurveThenChainConstructiveDualTwistEndpointAdjacentWorstG1NormalSolveIfG0Safe,
    PcurveThenWorstBoundaryRowGroupIfG0Safe,
}

impl SupportCellCompletionCandidate {
    pub(super) const ALL: [Self; 30] = [
        Self::PcurveOnly,
        Self::BoundaryOnly,
        Self::PcurveThenBoundary,
        Self::PcurveThenBoundaryIfG0Safe,
        Self::PcurveThenWorstBoundaryIfG0Safe,
        Self::PcurveThenTop3BoundariesIfG0Safe,
        Self::PcurveThenRankedBoundaryPrefixIfG0Safe,
        Self::PcurveThenWorstBoundaryNeighborhoodIfG0Safe,
        Self::PcurveThenWorstEndpointClusterIfG0Safe,
        Self::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointClusterThenThreeRowGroupsIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactBoundaryIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointRowGroupThenPartialBoundaryIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactAndPartialBoundaryIfG0Safe,
        Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactPartialTangentIfG0Safe,
        Self::PcurveThenChainCurveNetworkIfG0Safe,
        Self::PcurveThenChainNetworkCrossTangentRibbonIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistCellsIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistJointLegSolveIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistBetaEnvelopeIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistPairNormalSolveIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistActiveNormalSolveIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistWorstG1TangentIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistWorstG1NormalSolveIfG0Safe,
        Self::PcurveThenChainConstructiveDualTwistEndpointAdjacentWorstG1NormalSolveIfG0Safe,
        Self::PcurveThenChainConstructiveCellsIfG0Safe,
        Self::PcurveThenWorstBoundaryRowGroupIfG0Safe,
    ];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::PcurveOnly => "pcurve_only",
            Self::BoundaryOnly => "boundary_only",
            Self::PcurveThenBoundary => "pcurve_then_boundary",
            Self::PcurveThenBoundaryIfG0Safe => "pcurve_then_boundary_if_g0_safe",
            Self::PcurveThenWorstBoundaryIfG0Safe => "pcurve_then_worst_boundary_if_g0_safe",
            Self::PcurveThenTop3BoundariesIfG0Safe => "pcurve_then_top3_boundaries_if_g0_safe",
            Self::PcurveThenRankedBoundaryPrefixIfG0Safe => {
                "pcurve_then_ranked_boundary_prefix_if_g0_safe"
            }
            Self::PcurveThenWorstBoundaryNeighborhoodIfG0Safe => {
                "pcurve_then_worst_boundary_neighborhood_if_g0_safe"
            }
            Self::PcurveThenWorstEndpointClusterIfG0Safe => {
                "pcurve_then_worst_endpoint_cluster_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe => {
                "pcurve_then_ranked_boundary_prefix_then_worst_endpoint_cluster_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_cluster_then_row_group_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_cluster_then_two_row_groups_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointClusterThenThreeRowGroupsIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_cluster_then_three_row_groups_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_row_group_then_tangent_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactBoundaryIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_boundary_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointRowGroupThenPartialBoundaryIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_row_group_then_partial_boundary_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactAndPartialBoundaryIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_and_partial_boundary_if_g0_safe"
            }
            Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactPartialTangentIfG0Safe => {
                "pcurve_then_ranked_boundary_endpoint_row_group_then_exact_partial_tangent_if_g0_safe"
            }
            Self::PcurveThenChainCurveNetworkIfG0Safe => {
                "pcurve_then_chain_curve_network_if_g0_safe"
            }
            Self::PcurveThenChainNetworkCrossTangentRibbonIfG0Safe => {
                "pcurve_then_chain_network_cross_tangent_ribbon_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveCellsIfG0Safe => {
                "pcurve_then_chain_constructive_cells_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistCellsIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_cells_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistJointLegSolveIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_joint_leg_solve_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistBetaEnvelopeIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_beta_envelope_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistPairNormalSolveIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_pair_normal_solve_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistActiveNormalSolveIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_active_normal_solve_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistWorstG1TangentIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_worst_g1_tangent_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistWorstG1NormalSolveIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_worst_g1_normal_solve_if_g0_safe"
            }
            Self::PcurveThenChainConstructiveDualTwistEndpointAdjacentWorstG1NormalSolveIfG0Safe => {
                "pcurve_then_chain_constructive_dual_twist_endpoint_adjacent_worst_g1_normal_solve_if_g0_safe"
            }
            Self::PcurveThenWorstBoundaryRowGroupIfG0Safe => {
                "pcurve_then_worst_boundary_row_group_if_g0_safe"
            }
        }
    }

    pub(super) fn is_boundary_safe(self) -> bool {
        matches!(
            self,
            Self::PcurveThenBoundaryIfG0Safe
                | Self::PcurveThenWorstBoundaryIfG0Safe
                | Self::PcurveThenTop3BoundariesIfG0Safe
                | Self::PcurveThenRankedBoundaryPrefixIfG0Safe
                | Self::PcurveThenWorstBoundaryNeighborhoodIfG0Safe
                | Self::PcurveThenWorstEndpointClusterIfG0Safe
                | Self::PcurveThenRankedBoundaryPrefixThenWorstEndpointClusterIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointClusterThenRowGroupIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointClusterThenTwoRowGroupsIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointClusterThenThreeRowGroupsIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointRowGroupThenTangentIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactBoundaryIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointRowGroupThenPartialBoundaryIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactAndPartialBoundaryIfG0Safe
                | Self::PcurveThenRankedBoundaryEndpointRowGroupThenExactPartialTangentIfG0Safe
                | Self::PcurveThenChainCurveNetworkIfG0Safe
                | Self::PcurveThenChainNetworkCrossTangentRibbonIfG0Safe
                | Self::PcurveThenChainConstructiveCellsIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistCellsIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistJointLegSolveIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistBetaEnvelopeIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistPairNormalSolveIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistActiveNormalSolveIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistWorstG1TangentIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistWorstG1NormalSolveIfG0Safe
                | Self::PcurveThenChainConstructiveDualTwistEndpointAdjacentWorstG1NormalSolveIfG0Safe
                | Self::PcurveThenWorstBoundaryRowGroupIfG0Safe
        )
    }
}
