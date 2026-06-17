//! Chain-network authoring modes for support-cell row construction.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ChainNetworkMode {
    Curves,
    CurvesAndCrossTangentRibbons,
    CurvesAndConstructiveCells,
    /// Constructive cells emitted as dual-twist patches: each
    /// boundary's derived corner band coefficient lands in its
    /// own Chiyokura-Kimura candidate slot, so twist conflicts
    /// cease to exist instead of being compromised.
    CurvesAndConstructiveDualTwistCells,
    CurvesAndConstructiveDualTwistCellsWithJointLegSolve,
    CurvesAndConstructiveDualTwistCellsWithBetaEnvelope,
    CurvesAndConstructiveDualTwistCellsWithPairNormalSolve,
    CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve,
    CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve,
    CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve,
    CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve,
    CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve,
}

impl ChainNetworkMode {
    pub(super) fn constructive(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
        )
    }

    pub(super) fn dual_twist(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCells
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
        )
    }

    pub(super) fn joint_leg_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve
        )
    }

    pub(super) fn beta_envelope(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
        )
    }

    pub(super) fn pair_normal_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve
        )
    }

    pub(super) fn active_normal_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
        )
    }

    pub(super) fn curvature_solve(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
        )
    }

    pub(super) fn prefers_spanning_owner_row(self) -> bool {
        matches!(
            self,
            ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve
                | ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve
        )
    }

    pub(super) fn status_label(self, authored: bool) -> &'static str {
        match (self, authored) {
            (ChainNetworkMode::Curves, true) => "chain_curve_network_authored",
            (ChainNetworkMode::Curves, false) => "chain_curve_network_declined",
            (ChainNetworkMode::CurvesAndCrossTangentRibbons, true) => "chain_ribbon_network_authored",
            (ChainNetworkMode::CurvesAndCrossTangentRibbons, false) => "chain_ribbon_network_declined",
            (ChainNetworkMode::CurvesAndConstructiveCells, true) => "chain_constructive_cells_authored",
            (ChainNetworkMode::CurvesAndConstructiveCells, false) => "chain_constructive_cells_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCells, true) => "chain_constructive_dual_twist_cells_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCells, false) => "chain_constructive_dual_twist_cells_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve, true) => "chain_constructive_dual_twist_joint_leg_solve_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithJointLegSolve, false) => "chain_constructive_dual_twist_joint_leg_solve_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope, true) => "chain_constructive_dual_twist_beta_envelope_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelope, false) => "chain_constructive_dual_twist_beta_envelope_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve, true) => "chain_constructive_dual_twist_pair_normal_solve_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithPairNormalSolve, false) => "chain_constructive_dual_twist_pair_normal_solve_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve, true) => "chain_constructive_dual_twist_active_normal_solve_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithActiveNormalSolve, false) => "chain_constructive_dual_twist_active_normal_solve_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve, true) => "chain_constructive_dual_twist_beta_envelope_active_normal_solve_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve, false) => "chain_constructive_dual_twist_beta_envelope_active_normal_solve_declined",
            (ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve, true) => "chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_solve_authored",
            (ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalSolve, false) => "chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_solve_declined",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve, true) => "chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_authored",
            (ChainNetworkMode::CurvesAndConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve, false) => "chain_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_declined",
            (ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve, true) => "chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_authored",
            (ChainNetworkMode::CurvesAndOwnerRowConstructiveDualTwistCellsWithBetaEnvelopeActiveNormalCurvatureSolve, false) => "chain_owner_row_constructive_dual_twist_beta_envelope_active_normal_curvature_solve_declined",
        }
    }
}
