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
        )
    }
}
