//! Why a variational cap could not be produced.

use super::RimResidual;

/// A structured decline of [`super::fit_variational_cap`]. Every
/// variant names the offending side, stage or sample and the measured
/// quantity, so a caller can report it and a GUI can point at it.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum VariationalCapError {
    /// The opening has too few or too many sides.
    #[error("variational cap: {got} sides; the fit supports 3 to 12")]
    UnsupportedSideCount {
        /// Sides supplied.
        got: usize,
    },
    /// The degree is outside the supported range.
    #[error("variational cap: degree {degree} is unsupported (5 to 9)")]
    UnsupportedDegree {
        /// Degree requested.
        degree: usize,
    },
    /// Side `side`'s end and the next side's start do not meet.
    #[error(
        "variational cap: side {side} ends {gap:.3e} from the next side's start (tolerance {tolerance:.3e})"
    )]
    BoundaryNotClosed {
        /// Side whose end missed.
        side: usize,
        /// Measured gap.
        gap: f64,
        /// Tolerance it had to clear.
        tolerance: f64,
    },
    /// A side has (nearly) no length.
    #[error("variational cap: side {side} is degenerate (length {length:.3e})")]
    DegenerateSide {
        /// Offending side.
        side: usize,
        /// Measured chord length.
        length: f64,
    },
    /// The neighbour surface could not be evaluated at a rim sample.
    #[error(
        "variational cap: side {side}'s neighbour has no usable jet at parameter {at:.4}: {reason}"
    )]
    NeighbourJet {
        /// Offending side.
        side: usize,
        /// Side-curve parameter of the sample.
        at: f64,
        /// What failed.
        reason: String,
    },
    /// The hard constraints (exact spokes plus fixed rim rows) have no
    /// solution: the sides' 2-jets at their midpoints or the corners
    /// disagree.
    #[error("variational cap: hard constraints inconsistent at {stage} (residual {residual:.3e})")]
    Inconsistent {
        /// Where the inconsistency was measured.
        stage: &'static str,
        /// Normalised residual.
        residual: f64,
    },
    /// A linear system was numerically singular.
    #[error("variational cap: singular system at {stage}")]
    Singular {
        /// Which solve.
        stage: &'static str,
    },
    /// The targets were not met within the refinement budget.
    #[error(
        "variational cap: targets not met after {levels} refinement level(s) and {unknowns} unknowns; worst side {} at parameter {:.4}: g0 {:.3e}, g1 {:.3e}, g2 {:.3e}",
        worst.side, worst.parameter, worst.g0, worst.g1, worst.g2
    )]
    BudgetExhausted {
        /// The worst rim sample when the budget ran out.
        worst: RimResidual,
        /// Refinement levels used.
        levels: usize,
        /// Unknowns in the last solve.
        unknowns: usize,
    },
    /// A fitted sector folds or degenerates.
    #[error("variational cap: sector {sector} folds (min regularity {min_regularity:.3e})")]
    Folded {
        /// Offending sector.
        sector: usize,
        /// Smallest `|S_u × S_v|` relative to the sector mean.
        min_regularity: f64,
    },
}
