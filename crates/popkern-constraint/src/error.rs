use crate::constraint::ConstraintHandle;
use crate::elements::ElementHandle;
use serde::{Deserialize, Serialize};

/// Errors that can occur during constraint operations.
#[derive(Debug, Clone, thiserror::Error)]
pub enum ConstraintError {
    /// A handle references an element or constraint that does not exist.
    #[error("invalid handle")]
    InvalidHandle,

    /// The constraint does not have a modifiable dimensional value.
    #[error("constraint is not a dimensional type (Distance, Angle, or PointLineDistance)")]
    NotDimensional,

    /// `drag_to` was called without an active drag session.
    #[error("no drag is active — call begin_drag first")]
    NoDragActive,

    /// The Newton–Raphson solver failed to converge.
    #[error("solve failed: {0}")]
    SolveFailed(String),
}

/// The outcome of a solve attempt.
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// Whether and how the sketch is constrained.
    pub status: SolveStatus,
    /// Number of Newton–Raphson iterations performed.
    pub iterations: u32,
    /// L2 norm of the residual vector at termination.
    pub residual_norm: f64,
}

/// Constraint satisfaction state of the sketch after solving.
#[derive(Debug, Clone)]
pub enum SolveStatus {
    /// Every degree of freedom is constrained exactly.
    FullyConstrained,
    /// Some degrees of freedom remain free.
    UnderConstrained(Vec<FreedomInfo>),
    /// More constraints than degrees of freedom; conflicts detected.
    OverConstrained(Vec<ConstraintHandle>),
    /// Solver encountered a numerical failure.
    Failed(String),
}

/// Describes a free (unconstrained) degree of freedom on a sketch element.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FreedomInfo {
    /// The element that has remaining freedom.
    pub element: ElementHandle,
    /// Unit directions in which the element can move freely.
    /// Each tuple is a (dx, dy) unit vector.
    pub directions: Vec<(f64, f64)>,
}
