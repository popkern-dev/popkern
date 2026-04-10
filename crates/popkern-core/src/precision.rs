//! Precision tracking and diagnostics types.
//!
//! This module provides the type infrastructure for tracking, reporting, and
//! remediating precision degradation across BRep operations. Types are
//! organised into three layers:
//!
//! - **Entity addressing:** [`EntityHandle`], [`Containment`]
//! - **Diagnostics:** [`PrecisionDiagnostic`], [`PrecisionHealth`],
//!   [`PrecisionSummary`], [`OperationPrecisionReport`]
//! - **Solver feedback:** [`SolverConfig`], [`SolverDiagnostics`],
//!   [`ConvergenceOutcome`], [`FailureReason`]

use crate::topology::{EdgeHandle, FaceHandle, VertexHandle};
use std::collections::VecDeque;
use std::fmt;

// ---------------------------------------------------------------------------
// Entity addressing
// ---------------------------------------------------------------------------

/// Type-erased handle to a topological entity.
///
/// Wraps [`VertexHandle`], [`EdgeHandle`], or [`FaceHandle`] so that
/// diagnostics can reference any entity uniformly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EntityHandle {
    /// A vertex entity.
    Vertex(VertexHandle),
    /// An edge entity.
    Edge(EdgeHandle),
    /// A face entity.
    Face(FaceHandle),
}

impl From<VertexHandle> for EntityHandle {
    fn from(h: VertexHandle) -> Self {
        Self::Vertex(h)
    }
}

impl From<EdgeHandle> for EntityHandle {
    fn from(h: EdgeHandle) -> Self {
        Self::Edge(h)
    }
}

impl From<FaceHandle> for EntityHandle {
    fn from(h: FaceHandle) -> Self {
        Self::Face(h)
    }
}

impl fmt::Display for EntityHandle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Vertex(h) => write!(f, "Vertex({h:?})"),
            Self::Edge(h) => write!(f, "Edge({h:?})"),
            Self::Face(h) => write!(f, "Face({h:?})"),
        }
    }
}

/// Point-in-region classification.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Containment {
    /// Strictly inside the region.
    Inside,
    /// Strictly outside the region.
    Outside,
    /// On the boundary of the region.
    OnBoundary,
}

// ---------------------------------------------------------------------------
// Diagnostic severity and kinds
// ---------------------------------------------------------------------------

/// Severity level for precision diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiagnosticSeverity {
    /// Precision constraint violated -- operation may be invalid.
    Error,
    /// Precision is degraded but still within bounds.
    Warning,
    /// Informational (e.g. tracking budget consumption).
    Info,
}

/// Classification of a precision diagnostic.
#[derive(Clone, Debug)]
pub enum PrecisionDiagnosticKind {
    /// Measured imprecision exceeds the maximum allowed bound.
    ExceedsMaxBound,
    /// Value is below the numerical floor (near machine epsilon).
    BelowNumericalFloor,
    /// Point containment classification is inconsistent with geometry.
    ContainmentViolation,
    /// Two entities have overlapping imprecision zones.
    OverlappingImprecision,
    /// Precision has degraded significantly from a previous operation.
    SignificantDegradation,
    /// A numerical computation failed for the given reason.
    NumericalFailure {
        /// Human-readable description of the failure.
        reason: String,
    },
    /// Edge geometry is inconsistent with its vertex positions.
    EdgeInconsistency,
    /// Precision budget is being consumed; estimated operations remaining.
    PrecisionBudgetWarning {
        /// Estimated number of operations before the budget is exhausted.
        estimated_remaining: usize,
    },
}

// ---------------------------------------------------------------------------
// Remediation
// ---------------------------------------------------------------------------

/// Suggested action to recover from a precision issue.
#[derive(Clone, Debug)]
pub enum Remediation {
    /// Re-snap vertex position to its owning curve/surface.
    RecoverVertexPrecision(VertexHandle),
    /// Re-fit edge curve through its endpoint vertices.
    RecoverEdgePrecision(EdgeHandle),
    /// Re-check edge consistency (twins, orientation).
    RecoverEdgeConsistency(EdgeHandle),
    /// Free-form remediation description.
    Manual(String),
    /// No remediation available or needed.
    None,
}

impl fmt::Display for Remediation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RecoverVertexPrecision(h) => {
                write!(f, "recover vertex precision for {h:?}")
            }
            Self::RecoverEdgePrecision(h) => {
                write!(f, "recover edge precision for {h:?}")
            }
            Self::RecoverEdgeConsistency(h) => {
                write!(f, "recover edge consistency for {h:?}")
            }
            Self::Manual(msg) => write!(f, "{msg}"),
            Self::None => write!(f, "no remediation"),
        }
    }
}

// ---------------------------------------------------------------------------
// PrecisionDiagnostic
// ---------------------------------------------------------------------------

/// A single precision diagnostic record.
///
/// Captures what went wrong, where, and what to do about it.
#[derive(Clone, Debug)]
pub struct PrecisionDiagnostic {
    /// How severe the issue is.
    pub severity: DiagnosticSeverity,
    /// What kind of precision issue was detected.
    pub kind: PrecisionDiagnosticKind,
    /// Which entities are affected.
    pub entities: Vec<EntityHandle>,
    /// The measured imprecision value.
    pub measured: f64,
    /// The threshold that was compared against.
    pub threshold: f64,
    /// Suggested remediation action.
    pub remediation: Remediation,
}

// ---------------------------------------------------------------------------
// Health assessment
// ---------------------------------------------------------------------------

/// Overall precision health of a BRep or operation result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrecisionHealth {
    /// All entities are well within tolerance bounds.
    Good,
    /// Some entities are approaching tolerance bounds.
    Degraded,
    /// At least one entity exceeds tolerance bounds.
    Critical,
}

impl PrecisionHealth {
    /// Classify health from a budget ratio (remaining / total).
    ///
    /// - ratio >= 0.5 => Good
    /// - ratio >= 0.1 => Degraded
    /// - ratio < 0.1  => Critical
    pub fn from_budget_ratio(ratio: f64) -> Self {
        if ratio >= 0.5 {
            Self::Good
        } else if ratio >= 0.1 {
            Self::Degraded
        } else {
            Self::Critical
        }
    }
}

/// Aggregate precision summary for a BRep.
#[derive(Clone, Debug)]
pub struct PrecisionSummary {
    /// Normalised health score in `[0.0, 1.0]`.
    pub health_score: f64,
    /// Discrete health classification.
    pub health: PrecisionHealth,
    /// Entity with the worst (largest) imprecision, if any.
    pub worst_entity: Option<(EntityHandle, f64)>,
    /// Number of entities in Degraded state.
    pub degraded_count: usize,
    /// Number of entities in Critical state.
    pub critical_count: usize,
    /// Estimated operations remaining before budget exhaustion.
    pub estimated_remaining: Option<usize>,
}

// ---------------------------------------------------------------------------
// Operation report
// ---------------------------------------------------------------------------

/// Precision report produced after a geometric operation.
///
/// Contains per-entity before/after deltas and any diagnostics raised.
#[derive(Clone, Debug)]
pub struct OperationPrecisionReport {
    /// Per-entity (entity, before_imprecision, after_imprecision).
    pub entity_deltas: Vec<(EntityHandle, f64, f64)>,
    /// Maximum imprecision observed across all entities.
    pub max_imprecision: f64,
    /// Entities that exceeded their tolerance bound.
    pub exceeded_bound: Vec<EntityHandle>,
    /// Diagnostics raised during the operation.
    pub diagnostics: Vec<PrecisionDiagnostic>,
}

// ---------------------------------------------------------------------------
// Near-coincidence configuration
// ---------------------------------------------------------------------------

/// Configuration for near-coincidence detection.
///
/// Controls how aggressively the system flags entities that are close
/// to being coincident but not quite within tolerance.
#[derive(Clone, Copy, Debug)]
pub struct NearCoincidenceConfig {
    /// Multiplier applied to tolerance for the near-coincidence range.
    pub range_multiplier: f64,
    /// Maximum multiplier cap.
    pub max_multiplier: f64,
}

impl Default for NearCoincidenceConfig {
    fn default() -> Self {
        Self {
            range_multiplier: 3.0,
            max_multiplier: 10.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Degradation tracking
// ---------------------------------------------------------------------------

/// Tracks precision degradation across a sequence of operations.
///
/// Maintains a sliding window of (operation_name, max_imprecision) pairs
/// to estimate the rate of degradation and remaining budget.
pub struct DegradationTracker {
    history: VecDeque<(&'static str, f64)>,
    window_size: usize,
}

impl DegradationTracker {
    /// Create a new tracker with the given sliding window size.
    pub fn new(window_size: usize) -> Self {
        Self {
            history: VecDeque::with_capacity(window_size),
            window_size,
        }
    }

    /// Record the precision impact of an operation.
    pub fn record(&mut self, op_name: &'static str, report: &OperationPrecisionReport) {
        if self.history.len() >= self.window_size {
            self.history.pop_front();
        }
        self.history.push_back((op_name, report.max_imprecision));
    }

    /// Compute the average degradation rate (imprecision per operation).
    ///
    /// Returns 0.0 if fewer than 2 operations have been recorded.
    pub fn degradation_rate(&self) -> f64 {
        if self.history.len() < 2 {
            return 0.0;
        }
        let first = self.history.front().map(|(_, v)| *v).unwrap_or(0.0);
        let last = self.history.back().map(|(_, v)| *v).unwrap_or(0.0);
        let steps = (self.history.len() - 1) as f64;
        (last - first) / steps
    }

    /// Estimate how many more operations can be performed before
    /// `current_worst` reaches `max_bound`.
    ///
    /// Returns `None` if the degradation rate is non-positive (not degrading).
    pub fn estimated_operations_remaining(
        &self,
        current_worst: f64,
        max_bound: f64,
    ) -> Option<usize> {
        let rate = self.degradation_rate();
        if rate <= 0.0 {
            return None;
        }
        let remaining = (max_bound - current_worst) / rate;
        if remaining <= 0.0 {
            Some(0)
        } else {
            Some(remaining.ceil() as usize)
        }
    }
}

// ---------------------------------------------------------------------------
// Solver configuration and diagnostics
// ---------------------------------------------------------------------------

/// Configuration for iterative geometric solvers.
#[derive(Clone, Debug)]
pub struct SolverConfig {
    /// Maximum iterations for the primary solver.
    pub max_iterations: usize,
    /// Maximum iterations for fallback methods.
    pub max_fallback_iterations: usize,
    /// Whether to allow precision escalation (tightening tolerance mid-solve).
    pub allow_escalation: bool,
    /// Whether to allow fallback to alternative solver methods.
    pub allow_fallback: bool,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            max_iterations: 50,
            max_fallback_iterations: 100,
            allow_escalation: true,
            allow_fallback: true,
        }
    }
}

/// Diagnostics produced by an iterative solver run.
#[derive(Clone, Debug)]
pub struct SolverDiagnostics {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual at convergence (or termination).
    pub final_residual: f64,
    /// Precision actually achieved.
    pub achieved_precision: f64,
    /// How the solver terminated.
    pub convergence: ConvergenceOutcome,
}

/// How an iterative solver terminated.
#[derive(Clone, Debug)]
pub enum ConvergenceOutcome {
    /// Converged within the requested tolerance.
    Converged,
    /// Converged, but at a relaxed tolerance.
    ConvergedRelaxed {
        /// Originally requested precision.
        requested: f64,
        /// Actually achieved precision.
        achieved: f64,
    },
    /// Primary solver failed; a fallback method succeeded.
    FallbackMethod {
        /// Name of the fallback method used.
        method: &'static str,
    },
    /// Tolerance was tightened mid-solve (precision escalation).
    PrecisionEscalated,
    /// Solver failed to converge.
    Failed {
        /// Why the solver failed.
        reason: FailureReason,
    },
}

/// Why an iterative solver failed to converge.
#[derive(Clone, Debug)]
pub enum FailureReason {
    /// Hit the maximum iteration count.
    MaxIterations,
    /// Residual is increasing instead of decreasing.
    Divergence,
    /// Residual is oscillating without converging.
    Oscillation,
    /// System is ill-conditioned (near-singular Jacobian).
    IllConditioned {
        /// Estimated condition number of the system.
        condition_estimate: f64,
    },
    /// System is exactly singular.
    Singular,
}

#[cfg(test)]
#[path = "precision_tests.rs"]
mod precision_tests;
