//! Error types for the history/feature-tree subsystem.

use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors arising from tree-level operations (insert, suppress, reorder, etc.).
#[derive(Debug, Error)]
pub enum HistoryError {
    /// The requested feature handle does not exist in the tree.
    #[error("feature not found")]
    FeatureNotFound,
    /// The given timeline index is out of range.
    #[error("invalid timeline position: {0}")]
    InvalidPosition(usize),
    /// A dependency cycle was detected when inserting or reordering a feature.
    #[error("dependency cycle detected")]
    DependencyCycle,
    /// A dependency references a feature that comes after the insertion point.
    #[error("dependency {0} is after insertion point")]
    DependencyAfterInsert(String),
    /// The referenced body does not exist.
    #[error("body not found")]
    BodyNotFound,
    /// Serialization or deserialization of the tree failed.
    #[error("serialization error: {0}")]
    SerializeError(String),
}

/// Errors arising during replay (re-evaluating features against a BRep).
#[derive(Debug, Error)]
pub enum ReplayError {
    /// A specific feature failed during evaluation.
    #[error("feature '{label}' failed: {error}")]
    FeatureFailed {
        /// Human-readable label of the failed feature.
        label: String,
        /// Error message from the underlying operation.
        error: String,
    },
    /// A stable name reference could not be resolved in the current BRep.
    #[error("reference resolution failed: {0:?}")]
    ResolutionFailed(EntityName),
    /// A core geometry or topology error occurred.
    #[error("core error: {0}")]
    CoreError(#[from] popkern_core::error::CoreError),
}

/// The cached outcome of evaluating a single [`crate::feature::Feature`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FeatureResult {
    /// The feature evaluated successfully.  The [`HistoryMap`] records the
    /// stable names assigned to the topology produced by this feature.
    ///
    /// [`HistoryMap`]: popkern_core::history::HistoryMap
    Success {
        /// Stable-name tracking output from the operation.
        history: popkern_core::history::HistoryMap,
    },
    /// The feature was attempted but failed; the error message is stored so
    /// it can be displayed without re-evaluating.
    Failed {
        /// Human-readable description of the failure.
        error: String,
    },
    /// The feature is marked suppressed and was not evaluated.
    Suppressed,
}
