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
    /// The file was saved with a newer format version than this build supports.
    #[error("file version {file_version} not supported (max: {max_supported})")]
    UnsupportedVersion {
        /// Version number found in the file.
        file_version: u32,
        /// Highest version this build can read.
        max_supported: u32,
    },
    /// The file records an entity-naming scheme this build cannot read
    /// (newer than [`popkern_core::naming::NAMING_SCHEME_VERSION`], or one
    /// with no migration).
    #[error("naming scheme {file_scheme} not supported (this build mints scheme {supported})")]
    UnsupportedNamingScheme {
        /// Scheme version recorded in the file.
        file_scheme: u32,
        /// Scheme version this build mints and reads.
        supported: u32,
    },
    /// A feature or entity label is already used by another feature or
    /// entity of the same kind.
    #[error("duplicate label: {0}")]
    DuplicateLabel(String),
    /// A label must not be empty.
    #[error("label must not be empty")]
    EmptyLabel,
    /// No feature or entity carries the given label.
    #[error("label not found: {0}")]
    LabelNotFound(String),
    /// A cycle was detected while topologically sorting parameter expressions.
    #[error("expression cycle detected: {0:?}")]
    ExpressionCycle(Vec<crate::parameter::ParameterId>),
    /// A parameter was used with a unit whose [`DimensionKind`] does not match.
    ///
    /// [`DimensionKind`]: crate::unit::DimensionKind
    #[error("dimension mismatch: expected {expected:?}, got {got:?}")]
    DimensionMismatch {
        /// The parameter whose unit was violated.
        param: crate::parameter::ParameterId,
        /// The kind that was expected.
        expected: crate::unit::DimensionKind,
        /// The kind that was supplied.
        got: crate::unit::DimensionKind,
    },
    /// An expression references a parameter that could not be resolved.
    #[error("expression references unknown parameter {0:?}")]
    ExpressionUnresolved(crate::parameter::ParameterId),
    /// Attempted to remove a parameter that is still referenced by other
    /// parameters.
    #[error("parameter is referenced by: {0:?}")]
    ParameterInUse(Vec<crate::parameter::ParameterId>),
    /// Attempted to modify a read-only built-in parameter.
    #[error("cannot modify read-only parameter")]
    ReadOnlyParameter,
    /// A parameter name was already in use.
    #[error("duplicate parameter name: {0}")]
    DuplicateParameterName(String),
    /// A `suppress_if` expression must be dimensionless but had a different kind.
    #[error("suppress_if expression must be dimensionless, got {0:?}")]
    SuppressIfNotDimensionless(crate::unit::DimensionKind),
    /// The requested configuration does not exist.
    #[error("configuration not found")]
    ConfigNotFound,
    /// Adding the requested parent would create an inheritance cycle.
    #[error("configuration inheritance cycle detected")]
    InheritanceCycle,
    /// A configuration with this name already exists.
    #[error("duplicate configuration name: {0}")]
    DuplicateConfigName(String),
    /// The seed features list for a pattern is empty.
    #[error("seed features must not be empty")]
    EmptySeedFeatures,
    /// Seed features belong to more than one body.
    #[error("seed features must all be on the same body")]
    SeedBodyMismatch,
    /// Seed features do not occupy consecutive positions in the timeline.
    #[error("seed features must be contiguous in the timeline")]
    SeedNotContiguous,
    /// The first seed feature is not a Sketch.
    #[error("seed chain must start with a Sketch feature")]
    SeedMissingSketch,
    /// No seed feature has an [`OperationMode::Boss`] or [`OperationMode::Cut`] mode.
    #[error("seed chain must contain at least one feature with Boss or Cut mode")]
    SeedMissingMode,
    /// A seed feature is already used by another pattern.
    #[error("seed feature is already owned by another pattern")]
    SeedAlreadyOwned,
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
    /// A tree-level error (feature not found, invalid position, etc.).
    #[error("history error: {0}")]
    HistoryError(#[from] HistoryError),
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
