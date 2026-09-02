//! The scripting surface's error taxonomy (spec §10): one class per
//! failure kind, every class carrying the values a dialog needs; the
//! `Display` text is redundant with the fields, never the only carrier.

use popkern_history::{FeatureFailure, HistoryError, ReplayError};
use thiserror::Error;

/// What a selection was asked to produce.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionKind {
    /// Edges.
    Edges,
    /// Faces.
    Faces,
}

impl std::fmt::Display for SelectionKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            SelectionKind::Edges => "edges",
            SelectionKind::Faces => "faces",
        })
    }
}

/// Errors raised by the recorder. Record-time classes project
/// `popkern_history::HistoryError`; kernel declines arrive through the
/// replay report of the feature that was just appended.
#[derive(Debug, Error)]
pub enum ScriptError {
    /// The appended feature failed to replay; it has been removed
    /// again and the document is as it was before the call.
    #[error("feature '{label}' ({operation}) failed: {error}")]
    FeatureFailed {
        /// The label the feature would have had.
        label: String,
        /// Operation identifier (`Feature::operation_label`).
        operation: String,
        /// The structured kernel payload (spec §10.2 class): a blend
        /// decline carries `failures[]` with per-edge reasons and
        /// remediation, a boolean its stage, and so on. `Display` is
        /// the kernel's message. Boxed to keep the error small.
        error: Box<FeatureFailure>,
        /// The replay report up to and including the failure (spec
        /// §10.3 `report`), before the transactional rollback.
        report: Box<popkern_history::replay_report::ReplayReport>,
    },
    /// A selection matched nothing (and `allow_empty` was not set).
    #[error(
        "selection '{selection}' matched no {kind} (body has {available}); \
         call .allow_empty() if a no-op is intended"
    )]
    EmptySelection {
        /// Canonical spelling of the query.
        selection: String,
        /// Edges or faces.
        kind: SelectionKind,
        /// How many candidates existed before filtering.
        available: usize,
    },
    /// A single-entity argument resolved to several (or no) entities.
    #[error("selection '{selection}' must resolve to exactly {expected} {kind}, got {got}")]
    SelectionCardinality {
        /// Canonical spelling of the query.
        selection: String,
        /// Edges or faces.
        kind: SelectionKind,
        /// Required count.
        expected: usize,
        /// Resolved count.
        got: usize,
    },
    /// An argument was malformed before any kernel work ran.
    #[error("argument '{argument}': {reason}")]
    Argument {
        /// Which argument.
        argument: String,
        /// Why it was rejected.
        reason: String,
    },
    /// The model has no body to operate on yet.
    #[error("no active body: record a sketch and a solid feature first")]
    NoActiveBody,
    /// The body no longer exists (consumed by a boolean, or deleted).
    #[error("body {0} is not alive at the end of the timeline")]
    BodyNotAlive(u64),
    /// A gate failed: the tessellated body has open mesh edges.
    #[error("body is not watertight: {open_edges} open mesh edges")]
    NotWatertight {
        /// Number of single-use mesh edges.
        open_edges: usize,
    },
    /// The serving-boundary validation rejected the body.
    #[error("body failed validation: {detail}")]
    Validation {
        /// The kernel's validation message.
        detail: String,
    },
    /// File export or import failed.
    #[error("I/O failed for '{path}': {detail}")]
    Io {
        /// The file involved.
        path: String,
        /// The error text.
        detail: String,
    },
    /// A document/label/parameter/timeline error from the tree.
    #[error(transparent)]
    History(Box<HistoryError>),
    /// A replay error not attributable to the feature just recorded.
    #[error(transparent)]
    Replay(Box<ReplayError>),
    /// The `popkern-part` facade declined (gates, exports).
    #[error(transparent)]
    Part(Box<popkern_part::PartError>),
}

impl ScriptError {
    /// The structured kernel payload when this is a feature failure.
    pub fn feature_failure(&self) -> Option<&FeatureFailure> {
        match self {
            ScriptError::FeatureFailed { error, .. } => Some(error),
            _ => None,
        }
    }
}

impl From<HistoryError> for ScriptError {
    fn from(e: HistoryError) -> Self {
        ScriptError::History(Box::new(e))
    }
}

impl From<ReplayError> for ScriptError {
    fn from(e: ReplayError) -> Self {
        ScriptError::Replay(Box::new(e))
    }
}

impl From<popkern_part::PartError> for ScriptError {
    fn from(e: popkern_part::PartError) -> Self {
        ScriptError::Part(Box::new(e))
    }
}
