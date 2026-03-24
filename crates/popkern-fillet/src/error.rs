//! Error types for fillet and chamfer operations.

use popkern_core::topology::{EdgeHandle, VertexHandle};
use thiserror::Error;

/// Fatal error — the blend operation could not proceed.
#[derive(Debug, Error)]
pub enum BlendError {
    /// No valid edges to process.
    #[error("no valid edges: {0}")]
    NoValidEdges(String),

    /// Internal topology error.
    #[error("topology error: {0}")]
    Topology(#[from] popkern_core::error::CoreError),

    /// Some edges failed — partial result available.
    #[error("{failed} of {total} edges failed")]
    PartialFailure {
        /// Number of failed edges.
        failed: usize,
        /// Total edges attempted.
        total: usize,
        /// Partial result with successful blends applied.
        partial_result: Box<popkern_core::topology::BRep>,
        /// Per-edge failure details.
        failures: Vec<(EdgeHandle, BlendFailure)>,
    },
}

/// Per-edge failure reason.
#[derive(Debug, Clone)]
pub enum BlendFailure {
    /// Fillet radius exceeds available space.
    RadiusTooLarge {
        /// Maximum radius that would fit.
        max_radius: f64,
    },
    /// Surface offset computation failed.
    OffsetFailed {
        /// Description of the failure.
        reason: String,
    },
    /// Vertex blend could not be computed.
    VertexBlendFailed {
        /// The vertex that failed.
        vertex: VertexHandle,
        /// Description of the failure.
        reason: String,
    },
    /// Surface intersection failed during trim computation.
    IntersectionFailed {
        /// Description of the failure.
        reason: String,
    },
    /// Edge has no adjacent faces or topology is invalid.
    InvalidTopology {
        /// Description of the topology issue.
        reason: String,
    },
}

/// Non-fatal per-edge warning.
#[derive(Debug, Clone)]
pub enum BlendWarning {
    /// Edge was auto-propagated from a chain.
    Propagated,
    /// Radius was clamped to fit available space.
    RadiusClamped {
        /// Originally requested radius.
        requested: f64,
        /// Actual radius used after clamping.
        actual: f64,
    },
    /// Vertex mode was downgraded because the requested mode is not yet
    /// supported for this configuration.
    VertexModeFallback {
        /// The mode that was requested.
        requested: crate::types::VertexMode,
        /// The mode that was actually used.
        actual: crate::types::VertexMode,
    },
}
