//! Error types, warnings, and structured remediation hints for
//! fillet and chamfer operations.
//!
//! ## Design
//!
//! The taxonomy has four layers, mirroring `popkern-bool`:
//!
//! 1. [`BlendError`] — fatal errors that abort the whole operation.
//!    Returned by every public entry point.
//! 2. [`BlendFailure`] — per-edge failure reasons attached to the
//!    [`BlendError::PartialFailure`] variant. Multiple edges in a
//!    single batch can fail for different reasons.
//! 3. [`BlendWarning`] — non-fatal observations emitted alongside a
//!    successful result (e.g. a clamped radius, an overlap hint).
//! 4. [`Remediation`] — structured guidance the engine returns for
//!    any [`BlendError`] or [`BlendFailure`] variant, naming an
//!    [`RemediationAction`], an optional tolerance parameter to
//!    adjust, and a human-readable rationale.
//!
//! Every `BlendError` and every `BlendFailure` variant maps to a
//! non-trivial [`Remediation`] through
//! [`BlendError::remediation_hint`] / [`BlendFailure::remediation_hint`],
//! so both functions are total — no variant returns "I don't know".
//! This is the same §13.4-style invariant that `popkern-bool` relies
//! on for commercial-grade diagnostics.

use popkern_core::geometry::Point3;
use popkern_core::topology::{EdgeHandle, VertexHandle};
use thiserror::Error;

/// Fatal error — the blend operation could not proceed.
#[derive(Debug, Error)]
pub enum BlendError {
    /// A parameter supplied by the caller is invalid (NaN, out-of-range,
    /// unknown entity handle, etc.). Reported *before* any geometry
    /// work happens, so the input BRep is guaranteed unchanged.
    #[error("invalid parameter `{parameter}` = {value}: {reason}")]
    InvalidParameter {
        /// Which parameter failed validation (e.g. `"radius"`,
        /// `"chamfer.distance"`, `"edge"`).
        parameter: String,
        /// Rendered value the caller passed in.
        value: String,
        /// Why the value is invalid.
        reason: String,
    },

    /// No valid edges could be analysed — every edge in the spec list
    /// failed the pre-operation analysis phase.
    #[error("no valid edges: {0}")]
    NoValidEdges(String),

    /// Internal topology error bubbled up from `popkern-core`.
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

    /// The operation was cancelled via a
    /// [`crate::cancellation::CancellationToken`].
    ///
    /// The input BRep is guaranteed unchanged when this is
    /// returned — cancellation happens between pipeline phases and
    /// never leaves the model in a partially-modified state.
    #[error("operation cancelled by caller")]
    Cancelled,
}

/// Per-edge failure reason.
///
/// Every variant carries enough structure for a UI or logging layer
/// to render an actionable message without having to parse strings.
/// See [`Self::remediation_hint`] for the engine's built-in
/// remediation mapping.
#[derive(Debug, Clone)]
pub enum BlendFailure {
    /// Fillet radius exceeds available space — the rolling ball would
    /// self-intersect, escape the face, or collide with geometry on
    /// the opposite side of the blend.
    RadiusTooLarge {
        /// Requested radius.
        requested: f64,
        /// Maximum radius the engine estimates would fit, if one was
        /// determined. `None` if the safe radius could not be
        /// established (e.g. surface offset failed for other reasons).
        max_radius: Option<f64>,
        /// 3D location where the failure was detected, if available.
        location: Option<Point3>,
    },
    /// Surface offset computation failed.
    OffsetFailed {
        /// Description of the failure.
        reason: String,
        /// 3D location where the failure was detected, if available.
        location: Option<Point3>,
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
        /// 3D location where the failure was detected, if available.
        location: Option<Point3>,
    },
    /// Edge has no adjacent faces or topology is invalid.
    InvalidTopology {
        /// Description of the topology issue.
        reason: String,
    },
    /// The rolling-ball contact solver failed to converge on this
    /// edge — typically a near-tangent configuration or a singularity
    /// on the support surfaces.
    ContactSolverFailed {
        /// Description of the failure.
        reason: String,
        /// Spine parameter where convergence broke down, if known.
        spine_parameter: Option<f64>,
    },
    /// The edge is already part of a blend (already filleted or
    /// chamfered) and the engine is not configured for blend-on-blend.
    AlreadyBlended,
    /// The edge is a degenerate seam edge (same face on both sides),
    /// which cannot be filleted or chamfered without breaking the
    /// seam first.
    SeamEdge,
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
        /// The vertex where the fallback occurred, if identified.
        vertex: Option<VertexHandle>,
    },
    /// Two fillet trim curves overlap on a shared face.
    FilletOverlap {
        /// Index of the first fillet edge in the specs list.
        fillet_a: usize,
        /// Index of the second fillet edge in the specs list.
        fillet_b: usize,
    },
    /// The blend output surface may have self-intersection risk at
    /// the requested radius — the caller should consider reducing
    /// the radius or running a post-operation validation.
    SelfIntersectionRisk {
        /// The radius that triggered the warning.
        radius: f64,
        /// The estimated safe maximum radius.
        safe_radius: f64,
    },
}

// The Remediation types and the per-variant remediation-hint
// impls live in `crate::remediation` — broken out to keep each
// source file under the 500-line structural limit. Re-exported
// here so `crate::error::Remediation` etc. continue to work.
pub use crate::remediation::{Remediation, RemediationAction, RemediationContext};

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
