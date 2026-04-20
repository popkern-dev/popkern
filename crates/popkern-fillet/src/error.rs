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
        /// 3D location of the failing vertex, if available. Populated
        /// at construction so UI consumers do not have to resolve the
        /// `vertex` handle against the BRep.
        location: Option<Point3>,
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
        /// 3D location where the topology issue was detected, if
        /// available. Typically the midpoint of the offending edge.
        location: Option<Point3>,
    },
    /// The rolling-ball contact solver failed to converge on this
    /// edge — typically a near-tangent configuration or a singularity
    /// on the support surfaces.
    ContactSolverFailed {
        /// Description of the failure.
        reason: String,
        /// Spine parameter where convergence broke down, if known.
        spine_parameter: Option<f64>,
        /// 3D location along the spine where convergence broke down,
        /// if available. Typically the edge point at `spine_parameter`.
        location: Option<Point3>,
    },
    /// The edge is already part of a blend (already filleted or
    /// chamfered) and the engine is not configured for blend-on-blend.
    AlreadyBlended {
        /// 3D location of the already-blended edge, if available.
        /// Typically the midpoint of the edge.
        location: Option<Point3>,
    },
    /// The edge is a degenerate seam edge (same face on both sides),
    /// which cannot be filleted or chamfered without breaking the
    /// seam first.
    SeamEdge {
        /// 3D location of the seam edge, if available. Typically the
        /// midpoint of the edge.
        location: Option<Point3>,
    },
    /// The constructed fillet surface self-intersects (§13.2
    /// "Self-Intersection Detection"). Detected after
    /// `build_fillet_surface` returns but before the topology
    /// stitch-up phase by [`crate::surface_validate::validate_fillet_surface`].
    ///
    /// Typical causes:
    /// - a spindle torus (minor_radius > major_radius) from an
    ///   over-large fillet radius that still slipped past the
    ///   pre-offset curvature check;
    /// - a NURBS sweep whose control polygon has wrapped back on
    ///   itself (walker instability on near-tangent geometry);
    /// - a twist in the cross-section field (rotation > π between
    ///   adjacent UV samples).
    ///
    /// The remediation is typically to reduce the radius; see
    /// [`Self::remediation_hint`] for the full guidance.
    SelfIntersectingBlend {
        /// Structured description of the detected failure mode.
        reason: String,
        /// 3D location where the self-intersection was detected,
        /// if available (typically the edge midpoint).
        location: Option<Point3>,
    },
}

impl BlendFailure {
    /// Best-available 3D location of the failure, for UI consumers
    /// that want to 3D-pick every diagnostic.
    ///
    /// Returns `None` when the failure has no associated location
    /// (none of the construction sites had a usable point) or when
    /// the variant does not carry location data. Every current
    /// variant carries `Option<Point3>`, so the return is purely a
    /// reflection of whether the construction site populated it.
    pub fn location(&self) -> Option<Point3> {
        match self {
            BlendFailure::RadiusTooLarge { location, .. }
            | BlendFailure::OffsetFailed { location, .. }
            | BlendFailure::VertexBlendFailed { location, .. }
            | BlendFailure::IntersectionFailed { location, .. }
            | BlendFailure::InvalidTopology { location, .. }
            | BlendFailure::ContactSolverFailed { location, .. }
            | BlendFailure::SelfIntersectingBlend { location, .. }
            | BlendFailure::AlreadyBlended { location }
            | BlendFailure::SeamEdge { location } => *location,
        }
    }
}

// `BlendWarning` moved to `crate::blend_warnings` to keep this file
// under the 500-LoC structural limit and to give Phase I (polar K-P)
// headroom for its fallback-ladder variants. Re-exported below.
pub use crate::blend_warnings::BlendWarning;

// The Remediation types and the per-variant remediation-hint
// impls live in `crate::remediation` — broken out to keep each
// source file under the 500-line structural limit. Re-exported
// here so `crate::error::Remediation` etc. continue to work.
pub use crate::remediation::{Remediation, RemediationAction, RemediationContext};

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
