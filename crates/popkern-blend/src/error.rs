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
    /// Two requested fillets' removal strips overlap mid-face on a
    /// shared support face (their edges share no vertex, so this is
    /// not a corner blend): completing would require blend×blend
    /// intersection surgery the engine does not construct yet, and
    /// proceeding produced silent wrong-shape output (measured:
    /// r=0.6+0.6 across a width-1.0 face applied ONE band trimmed at
    /// the other's rail and dropped the second fillet with no
    /// diagnostic). The operation declines transactionally — the
    /// input BRep is unchanged.
    #[error(
        "fillet strips interfere mid-face: fillet {fillet_a} and fillet {fillet_b} overlap on a \
         shared support face; blend-blend intersection is not yet constructed — operation \
         declined, input unchanged"
    )]
    FilletInterference {
        /// Spec index of the first fillet.
        fillet_a: usize,
        /// Spec index of the second fillet.
        fillet_b: usize,
    },

    /// A face–face blend attach (TrimBoth) declined during trim
    /// surgery. The operation works on a private copy, so the input
    /// BRep is unchanged. `mechanism` names the surgery stage that
    /// declined (e.g. `"end-arc support"`, `"trim side"`,
    /// `"doomed-region flood"`), `detail` carries the measured
    /// evidence (residuals, handles, counts).
    #[error("face blend attach declined at {mechanism}: {detail} — input unchanged")]
    AttachDeclined {
        /// Surgery stage that declined.
        mechanism: String,
        /// Measured evidence for the decline.
        detail: String,
    },

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

impl BlendError {
    /// Build the transactional [`BlendError::PartialFailure`] with the
    /// count invariant enforced by construction: `failures` is
    /// deduplicated per edge (first report wins — a trim stage that
    /// defers both sides of one fillet edge must not report the edge
    /// twice; measured: "24 of 12 edges failed" on the oversized cube
    /// fillet, PK-380), and `failed` is the distinct-edge count,
    /// clamped to `total` so `failed <= total` holds for every
    /// consumer.
    pub fn partial_failure(
        total: usize,
        partial_result: popkern_core::topology::BRep,
        failures: Vec<(EdgeHandle, BlendFailure)>,
    ) -> BlendError {
        let mut seen: Vec<EdgeHandle> = Vec::with_capacity(failures.len());
        let failures: Vec<(EdgeHandle, BlendFailure)> = failures
            .into_iter()
            .filter(|(edge, _)| {
                if seen.contains(edge) {
                    false
                } else {
                    seen.push(*edge);
                    true
                }
            })
            .collect();
        BlendError::PartialFailure {
            failed: failures.len().min(total.max(1)),
            total,
            partial_result: Box::new(partial_result),
            failures,
        }
    }
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
    /// The served blend territory is geometrically inconsistent: an
    /// edge of a face the operation created or rewrote leaves its
    /// faces beyond its own band (check 9), or declares a band above
    /// the model's Precision ceiling. The construction produced
    /// geometry that is structurally sound but not a valid solid, so
    /// the operation declines instead of serving it (PK-530: strip
    /// ends on oblique faces once served 8.4e-2 off).
    InconsistentGeometry {
        /// The measured violation: edge, deviation or band, count.
        reason: String,
        /// Midpoint of the worst edge, if available.
        location: Option<Point3>,
    },
    /// A strict corner-cap preference
    /// ([`crate::CornerCapPreference::strict`]) named a construction
    /// that could not serve this corner, so the fillet declines
    /// instead of falling back (PK-603).
    CornerCapUnavailable {
        /// The corner vertex.
        vertex: VertexHandle,
        /// The construction the preference required.
        construction: crate::CornerConstruction,
        /// Why that construction declined (carve, postcondition, or
        /// certified continuity below the request).
        reason: String,
        /// 3D location of the corner, if available.
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
            | BlendFailure::InconsistentGeometry { location, .. }
            | BlendFailure::CornerCapUnavailable { location, .. }
            | BlendFailure::AlreadyBlended { location }
            | BlendFailure::SeamEdge { location } => *location,
        }
    }
}

// `BlendWarning` moved to `crate::blend_warnings` to keep this file
// under the 500-LoC structural limit and to give the polar-cap campaign (polar K-P)
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
