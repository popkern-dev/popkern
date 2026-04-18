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
    /// A fillet spec was soft-skipped from the batch by the
    /// §13.1(5) collapsed-edge pre-pass because the target edge's
    /// two endpoint vertices were separated by less than the
    /// geometric tolerance. The rest of the batch proceeds
    /// normally; the caller should check `warnings` in addition
    /// to `failures` to understand the full set of spec outcomes.
    CollapsedEdgeSkipped {
        /// The edge that was dropped.
        edge: EdgeHandle,
        /// Measured 3D length between the edge's endpoint vertices.
        length: f64,
    },
    /// A sliver face was detected in the fillet output topology but
    /// could not be absorbed into a neighbour (no suitable target
    /// or the topology surgery would be invalid). The output BRep
    /// contains the sliver unchanged; downstream consumers may
    /// choose to ignore it, run a separate healing pass, or fail.
    ///
    /// §13.1(2) "remove, stitch neighbors" — commercial kernels
    /// ship this warning for the residual cases where absorption
    /// is impossible. See `construct::sliver_absorb` for the
    /// classifier details.
    SliverRetained {
        /// The residual sliver face.
        face: popkern_core::topology::FaceHandle,
        /// Short description of the sliver category that was
        /// classified but not absorbed.
        kind: &'static str,
    },
    /// §13.4 soft-supports mode was active for this fillet
    /// operation: the pipeline's internal tolerance was scaled up
    /// from the caller's value to accommodate low-precision input
    /// geometry. Downstream consumers can distinguish soft-mode
    /// results from strict ones by checking for this warning.
    ///
    /// Emitted exactly once per fillet call when
    /// [`crate::types::FilletOptions::soft_supports`] is `true`.
    SoftSupportTolerancesApplied {
        /// The multiplicative factor applied to `tolerance.linear()`
        /// when producing the relaxed tolerance. Defaults to 10,
        /// set via
        /// [`crate::types::FilletOptions::soft_support_factor`].
        factor: f64,
    },
    /// §8.3 cap was produced at a blend endpoint. Carries the
    /// requested kind and a flag describing whether topology was
    /// modified to insert the cap face (on an open-shell or
    /// free-boundary endpoint) or whether the cap is represented
    /// implicitly by an existing adjacent solid face (the common
    /// case on closed bodies, where the adjacent third face at a
    /// corner already serves as the natural "flat" cap).
    ///
    /// Emitted once per `(edge, endpoint)` pair whenever
    /// [`crate::types::FilletOptions::cap_kind`] is non-`Open`.
    /// Callers can inspect the per-edge `EdgeGeometry::start_cap` /
    /// `end_cap` fields for the raw cap geometry (surface +
    /// boundary polyline) regardless of the warning.
    CapApplied {
        /// The filleted edge the cap terminates.
        edge: EdgeHandle,
        /// `true` if the cap is at the spine's start endpoint
        /// (`t = 0`); `false` if at the end (`t = 1`).
        at_start: bool,
        /// The cap kind the caller requested.
        kind: crate::types::CapKind,
        /// `true` if cap topology was inserted into the output
        /// BRep; `false` if the endpoint was already closed by an
        /// adjacent face and no insertion was needed.
        topology_inserted: bool,
    },
    /// §7 continuity fallback: the caller requested a higher
    /// continuity level than the engine was able to deliver for
    /// this edge. The output is still a valid fillet — it just
    /// sits at a lower level than requested.
    ///
    /// **When this fires.** The primary trigger is a G2 request
    /// where either (a) the spine walker couldn't produce contact
    /// solutions for the surface pair (e.g. degenerate normals,
    /// analytic plane-plane path that bypasses the walker), or
    /// (b) `g2_sweep_from_contact_solutions` errored on at least
    /// one spine sample (`G2Error::ZeroChord`,
    /// `DegenerateNormal`, or `TangentDegenerate`).
    ///
    /// **UI guidance.** Display this alongside
    /// [`crate::quality::BlendQuality::achieved_continuity`] so the
    /// caller sees both "you asked for G2" and "you got G1 because
    /// the cross-section at sample 3 had a zero chord." Consumer
    /// code that requires strict G2 should treat this as a failure
    /// and reject the output.
    ContinuityDowngraded {
        /// The filleted edge whose blend surface was downgraded.
        edge: EdgeHandle,
        /// What the caller asked for (from
        /// [`crate::types::FilletOptions::continuity`]).
        requested: crate::types::ContinuityLevel,
        /// What the engine actually delivered on this edge's blend
        /// face. Always lower than `requested`.
        achieved: crate::types::ContinuityLevel,
        /// Human-readable reason for the downgrade. Comes from the
        /// underlying walker / G2-sweep error so UI consumers can
        /// surface the specific geometric condition that forced
        /// the fallback.
        reason: String,
    },
    /// §5.3 corner-patch continuity fallback at a high-valence
    /// vertex. Fires when either (a) the construction reported an
    /// achieved continuity lower than `requested`, or (b) the
    /// corner has five or more incident fillet strips
    /// (`n_sides >= 5`), where even the curvature-matched polar
    /// NURBS patch degrades in interior quality because a single
    /// degree-2 radial cross-section can't capture the full
    /// curvature variation around N boundaries.
    ///
    /// **UI guidance.** Display alongside the
    /// [`crate::CornerPatchInfo`] for the same vertex so the caller
    /// sees `(method = Gregory, n_sides = 6, achieved = G1)`
    /// paired with "corner continuity downgraded because N ≥ 5."
    /// Consumer code that requires strict Class-A at N≥5 should
    /// treat this warning as a failure and reject the corner
    /// (Plate-minimisation / true rational Gregory for N≥5 is
    /// planned for §5.3 C6+).
    ///
    /// Keyed on a representative incident edge on the outer
    /// `Vec<(EdgeHandle, BlendWarning)>` to keep the warning
    /// plumbing uniform with other vertex-scoped warnings
    /// (`VertexModeFallback`); `vertex` inside the body carries
    /// the authoritative handle.
    CornerContinuityDowngraded {
        /// The vertex whose corner patch was downgraded.
        vertex: VertexHandle,
        /// Valence of the corner — number of incident filleted
        /// edges meeting at this vertex.
        n_sides: usize,
        /// The continuity level the caller asked for (from
        /// [`crate::types::FilletOptions::continuity`]).
        requested: crate::types::ContinuityLevel,
        /// The continuity level the corner construction actually
        /// delivered.
        achieved: crate::types::ContinuityLevel,
        /// Human-readable reason (e.g. "N >= 5, polar NURBS
        /// interior may degrade").
        reason: String,
    },
    /// §5.3 item #6 — a face-edge bridge in the vertex-blend
    /// boundary construction fell back to a straight-line cubic
    /// Bezier because the generalised surface-hugging bridge
    /// couldn't be built on the face's supporting surface (e.g.,
    /// degenerate UV projection at a singular NURBS point).
    ///
    /// The straight-line bridge still interpolates the arc
    /// endpoints correctly and the downstream Gregory fit
    /// re-samples the bridge — so the blend surface is valid
    /// but may sit slightly *off* the supporting face near the
    /// bridge, rather than hugging it. Surfaces downstream of
    /// the blend (mesh tessellation, boolean consumers) may
    /// show a visible seam where the bridge leaves the face.
    ///
    /// Callers can surface this warning in a UI as "curved-
    /// support blend fell back to planar approximation" — in
    /// the common case (planar supports), this variant is never
    /// emitted.
    SurfaceBridgeFallback {
        /// The vertex whose blend boundary triggered the fallback.
        vertex: VertexHandle,
        /// Short description of which surface type the fallback
        /// occurred on (e.g., `"Cylinder"`, `"NURBS"`) — for
        /// logging / UI display.
        surface_kind: &'static str,
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
