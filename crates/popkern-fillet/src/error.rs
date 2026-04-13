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

// ---------------------------------------------------------------------------
// Remediation
// ---------------------------------------------------------------------------

/// A structured suggestion for how to recover from a [`BlendError`]
/// or [`BlendFailure`]. The engine returns one for every error
/// variant, so callers never see "unknown failure" — they always get
/// a concrete next step.
#[derive(Debug, Clone)]
pub struct Remediation {
    /// What class of action to take.
    pub action: RemediationAction,
    /// Which tolerance parameter to adjust, if applicable.
    pub parameter: Option<String>,
    /// Suggested new numeric value for the parameter, if applicable.
    pub value: Option<f64>,
    /// Why this action is expected to help.
    pub rationale: String,
}

impl Remediation {
    /// Build a remediation with just an action and rationale (no
    /// tolerance-parameter metadata).
    pub fn new(action: RemediationAction, rationale: impl Into<String>) -> Self {
        Self {
            action,
            parameter: None,
            value: None,
            rationale: rationale.into(),
        }
    }

    /// Build a tolerance-adjustment remediation that proposes
    /// `new_value` for `parameter`.
    pub fn adjust_tolerance(
        action: RemediationAction,
        parameter: impl Into<String>,
        new_value: f64,
        rationale: impl Into<String>,
    ) -> Self {
        Self {
            action,
            parameter: Some(parameter.into()),
            value: Some(new_value),
            rationale: rationale.into(),
        }
    }

    /// Build a radius-reduction remediation that proposes a specific
    /// `new_radius`.
    pub fn reduce_radius(new_radius: f64, rationale: impl Into<String>) -> Self {
        Self {
            action: RemediationAction::ReduceRadius,
            parameter: Some("radius".to_string()),
            value: Some(new_radius),
            rationale: rationale.into(),
        }
    }
}

/// Context a caller can supply to
/// [`BlendError::remediation_hint_with_context`] so the engine can
/// tailor its suggestions. Default-constructed context yields the
/// same hints as [`BlendError::remediation_hint`].
#[derive(Debug, Clone, Copy, Default)]
pub struct RemediationContext<'a> {
    /// The tolerance the failing operation ran with. When set,
    /// tolerance-adjusting remediations include a specific new value
    /// derived from the current tolerance (e.g. `linear × 10`).
    pub tolerance: Option<&'a popkern_core::tolerance::Tolerance>,
    /// When set, the engine knows the caller has already tried a
    /// workaround attempt (e.g. retried with a perturbed input or a
    /// reduced radius). Remediations that would re-suggest the same
    /// workaround are suppressed so callers don't loop.
    pub workaround_attempted: bool,
}

/// Possible remediation actions the engine suggests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RemediationAction {
    /// The caller passed an invalid parameter — fix the input and retry.
    FixInput,
    /// Increase a tolerance parameter (typically `linear` or
    /// `angular`) to widen the solver's convergence envelope.
    IncreaseTolerance,
    /// Reduce the fillet radius — the requested value exceeds the
    /// local geometry's feasible range.
    ReduceRadius,
    /// Heal the input geometry (close gaps, remove duplicates, fix
    /// degenerate faces) before retrying.
    HealInput,
    /// Split the operation into smaller pieces, or run a pre-check
    /// to identify problematic edges before retrying.
    Precheck,
    /// Break the offending topology (e.g. split a seam edge so the
    /// fillet has two distinct support faces) before retrying.
    BreakTopology,
    /// Retry with a different vertex-blend mode or profile variant.
    UseAlternativeMode,
    /// The operation is not applicable to this input (e.g.
    /// already-blended edge, seam edge, disjoint geometry). No
    /// engine-level action will help — the caller must revise the spec.
    NoAction,
}

// ---------------------------------------------------------------------------
// BlendError → Remediation
// ---------------------------------------------------------------------------

impl BlendError {
    /// Best-fit remediation hint for this error variant.
    ///
    /// Returns a structured [`Remediation`] derived purely from the
    /// error kind — no access to tolerance values is needed. Use
    /// [`Self::remediation_hint_with_context`] for a version that
    /// populates concrete tolerance values from the caller's
    /// [`popkern_core::tolerance::Tolerance`].
    ///
    /// Every variant returns a non-empty rationale — the function
    /// is total.
    pub fn remediation_hint(&self) -> Remediation {
        self.remediation_hint_with_context(&RemediationContext::default())
    }

    /// Best-fit remediation hint with caller-supplied context.
    ///
    /// When `ctx.tolerance` is `Some`, tolerance-adjusting
    /// remediations include a concrete suggested `value` (e.g. the
    /// current `linear` tolerance multiplied by 10). When
    /// `ctx.workaround_attempted` is `true`, the hint mapping
    /// avoids re-suggesting the same workaround so callers don't
    /// loop.
    pub fn remediation_hint_with_context(&self, ctx: &RemediationContext) -> Remediation {
        match self {
            BlendError::InvalidParameter {
                parameter, reason, ..
            } => Remediation::new(
                RemediationAction::FixInput,
                format!(
                    "Parameter `{parameter}` is invalid: {reason}. \
                     Fix the caller-supplied value and retry."
                ),
            ),
            BlendError::NoValidEdges(detail) => Remediation::new(
                RemediationAction::Precheck,
                format!(
                    "No edges in the spec list passed the pre-operation \
                     analysis phase ({detail}). Run `precheck_fillet()` \
                     on each edge to identify specific problems before \
                     building the spec list."
                ),
            ),
            BlendError::Topology(_) => Remediation::new(
                RemediationAction::HealInput,
                "A core topology operation failed. This usually means \
                 the input BRep has a pre-existing validity problem. \
                 Run the healing pipeline on the input and re-validate \
                 before retrying the fillet.",
            ),
            BlendError::PartialFailure { failed, total, .. } => Remediation::new(
                RemediationAction::Precheck,
                format!(
                    "{failed} of {total} edges failed individually. \
                     Inspect the `failures` list — each entry has its \
                     own remediation hint via \
                     `BlendFailure::remediation_hint()`. Consider \
                     running `precheck_fillet()` before the next attempt."
                ),
            ),
        }
        .fill_tolerance_value(ctx)
    }
}

impl Remediation {
    /// When the remediation action is `IncreaseTolerance` and the
    /// context provides a tolerance, fill in a concrete suggested
    /// value. A no-op otherwise.
    fn fill_tolerance_value(mut self, ctx: &RemediationContext) -> Self {
        if self.action == RemediationAction::IncreaseTolerance
            && self.value.is_none()
            && let Some(tol) = ctx.tolerance
            && let Some(param) = self.parameter.as_deref()
        {
            self.value = Some(match param {
                "linear" => tol.linear() * 10.0,
                "angular" => tol.angular() * 10.0,
                _ => tol.linear() * 10.0,
            });
        }
        self
    }
}

// ---------------------------------------------------------------------------
// BlendFailure → Remediation
// ---------------------------------------------------------------------------

impl BlendFailure {
    /// Best-fit remediation hint for this per-edge failure.
    ///
    /// Returns a structured [`Remediation`]. The function is total —
    /// every variant yields a non-empty rationale.
    pub fn remediation_hint(&self) -> Remediation {
        self.remediation_hint_with_context(&RemediationContext::default())
    }

    /// Best-fit remediation with caller-supplied context.
    pub fn remediation_hint_with_context(&self, ctx: &RemediationContext) -> Remediation {
        match self {
            BlendFailure::RadiusTooLarge {
                requested,
                max_radius,
                ..
            } => {
                if let Some(max) = max_radius {
                    // Suggest a radius that is 95% of the safe maximum
                    // — leaving a small margin for numerical headroom.
                    let suggested = max * 0.95;
                    Remediation::reduce_radius(
                        suggested,
                        format!(
                            "The requested radius {requested} exceeds the \
                             rolling-ball feasibility bound on this edge \
                             (max ≈ {max}). Retry with a radius no \
                             greater than {suggested:.6} — this leaves a \
                             5% safety margin to avoid self-intersection."
                        ),
                    )
                } else {
                    Remediation::new(
                        RemediationAction::ReduceRadius,
                        format!(
                            "The requested radius {requested} exceeds the \
                             feasibility bound on this edge, but the \
                             engine could not determine a specific safe \
                             maximum. Reduce the radius and retry, or \
                             call `feasible_radius_bounds()` to compute \
                             the precise range."
                        ),
                    )
                }
            }
            BlendFailure::OffsetFailed { reason, .. } => {
                if ctx.workaround_attempted {
                    Remediation::new(
                        RemediationAction::HealInput,
                        format!(
                            "Surface offset failed again after a prior \
                             workaround ({reason}). The support surface \
                             is likely degenerate — run the healing \
                             pipeline or substitute a cleaner surface."
                        ),
                    )
                } else {
                    Remediation::new(
                        RemediationAction::IncreaseTolerance,
                        format!(
                            "Surface offset failed ({reason}). Offset \
                             routines are sensitive to near-tangent or \
                             near-coincident surfaces. Try increasing \
                             `linear` tolerance, or heal the support \
                             surface if it has sliver faces."
                        ),
                    )
                    .with_parameter("linear")
                }
            }
            BlendFailure::VertexBlendFailed { reason, .. } => Remediation::new(
                RemediationAction::UseAlternativeMode,
                format!(
                    "The vertex blend solver could not converge \
                     ({reason}). Retry with `VertexMode::Mitred` — the \
                     mitred variant is more robust at the cost of \
                     losing G1 continuity at the vertex. If the \
                     mitred mode also fails, this vertex configuration \
                     is likely beyond the current engine's capability."
                ),
            ),
            BlendFailure::IntersectionFailed { reason, .. } => Remediation::new(
                RemediationAction::IncreaseTolerance,
                format!(
                    "Surface-surface intersection failed during trim \
                     computation ({reason}). The trim curve solver is \
                     sensitive to tangent and near-tangent pairs. \
                     Try increasing `linear` tolerance, reducing the \
                     radius, or healing the input to eliminate \
                     sliver faces adjacent to the blend edge."
                ),
            )
            .with_parameter("linear"),
            BlendFailure::InvalidTopology { reason } => Remediation::new(
                RemediationAction::HealInput,
                format!(
                    "The edge's surrounding topology is not usable by \
                     the fillet engine ({reason}). Heal the input BRep \
                     (close gaps, fix non-manifold edges, split seams) \
                     before retrying."
                ),
            ),
            BlendFailure::ContactSolverFailed {
                reason,
                spine_parameter,
            } => {
                let where_ = spine_parameter
                    .map(|t| format!(" at spine t = {t:.3}"))
                    .unwrap_or_default();
                Remediation::new(
                    RemediationAction::ReduceRadius,
                    format!(
                        "The rolling-ball contact solver failed to \
                         converge{where_} ({reason}). This is usually a \
                         sign that the requested radius is too close to \
                         the edge's local feasibility limit. Reduce the \
                         radius, or split the chain at this spine \
                         parameter and blend the two sub-chains \
                         separately."
                    ),
                )
            }
            BlendFailure::AlreadyBlended => Remediation::new(
                RemediationAction::NoAction,
                "The edge is already part of an existing blend — the \
                 current engine does not support blend-on-blend \
                 construction. Remove the existing blend with \
                 `unfillet()` before re-applying, or target a \
                 different edge.",
            ),
            BlendFailure::SeamEdge => Remediation::new(
                RemediationAction::BreakTopology,
                "The edge is a seam (same face on both sides, e.g. the \
                 U=0 seam on a cylinder). Seams cannot be filleted \
                 directly — split the face along the seam into two \
                 distinct faces before applying the fillet.",
            ),
        }
    }
}

impl Remediation {
    /// Replace the `parameter` field. Private chaining helper.
    fn with_parameter(mut self, parameter: &str) -> Self {
        self.parameter = Some(parameter.to_string());
        self
    }
}

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
