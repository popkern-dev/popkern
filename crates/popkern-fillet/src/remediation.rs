//! Structured remediation suggestions for [`crate::BlendError`] and
//! [`crate::BlendFailure`].
//!
//! Every error variant the engine produces maps to a non-trivial
//! [`Remediation`] via the total functions
//! [`crate::BlendError::remediation_hint`] and
//! [`crate::BlendFailure::remediation_hint`]. This module holds
//! the `Remediation`, `RemediationContext`, and `RemediationAction`
//! types plus the per-variant impl mapping.
//!
//! Broken out from `error.rs` to keep each source file under the
//! 500-line structural cap.

use crate::error::{BlendError, BlendFailure};

/// A structured suggestion for how to recover from a [`BlendError`]
/// or [`BlendFailure`]. The engine returns one for every error
/// variant, so callers never see "unknown failure" — they always
/// get a concrete next step.
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

    /// When the remediation action is `IncreaseTolerance` and the
    /// context provides a tolerance, fill in a concrete suggested
    /// value. A no-op otherwise.
    pub(crate) fn fill_tolerance_value(mut self, ctx: &RemediationContext) -> Self {
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

    /// Replace the `parameter` field. Private chaining helper.
    pub(crate) fn with_parameter(mut self, parameter: &str) -> Self {
        self.parameter = Some(parameter.to_string());
        self
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
            BlendError::Cancelled => Remediation::new(
                RemediationAction::NoAction,
                "The operation was cancelled by the caller via a \
                 cancellation token. No action needed — the input \
                 model is unchanged.",
            ),
        }
        .fill_tolerance_value(ctx)
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
