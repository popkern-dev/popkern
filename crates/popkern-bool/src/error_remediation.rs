//! Remediation-hint mapping for `BoolError` variants.
//!
//! Split out of `error.rs` to keep the parent under the 500-line cap.

use super::error::{BoolError, Remediation, RemediationAction, RemediationContext};

impl BoolError {
    /// Best-fit remediation hint for this error variant.
    ///
    /// Returns a structured [`Remediation`] suggestion derived
    /// purely from the error kind — no access to tolerances or
    /// input geometry is needed. Callers can use this as a default
    /// when the failure site didn't populate its own entry-level
    /// `Remediation`, or to build user-facing error messages.
    ///
    /// §13.4 acceptance criterion: every `BoolError` variant
    /// returns a non-`None` hint, so `remediation_hint()` is a
    /// total function.
    ///
    /// For richer hints that reference a specific tolerance value
    /// or that avoid looping on already-attempted remediations,
    /// use [`BoolError::remediation_hint_with_context`].
    pub fn remediation_hint(&self) -> Remediation {
        self.remediation_hint_with_context(&RemediationContext::default())
    }

    /// Best-fit remediation hint with caller-supplied context.
    ///
    /// When `ctx.tolerance` is `Some`, tolerance-adjusting
    /// remediations include a concrete suggested `value` (e.g. the
    /// current `linear` tolerance multiplied by 10) in addition to
    /// the action and rationale. When `ctx.perturbation_attempted`
    /// is `true`, the hint mapping replaces
    /// `RetryWithPerturbation` with the next-best alternative to
    /// avoid looping — §13.4 "cross-reference with the existing
    /// perturbation fallback so suggestions don't loop".
    pub fn remediation_hint_with_context(&self, ctx: &RemediationContext) -> Remediation {
        match self {
            BoolError::NoIntersection => Remediation::new(
                RemediationAction::NoAction,
                "The two solids do not intersect. No boolean result exists \
                 for this operation. If you expected an intersection, verify \
                 that the bodies are positioned as intended — a coordinate \
                 error in one operand is the most common cause.",
            ),
            BoolError::IntersectionFailed(_) => {
                if ctx.perturbation_attempted {
                    // Perturbation already tried and failed — suggest
                    // tolerance widening instead, with a concrete value
                    // if the caller provided one.
                    tolerance_bump_hint(
                        ctx,
                        "linear",
                        10.0,
                        "Surface intersection solver failed even after \
                         geometric perturbation. Try increasing `linear` \
                         tolerance to widen the solver's convergence \
                         envelope, or heal the input to eliminate \
                         near-tangent face pairs.",
                    )
                } else {
                    Remediation::new(
                        RemediationAction::RetryWithPerturbation,
                        "Surface intersection solver failed. Try \
                         `boolean_with_perturbation()` to slightly offset \
                         the geometry and retry, or increase `linear` \
                         tolerance to widen the solver's convergence \
                         envelope.",
                    )
                }
            }
            BoolError::SplitFailed { .. } => tolerance_bump_hint(
                ctx,
                "linear",
                10.0,
                "Face splitting failed — typically caused by an \
                 intersection curve that lies too close to an existing \
                 edge or vertex. Try increasing `linear` tolerance, or \
                 heal the input to snap near-coincident features before \
                 retrying.",
            ),
            BoolError::ClassificationFailed { .. } => {
                if ctx.perturbation_attempted {
                    tolerance_bump_hint(
                        ctx,
                        "coincidence",
                        2.0,
                        "Point-in-solid classification failed even after \
                         perturbation. Try increasing `coincidence` \
                         tolerance so the multi-ray sampler treats \
                         near-boundary points deterministically, or heal \
                         the input to eliminate near-tangent face pairs.",
                    )
                } else {
                    tolerance_bump_hint(
                        ctx,
                        "coincidence",
                        2.0,
                        "Point-in-solid classification failed for a face \
                         fragment. Try increasing `coincidence` tolerance \
                         so the multi-ray sampler treats near-boundary \
                         points deterministically, or heal the input to \
                         eliminate near-tangent face pairs.",
                    )
                }
            }
            BoolError::MergeFailed { .. } => Remediation::new(
                RemediationAction::HealInput,
                "Assembling the result solid from selected face \
                 fragments failed. This usually indicates non-manifold \
                 topology in the intersection graph. Heal the input \
                 geometry (close gaps, remove duplicate vertices) before \
                 retrying.",
            ),
            BoolError::ValidationFailed(_) => Remediation::new(
                RemediationAction::HealInput,
                "Post-merge validation flagged the result solid as \
                 invalid. The boolean pipeline may have produced slivers, \
                 non-manifold edges, or open loops. Try healing the \
                 result, or retry with tighter tolerances so the merge \
                 phase produces cleaner topology.",
            ),
            BoolError::DegenerateInput(_) => Remediation::new(
                RemediationAction::HealInput,
                "One or both inputs are geometrically degenerate (zero \
                 volume, duplicate vertices, or zero-area faces). Run \
                 the healing pipeline on the inputs before retrying the \
                 boolean.",
            ),
            BoolError::CoincidentFaces { .. } => {
                if ctx.perturbation_attempted {
                    // Perturbation already tried; fall back to the
                    // coincidence-mode suggestion.
                    Remediation::new(
                        RemediationAction::UseCoincidenceMode,
                        "Two faces are coincident and perturbation did \
                         not break the ambiguity. Set \
                         `BoolTolerance::mode` to `AssumedFull` or \
                         `ToleranceEscalation` so the coincidence engine \
                         handles the overlap explicitly.",
                    )
                } else {
                    Remediation::new(
                        RemediationAction::UseCoincidenceMode,
                        "Two faces of the input solids are geometrically \
                         coincident, creating an ambiguous boundary. Set \
                         the `BoolTolerance::mode` to `AssumedFull` or \
                         `ToleranceEscalation` so the coincidence engine \
                         handles the overlap explicitly instead of the \
                         transversal solver.",
                    )
                }
            }
            BoolError::CoreError(_) => Remediation::new(
                RemediationAction::HealInput,
                "A core topology or geometry operation failed. This is \
                 usually a sign of invalid input; running the healing \
                 pipeline and re-running validation before retrying the \
                 boolean is the safest first step.",
            ),
            BoolError::Cancelled => Remediation::new(
                RemediationAction::NoAction,
                "The operation was cancelled by the caller via a \
                 cancellation token. No action needed — the input \
                 model is unchanged.",
            ),
            BoolError::OpenShellNotSupported { .. } => Remediation::new(
                RemediationAction::Precheck,
                "At least one input is an open shell / sheet body. \
                 Route through `boolean_with_sheets()` instead of the \
                 standard `boolean()` entry point — the sheet-aware \
                 dispatcher handles the operation matrix for open \
                 shells.",
            ),
            BoolError::InterferenceFailed(_) => {
                if ctx.perturbation_attempted {
                    tolerance_bump_hint(
                        ctx,
                        "linear",
                        10.0,
                        "Interference detection failed after perturbation. \
                         Try increasing `linear` tolerance so the \
                         BVH-backed query treats near-contact pairs as \
                         overlapping.",
                    )
                } else {
                    Remediation::new(
                        RemediationAction::RetryWithPerturbation,
                        "Interference detection could not determine whether \
                         the bodies overlap. This is usually a near-tangent \
                         or near-coincident configuration. Try \
                         `boolean_with_perturbation()` or increase `linear` \
                         tolerance.",
                    )
                }
            }
            BoolError::SystemError(_) => Remediation::new(
                RemediationAction::Precheck,
                "A system-level resource failure occurred (out of \
                 memory, per-stage iteration limit exceeded). Simplify \
                 the input (fewer faces, decimated NURBS surfaces, \
                 smaller models) and retry.",
            ),
        }
    }
}

/// Build an `IncreaseTolerance` remediation, populating the concrete
/// suggested value from `ctx.tolerance` when available. The factor
/// controls how much the current value is scaled up — 10× is the
/// default for linear failures, 2× for coincidence.
fn tolerance_bump_hint(
    ctx: &RemediationContext,
    parameter: &str,
    factor: f64,
    rationale: &str,
) -> Remediation {
    let suggested = ctx.tolerance.map(|t| match parameter {
        "linear" => t.linear() * factor,
        "angular" => t.angular() * factor,
        "parametric" => t.parametric_base() * factor,
        "coincidence" => t.linear() * 10.0 * factor,
        _ => t.linear() * factor,
    });
    Remediation {
        action: RemediationAction::IncreaseTolerance,
        parameter: Some(parameter.to_string()),
        value: suggested,
        rationale: rationale.to_string(),
    }
}
