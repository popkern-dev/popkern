//! Error types for Phase 6c-A2 apex construction.
//!
//! Per ADR §4.2 and §5.4: every variant carries enough context to
//! be **user-actionable** when surfaced through a CAD UI.  Display
//! impls produce end-user-readable strings (no developer jargon).

use popkern_core::error::ValidationError;
use popkern_core::topology::VertexHandle;

/// Errors from [`plan_apex`](super::plan_apex).
///
/// Returned when the planner cannot produce a valid [`ApexPlan`]
/// from the inputs.  Each variant identifies a specific failure
/// class and includes context for user-actionable messaging.
#[derive(thiserror::Error, Debug)]
pub enum ApexPlanError {
    /// The apex vertex has fewer than 3 incident fillets.
    /// Sphere-fit apex resolution requires at least 3 face planes
    /// meeting at the apex to pin down the cap centre.
    #[error("apex {apex:?} has only {n} fillets; sphere-fit apex resolution requires at least 3")]
    InsufficientFillets {
        /// The apex vertex.
        apex: VertexHandle,
        /// Number of incident fillets actually present.
        n: usize,
    },

    /// The cap-centre linear system is singular: the three (or more)
    /// face planes through the apex are nearly coplanar, so no unique
    /// cap centre exists.  The condition number quantifies how singular.
    ///
    /// User-actionable: "Faces nearly coplanar; consider relocating
    /// fillet edges or accepting a lower-rigor fallback."
    #[error(
        "cap-centre system is singular at apex {apex:?} (condition number {condition_number:.3e}); \
         faces nearly coplanar"
    )]
    SingularCapCentre {
        /// The apex vertex.
        apex: VertexHandle,
        /// Reciprocal of the smallest singular value, indicating
        /// proximity to coplanarity.
        condition_number: f64,
    },

    /// A fillet's nominal radius exceeds the maximum compatible with
    /// this corner's extent.  The cap sphere would protrude past the
    /// adjacent face boundaries.
    ///
    /// User-actionable: "Reduce fillet radius below `max_for_corner`."
    #[error(
        "fillet radius {radius:.4} m exceeds maximum {max_for_corner:.4} m for corner at {apex:?}; \
         try a smaller radius"
    )]
    RadiusExceedsCornerExtent {
        /// The apex vertex.
        apex: VertexHandle,
        /// The fillet radius the user requested.
        radius: f64,
        /// Maximum compatible radius for this corner.
        max_for_corner: f64,
    },

    /// The cap sphere (computed from the cap-centre solution and the
    /// fillet's radius) does not intersect the fillet's cylinder.
    /// Geometrically, the rolling sphere can't reach the fillet.
    ///
    /// User-actionable: "Fillet `fillet_idx` is too small for this
    /// corner geometry.  Increase its radius."
    #[error(
        "cap sphere does not intersect fillet {fillet_idx}'s cylinder; sphere-fit termination undefined"
    )]
    SphereCylinderNoIntersection {
        /// The fillet index (0-based) that fails to terminate.
        fillet_idx: usize,
    },

    /// Generic invalid input: NaN/Inf/zero radii, malformed
    /// `FilletSpec`, or an apex vertex that doesn't exist in the BRep.
    /// `field` identifies which input is invalid; `value` shows what
    /// was supplied (for debug; UIs may not surface this directly).
    #[error("invalid input for {field}: {value}")]
    InvalidInput {
        /// Which input field failed validation.
        field: &'static str,
        /// String form of the offending value (for diagnostic).
        value: String,
    },

    /// The fillet uses a surface kind A2 v1 doesn't support
    /// (currently: only Cylinder fillets are supported).
    /// User-actionable: "Use a cylinder fillet, or wait for the
    /// corresponding A2 extension."
    #[error(
        "fillet {fillet_idx} uses surface kind {kind} which A2 v1 does not yet support; \
         use a cylinder fillet"
    )]
    UnsupportedFilletKind {
        /// The fillet index (0-based).
        fillet_idx: usize,
        /// The surface kind name (e.g. "Torus", "Cone").
        kind: &'static str,
    },

    /// Operation cancelled via the supplied [`CancelToken`](crate::cancellation::CancelToken).
    #[error("apex planning cancelled")]
    Cancelled,
}

/// Errors from [`construct_apex`](super::construct_apex).
///
/// Returned when applying an [`ApexPlan`] to the BRep fails.
/// On any of these errors, the BRep is rolled back to its
/// pre-construct state.
#[derive(thiserror::Error, Debug)]
pub enum ApexBuildError {
    /// An Euler operator failed at a specific step.  `op` names the
    /// operator (e.g. "mev", "mef"); `cause` describes why.
    #[error("Euler op {op} failed at step {step}: {cause}")]
    EulerOpFailed {
        /// Which Euler operator failed.
        op: &'static str,
        /// Index in the build step sequence.
        step: usize,
        /// Why it failed.
        cause: String,
    },

    /// Post-build validation reported an error.  The BRep is rolled
    /// back to its pre-construct state.  `error` is the underlying
    /// validation failure.
    #[error("post-build validation failed: {error}")]
    ValidationFailed {
        /// The validation error.
        error: ValidationError,
    },

    /// The plan is inconsistent with the current BRep state (e.g.
    /// vertex handles in the plan don't exist, or geometry has
    /// drifted since planning).
    #[error("plan inconsistent with BRep state: {reason}")]
    PlanMismatch {
        /// What went wrong.
        reason: String,
    },

    /// Operation cancelled via the supplied [`CancelToken`](crate::cancellation::CancelToken).
    #[error("apex construction cancelled")]
    Cancelled,
}

/// Combined error type for [`resolve_apex`](super::resolve_apex)
/// (combined plan + construct).  Either planning failed
/// ([`ApexPlanError`]) or construction failed ([`ApexBuildError`]).
#[derive(thiserror::Error, Debug)]
pub enum ApexError {
    /// Planning failed.
    #[error(transparent)]
    Plan(#[from] ApexPlanError),
    /// Construction failed.
    #[error(transparent)]
    Build(#[from] ApexBuildError),
}
