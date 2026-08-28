//! Error types for direct construction.
//!
//! Every variant names the specific corner / spoke / edge /
//! tolerance that triggered it, and in most cases a suggested
//! remediation. This is the UX-quality gate from the direct-construction
//! plan's "Actionable errors" requirement — no
//! `Unsupported(String)` catch-alls.

use crate::geometry::gregory::direct_types::CornerKind;

/// Reasons a direct-construction build can fail.
///
/// The caller chooses the remediation:
/// - For `SelfCheckFailedBug`: file a report, retry is unlikely
///   to help (deterministic math bug on this input).
/// - For everything else: fall back to the sampled-damped path
///   and emit a remediation warning via the
///   `CornerRemediationDirectConstructionFailed` machinery.
#[derive(Debug, Clone, thiserror::Error)]
pub enum DirectConstructionError {
    /// The boundary tangent magnitude at `t = t` on spoke
    /// `spoke_idx` is below `tolerance · characteristic_length`.
    /// The spoke curve is degenerately parameterised.
    #[error(
        "boundary tangent near-zero at spoke {spoke_idx} t={t}: \
         |T|={magnitude} (threshold {tolerance})"
    )]
    BoundaryTangentDegenerate {
        /// Index of the offending spoke in the ring.
        spoke_idx: usize,
        /// Parameter value at which the tangent was evaluated.
        t: f64,
        /// Tangent magnitude measured.
        magnitude: f64,
        /// Threshold that was violated.
        tolerance: f64,
    },

    /// The corner-tangent-plane matrix has a singular value
    /// below the threshold — the N incident tangents are
    /// nearly colinear, so the tangent plane at this corner is
    /// not well defined.
    #[error(
        "corner {corner_idx} ({kind:?}): tangent-plane singular \
         (smallest SV {smallest_sv} < {tolerance}); \
         try reducing radius or excluding edge {edge_hint}"
    )]
    CornerTangentPlaneSingular {
        /// Index of the offending corner.
        corner_idx: usize,
        /// Which kind of corner (Central or OuterHex).
        kind: CornerKind,
        /// Smallest singular value of the tangent-plane matrix.
        smallest_sv: f64,
        /// Threshold that was violated.
        tolerance: f64,
        /// Hint naming the edge most likely responsible.
        edge_hint: String,
    },

    /// An input edge's radius is too small relative to the
    /// characteristic length. Direct construction requires
    /// `radius / characteristic_length ≥ min_radius_ratio`.
    #[error(
        "radius too small at edge {edge_id}: {radius} < {min_allowed} \
         (characteristic length {char_len})"
    )]
    RadiusTooSmall {
        /// Identifier of the offending edge (for user-facing messages).
        edge_id: String,
        /// Input radius value at that edge.
        radius: f64,
        /// Minimum allowed radius (absolute units).
        min_allowed: f64,
        /// Characteristic length used for the ratio test.
        char_len: f64,
    },

    /// An input edge's radius exceeds the maximum ratio.
    #[error("radius too large at edge {edge_id}: {radius} > {max_allowed}")]
    RadiusTooLarge {
        /// Identifier of the offending edge.
        edge_id: String,
        /// Input radius value.
        radius: f64,
        /// Maximum allowed radius (absolute units).
        max_allowed: f64,
    },

    /// The max/min radius ratio at a vertex exceeds the
    /// asymmetry ceiling. Direct formulas are poorly
    /// conditioned; caller should fall back to sampled-damped.
    #[error(
        "radius asymmetry {ratio} exceeds {max_allowed} at vertex — \
         fall back to sampled-damped path"
    )]
    RadiusAsymmetryTooHigh {
        /// Measured max/min radius ratio.
        ratio: f64,
        /// Ceiling from the options struct.
        max_allowed: f64,
    },

    /// Zero-radius edge in the blend set. A zero-radius edge is
    /// a sharp edge and cannot be filleted; the caller should
    /// drop it from the blend request.
    #[error(
        "zero-radius edge {edge_id} in the blend set — \
         direct construction requires every edge to have a positive radius"
    )]
    RadiusZero {
        /// Identifier of the zero-radius edge.
        edge_id: String,
    },

    /// Two spokes coincide — the hex topology is degenerate.
    #[error("degenerate hex topology at vertex: spokes {spoke_a} and {spoke_b} coincide")]
    DegenerateHexTopology {
        /// First coinciding spoke index.
        spoke_a: usize,
        /// Second coinciding spoke index.
        spoke_b: usize,
    },

    /// Valence ≥ 5 with the narrow-pie-slice fold. Out of
    /// scope for the Gregory direct construction; caller falls back to polar-NURBS.
    #[error(
        "N={n} requires a boundary-level fold fix (narrow pie slices fold); \
         falling back to polar-NURBS"
    )]
    ValenceRequiresPolar {
        /// Valence of the offending vertex.
        n: usize,
    },

    /// Post-build self-check detected residual above threshold.
    /// This indicates a **math bug in direct construction**, not
    /// a user error. Every trace of it in a log is a bug report.
    #[error(
        "construction emitted CPs but post-solve self-check failed: \
         dense G1 min_cos {measured} < {required} — LIKELY A MATH BUG, \
         please report with the offending fixture"
    )]
    SelfCheckFailedBug {
        /// Dense min_cos actually measured post-build.
        measured: f64,
        /// Required min_cos per `self_check_g1_tol`.
        required: f64,
    },

    /// Caller cancelled the build via the cancellation flag.
    #[error("cancelled by caller")]
    Cancelled,

    /// Input has `n_sub_patches < 3` — direct construction only
    /// applies to proper N-gon vertex blends (N ≥ 3). Common
    /// cause: the decomposition builder produced 0 or 2
    /// sub-patches for a degenerate input.
    #[error("too few sub-patches for direct construction: got {got}, need ≥ 3")]
    TooFewSubPatches {
        /// Count actually received.
        got: usize,
    },

    /// The upstream three-quad decomposition couldn't build the
    /// input the caller provided. Detail string should name the
    /// specific topological or geometric issue.
    #[error("upstream decomposition failed: {detail}")]
    DecompositionFailed {
        /// Specific failure reason from the decomposition layer.
        detail: String,
    },
}

impl DirectConstructionError {
    /// True if this is a "fall back to sampled-damped" error.
    /// `SelfCheckFailedBug` and `Cancelled` do NOT trigger
    /// fallback (one is a bug, the other is user intent).
    pub fn should_fall_back_to_sampled(&self) -> bool {
        matches!(
            self,
            Self::BoundaryTangentDegenerate { .. }
                | Self::CornerTangentPlaneSingular { .. }
                | Self::RadiusAsymmetryTooHigh { .. }
                | Self::ValenceRequiresPolar { .. }
                | Self::DecompositionFailed { .. }
        )
    }

    /// Short machine-readable tag for UI / telemetry.
    pub fn tag(&self) -> &'static str {
        match self {
            Self::BoundaryTangentDegenerate { .. } => "boundary_tangent_degenerate",
            Self::CornerTangentPlaneSingular { .. } => "corner_tangent_plane_singular",
            Self::RadiusTooSmall { .. } => "radius_too_small",
            Self::RadiusTooLarge { .. } => "radius_too_large",
            Self::RadiusAsymmetryTooHigh { .. } => "radius_asymmetry_too_high",
            Self::RadiusZero { .. } => "radius_zero",
            Self::DegenerateHexTopology { .. } => "degenerate_hex_topology",
            Self::ValenceRequiresPolar { .. } => "valence_requires_polar",
            Self::SelfCheckFailedBug { .. } => "self_check_failed_bug",
            Self::Cancelled => "cancelled",
            Self::TooFewSubPatches { .. } => "too_few_sub_patches",
            Self::DecompositionFailed { .. } => "decomposition_failed",
        }
    }
}

#[cfg(test)]
#[path = "direct_errors_tests.rs"]
mod tests;
