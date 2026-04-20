//! Result and diagnostic types for direct-construction (Phase G).
//!
//! These types are the programmatic API surface for
//! visualisation, UI preview, logging, and per-corner quality
//! badging — "UI-ready diagnostics" requirement from the Phase G
//! plan.

use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::healing::quality::ContinuityLevel;

/// Runtime validation intensity for [`super::direct_options::DirectConstructionOptions`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationLevel {
    /// No self-check. Fastest path. For trusted batch pipelines
    /// where you've already validated the inputs.
    None,
    /// Default. Post-build G1 identity residual self-check.
    /// Fires [`super::direct_errors::DirectConstructionError::SelfCheckFailedBug`]
    /// on residual above `self_check_g1_tol` — this is a bug
    /// detector, not a user-facing error path.
    Standard,
    /// Standard plus dense residual sampling suitable for test
    /// harnesses and UI preview. Populates every `residual`
    /// field in the diagnostics; slower.
    Full,
}

/// Which sort of corner a [`CornerDiagnostic`] describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CornerKind {
    /// The common central vertex `C` where all `N` sub-patches
    /// meet (origin of every spoke).
    Central,
    /// An outer hex corner `V_{2k+2}` shared between two
    /// adjacent sub-patches.
    OuterHex,
}

/// Which published formula was used at this corner or spoke.
/// Exposed so a GUI can badge "C-K 1983 baseline" vs
/// "C-K tangent-corner variant" vs the variable-radius form, and
/// so telemetry can track which formula families the kernel
/// exercises on real workloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormulaFamily {
    /// Chiyokura-Kimura 1983 baseline G1 construction.
    ChiyokuraKimuraStandard,
    /// C-K with the tangent-corner variant (coplanar adjacent
    /// faces at the corner).
    ChiyokuraKimuraTangentCorner,
    /// Kato-Kimura-Chiyokura 1991 vertex-blend-specific form —
    /// handles variable radius cleanly.
    KatoKimuraChiyokura1991,
    /// Shirman-Séquin N-way vertex enclosure for the central
    /// vertex corner.
    ShirmanSequinEnclosure,
    /// Peters 2000 two-stripe row-1 extension (extra DOFs
    /// beyond classical C-K; enables exact G1 on asymmetric
    /// geometry).
    PetersTwoStripeRow1,
}

/// Per-corner diagnostic data — one entry per sub-patch corner.
#[derive(Debug, Clone)]
pub struct CornerDiagnostic {
    /// Which kind of corner.
    pub corner_kind: CornerKind,
    /// Which published formula was used here.
    pub formula_family: FormulaFamily,
    /// In `[0.0, 1.0]`. `1.0` = perfectly well-conditioned
    /// (symmetric, radii all similar); `0.0` = on the edge of a
    /// degenerate case. A future GUI badges corners below ~0.5.
    pub near_degeneracy_score: f64,
    /// Which tolerance bound was binding on the construction here.
    pub tolerance_used: f64,
}

/// Per-spoke diagnostic data — one entry per shared spoke (= N
/// per N-sided vertex blend).
#[derive(Debug, Clone)]
pub struct SpokeDiagnostic {
    /// Formula family used for this spoke's row-1 CPs.
    pub formula_family: FormulaFamily,
    /// Dense `1 − |n_A · n_B|` max over [0, 1], measured at
    /// `ValidationLevel::Full`. `NaN` if not measured.
    pub g1_residual: f64,
    /// Dense second-fundamental-form residual for G2 target.
    /// `None` if G1-only target. `NaN` if not measured.
    pub g2_residual: Option<f64>,
    /// Degree of the `α(t)` blending polynomial chosen.
    pub alpha_polynomial_degree: usize,
    /// Degree of the `β(t)` blending polynomial chosen.
    pub beta_polynomial_degree: usize,
    /// `[0.0, 1.0]` conditioning score for this spoke.
    pub near_degeneracy_score: f64,
}

/// The full output of a successful direct-construction build.
#[derive(Debug, Clone)]
pub struct DirectConstructionResult {
    /// The `N` constructed two-stripe sub-patches, in spoke
    /// order. Callers install these onto the BRep's sub-faces
    /// at the three-quad poke positions.
    pub sub_patches: Vec<RectangularGregoryTwoStripe>,
    /// Continuity level actually delivered. Usually matches
    /// `options.target_continuity` but may degrade (e.g. G2→G1)
    /// on tangent corners where G2 is geometrically
    /// ill-defined.
    pub achieved_continuity: ContinuityLevel,

    /// Per-corner diagnostics in `sub_patches` order (one per
    /// sub-patch per corner kind — 2·N entries total for a
    /// three-quad N-gon).
    pub corner_diagnostics: Vec<CornerDiagnostic>,
    /// Per-spoke diagnostics in canonical spoke order.
    pub spoke_diagnostics: Vec<SpokeDiagnostic>,

    /// Dense-sample witness: `min_cos ≥ 0.9999` for commercial
    /// grade. Populated at `Standard` or `Full` validation.
    /// `NaN` at `ValidationLevel::None`.
    pub dense_min_cos_g1: f64,
    /// Dense-sample G2 witness. `None` when target is G1.
    /// `NaN` when target is G2 but validation is `None`.
    pub dense_min_cos_g2: Option<f64>,
    /// Max absolute position gap at any shared boundary sample.
    /// Expected `≤ 1e-12` for exact construction.
    pub max_position_gap_g0: f64,
}

/// Warnings emitted by [`super::direct_build::can_build`] that
/// don't invalidate the feasibility decision but are useful for
/// UI preview.
#[derive(Debug, Clone)]
pub enum FeasibilityWarning {
    /// Build would proceed but with tangent-corner variant at
    /// the named corner.
    TangentCornerVariant {
        /// Index of the corner involved.
        corner_idx: usize,
    },
    /// Asymmetry is below the rejection ceiling but close to it.
    AsymmetryNearCeiling {
        /// Measured max/min radius ratio.
        ratio: f64,
        /// Ceiling from the options struct.
        ceiling: f64,
    },
    /// N ≥ 5 — direct construction will fall through to polar.
    ValenceRequiresPolar {
        /// Valence that triggered the polar fallback.
        n: usize,
    },
    /// A specific corner has a low conditioning score.
    CornerNearDegenerate {
        /// Index of the weak corner.
        corner_idx: usize,
        /// Measured conditioning score in `[0.0, 1.0]`.
        score: f64,
    },
}

/// Output of the [`super::direct_build::can_build`] feasibility
/// query. Cheap — no CP emission.
#[derive(Debug, Clone)]
pub struct FeasibilityReport {
    /// Best continuity direct construction can deliver for this
    /// input. May be less than the requested target (e.g. G2
    /// requested but only G1 reachable due to tangent corners).
    pub achievable_continuity: ContinuityLevel,
    /// Worst conditioning score over all corners, in
    /// `[0.0, 1.0]`. 1.0 = clean, 0.0 = on the edge of a
    /// rejection.
    pub near_degeneracy_score: f64,
    /// Non-fatal warnings (tangent-corner variants, near-
    /// degenerate corners, etc.).
    pub warnings: Vec<FeasibilityWarning>,
    /// Rough estimate for UI budget forecasting. Not a hard
    /// promise; measured by the benchmark suite on reference
    /// fixtures.
    pub estimated_construction_us: u64,
}

#[cfg(test)]
#[path = "direct_types_tests.rs"]
mod tests;
