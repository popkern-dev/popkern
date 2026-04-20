//! Options and tolerance knobs for the direct boundary-based
//! vertex-blend construction (Phase G).
//!
//! Per `docs/superpowers/plans/2026-04-19-vertex-blend-direct-construction.md`:
//! every numerical threshold is an explicit parameter here, never
//! a global `const EPSILON`. The struct is `Clone` so a caller
//! can tweak a single knob for a specific build without rebuilding
//! the whole thing.

use crate::geometry::gregory::direct_types::ValidationLevel;
use crate::healing::quality::ContinuityLevel;

/// Configuration for a direct-construction build.
///
/// Defaults are chosen conservatively: `ValidationLevel::Standard`
/// catches math bugs at construction time, tolerances are tight
/// enough to assert machine precision while loose enough to tolerate
/// ordinary accumulated f64 error from the upstream biseptic seed.
///
/// Override per-call — never a global mutable state.
#[derive(Debug, Clone)]
pub struct DirectConstructionOptions {
    /// Characteristic scale of the input BRep (typically the
    /// bounding-box diagonal). Several thresholds are interpreted
    /// as a fraction of this scale, so constructions on a 1 mm
    /// part and a 1 km part both pass the same tolerance tests.
    /// **Required** — no sensible default.
    pub characteristic_length: f64,

    /// Angular tolerance (radians) below which two boundary
    /// tangents at a corner are considered coplanar — dispatches
    /// the tangent-corner formula variant (§Chiyokura-Kimura 1983
    /// §3.4 degenerate case).
    pub tangent_corner_rad: f64,

    /// Minimum boundary-tangent magnitude relative to
    /// `characteristic_length`. Below this the boundary is
    /// degenerate and construction emits
    /// [`super::direct_errors::DirectConstructionError::BoundaryTangentDegenerate`].
    pub boundary_tangent_eps: f64,

    /// Smallest singular value of the corner-tangent-plane matrix
    /// below which the corner construction is ill-conditioned.
    pub corner_singularity_eps: f64,

    /// Threshold on the post-build G1 identity residual for the
    /// Standard-level self-check. Exceeding this fires
    /// [`super::direct_errors::DirectConstructionError::SelfCheckFailedBug`].
    /// This is a _math-bug detector_, not a user-input sanity
    /// check: users should never see it.
    pub self_check_g1_tol: f64,

    /// Threshold on the post-build G2 residual. Looser than G1
    /// because second derivatives are more numerically delicate.
    pub self_check_g2_tol: f64,

    /// Minimum ratio `radius / characteristic_length`. Below this
    /// the blend is degenerate; rejected with `RadiusTooSmall`.
    pub min_radius_ratio: f64,

    /// Maximum ratio `radius / characteristic_length`. Above this
    /// the blend is geometrically unreasonable; rejected with
    /// `RadiusTooLarge`.
    pub max_radius_ratio: f64,

    /// Maximum ratio `max_radius / min_radius` across a vertex's
    /// incident edges. Beyond this the direct formulas are
    /// numerically poorly conditioned and the caller should fall
    /// back to the sampled-damped path.
    pub max_radius_asymmetry: f64,

    /// Runtime validation intensity. `Standard` is the default
    /// and runs a post-build G1 self-check; `Full` adds dense
    /// residual sampling suitable for test harnesses and UI
    /// preview; `None` skips all checks.
    pub validation_level: ValidationLevel,

    /// Target continuity — `G1` for mechanical CAD, `G2` for
    /// Class-A automotive surfacing. `G2` requires the
    /// second-derivative formulas in [`super::direct_row2_g2`]
    /// and produces a bit more work per corner.
    pub target_continuity: ContinuityLevel,
}

impl DirectConstructionOptions {
    /// Standard-default options for a given characteristic length.
    /// All tolerance knobs take the values recommended in the
    /// Phase G plan.
    pub fn with_characteristic_length(characteristic_length: f64) -> Self {
        Self {
            characteristic_length,
            tangent_corner_rad: 1e-6,
            boundary_tangent_eps: 1e-10,
            corner_singularity_eps: 1e-8,
            self_check_g1_tol: 1e-9,
            self_check_g2_tol: 1e-7,
            min_radius_ratio: 1e-4,
            max_radius_ratio: 10.0,
            max_radius_asymmetry: 50.0,
            validation_level: ValidationLevel::Standard,
            target_continuity: ContinuityLevel::G1,
        }
    }

    /// Clone with a different target continuity.
    pub fn with_continuity(mut self, target: ContinuityLevel) -> Self {
        self.target_continuity = target;
        self
    }

    /// Clone with a different validation level.
    pub fn with_validation(mut self, level: ValidationLevel) -> Self {
        self.validation_level = level;
        self
    }

    /// Minimum acceptable edge radius in absolute units (not a
    /// ratio). Convenience for callers that already know their
    /// minimum radius rather than the characteristic length.
    pub fn min_absolute_radius(&self) -> f64 {
        self.min_radius_ratio * self.characteristic_length
    }

    /// Maximum acceptable edge radius in absolute units.
    pub fn max_absolute_radius(&self) -> f64 {
        self.max_radius_ratio * self.characteristic_length
    }
}

#[cfg(test)]
#[path = "direct_options_tests.rs"]
mod tests;
