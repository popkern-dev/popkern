//! Plan telemetry for Phase 6c-A2 apex construction.
//!
//! Captures the planner's reasoning and numeric inputs so the
//! caller (especially UIs) can:
//!
//! * Display human-readable explanations of cap-kind decisions.
//! * Show diagnostic numerics (cap centre, sphere radius, condition
//!   numbers, T positions).
//! * Identify which fallback rung was used and why.
//! * Serialise to JSON for offline inspection or remote debugging.
//!
//! Replaces the legacy 5+ env-flag tracing approach with a
//! structured query interface.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Construction rung enumeration.  Ordered from preferred (Direct)
/// to last-resort (Open).  The planner attempts each rung in order
/// and emits a [`PlanTelemetry`] recording which rung succeeded.
///
/// Each rung corresponds to a [`BlendWarning::ContinuityDowngraded`]
/// emission when a more-preferred rung was attempted but failed
/// (legacy fillet pipeline pattern; preserved here).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum ConstructionRung {
    /// Direct sphere-fit apex construction — the A2 default.
    /// All fillets meet at sphere-fit terminations; no compromise.
    Direct,

    /// Setback rung: terminate fillets short of the apex by
    /// `setback_distance`, leaving a small gap.  Used when Direct
    /// fails sphere-fit math (e.g. sphere doesn't intersect
    /// cylinder, cap centre singular).  Lower-quality but still
    /// well-defined.
    Setback {
        /// Distance to set each fillet back from its computed
        /// termination, in metres.
        setback_distance: f64,
    },

    /// Open rung: produce no apex cap at all.  The fillets just
    /// terminate at their nominal endpoints, leaving an exposed
    /// gap.  Last resort when Setback also fails.
    Open,
}

/// Plan telemetry — explains why the planner chose what it did.
///
/// Every [`ApexPlan`](super::ApexPlan) carries a `PlanTelemetry`
/// describing the decision process.  UIs render this for users
/// to understand the fillet result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanTelemetry {
    /// Human-readable reason for the cap-kind choice.
    /// E.g. "uniform radii within 0.1% — Sphere cap chosen" or
    /// "variable radii (0.04..0.07) — Gregory cap chosen".
    pub cap_kind_reason: String,

    /// Construction rung used by this plan (Direct, Setback, Open).
    pub rung: ConstructionRung,

    /// Reason for any fallback (only populated when `rung != Direct`).
    /// E.g. "sphere does not intersect fillet[2]'s cylinder" for
    /// Setback after a SphereCylinderNoIntersection.
    pub fallback_reason: Option<String>,

    /// Numeric diagnostic values keyed by name.  Examples:
    /// * `"cap_centre_x"`, `"cap_centre_y"`, `"cap_centre_z"`
    /// * `"sphere_radius"`
    /// * `"max_bridge_collapse"` — the largest |T_a[k+1] − T_b[k]|;
    ///   used for Sphere/Gregory dispatch.
    /// * `"condition_number"` — cap-centre system conditioning.
    pub diagnostic_values: BTreeMap<String, f64>,
}

impl PlanTelemetry {
    /// Create an empty telemetry block, defaulting to Direct rung.
    /// Callers populate fields as planning progresses.
    pub fn new() -> Self {
        Self {
            cap_kind_reason: String::new(),
            rung: ConstructionRung::Direct,
            fallback_reason: None,
            diagnostic_values: BTreeMap::new(),
        }
    }

    /// Record a numeric diagnostic value keyed by name.
    pub fn set_value(&mut self, key: impl Into<String>, value: f64) {
        self.diagnostic_values.insert(key.into(), value);
    }
}

impl Default for PlanTelemetry {
    fn default() -> Self {
        Self::new()
    }
}
