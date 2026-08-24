//! Typed declines for surface offsetting.
//!
//! An offset decline is a property of the *request* (surface class ×
//! distance × direction), not of one caller's health, so the reasons are
//! typed here where the geometry is decided. Callers translate them into
//! their own error vocabulary (e.g. `ShellError::OffsetTooLarge`)
//! without re-deriving the geometry.

/// Why a surface offset could not be produced.
#[derive(Debug, Clone, thiserror::Error)]
pub enum OffsetError {
    /// A radius-form surface (cylinder/sphere/torus minor) was driven
    /// through zero: the offset surface does not exist at the requested
    /// distance. Approximating it would serve inverted geometry.
    #[error(
        "{class} offset radius would be {new_radius:.6} (≤ 0); \
         max safe offset {max_safe:.6}"
    )]
    RadiusThroughZero {
        /// Surface class name ("cylinder", "sphere", "torus").
        class: &'static str,
        /// The (non-positive) radius the offset would have produced.
        new_radius: f64,
        /// Largest offset distance in this direction that keeps the
        /// radius positive (0.99 of the exact bound, matching the
        /// collapse-gate serve convention).
        max_safe: f64,
    },

    /// A torus offset grew the minor radius to or past the major
    /// radius: the result is a self-intersecting spindle torus. Routed
    /// to a decline until global trimming (PK-327) can serve it.
    #[error(
        "torus offset minor radius {new_minor:.6} ≥ major radius \
         {major:.6} (spindle torus); max safe offset {max_safe:.6}"
    )]
    SpindleTorus {
        /// The minor radius the offset would have produced.
        new_minor: f64,
        /// The (unchanged) major radius it would have met or crossed.
        major: f64,
        /// Largest offset distance in this direction that keeps
        /// `minor < major` (0.99 of the exact bound).
        max_safe: f64,
    },

    /// The surface's defining axis or normal is degenerate (zero
    /// length, or a cone half-angle outside `(0, π/2)`).
    #[error("{class} has degenerate defining geometry")]
    DegenerateGeometry {
        /// Surface class name.
        class: &'static str,
    },

    /// No closed-form offset arm exists for this surface class.
    /// Produced only by the compatibility wrapper — dispatch code uses
    /// [`super::AnalyticOffset::Unsupported`] to route to approximation.
    #[error("analytic offset not supported for {class} surfaces")]
    Unsupported {
        /// Surface class name.
        class: &'static str,
    },

    /// The offset *map* folds inside the sampled face window (distance
    /// exceeds the local curvature radius), so a sampled approximation
    /// would fit inverted geometry. Only raised on the approximation
    /// path — classes with exact offsets handle folds downstream via
    /// trimming and the collapse gate.
    #[error(
        "offset map folds within the face window (distance {requested:.6} \
         exceeds the local curvature radius); max safe offset {max_safe:.6}"
    )]
    OffsetMapFolds {
        /// The requested offset distance.
        requested: f64,
        /// Largest fold-free offset distance in this direction.
        max_safe: f64,
    },

    /// The least-squares fit machinery failed (singular normal
    /// equations — degenerate sampling).
    #[error("NURBS offset fit failed: {reason}")]
    FitFailed {
        /// Why the fit could not be produced.
        reason: String,
    },
}

impl OffsetError {
    /// The largest safe offset distance carried by this decline, when
    /// the reason bounds one (radius-through-zero, spindle, fold).
    pub fn max_safe(&self) -> Option<f64> {
        match self {
            OffsetError::RadiusThroughZero { max_safe, .. }
            | OffsetError::SpindleTorus { max_safe, .. }
            | OffsetError::OffsetMapFolds { max_safe, .. } => Some(*max_safe),
            _ => None,
        }
    }
}
