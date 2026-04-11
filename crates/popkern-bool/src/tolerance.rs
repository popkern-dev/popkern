//! Extended tolerance parameters for boolean operations.
//!
//! [`BoolTolerance`] wraps the core [`Tolerance`] and adds parameters specific
//! to the boolean pipeline: snapping, coincidence detection, sliver prevention,
//! and micro-edge elimination.

use popkern_core::tolerance::Tolerance;

/// Extended tolerance parameters for the boolean pipeline.
///
/// Maintains the ordering invariant:
/// `linear <= snap <= sliver <= coincidence`
///
/// Violating this ordering creates inconsistencies — e.g., the engine
/// might snap to prevent slivers but the coincidence threshold is narrower
/// than the snap range, causing snapped geometry to be misclassified.
///
/// # Example
///
/// ```
/// use popkern_core::tolerance::Tolerance;
/// use popkern_bool::tolerance::BoolTolerance;
///
/// let tol = BoolTolerance::from(Tolerance::default());
/// assert!(tol.snap >= tol.core.linear());
/// assert!(tol.sliver >= tol.snap);
/// assert!(tol.coincidence >= tol.sliver);
/// ```
#[derive(Debug, Clone)]
pub struct BoolTolerance {
    /// The underlying core tolerance (linear, angular, parametric).
    pub core: Tolerance,
    /// Snapping tolerance: intersection endpoints within this distance of
    /// an existing vertex are merged instead of creating a new vertex.
    /// Default: 2x linear tolerance.
    pub snap: f64,
    /// Coincidence tolerance: faces whose surfaces are everywhere within
    /// this distance are treated as coincident. Default: 10x linear tolerance.
    pub coincidence: f64,
    /// Sliver tolerance: minimum face width / minimum edge length below
    /// which faces are classified as slivers. Default: 5x linear tolerance.
    pub sliver: f64,
    /// Micro-edge tolerance: edges shorter than this are collapsed.
    /// Default: 1x linear tolerance.
    pub micro_edge: f64,
}

impl BoolTolerance {
    /// Create a new `BoolTolerance` with explicit parameters.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if the ordering invariant is violated:
    /// `core.linear <= snap <= sliver <= coincidence`.
    pub fn new(core: Tolerance, snap: f64, coincidence: f64, sliver: f64, micro_edge: f64) -> Self {
        debug_assert!(
            core.linear() <= snap && snap <= sliver && sliver <= coincidence,
            "BoolTolerance ordering violated: linear={} <= snap={} <= sliver={} <= coincidence={}",
            core.linear(),
            snap,
            sliver,
            coincidence
        );
        Self {
            core,
            snap,
            coincidence,
            sliver,
            micro_edge,
        }
    }

    /// Create with custom multipliers applied to the linear tolerance.
    ///
    /// This is useful when you want to scale all boolean tolerances
    /// proportionally to the base precision.
    pub fn with_multipliers(
        core: Tolerance,
        snap_mult: f64,
        coincidence_mult: f64,
        sliver_mult: f64,
        micro_edge_mult: f64,
    ) -> Self {
        Self::new(
            core,
            core.linear() * snap_mult,
            core.linear() * coincidence_mult,
            core.linear() * sliver_mult,
            core.linear() * micro_edge_mult,
        )
    }
}

impl From<Tolerance> for BoolTolerance {
    /// Create a `BoolTolerance` from a core `Tolerance` using sensible defaults:
    /// - `snap` = 2x linear
    /// - `sliver` = 5x linear
    /// - `coincidence` = 10x linear
    /// - `micro_edge` = 1x linear
    fn from(core: Tolerance) -> Self {
        let linear = core.linear();
        Self {
            core,
            snap: linear * 2.0,
            coincidence: linear * 10.0,
            sliver: linear * 5.0,
            micro_edge: linear,
        }
    }
}

impl From<&Tolerance> for BoolTolerance {
    fn from(core: &Tolerance) -> Self {
        Self::from(*core)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_ordering_invariant() {
        let bt = BoolTolerance::from(Tolerance::default());
        assert!(bt.core.linear() <= bt.snap);
        assert!(bt.snap <= bt.sliver);
        assert!(bt.sliver <= bt.coincidence);
    }

    #[test]
    fn custom_multipliers() {
        let core = Tolerance::default().with_linear(1e-4);
        let bt = BoolTolerance::with_multipliers(core, 3.0, 20.0, 8.0, 0.5);
        assert!((bt.snap - core.linear() * 3.0).abs() < 1e-15);
        assert!((bt.coincidence - core.linear() * 20.0).abs() < 1e-15);
        assert!((bt.sliver - core.linear() * 8.0).abs() < 1e-15);
        assert!((bt.micro_edge - core.linear() * 0.5).abs() < 1e-15);
    }

    #[test]
    fn from_ref() {
        let core = Tolerance::default();
        let bt = BoolTolerance::from(&core);
        assert!((bt.core.linear() - core.linear()).abs() < 1e-15);
    }
}
