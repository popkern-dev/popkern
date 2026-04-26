//! Error variants emitted by the factory and by `EdgeTriple::verify`.

use crate::geometry::Point3;
use crate::geometry::edge_triple::inputs::SurfaceSide;

/// Failure modes from
/// [`crate::geometry::edge_triple::compute_edge_geometry`]. The
/// factory never silently falls back; callers read the variant and
/// either tolerate (keep existing geometry) or propagate.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum EdgeTripleError {
    /// Surfaces do not intersect anywhere within the search region.
    #[error("surfaces do not intersect")]
    NoIntersection,
    /// An anchor point is more than `tol.linear()` away from one of
    /// the surfaces.
    #[error("anchor on surface {which:?} is {distance} off the surface")]
    AnchorOffSurface {
        /// Which surface the anchor is too far from.
        which: SurfaceSide,
        /// Achieved distance from anchor to the surface.
        distance: f64,
    },
    /// Marching solver hit `max_marching_iterations` without
    /// converging to within `tol.linear()`.
    #[error(
        "marching solver did not converge: residual {residual_at_abort} after {iterations} iterations"
    )]
    NonConvergent {
        /// Residual at the point the solver gave up.
        residual_at_abort: f64,
        /// Number of iterations consumed.
        iterations: u32,
    },
    /// Multiple valid branches and no hint provided. Caller must
    /// disambiguate.
    #[error("ambiguous intersection: {count} candidate branches; supply a hint")]
    AmbiguousBranch {
        /// Number of equally-good candidate branches.
        count: usize,
    },
    /// Cancellation token was triggered during a marching solve.
    #[error("cancelled by caller")]
    Cancelled,
    /// Surface pair is degenerate (coincident, zero-area, …).
    #[error("degenerate surface pair: {reason}")]
    DegenerateSurfacePair {
        /// Static-string reason for downstream diagnostic display.
        reason: &'static str,
    },
    /// The two surfaces are *identical* (coincident — same algebraic
    /// surface modulo scaling). The "intersection" is the entire
    /// surface, which is not a 1-D edge. Detected exactly via the
    /// symbolic Q[√d] path (`A = c·B` for some non-zero rational
    /// `c`). CAD scenario: duplicated surface in a BRep model.
    #[error("surfaces are coincident: {reason}")]
    CoincidentSurfaces {
        /// Static-string reason for downstream diagnostic display.
        reason: &'static str,
    },
    /// Stage 1 placeholder. Removed once dispatch stages land.
    #[error("dispatch path not yet implemented for this surface pair")]
    Unimplemented,
}

/// Failure mode from
/// [`crate::geometry::edge_triple::EdgeTriple::verify`].
#[derive(Debug, thiserror::Error)]
pub enum EdgeTripleVerifyError {
    /// G2 broken at parameter `t` with residual exceeding tolerance.
    #[error("positional invariant broken at t={t}: residual {residual}")]
    PositionalInvariantBroken {
        /// Curve parameter at which G2 was violated.
        t: f64,
        /// Achieved positional residual at `t`.
        residual: f64,
    },
    /// G3/G5 broken: `t_end < t_start`.
    #[error("direction reversed: t_end < t_start")]
    DirectionReversed,
    /// G4 broken: a pcurve CP lies outside the surface's UV box.
    #[error("pcurve on {side:?} has CP outside the surface UV box: {cp:?}")]
    PcurveOutsideUvBox {
        /// Which side's pcurve has the offending CP.
        side: SurfaceSide,
        /// The offending CP coordinate.
        cp: Point3,
    },
}
