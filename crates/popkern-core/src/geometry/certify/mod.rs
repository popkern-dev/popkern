//! Certified (interval) enclosures of surface and pcurve jets, the
//! substrate of certified seam continuity (PK-499; design in
//! `docs/architecture/certified-seam-continuity.md`).
//!
//! Everything here computes rigorous enclosures: outward-rounded
//! [`Iv`] arithmetic ([`interval`]), forward-mode interval derivatives
//! ([`Taylor`] series), Bernstein-hull enclosures of rational spline spans
//! and closed forms for the analytic kinds. A kind that cannot be
//! enclosed rigorously is declined with a typed [`CertifyDecline`],
//! never approximated.

pub mod analytic_jet;
pub mod bernstein;
mod boundary_surrogate;
pub mod interval;
pub mod interval_vec;
pub mod nurbs_jet;
mod pcurve_affine;
pub mod pcurve_jet;
pub mod seam;
mod seam_interval;
pub mod seam_projection;
mod seam_refine;
pub mod seam_state;
pub mod surface_jet;
pub mod taylor;

pub use interval::Iv;
pub use interval_vec::{IvVec2, IvVec3};
pub use pcurve_jet::{PcurveJet, PcurveJetSource};
pub use seam::{
    CertificateStatus, CertifyTargets, InconclusiveReason, SeamCertificate, SeamCertifyOptions,
    SeamCorrespondence, certify_seam,
};
pub use seam_state::SeamSideSource;
pub use surface_jet::{SurfaceJet, SurfaceJetSource};
pub use taylor::{Taylor, TaylorVec3};

/// A geometry form the certifier has no exact enclosure for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum UnsupportedGeometry {
    /// Dual-twist bi-nonic Gregory: its NURBS form is a midpoint
    /// containment, not exact.
    DualTwistGregory,
    /// Triangular Gregory: its exact form uses another chart than the
    /// patch evaluator.
    TriangularGregory,
    /// Rational-quadratic arc pcurve spanning half a turn or more
    /// (atan2 no longer reduces to atan).
    RationalQuadArcHalfTurn,
}

/// Why an enclosure could not be produced. Every decline is an
/// honest "cannot certify", never a pass.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, thiserror::Error, serde::Serialize, serde::Deserialize,
)]
pub enum CertifyDecline {
    /// The surface kind has no exact enclosure in its own chart.
    #[error("unsupported surface: {0:?}")]
    UnsupportedSurface(UnsupportedGeometry),
    /// The pcurve kind has no exact enclosure.
    #[error("unsupported pcurve: {0:?}")]
    UnsupportedPcurve(UnsupportedGeometry),
    /// A pcurve with an empty parameter range.
    #[error("degenerate pcurve parameter range")]
    DegeneratePcurve,
    /// Knot, weight or control-point counts disagree.
    #[error("malformed NURBS data")]
    MalformedNurbs,
    /// A weight is zero, negative or not finite.
    #[error("non-positive rational weight")]
    NonPositiveWeight,
    /// The rational denominator's enclosure contains zero (a Gregory
    /// corner or a vanishing weight function).
    #[error("rational denominator enclosure contains zero")]
    DenominatorContainsZero,
    /// A parameter box straddles an interior knot; the caller splits
    /// there.
    #[error("parameter box straddles a knot")]
    SpanStraddle,
    /// An enclosure overflowed or became undefined.
    #[error("non-finite enclosure")]
    NonFinite,
}

#[cfg(test)]
#[path = "certify_substrate_tests.rs"]
mod certify_substrate_tests;

#[cfg(test)]
#[path = "certify_seam_tests.rs"]
mod certify_seam_tests;

#[cfg(test)]
#[path = "certify_projection_tests.rs"]
mod certify_projection_tests;
