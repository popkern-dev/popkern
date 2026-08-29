//! Blend-option enums split out of `types.rs` (500-line cap):
//! the cross-section profile (PK-189) and the strip-interference
//! resolution policy (PK-185).

/// Cross-section profile of a blend (requirements 2.2, PK-189).
///
/// The rolling-ball contact curves and end tangents are identical
/// for every profile — only the section INTERIOR differs — so
/// profiles compose with variable radius and with the corner-cap
/// dispatch exactly like the circular sweep does.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BlendProfile {
    /// Circular rolling-ball section (default; the historical
    /// behaviour, analytic torus/cylinder fast paths eligible).
    #[default]
    Circular,
    /// Conic (rational quadratic) section with an explicit shoulder
    /// factor. `rho` is the classic conic fullness parameter in
    /// (0, 1): the section's shoulder point sits at
    /// `M + rho·(A − M)` where `M` is the contact-chord midpoint
    /// and `A` the tangent-intersection apex. `rho = 0.5` is
    /// parabolic; smaller is flatter, larger is fuller. At each
    /// station the CIRCULAR section corresponds to
    /// `rho = cos(θ/2) / (1 + cos(θ/2))` for section opening angle
    /// `θ`, and supplying that value reproduces the circular blend
    /// to roundoff — the conic path subsumes the arc, it does not
    /// fork it.
    ///
    /// Conic sections are exact rational-quadratic rows — there is
    /// no approximation or fit. Stations whose opening angle
    /// exceeds the single-segment rational range decline with a
    /// structured error (no silent downgrade to circular).
    Conic {
        /// Shoulder fullness in (0, 1); see the variant docs.
        rho: f64,
    },
    /// Curvature-continuous degree-5 section — the existing G2
    /// sweep. Selecting this is equivalent to requesting
    /// [`crate::types::ContinuityLevel::G2`].
    CurvatureContinuous,
}

impl BlendProfile {
    /// Rational-quadratic middle weight for a conic profile:
    /// `w = rho / (1 − rho)`.
    pub fn conic_weight(rho: f64) -> f64 {
        rho / (1.0 - rho)
    }

    /// The rho value at which a conic section reproduces the
    /// circular arc for a station with opening angle
    /// `theta` (radians): `cos(θ/2) / (1 + cos(θ/2))`.
    pub fn circular_rho(theta: f64) -> f64 {
        let w = (0.5 * theta).cos();
        w / (1.0 + w)
    }
}

/// PK-185 — resolution policy for non-completable mid-face strip
/// interference. See [`crate::types::FilletOptions::interference_policy`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InterferencePolicy {
    /// Decline the whole operation transactionally (default — the
    /// historical behaviour).
    #[default]
    Decline,
    /// Drop both members of each interfering pair and serve the
    /// remaining fillets; every dropped pair is reported with a
    /// [`BlendWarning::FilletPairDropped`](crate::error::BlendWarning::FilletPairDropped)
    /// carrying the reason. When nothing remains after dropping,
    /// the operation still declines (there is no partial result to
    /// serve).
    ServeNonInterfering,
}
