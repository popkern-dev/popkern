//! Heaviside projections of the filtered density x̃ to the physical
//! density x̄, identical to the density filter as β → 0 and a step as β
//! grows: they remove the grey band a density filter leaves (PH-114:
//! 91% of a coarse cantilever's entries grey).
//!
//! - Guest et al. 2004 as top88 implements it (Andreassen 2011 §5.4,
//!   eqns 21–23): x̄ = 1 − e^(−βx̃) + x̃ e^(−β), pushing every x̃ above 0
//!   towards 1. It raises the volume each time β grows.
//! - The threshold projection of Wang, Lazarov & Sigmund 2011 (eq 9):
//!   x̄ = (tanh βη + tanh β(x̃ − η)) / (tanh βη + tanh β(1 − η)),
//!   projecting about η; at η = ½ it nearly keeps the volume as β grows.
//!
//! Under the optimality criteria update neither settled (measured
//! 2026-10-13: the threshold form two-cycled at the full move limit at
//! any β; Guest's stalled on 2 of 5 grid offsets of a thin cantilever);
//! MMA carries them (PH-122).

use popkern_math::PortableMath;

/// Which projection.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub enum ProjectionKind {
    /// Guest's exponential projection towards solid.
    Guest,
    /// The tanh projection about a threshold η in (0, 1).
    Threshold(f64),
}

/// A projection at sharpness β (> 0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Projection {
    /// Sharpness.
    pub beta: f64,
    /// The projection's form.
    pub kind: ProjectionKind,
}

impl Projection {
    /// x̄ of `t`.
    pub fn apply(&self, t: f64) -> f64 {
        self.curve().apply(t)
    }

    /// dx̄/dx̃ at `t`.
    pub fn slope(&self, t: f64) -> f64 {
        self.curve().slope(t)
    }

    /// The projection with its terms that do not depend on x̃ evaluated
    /// once, for a sweep over a design's entries: the same arithmetic,
    /// so the same bits as [`Projection::apply`], at a third of the
    /// transcendentals for the threshold form (the update's volume
    /// closure applies it ~60 times per update, PK-1223).
    pub fn curve(&self) -> ProjectionCurve {
        let b = self.beta;
        match self.kind {
            ProjectionKind::Guest => ProjectionCurve::Guest {
                beta: b,
                tail: (-b).pexp(),
            },
            ProjectionKind::Threshold(eta) => {
                let low = (b * eta).ptanh();
                ProjectionCurve::Threshold {
                    beta: b,
                    eta,
                    low,
                    span: low + (b * (1.0 - eta)).ptanh(),
                }
            }
        }
    }
}

/// A [`Projection`] with its constant terms evaluated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProjectionCurve {
    /// Guest's form; `tail` is e^(−β).
    Guest {
        /// Sharpness.
        beta: f64,
        /// e^(−β).
        tail: f64,
    },
    /// The threshold form; `low` is tanh βη, `span` tanh βη + tanh β(1 − η).
    Threshold {
        /// Sharpness.
        beta: f64,
        /// Threshold.
        eta: f64,
        /// tanh βη.
        low: f64,
        /// tanh βη + tanh β(1 − η).
        span: f64,
    },
}

impl ProjectionCurve {
    /// x̄ of `t`.
    pub fn apply(&self, t: f64) -> f64 {
        match *self {
            ProjectionCurve::Guest { beta, tail } => 1.0 - (-beta * t).pexp() + t * tail,
            ProjectionCurve::Threshold {
                beta,
                eta,
                low,
                span,
            } => (low + (beta * (t - eta)).ptanh()) / span,
        }
    }

    /// dx̄/dx̃ at `t`.
    pub fn slope(&self, t: f64) -> f64 {
        match *self {
            ProjectionCurve::Guest { beta, tail } => beta * (-beta * t).pexp() + tail,
            ProjectionCurve::Threshold {
                beta, eta, span, ..
            } => {
                let th = (beta * (t - eta)).ptanh();
                beta * (1.0 - th * th) / span
            }
        }
    }
}
