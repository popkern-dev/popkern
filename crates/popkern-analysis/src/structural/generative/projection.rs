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
        let b = self.beta;
        match self.kind {
            ProjectionKind::Guest => 1.0 - (-b * t).pexp() + t * (-b).pexp(),
            ProjectionKind::Threshold(eta) => {
                let low = (b * eta).ptanh();
                (low + (b * (t - eta)).ptanh()) / (low + (b * (1.0 - eta)).ptanh())
            }
        }
    }

    /// dx̄/dx̃ at `t`.
    pub fn slope(&self, t: f64) -> f64 {
        let b = self.beta;
        match self.kind {
            ProjectionKind::Guest => b * (-b * t).pexp() + (-b).pexp(),
            ProjectionKind::Threshold(eta) => {
                let th = (b * (t - eta)).ptanh();
                b * (1.0 - th * th) / ((b * eta).ptanh() + (b * (1.0 - eta)).ptanh())
            }
        }
    }
}
