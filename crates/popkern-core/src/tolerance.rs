//! Structured tolerance type for geometric comparisons.
//!
//! [`Tolerance`] encapsulates six precision thresholds with a strict ordering
//! invariant: `numerical_floor < solver < linear`. Fields are private; access
//! via getter methods and construct via [`Default`], builder methods,
//! [`Tolerance::from_model_scale`], or [`Tolerance::preset`].

use crate::geometry::Point3;
use serde::{Deserialize, Serialize};

/// Named tolerance presets for common workflows.
///
/// Use with [`Tolerance::preset`] to create a tolerance scaled to model size.
///
/// | Preset | linear | solver | approximation |
/// |---|---|---|---|
/// | Draft | 1e-4 | 1e-8 | 1e-3 |
/// | Standard | 1e-6 | 1e-10 | 1e-5 |
/// | Fine | 1e-8 | 1e-12 | 1e-7 |
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TolerancePreset {
    /// Coarse -- fast operations, lower quality. For draft/preview work.
    Draft,
    /// Balanced precision and performance. Default for most workflows.
    Standard,
    /// High precision, slower operations. For final production geometry.
    Fine,
}

/// Tolerance values for geometric comparisons.
///
/// Six precision thresholds organised from tightest to loosest:
/// - `numerical_floor` (1e-14) -- near machine epsilon, for zero checks
/// - `solver` (1e-10) -- convergence criterion for iterative solvers
/// - `linear` (1e-6) -- physical distance tolerance (point coincidence)
/// - `angular` (1e-4) -- angle tolerance in radians
/// - `parametric_base` (1e-8) -- floor for parameter-space comparisons
/// - `approximation` (1e-5) -- display / tessellation adequacy
///
/// **Ordering invariant:** `numerical_floor < solver < linear`.
///
/// Fields are private. Use getter methods to read values, and builder
/// methods or factories to construct.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Tolerance {
    /// Near machine epsilon. Fixed at 1e-14 -- invariant under unit scaling.
    #[serde(default = "default_numerical_floor")]
    numerical_floor: f64,

    /// Solver convergence threshold. Default: 1e-10.
    #[serde(default = "default_solver")]
    solver: f64,

    /// Distance tolerance for point coincidence and edge length. Default: 1e-6.
    linear: f64,

    /// Angle tolerance in radians. Default: 1e-4.
    angular: f64,

    /// Base parametric tolerance (floor). Use `parametric_at()` for
    /// derivative-scaled comparisons. Default: 1e-8.
    #[serde(alias = "parametric")]
    parametric_base: f64,

    /// Display / tessellation adequacy threshold. Default: 1e-5.
    #[serde(default = "default_approximation")]
    approximation: f64,
}

fn default_numerical_floor() -> f64 {
    1e-14
}
fn default_solver() -> f64 {
    1e-10
}
fn default_approximation() -> f64 {
    1e-5
}

impl Default for Tolerance {
    fn default() -> Self {
        Self {
            numerical_floor: 1e-14,
            solver: 1e-10,
            linear: 1e-6,
            angular: 1e-4,
            parametric_base: 1e-8,
            approximation: 1e-5,
        }
    }
}

// ---------------------------------------------------------------------------
// Getters
// ---------------------------------------------------------------------------

impl Tolerance {
    /// Near-machine-epsilon floor for zero checks.
    pub fn numerical_floor(&self) -> f64 {
        self.numerical_floor
    }

    /// Solver convergence threshold.
    pub fn solver(&self) -> f64 {
        self.solver
    }

    /// Distance tolerance for point coincidence.
    pub fn linear(&self) -> f64 {
        self.linear
    }

    /// Angle tolerance in radians.
    pub fn angular(&self) -> f64 {
        self.angular
    }

    /// Base parametric tolerance (floor). Prefer `parametric_at()` when
    /// derivative magnitude is available.
    pub fn parametric_base(&self) -> f64 {
        self.parametric_base
    }

    /// Display / tessellation adequacy threshold.
    pub fn approximation(&self) -> f64 {
        self.approximation
    }
}

// ---------------------------------------------------------------------------
// Builders (backwards-compatible + new)
// ---------------------------------------------------------------------------

impl Tolerance {
    /// Set linear tolerance.
    pub fn with_linear(mut self, linear: f64) -> Self {
        self.linear = linear;
        self
    }

    /// Set angular tolerance.
    pub fn with_angular(mut self, angular: f64) -> Self {
        self.angular = angular;
        self
    }

    /// Set parametric base tolerance.
    ///
    /// This replaces the old `with_parametric` name but the behaviour is
    /// identical -- sets the `parametric_base` field.
    pub fn with_parametric(mut self, parametric: f64) -> Self {
        self.parametric_base = parametric;
        self
    }

    /// Set solver convergence threshold.
    pub fn with_solver(mut self, solver: f64) -> Self {
        self.solver = solver;
        self
    }

    /// Set numerical floor.
    pub fn with_numerical_floor(mut self, floor: f64) -> Self {
        self.numerical_floor = floor;
        self
    }

    /// Set approximation (display/tessellation) threshold.
    pub fn with_approximation(mut self, approx: f64) -> Self {
        self.approximation = approx;
        self
    }
}

// ---------------------------------------------------------------------------
// Comparison methods
// ---------------------------------------------------------------------------

impl Tolerance {
    /// Returns true if `v` is effectively zero (below numerical floor).
    pub fn nearly_zero(&self, v: f64) -> bool {
        v.abs() < self.numerical_floor
    }

    /// Returns true if a solver residual indicates convergence.
    pub fn solver_converged(&self, residual: f64) -> bool {
        residual < self.solver
    }

    /// Returns true if two points are within linear tolerance.
    pub fn points_coincident(&self, a: &Point3, b: &Point3) -> bool {
        a.distance_to(b) < self.linear
    }

    /// Returns true if two angles (in radians) are within angular tolerance.
    pub fn angles_equal(&self, a: f64, b: f64) -> bool {
        (a - b).abs() < self.angular
    }

    /// Returns true if two distances are within linear tolerance.
    pub fn distances_equal(&self, d1: f64, d2: f64) -> bool {
        (d1 - d2).abs() < self.linear
    }

    /// Returns true if two parameter values are within parametric tolerance,
    /// scaled by derivative magnitude.
    pub fn params_equal_at(&self, t1: f64, t2: f64, deriv_mag: f64) -> bool {
        (t1 - t2).abs() < self.parametric_at(deriv_mag)
    }

    /// Returns true if two parameter values are within the base parametric
    /// tolerance.
    #[deprecated(note = "Use params_equal_at() with derivative magnitude, \
        or compare against parametric_base() directly")]
    pub fn params_equal(&self, t1: f64, t2: f64) -> bool {
        (t1 - t2).abs() < self.parametric_base
    }

    /// Returns true if an approximation error is adequate for display.
    pub fn approx_adequate(&self, error: f64) -> bool {
        error < self.approximation
    }
}

// ---------------------------------------------------------------------------
// Derivative-scaled parametric tolerance
// ---------------------------------------------------------------------------

impl Tolerance {
    /// Parametric tolerance scaled by local derivative magnitude.
    ///
    /// Returns `(linear / derivative_mag).max(parametric_base)`. When the
    /// derivative magnitude is degenerate (below `numerical_floor`), returns
    /// `parametric_base` directly.
    pub fn parametric_at(&self, derivative_mag: f64) -> f64 {
        if derivative_mag < self.numerical_floor {
            return self.parametric_base;
        }
        (self.linear / derivative_mag).max(self.parametric_base)
    }
}

// ---------------------------------------------------------------------------
// Factories
// ---------------------------------------------------------------------------

impl Tolerance {
    /// Derive tolerance from model bounding-box diagonal.
    ///
    /// Linear tolerance is 1 ppm of the diagonal, with a floor of 1e-12 for
    /// tiny models. Other fields scale proportionally.
    pub fn from_model_scale(diagonal: f64) -> Self {
        let linear = (diagonal * 1e-6).max(1e-12);
        Self {
            numerical_floor: 1e-14,
            solver: linear * 1e-4,
            linear,
            angular: 1e-4,
            parametric_base: linear * 0.01,
            approximation: linear * 10.0,
        }
    }

    /// Create from a named preset, optionally scaled to model size.
    ///
    /// When `model_diagonal` is `Some`, the preset's linear tolerance is
    /// further scaled by the model's size relative to Standard (1 m diagonal).
    pub fn preset(preset: TolerancePreset, model_diagonal: Option<f64>) -> Self {
        let base = match preset {
            TolerancePreset::Draft => Self {
                numerical_floor: 1e-14,
                solver: 1e-8,
                linear: 1e-4,
                angular: 1e-4,
                parametric_base: 1e-8,
                approximation: 1e-3,
            },
            TolerancePreset::Standard => Self::default(),
            TolerancePreset::Fine => Self {
                numerical_floor: 1e-14,
                solver: 1e-12,
                linear: 1e-8,
                angular: 1e-4,
                parametric_base: 1e-8,
                approximation: 1e-7,
            },
        };
        match model_diagonal {
            Some(d) => {
                let scale_factor = d; // relative to 1m reference
                base.scale(scale_factor)
            }
            None => base,
        }
    }

    /// Scale distance-based precision values for unit conversion.
    ///
    /// Three fields are invariant under unit scaling:
    /// - `numerical_floor` -- property of IEEE 754, not model
    /// - `angular` -- radians are dimensionless
    /// - `parametric_base` -- parameter space is dimensionless
    pub fn scale(&self, factor: f64) -> Self {
        Self {
            numerical_floor: self.numerical_floor,
            solver: self.solver * factor,
            linear: self.linear * factor,
            angular: self.angular,
            parametric_base: self.parametric_base,
            approximation: self.approximation * factor,
        }
    }
}

#[cfg(test)]
#[path = "tolerance_tests.rs"]
mod tolerance_tests;
