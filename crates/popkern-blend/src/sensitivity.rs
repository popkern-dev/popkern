//! §14.8 quality sensitivity analysis.
//!
//! The spec says:
//!
//! > Feasible range [r_min, r_max], full curvature analysis,
//! > sensitivity **d(Q)/d(r)** and d(Q)/d(ρ) for gradient-based
//! > optimization, highlight-line preview for specified light
//! > direction.
//!
//! This module computes the radius sensitivity **d(Q)/d(r)** by
//! finite differences: run the fillet at two nearby radii, compute
//! [`crate::precheck::BlendQuality`] for each result, and return the
//! forward-difference quotient.
//!
//! **Scope.** Radius sensitivity is the only variational parameter
//! the rolling-ball engine exposes today. The ρ parameter from
//! §2.2 conic profiles is not implemented yet (see
//! `fillets-vs-requirements.md`); when it lands, this module will
//! grow a companion `rho_sensitivity` function.
//!
//! **Cost.** Each call runs the fillet pipeline **twice** (at `r`
//! and `r + ε`), so it is materially more expensive than a normal
//! fillet call. Cache the result when presenting a "quality vs
//! radius" slider to the user.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::BlendError;
use crate::precheck::{self, BlendQuality, ComputeQualityInput};
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{fillet_edges, geometry_compute_feasibility};

/// Report from [`quality_sensitivity`].
///
/// Contains the two underlying [`BlendQuality`] samples the
/// derivative was computed from (useful for logging / UI display)
/// and three derivative views: composite, `r_min`, and surface area.
#[derive(Debug, Clone, Copy)]
pub struct SensitivityReport {
    /// The radius `r` the baseline quality was computed at.
    pub radius_base: f64,
    /// The step size `ε` used in the forward-difference estimate.
    pub epsilon: f64,
    /// Baseline quality at `r`.
    pub quality_base: BlendQuality,
    /// Quality at `r + ε`.
    pub quality_step: BlendQuality,
    /// Derivative of the composite quality score: `(Q(r+ε) − Q(r)) / ε`.
    /// Units: 1/length. Positive → quality improves as radius grows.
    pub d_composite_d_r: f64,
    /// Derivative of the minimum radius-of-curvature metric.
    /// A healthy fillet has `d(R_min)/dr ≈ 1` (the blend's minimum
    /// curvature radius scales linearly with the requested radius).
    /// Values far from 1 flag constrained geometry where growing
    /// the radius doesn't translate proportionally into surface
    /// curvature room.
    pub d_r_min_d_r: f64,
    /// Derivative of the blend surface area. Positive on a healthy
    /// convex fillet (bigger radius = more surface); near zero or
    /// negative flags the radius hitting a size constraint.
    pub d_area_d_r: f64,
}

impl SensitivityReport {
    /// True if the composite score responds positively to radius
    /// increase — the "room to grow" heuristic for a UI slider.
    pub fn grows_with_radius(&self) -> bool {
        self.d_composite_d_r > 0.0
    }

    /// Radius value where the linear extrapolation of the composite
    /// score from the base point would hit 1.0 (ideal quality), or
    /// `None` if the derivative is non-positive. Rough "which
    /// radius would be perfect" hint for a UI suggestion — not a
    /// guarantee, because the real quality curve is non-linear.
    pub fn radius_for_perfect_score(&self) -> Option<f64> {
        if self.d_composite_d_r <= 0.0 {
            return None;
        }
        let headroom = 1.0 - self.quality_base.composite_score;
        if headroom <= 0.0 {
            return Some(self.radius_base);
        }
        Some(self.radius_base + headroom / self.d_composite_d_r)
    }
}

/// Compute the forward-difference quality sensitivity `d(Q)/d(r)`
/// at a given radius.
///
/// Algorithm:
/// 1. Run [`fillet_edges`] at `r = spec.radius`.
/// 2. Run [`fillet_edges`] at `r + epsilon`.
/// 3. Compute [`BlendQuality`] for each result using
///    [`precheck::fillet_quality`] against the input BRep's pre-face
///    set as the baseline.
/// 4. Return `(Q(r+ε) − Q(r)) / ε` for composite, `R_min`, and
///    surface area, packaged as a [`SensitivityReport`].
///
/// Both fillet operations must succeed. If either one fails the
/// underlying [`BlendError`] is bubbled up — the caller typically
/// wraps this in a "try to reduce `epsilon`" retry loop.
///
/// **Choosing `epsilon`.** A sensible default is `radius * 1e-3`:
/// small enough that the finite-difference error is negligible but
/// large enough to escape numerical noise in the quality pipeline.
/// For radii comparable to the tolerance floor, fall back to a
/// fixed absolute step (e.g. `tolerance.linear() * 10`).
///
/// **Scope limitation.** The spec's companion `d(Q)/d(ρ)` for conic
/// profiles (§2.2) is not implemented — the rolling-ball engine has
/// no ρ parameter. This function only produces the `dr` derivative.
pub fn quality_sensitivity(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    epsilon: f64,
) -> Result<SensitivityReport, BlendError> {
    if specs.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "specs".into(),
            value: "<empty>".into(),
            reason: "sensitivity analysis requires at least one fillet spec".into(),
        });
    }
    if !epsilon.is_finite() || epsilon <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "epsilon".into(),
            value: format!("{epsilon}"),
            reason: "epsilon must be a positive finite step size".into(),
        });
    }

    let radius_base = specs
        .iter()
        .map(|s| geometry_compute_feasibility::largest_requested_radius(&s.radius))
        .fold(0.0_f64, f64::max);

    // Build the step spec set — every spec's constant / variable
    // profile is scaled up by `epsilon` in absolute terms.
    let step_specs: Vec<FilletSpec> = specs.iter().map(|s| shifted_spec(s, epsilon)).collect();

    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    let base_result = fillet_edges(brep, specs, options, tolerance, None)?;
    let step_result = fillet_edges(brep, &step_specs, options, tolerance, None)?;

    let quality_base = precheck::fillet_quality(&ComputeQualityInput {
        post_brep: &base_result.brep,
        pre_faces: &pre_faces,
        requested_radius: radius_base,
    });
    let quality_step = precheck::fillet_quality(&ComputeQualityInput {
        post_brep: &step_result.brep,
        pre_faces: &pre_faces,
        requested_radius: radius_base + epsilon,
    });

    let d_composite_d_r = (quality_step.composite_score - quality_base.composite_score) / epsilon;
    let d_r_min_d_r = (quality_step.r_min - quality_base.r_min) / epsilon;
    let d_area_d_r = (quality_step.surface_area - quality_base.surface_area) / epsilon;

    Ok(SensitivityReport {
        radius_base,
        epsilon,
        quality_base,
        quality_step,
        d_composite_d_r,
        d_r_min_d_r,
        d_area_d_r,
    })
}

/// Return a new [`FilletSpec`] with every radius value inside the
/// profile shifted by `epsilon`.
fn shifted_spec(s: &FilletSpec, epsilon: f64) -> FilletSpec {
    let radius = match &s.radius {
        FilletRadius::Constant(r) => FilletRadius::Constant(r + epsilon),
        FilletRadius::Variable { start, end } => FilletRadius::Variable {
            start: start + epsilon,
            end: end + epsilon,
        },
        FilletRadius::ByPoints(points) => {
            FilletRadius::ByPoints(points.iter().map(|(t, r)| (*t, r + epsilon)).collect())
        }
    };
    FilletSpec {
        edge: s.edge,
        radius,
    }
}

#[cfg(test)]
#[path = "sensitivity_tests.rs"]
mod tests;
