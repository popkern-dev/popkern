//! Normal-curvature evaluation along an ambient tangent direction.
//!
//! Promoted from the commercial-gate test helpers so that library
//! code (support-cell shared-trim quality, ring continuity gates)
//! can measure G2 continuity through one shared implementation.
//!
//! The normal curvature at a surface point along a tangent
//! direction `t` is `κ_n = II(t, t) / I(t, t)` where `I` and `II`
//! are the first and second fundamental forms. Two faces joined
//! along a trim are G2 along that trim when their normal
//! curvatures agree for every shared tangent direction.

use thiserror::Error;

use crate::geometry::Vec3;
use crate::geometry::gregory::GregoryPatch;

/// Why a normal-curvature sample could not be produced.
///
/// Callers measuring G2 quality must treat any error as
/// "unmeasured" and fail closed (not as a passing sample).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum NormalCurvatureError {
    /// The patch variant does not expose second partial
    /// derivatives, so the second fundamental form cannot be
    /// evaluated.
    #[error("gregory patch variant {variant} has no second-partial evaluation")]
    UnsupportedVariant {
        /// Name of the unsupported [`GregoryPatch`] variant.
        variant: &'static str,
    },
    /// The first fundamental form is singular at the sample
    /// (`E·G − F² ≈ 0`), e.g. at a collapsed parameterization.
    #[error("first fundamental form is degenerate at the sample")]
    DegenerateFirstFundamentalForm,
    /// The requested ambient direction projects to a (near-)zero
    /// tangent vector, so no direction exists to measure along.
    #[error("ambient direction projects to a degenerate surface tangent")]
    DegenerateDirection,
    /// The surface normal is undefined at the sample
    /// (`∂u × ∂v ≈ 0`).
    #[error("surface normal is degenerate at the sample")]
    DegenerateNormal,
}

/// Decomposition of one normal-curvature measurement.
///
/// `projected_u` and `projected_v` are the surface-coordinate
/// weights of the requested ambient direction after tangent-plane
/// projection and unit normalization in the first fundamental
/// form. `second_form_l/m/n` are the second fundamental form
/// coefficients, so `curvature` is
/// `L·u² + 2M·u·v + N·v²`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalCurvatureComponents {
    /// Normal curvature along the projected unit tangent.
    pub curvature: f64,
    /// Unit tangent coefficient along `∂S/∂u`.
    pub projected_u: f64,
    /// Unit tangent coefficient along `∂S/∂v`.
    pub projected_v: f64,
    /// Patch normal used to sign the second fundamental form.
    pub normal: Vec3,
    /// Ambient vector `S_uu·u² + 2S_uv·u·v + S_vv·v²`.
    pub directional_second: Vec3,
    /// Second fundamental form coefficient `n · ∂²S/∂u²`.
    pub second_form_l: f64,
    /// Second fundamental form coefficient `n · ∂²S/∂u∂v`.
    pub second_form_m: f64,
    /// Second fundamental form coefficient `n · ∂²S/∂v²`.
    pub second_form_n: f64,
}

/// Normal curvature of `patch` at `(u, v)` along the ambient
/// direction `direction`.
///
/// `direction` is projected onto the tangent plane via the dual
/// basis of the first fundamental form, then normalized in the
/// surface metric, so only its tangential component matters and
/// its magnitude is irrelevant. The sign convention follows the
/// patch's geometric normal `∂u × ∂v` — callers comparing faces
/// of opposite orientation must flip the sign for the flipped
/// face.
///
/// `degeneracy_tolerance` bounds the parametric degeneracy
/// guards (first-form determinant and projected-direction norm
/// squared). It is a dimensionless parametric guard, not a model
/// tolerance; gate helpers use `1e-14`, which rejects only
/// samples where the result would be numerically meaningless.
pub fn normal_curvature_along(
    patch: &GregoryPatch,
    u: f64,
    v: f64,
    direction: Vec3,
    degeneracy_tolerance: f64,
) -> Result<f64, NormalCurvatureError> {
    normal_curvature_components_along(patch, u, v, direction, degeneracy_tolerance)
        .map(|components| components.curvature)
}

/// Normal curvature and second-fundamental-form decomposition of
/// `patch` at `(u, v)` along the ambient direction `direction`.
///
/// This uses the same projection, normalization, degeneracy
/// guards, and sign convention as [`normal_curvature_along`].
pub fn normal_curvature_components_along(
    patch: &GregoryPatch,
    u: f64,
    v: f64,
    direction: Vec3,
    degeneracy_tolerance: f64,
) -> Result<NormalCurvatureComponents, NormalCurvatureError> {
    let (suu, suv, svv) = second_partials(patch, u, v)?;
    let su = patch.partial_u(u, v);
    let sv = patch.partial_v(u, v);
    let normal = su
        .cross(&sv)
        .normalized()
        .ok_or(NormalCurvatureError::DegenerateNormal)?;
    let e = su.dot(&su);
    let f = su.dot(&sv);
    let g = sv.dot(&sv);
    let det = e * g - f * f;
    if det.abs() < degeneracy_tolerance {
        return Err(NormalCurvatureError::DegenerateFirstFundamentalForm);
    }
    let su_dot_t = su.dot(&direction);
    let sv_dot_t = sv.dot(&direction);
    let a = (g * su_dot_t - f * sv_dot_t) / det;
    let b = (e * sv_dot_t - f * su_dot_t) / det;
    let norm_sq = e * a * a + 2.0 * f * a * b + g * b * b;
    if norm_sq < degeneracy_tolerance {
        return Err(NormalCurvatureError::DegenerateDirection);
    }
    let a = a / norm_sq.sqrt();
    let b = b / norm_sq.sqrt();
    let l = normal.dot(&suu);
    let m = normal.dot(&suv);
    let n = normal.dot(&svv);
    let directional_second = suu * (a * a) + suv * (2.0 * a * b) + svv * (b * b);
    Ok(NormalCurvatureComponents {
        curvature: normal.dot(&directional_second),
        projected_u: a,
        projected_v: b,
        normal,
        directional_second,
        second_form_l: l,
        second_form_m: m,
        second_form_n: n,
    })
}

/// Second partial derivatives `(∂²S/∂u², ∂²S/∂u∂v, ∂²S/∂v²)` for
/// the patch variants that expose them.
fn second_partials(
    patch: &GregoryPatch,
    u: f64,
    v: f64,
) -> Result<(Vec3, Vec3, Vec3), NormalCurvatureError> {
    match patch {
        GregoryPatch::RectangularBiNonic(r) => Ok((
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        )),
        GregoryPatch::RectangularBiNonicDualTwist(r) => Ok((
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        )),
        GregoryPatch::RectangularTwoStripe(r) => Ok((
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        )),
        GregoryPatch::RectangularBiseptic(r) => Ok((
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        )),
        GregoryPatch::RectangularBiquinticFull(r) => Ok((
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        )),
        GregoryPatch::Rectangular(_) => Err(NormalCurvatureError::UnsupportedVariant {
            variant: "Rectangular",
        }),
        GregoryPatch::RectangularDeg5(_) => Err(NormalCurvatureError::UnsupportedVariant {
            variant: "RectangularDeg5",
        }),
        GregoryPatch::Triangular(_) => Err(NormalCurvatureError::UnsupportedVariant {
            variant: "Triangular",
        }),
    }
}

#[cfg(test)]
#[path = "normal_curvature_tests.rs"]
mod tests;
