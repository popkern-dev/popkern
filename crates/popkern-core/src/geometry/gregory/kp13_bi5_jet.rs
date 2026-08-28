//! Taylor jets at a corner, plus chain rule for
//! composition.
//!
//! KP13 §3 eq (4) minimises `Σ_k Σ_{t ∈ {m_k, v_k}} ‖H_t(g△^k ∘ χ^k)
//! − H_t(b^k)‖²`. The `H_t` jet of a composite is most cleanly
//! expressed in **Taylor form** (value + first + second partials at
//! `t`), then converted to whatever BB representation the LSQ
//! residual uses.
//!
//! This module provides:
//!
//! * [`TaylorJet2`] — a 6-component jet (value + `∂u`, `∂v`, `∂uu`,
//!   `∂uv`, `∂vv`) of a scalar function of two variables.
//! * [`compose_chain_rule`] — given the outer function's 2-D Hessian
//!   data (evaluated at `χ(t)`) and the inner map's jet (`ChiJet`),
//!   computes the composite jet at `t` via the second-order chain
//!   rule.
//!
//! # Scope
//!
//! This handles the *geometry* side — the composite `(g△^k ∘ χ^k)(t)`
//! at a single corner — for a fixed `g△^k`. The LSQ (I.11e) will
//! call this for each free BB-coef of `g△^k` with that coef set to
//! 1 (and others to 0), giving the per-coef contribution to the
//! jet. That yields the linear map from `g△^k` coefs to jet
//! components needed by the LSQ design matrix.

#![allow(dead_code)]

use super::bernstein_tri5::{
    TRI5_CP_COUNT, bernstein_tri5_all_extended, bernstein_tri5_partial_u_extended,
    bernstein_tri5_partial_v_extended, bernstein_tri5_second_partial_uu_extended,
    bernstein_tri5_second_partial_uv_extended, bernstein_tri5_second_partial_vv_extended,
};
use super::kp13_bi5_chi::ChiJet;
use super::kp13_bi5_gtri5::Gtri5Sector;
use crate::geometry::Point3;

/// Second-order Taylor jet of a scalar-triple-valued function `f`
/// at a point. All fields evaluated at the same point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TaylorJet2 {
    /// `f(t)`.
    pub value: Point3,
    /// `∂f / ∂u` at `t`.
    pub d_u: Point3,
    /// `∂f / ∂v`.
    pub d_v: Point3,
    /// `∂²f / ∂u²`.
    pub d_uu: Point3,
    /// `∂²f / ∂u∂v`.
    pub d_uv: Point3,
    /// `∂²f / ∂v²`.
    pub d_vv: Point3,
}

impl TaylorJet2 {
    pub fn zero() -> Self {
        let z = Point3::new(0.0, 0.0, 0.0);
        Self {
            value: z,
            d_u: z,
            d_v: z,
            d_uu: z,
            d_uv: z,
            d_vv: z,
        }
    }

    /// Component-wise difference — useful for LSQ residuals.
    pub fn sub(&self, rhs: &Self) -> Self {
        Self {
            value: sub(self.value, rhs.value),
            d_u: sub(self.d_u, rhs.d_u),
            d_v: sub(self.d_v, rhs.d_v),
            d_uu: sub(self.d_uu, rhs.d_uu),
            d_uv: sub(self.d_uv, rhs.d_uv),
            d_vv: sub(self.d_vv, rhs.d_vv),
        }
    }

    /// `‖.‖²` as the sum of squared Euclidean distances of each
    /// component.
    pub fn sq_norm(&self) -> f64 {
        sq(self.value) + sq(self.d_u) + sq(self.d_v) + sq(self.d_uu) + sq(self.d_uv) + sq(self.d_vv)
    }
}

/// Second-order Hessian data of the outer function `f` at a point
/// in its domain — the operands for the chain rule that produces
/// `(f ∘ χ)`'s jet.
///
/// `f` here is the one whose gradient/Hessian lives in the
/// `χ`-target domain (Ω), i.e., `f(x, y)`. In the KP13 setting,
/// `f := g△^k : Ω → R³`.
#[derive(Debug, Clone, Copy)]
pub struct OuterHessian {
    /// `f(x, y)` at `χ(t)`.
    pub value: Point3,
    /// `∂f / ∂x` at `χ(t)`.
    pub d_x: Point3,
    /// `∂f / ∂y` at `χ(t)`.
    pub d_y: Point3,
    /// `∂²f / ∂x²`.
    pub d_xx: Point3,
    /// `∂²f / ∂x∂y`.
    pub d_xy: Point3,
    /// `∂²f / ∂y²`.
    pub d_yy: Point3,
}

/// Second-order chain rule `(f ∘ χ)` at a corner of `χ`'s domain.
///
/// Denote `x = χ_x(u, v)`, `y = χ_y(u, v)`. Then:
///
/// ```text
///   (f ∘ χ)(t)    = f(x, y)
///   ∂u(f ∘ χ)    = f_x · x_u + f_y · y_u
///   ∂v(f ∘ χ)    = f_x · x_v + f_y · y_v
///   ∂uu(f ∘ χ)   = f_xx · x_u² + 2·f_xy · x_u·y_u + f_yy · y_u²
///                + f_x · x_uu + f_y · y_uu
///   ∂uv(f ∘ χ)   = f_xx · x_u·x_v + f_xy · (x_u·y_v + y_u·x_v)
///                + f_yy · y_u·y_v + f_x · x_uv + f_y · y_uv
///   ∂vv(f ∘ χ)   = f_xx · x_v² + 2·f_xy · x_v·y_v + f_yy · y_v²
///                + f_x · x_vv + f_y · y_vv
/// ```
pub fn compose_chain_rule(outer: &OuterHessian, chi: &ChiJet) -> TaylorJet2 {
    let (xu, yu) = chi.d_u;
    let (xv, yv) = chi.d_v;
    let (xuu, yuu) = chi.d_uu;
    let (xuv, yuv) = chi.d_uv;
    let (xvv, yvv) = chi.d_vv;

    let d_u = linear2(outer.d_x, xu, outer.d_y, yu);
    let d_v = linear2(outer.d_x, xv, outer.d_y, yv);

    let d_uu = {
        let quad = linear3(
            outer.d_xx,
            xu * xu,
            outer.d_xy,
            2.0 * xu * yu,
            outer.d_yy,
            yu * yu,
        );
        let linpart = linear2(outer.d_x, xuu, outer.d_y, yuu);
        add(quad, linpart)
    };

    let d_uv = {
        let quad = linear3(
            outer.d_xx,
            xu * xv,
            outer.d_xy,
            xu * yv + yu * xv,
            outer.d_yy,
            yu * yv,
        );
        let linpart = linear2(outer.d_x, xuv, outer.d_y, yuv);
        add(quad, linpart)
    };

    let d_vv = {
        let quad = linear3(
            outer.d_xx,
            xv * xv,
            outer.d_xy,
            2.0 * xv * yv,
            outer.d_yy,
            yv * yv,
        );
        let linpart = linear2(outer.d_x, xvv, outer.d_y, yvv);
        add(quad, linpart)
    };

    TaylorJet2 {
        value: outer.value,
        d_u,
        d_v,
        d_uu,
        d_uv,
        d_vv,
    }
}

/// Evaluate the tri-5 sector polynomial at natural-simplex
/// coordinates `(u, v)` and package all values needed for the
/// chain rule into an [`OuterHessian`].
///
/// Uses `*_extended` variants so evaluation works for `(u, v)`
/// outside the natural simplex (i.e., `u + v > 1`) — needed for
/// the `v_k` corner of the χ-domain.
///
/// # Domain reconciliation
///
/// In the paper, `g△^k` is a function on Ω (a 2-D ambient region).
/// Its tri-5 BB representation parameterises the wedge △_k via
/// barycentric `(u, v, w = 1 − u − v)`. When we compose with `χ`
/// using the bilinear approximation `χ^k(u_χ, v_χ) := u_χ · m_k +
/// v_χ · m_{k+1}`, the natural-simplex `(u, v)` **equals** the
/// χ-domain `(u_χ, v_χ)` — same parameter. So this helper takes
/// `(u, v)` directly without any coordinate transform.
///
/// For a full Catmull-Clark `χ`, callers would first invert
/// `χ^k` to get simplex `(u, v)` from a χ-domain point. That's
/// future work.
pub fn outer_hessian_at_natural(sector: &Gtri5Sector, u: f64, v: f64) -> OuterHessian {
    let basis = bernstein_tri5_all_extended(u, v);
    let partial_u = bernstein_tri5_partial_u_extended(u, v);
    let partial_v = bernstein_tri5_partial_v_extended(u, v);
    let partial_uu = bernstein_tri5_second_partial_uu_extended(u, v);
    let partial_uv = bernstein_tri5_second_partial_uv_extended(u, v);
    let partial_vv = bernstein_tri5_second_partial_vv_extended(u, v);
    OuterHessian {
        value: combine(&basis, &sector.cps),
        d_x: combine(&partial_u, &sector.cps),
        d_y: combine(&partial_v, &sector.cps),
        d_xx: combine(&partial_uu, &sector.cps),
        d_xy: combine(&partial_uv, &sector.cps),
        d_yy: combine(&partial_vv, &sector.cps),
    }
}

#[inline]
fn combine(basis: &[f64; TRI5_CP_COUNT], cps: &[Point3; TRI5_CP_COUNT]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..TRI5_CP_COUNT {
        out.x += basis[i] * cps[i].x;
        out.y += basis[i] * cps[i].y;
        out.z += basis[i] * cps[i].z;
    }
    out
}

// ---- scalar helpers ----

#[inline]
fn sub(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

#[inline]
fn add(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

#[inline]
fn linear2(a: Point3, sa: f64, b: Point3, sb: f64) -> Point3 {
    Point3::new(
        sa * a.x + sb * b.x,
        sa * a.y + sb * b.y,
        sa * a.z + sb * b.z,
    )
}

#[inline]
fn linear3(a: Point3, sa: f64, b: Point3, sb: f64, c: Point3, sc: f64) -> Point3 {
    Point3::new(
        sa * a.x + sb * b.x + sc * c.x,
        sa * a.y + sb * b.y + sc * c.y,
        sa * a.z + sb * b.z + sc * c.z,
    )
}

#[inline]
fn sq(p: Point3) -> f64 {
    p.x * p.x + p.y * p.y + p.z * p.z
}

#[cfg(test)]
#[path = "kp13_bi5_jet_tests.rs"]
mod tests;
