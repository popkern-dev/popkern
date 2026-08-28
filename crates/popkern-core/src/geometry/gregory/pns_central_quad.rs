//! PnS cap — Central quadratic expansion via `a ∘ σ` chain rule.
//!
//! The cap's 6 order-`≤ 2` BB coefs per sector (`p^k_{ij}` with
//! `i + j ≤ 2`) are determined by composing an **undetermined
//! quadratic map** `a: R² → R³` with the sector's σ
//! reparameterization (see `pns_sigma`). Paper `14bi234.pdf` §4:
//!
//! > "The partial derivatives `∂^i_s ∂^j_t (a ∘ σ)`, `i + j ≤ 2`
//! > are then computed for a still undetermined quadratic map `a`
//! > (for this the composition involves only `σ_{ij}`, `i + j ≤ 2`).
//! > Then the quadratic expansion in terms of the partial
//! > derivatives is expressed in the BB-coefficients `p^k_{ij}`,
//! > `i + j ≤ 2` (circles in Fig. 8a) and this satisfies (2)."
//!
//! This module implements that composition explicitly. The
//! quadratic `a` is shared across all N sectors (gives a consistent
//! `G²` expansion at the pole); each sector's contribution uses a
//! rotated σ jet.
//!
//! # Convention for `a`
//!
//! ```text
//!   a(x, y) = a_0 + a_x·x + a_y·y + (1/2)·a_xx·x²
//!                                  + a_xy·x·y
//!                                  + (1/2)·a_yy·y²
//! ```
//!
//! So `∂²a/∂x² = a_xx`, `∂²a/∂y² = a_yy`, `∂²a/∂x∂y = a_xy`.
//! `a_0` is pinned to the **pole** (central point from eq 5).
//! The remaining `(a_x, a_y, a_xx, a_xy, a_yy)` are 5 `Point3`
//! quantities = 15 free scalars, determined jointly with the cap's
//! other free DoFs by the F_κ minimisation.
//!
//! # Sector-local σ via rotation
//!
//! Each sector's σ is a rotation of the canonical (sector-0) σ:
//! `σ^k(u, v) = R_{2π·k/n} · σ(u, v)`. In jet form this means each
//! R²-valued partial derivative is rotated by `2π·k/n` in image
//! space. [`rotate_sigma_jet`] handles that.

#![allow(dead_code)]

use super::pns_sigma::SigmaCenterJet;
use crate::geometry::Point3;

/// Jet of the undetermined quadratic map `a: R² → R³` used to build
/// the cap's central quadratic expansion.
#[derive(Debug, Clone, Copy)]
pub struct QuadMapJet {
    /// `a(0, 0)` — pinned to the pole (central point from eq 5).
    pub a_0: Point3,
    /// `∂a/∂x` at `(0, 0)`.
    pub a_x: Point3,
    /// `∂a/∂y` at `(0, 0)`.
    pub a_y: Point3,
    /// `∂²a/∂x²` at `(0, 0)`.
    pub a_xx: Point3,
    /// `∂²a/∂x∂y` at `(0, 0)`.
    pub a_xy: Point3,
    /// `∂²a/∂y²` at `(0, 0)`.
    pub a_yy: Point3,
}

/// Six order-`≤ 2` BB coefs of one bi-4 sector patch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CentralSixBB {
    /// `p_{0, 0}` — shared pole.
    pub p00: Point3,
    /// `p_{1, 0}`.
    pub p10: Point3,
    /// `p_{0, 1}`.
    pub p01: Point3,
    /// `p_{2, 0}`.
    pub p20: Point3,
    /// `p_{1, 1}`.
    pub p11: Point3,
    /// `p_{0, 2}`.
    pub p02: Point3,
}

/// Rotate σ's image-space R² jet by angle `theta` (radians).
///
/// Used to map sector-0's canonical σ jet to any sector `k`'s jet
/// via `θ = 2π·k/n`.
pub fn rotate_sigma_jet(jet: &SigmaCenterJet, theta: f64) -> SigmaCenterJet {
    let (sin_t, cos_t) = theta.sin_cos();
    let rot =
        |v: [f64; 2]| -> [f64; 2] { [cos_t * v[0] - sin_t * v[1], sin_t * v[0] + cos_t * v[1]] };
    SigmaCenterJet {
        s00: rot(jet.s00),
        ds_u: rot(jet.ds_u),
        ds_v: rot(jet.ds_v),
        ds_uu: rot(jet.ds_uu),
        ds_vv: rot(jet.ds_vv),
        ds_uv: rot(jet.ds_uv),
    }
}

/// Chain-rule composition: compute the 6 central BB coefs
/// `p_{ij}` (`i + j ≤ 2`) of one bi-4 sector patch from the
/// quadratic `a` and the (possibly sector-rotated) σ jet.
///
/// # Derivation
///
/// For `(a ∘ σ)` with `σ(0, 0) = (0, 0)`, the order-`≤ 2` Taylor
/// coefficients at `(0, 0)` are:
/// ```text
///   (a∘σ)(0, 0) = a_0
///   ∂_u (a∘σ)   = a_x · s^x_u + a_y · s^y_u
///   ∂_v (a∘σ)   = a_x · s^x_v + a_y · s^y_v
///   ∂²_u (a∘σ)  = a_xx·(s^x_u)² + 2·a_xy·s^x_u·s^y_u + a_yy·(s^y_u)²
///                + a_x · s^x_uu + a_y · s^y_uu
///   ∂²_v (a∘σ)  = a_xx·(s^x_v)² + 2·a_xy·s^x_v·s^y_v + a_yy·(s^y_v)²
///                + a_x · s^x_vv + a_y · s^y_vv
///   ∂²_{uv}(a∘σ)= a_xx·s^x_u·s^x_v + a_xy·(s^x_u·s^y_v + s^x_v·s^y_u)
///                + a_yy·s^y_u·s^y_v + a_x·s^x_uv + a_y·s^y_uv
/// ```
///
/// These Taylor coefficients convert to bi-4 BB coefs via the
/// standard formulas:
/// ```text
///   p_{00} = (a∘σ)(0, 0)
///   p_{10} = p_{00} + (1/4) ∂_u (a∘σ)
///   p_{01} = p_{00} + (1/4) ∂_v (a∘σ)
///   p_{20} = p_{00} + (1/2) ∂_u (a∘σ) + (1/12) ∂²_u (a∘σ)
///   p_{02} = p_{00} + (1/2) ∂_v (a∘σ) + (1/12) ∂²_v (a∘σ)
///   p_{11} = p_{00} + (1/4) ∂_u (a∘σ) + (1/4) ∂_v (a∘σ)
///          + (1/16) ∂²_{uv} (a∘σ)
/// ```
pub fn central_bb_from_a_and_sigma(a: &QuadMapJet, sigma: &SigmaCenterJet) -> CentralSixBB {
    // σ jet values (R² vectors).
    let sx_u = sigma.ds_u[0];
    let sy_u = sigma.ds_u[1];
    let sx_v = sigma.ds_v[0];
    let sy_v = sigma.ds_v[1];
    let sx_uu = sigma.ds_uu[0];
    let sy_uu = sigma.ds_uu[1];
    let sx_vv = sigma.ds_vv[0];
    let sy_vv = sigma.ds_vv[1];
    let sx_uv = sigma.ds_uv[0];
    let sy_uv = sigma.ds_uv[1];

    // (a∘σ)'s Taylor coefs at (0, 0).
    let t_u = add3(scale3(sx_u, a.a_x), scale3(sy_u, a.a_y));
    let t_v = add3(scale3(sx_v, a.a_x), scale3(sy_v, a.a_y));

    let t_uu = sum5(
        scale3(sx_u * sx_u, a.a_xx),
        scale3(2.0 * sx_u * sy_u, a.a_xy),
        scale3(sy_u * sy_u, a.a_yy),
        scale3(sx_uu, a.a_x),
        scale3(sy_uu, a.a_y),
    );
    let t_vv = sum5(
        scale3(sx_v * sx_v, a.a_xx),
        scale3(2.0 * sx_v * sy_v, a.a_xy),
        scale3(sy_v * sy_v, a.a_yy),
        scale3(sx_vv, a.a_x),
        scale3(sy_vv, a.a_y),
    );
    let t_uv = sum5(
        scale3(sx_u * sx_v, a.a_xx),
        scale3(sx_u * sy_v + sx_v * sy_u, a.a_xy),
        scale3(sy_u * sy_v, a.a_yy),
        scale3(sx_uv, a.a_x),
        scale3(sy_uv, a.a_y),
    );

    // Convert to bi-4 BB coefs.
    let p00 = a.a_0;
    let p10 = add3(p00, scale3(0.25, t_u));
    let p01 = add3(p00, scale3(0.25, t_v));
    let p20 = sum3(p00, scale3(0.5, t_u), scale3(1.0 / 12.0, t_uu));
    let p02 = sum3(p00, scale3(0.5, t_v), scale3(1.0 / 12.0, t_vv));
    let p11 = sum4(
        p00,
        scale3(0.25, t_u),
        scale3(0.25, t_v),
        scale3(1.0 / 16.0, t_uv),
    );

    CentralSixBB {
        p00,
        p10,
        p01,
        p20,
        p11,
        p02,
    }
}

/// Convenience: build all N sector central BB sextuples from a
/// single quadratic `a` and canonical σ jet by rotating σ per sector.
pub fn central_bb_all_sectors(
    a: &QuadMapJet,
    sigma_canonical: &SigmaCenterJet,
    n: usize,
) -> Vec<CentralSixBB> {
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let sigma_k = rotate_sigma_jet(sigma_canonical, theta);
        out.push(central_bb_from_a_and_sigma(a, &sigma_k));
    }
    out
}

// ── Small Point3 arithmetic helpers ──────────────────────────────

#[inline]
fn scale3(s: f64, p: Point3) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

#[inline]
fn add3(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

#[inline]
fn sum3(a: Point3, b: Point3, c: Point3) -> Point3 {
    Point3::new(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z)
}

#[inline]
fn sum4(a: Point3, b: Point3, c: Point3, d: Point3) -> Point3 {
    Point3::new(
        a.x + b.x + c.x + d.x,
        a.y + b.y + c.y + d.y,
        a.z + b.z + c.z + d.z,
    )
}

#[inline]
fn sum5(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Point3 {
    Point3::new(
        a.x + b.x + c.x + d.x + e.x,
        a.y + b.y + c.y + d.y + e.y,
        a.z + b.z + c.z + d.z + e.z,
    )
}

#[cfg(test)]
#[path = "pns_central_quad_tests.rs"]
mod tests;
