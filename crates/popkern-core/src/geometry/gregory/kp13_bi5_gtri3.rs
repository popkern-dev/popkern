//! KP13 § 5 eq (19), (20) — tri-3 `C²` guide `g̃^△` and its star-
//! jet extractor, used by I.14.e to set `H_o p^{o,k}`.
//!
//! Paper § 5 eq (19) defines the free-param set
//!
//! ```text
//!   G_3 := {g̃^1_{ijm} | i > 0} ∪ {g̃^k_{030}}_{k=0..n-1}
//! ```
//!
//! for a `C²` piecewise polynomial `g̃^△` of total degree 3.
//! Six CPs from sector 1 (those with first paper-index `a > 0` in
//! the `(a, b, c)` BB labelling, `a + b + c = 3`) plus one `m_k`
//! mid-edge CP per sector = `n + 6` scalars total.
//!
//! Paper eq (20) states `H_o p^{o,k} := H_o(g̃^△^k ∘ σ^o)`, so the
//! `O` piece's star jet is the tri-3 star jet of `g̃^△^k` after
//! composition with `σ^o`.
//!
//! # Tri-3 star jet — closed form
//!
//! For a tri-3 patch with BB coefs at `(a, b, c)` with `a+b+c=3`,
//! the 2nd-order Taylor jet at the star corner `(u, v) = (0, 0)`
//! (simplex barycentric) depends on exactly six CPs — the same
//! six CPs that are `a > 0` in the paper label, i.e., sector 1's
//! `G_3` contribution:
//!
//! ```text
//!   value  = cps[(3, 0, 0)]                                    = c_300
//!   ∂_u    = 3 · (cps[(2, 1, 0)] − cps[(3, 0, 0)])              = 3(c_210 − c_300)
//!   ∂_v    = 3 · (cps[(2, 0, 1)] − cps[(3, 0, 0)])              = 3(c_201 − c_300)
//!   ∂²_u   = 6 · (cps[(1, 2, 0)] − 2·cps[(2, 1, 0)] + cps[(3, 0, 0)])
//!   ∂²_v   = 6 · (cps[(1, 0, 2)] − 2·cps[(2, 0, 1)] + cps[(3, 0, 0)])
//!   ∂_u∂_v = 6 · (cps[(1, 1, 1)] − cps[(2, 1, 0)] − cps[(2, 0, 1)] + cps[(3, 0, 0)])
//! ```
//!
//! # Inter-sector `C²` at the star via 2-D rotation
//!
//! At the extraordinary vertex, adjacent sectors meet with their
//! 2-D local frames rotated by `θ_k = 2π k / n`. `C²` of
//! `g̃^△` at the star means the 3-D **star jet** is invariant in a
//! suitable rotational sense: value is shared, and the 1st- and
//! 2nd-order partials transform as a 2-D vector and a symmetric
//! 2-D tensor respectively under the sector-to-sector rotation.
//!
//! So sector `k`'s star jet = rotate(sector 0's star jet, θ_k),
//! with the rotation applied to the 2-D derivative *domain*:
//!
//! ```text
//!   ∂_u' = cos θ · ∂_u + sin θ · ∂_v
//!   ∂_v' = -sin θ · ∂_u + cos θ · ∂_v
//!   ∂²_u' = c² ∂²_u + 2cs ∂_u∂_v + s² ∂²_v
//!   ∂²_v' = s² ∂²_u − 2cs ∂_u∂_v + c² ∂²_v
//!   ∂_u∂_v' = cs (∂²_u − ∂²_v) + (c² − s²) ∂_u∂_v
//! ```
//!
//! where `c = cos θ`, `s = sin θ`.
//!
//! # Scope
//!
//! This module provides the narrow interface needed by I.14.e:
//! - [`Gtri3SixCps`] — sector 0's six "a > 0" CPs of `g̃^△`.
//! - [`tri3_star_jet`] — analytical star-jet extraction.
//! - [`rotate_star_jet_to_sector`] — 2-D rotation for sector `k`.
//!
//! The full `g̃^△` struct with all `10n` CPs per sector + inter-
//! sector `C²` propagation is a paper-faithful follow-up that
//! would thread through the mid-edge `g̃^k_{030}` coefs too.

#![allow(dead_code)]

use popkern_math::PortableMath;
use std::f64::consts::PI;

use super::kp13_bi5_jet::TaylorJet2;
use crate::geometry::Point3;

/// Sector 0's six `a > 0` CPs of `g̃^△`, laid out in the canonical
/// order used by this module:
///
/// - `cps[0]` = paper `(3, 0, 0)` — **star** (shared).
/// - `cps[1]` = paper `(2, 0, 1)` — toward `m_{k+1}` axis.
/// - `cps[2]` = paper `(1, 0, 2)` — 2-ring on `m_{k+1}` axis.
/// - `cps[3]` = paper `(2, 1, 0)` — toward `m_k` axis.
/// - `cps[4]` = paper `(1, 1, 1)` — star-interior.
/// - `cps[5]` = paper `(1, 2, 0)` — 2-ring on `m_k` axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gtri3SixCps {
    pub cps: [Point3; 6],
}

impl Gtri3SixCps {
    pub fn zero() -> Self {
        Self {
            cps: [Point3::new(0.0, 0.0, 0.0); 6],
        }
    }
}

/// Star jet of a tri-3 patch from its six `a > 0` BB coefs. Uses
/// the standard tri-3 partial-derivative formulas at the simplex
/// `(u, v) = (0, 0)` corner (`w = 1`).
pub fn tri3_star_jet(s: &Gtri3SixCps) -> TaylorJet2 {
    let c300 = s.cps[0];
    let c201 = s.cps[1];
    let c102 = s.cps[2];
    let c210 = s.cps[3];
    let c111 = s.cps[4];
    let c120 = s.cps[5];

    TaylorJet2 {
        value: c300,
        d_u: scale(sub(c210, c300), 3.0),
        d_v: scale(sub(c201, c300), 3.0),
        d_uu: scale(
            Point3::new(
                c120.x - 2.0 * c210.x + c300.x,
                c120.y - 2.0 * c210.y + c300.y,
                c120.z - 2.0 * c210.z + c300.z,
            ),
            6.0,
        ),
        d_vv: scale(
            Point3::new(
                c102.x - 2.0 * c201.x + c300.x,
                c102.y - 2.0 * c201.y + c300.y,
                c102.z - 2.0 * c201.z + c300.z,
            ),
            6.0,
        ),
        d_uv: scale(
            Point3::new(
                c111.x - c210.x - c201.x + c300.x,
                c111.y - c210.y - c201.y + c300.y,
                c111.z - c210.z - c201.z + c300.z,
            ),
            6.0,
        ),
    }
}

/// Rotate a star jet by `θ = 2π k / n` to obtain sector `k`'s
/// star jet from sector 0's. The 1st-order partials rotate as a
/// 2-D vector; the 2nd-order partials rotate as a symmetric 2-D
/// tensor. The scalar value is invariant.
pub fn rotate_star_jet_to_sector(jet: &TaylorJet2, k: usize, n: usize) -> TaylorJet2 {
    let theta = 2.0 * PI * k as f64 / n as f64;
    let c = theta.pcos();
    let s = theta.psin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;

    let d_u_rot = Point3::new(
        c * jet.d_u.x + s * jet.d_v.x,
        c * jet.d_u.y + s * jet.d_v.y,
        c * jet.d_u.z + s * jet.d_v.z,
    );
    let d_v_rot = Point3::new(
        -s * jet.d_u.x + c * jet.d_v.x,
        -s * jet.d_u.y + c * jet.d_v.y,
        -s * jet.d_u.z + c * jet.d_v.z,
    );
    let d_uu_rot = Point3::new(
        c2 * jet.d_uu.x + 2.0 * cs * jet.d_uv.x + s2 * jet.d_vv.x,
        c2 * jet.d_uu.y + 2.0 * cs * jet.d_uv.y + s2 * jet.d_vv.y,
        c2 * jet.d_uu.z + 2.0 * cs * jet.d_uv.z + s2 * jet.d_vv.z,
    );
    let d_vv_rot = Point3::new(
        s2 * jet.d_uu.x - 2.0 * cs * jet.d_uv.x + c2 * jet.d_vv.x,
        s2 * jet.d_uu.y - 2.0 * cs * jet.d_uv.y + c2 * jet.d_vv.y,
        s2 * jet.d_uu.z - 2.0 * cs * jet.d_uv.z + c2 * jet.d_vv.z,
    );
    let d_uv_rot = Point3::new(
        cs * (jet.d_uu.x - jet.d_vv.x) + (c2 - s2) * jet.d_uv.x,
        cs * (jet.d_uu.y - jet.d_vv.y) + (c2 - s2) * jet.d_uv.y,
        cs * (jet.d_uu.z - jet.d_vv.z) + (c2 - s2) * jet.d_uv.z,
    );

    TaylorJet2 {
        value: jet.value,
        d_u: d_u_rot,
        d_v: d_v_rot,
        d_uu: d_uu_rot,
        d_uv: d_uv_rot,
        d_vv: d_vv_rot,
    }
}

/// Scale a [`TaylorJet2`]'s derivative components (but not value)
/// by `(s, s, s², s², s²)`. Used to convert a sector-local star
/// jet into its piece-local equivalent (`s = 1/2` for the O piece
/// covering the `[0, 1/2]²` quadrant).
pub fn scale_jet_for_piece_local(jet: &TaylorJet2, s: f64) -> TaylorJet2 {
    let s2 = s * s;
    TaylorJet2 {
        value: jet.value,
        d_u: scale(jet.d_u, s),
        d_v: scale(jet.d_v, s),
        d_uu: scale(jet.d_uu, s2),
        d_uv: scale(jet.d_uv, s2),
        d_vv: scale(jet.d_vv, s2),
    }
}

#[inline]
fn sub(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

#[inline]
fn scale(p: Point3, s: f64) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

#[cfg(test)]
#[path = "kp13_bi5_gtri3_tests.rs"]
mod tests;
