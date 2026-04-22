//! KP13 § 5 I.14.c — eq (17) `g → a` coupling for the `O` piece.
//!
//! Paper eq (17):
//!
//! ```text
//!   H_t a^{o,k} := H_t (C(g)^k · (1/λ) · L⁻¹ ∘ σ^o)
//! ```
//!
//! — sets the 2nd-order BB-jet of the `O` piece at each of its
//! four corners `t` to the corresponding jet of the composite
//! `C(g)^k · (1/λ) · L⁻¹ ∘ σ^o`. Here:
//!
//! - `C(g)^k` is the guide after one step of Catmull-Clark
//!   subdivision of the cap-net (paper § 5: "to smooth out abrupt
//!   transitions in curvature before computing (4) and (7)").
//! - `1/λ` compensates the CC contraction so the subdivided guide
//!   covers the same parametric region as the original.
//! - `L⁻¹ ∘ σ^o` reparameterises the `O` piece's parameter domain
//!   into the guide's parameter domain via the linear shear `L`
//!   (`[0, 1]² → parallelogram`) and the sector `σ` map.
//!
//! # Scope of this commit (MVP with documented simplifications)
//!
//! ## `C(g) = g`  (no CC subdivision)
//!
//! The paper's CC subdivision step is a **quality** measure for
//! smoothing abrupt curvature transitions in the input cap-net
//! before the guide is built. It is **not required for
//! correctness** — a cap built from the un-subdivided guide
//! still satisfies all the C²/G² constraints downstream; it just
//! tracks the raw tensor-border data more directly.
//!
//! Implementing `C(g)` requires the Catmull-Clark subdivision
//! operator on the KP13 cap-net (star + outer-ring vertices),
//! which is a substantial separate piece. Deferred.
//!
//! ## `σ^o = identity` in piece-local coordinates
//!
//! The paper's `σ^o` is the "top" (star-near) `σ` reparameterisation
//! from § 4.3, derived by minimising `F_3(σ)` subject to boundary
//! constraints (implemented in
//! [`super::kp13_bi5_sigma_solve::solve_sigma_f3_minimal`]).
//! Threading the full `σ^o` through the chain rule is mechanical
//! additive work; for the MVP we use `σ^o(u, v) = (u, v)`.
//!
//! ## `L⁻¹` folded into the piece-local ↔ sector-local coord change
//!
//! The `O` piece covers sector-local `(u_s, v_s) ∈ [0, 1/2]²`; in
//! piece-local coordinates `(u_p, v_p) ∈ [0, 1]²` with `(u_s, v_s)
//! = (u_p/2, v_p/2)`. The linear shear `L` pulls sector-local
//! coordinates into the 2-D ambient parallelogram. Under identity
//! `σ^o`, the paper's `L⁻¹ ∘ σ^o` composed with the piece-local
//! coord change is just the piece-local → sector-local pull-back
//! at scale `1/λ` — which is the only paper-prescribed correction
//! we apply here.
//!
//! # What the simplified formula reduces to
//!
//! For piece-local corner `(u_p, v_p)` the MVP evaluates the guide
//! at position `p = (u_p / (2λ), v_p / (2λ))` in the guide's
//! parameter domain, with derivatives chain-ruled by the
//! `1 / (2λ)` factor in each variable:
//!
//! ```text
//!   a^{o,k}(u_p, v_p)      = g^k(u_p/(2λ), v_p/(2λ))
//!   ∂u a^{o,k}            = (1/(2λ)) ∂u g^k
//!   ∂²u a^{o,k}           = (1/(2λ))² ∂²u g^k        (etc.)
//! ```
//!
//! The 2nd-order jet at each of the four corners is then inverted
//! onto the piece's six corner-block CPs via the same analytical
//! formulas used in [`super::kp13_bi5_cap_outer`].

#![allow(dead_code)]

use super::kp13_bi5_cap::{Bi5Cap, CapPieceId};
use super::kp13_bi5_cap_outer::{
    set_corner_jet_00, set_corner_jet_01, set_corner_jet_10, set_corner_jet_11,
};
use super::kp13_bi5_chi::cc_subdominant_eigenvalue;
use super::kp13_bi5_guide::Bi5Guide;
use super::kp13_bi5_jet::TaylorJet2;
use crate::geometry::Point3;

/// Populate the `O` piece of every sector of `cap` from `guide`
/// via the simplified eq (17) (see module header).
///
/// For each sector `k` and each of the four piece-local corners
/// `t ∈ {(0,0), (1,0), (0,1), (1,1)}`:
///
/// 1. Compute the guide's 2nd-order jet at the scaled position
///    `(t.u / (2λ), t.v / (2λ))`.
/// 2. Chain-rule the derivatives by the `1 / (2λ)` scaling.
/// 3. Invert the resulting piece-local jet onto the 6 bi-5 CPs
///    controlling the corresponding corner block.
///
/// Writes exactly `4 × 6 × n = 24n` CPs total; the other 12 CPs
/// per `O` piece (at positions `(2, 1), (1, 2), (2, 2), (3, 1),
/// (1, 3), (3, 2), (2, 3), (3, 3), (4, 1), (1, 4), (4, 2), (2, 4)`
/// — the non-corner-block interior) stay at their pre-call values
/// and will be populated by I.14.d (C² internal connection).
pub fn populate_o_piece_from_guide_simplified(cap: &mut Bi5Cap, guide: &Bi5Guide) {
    assert_eq!(
        cap.n, guide.n,
        "populate_o_piece_from_guide_simplified: cap and guide n mismatch"
    );
    let n = cap.n;
    let lambda = cc_subdominant_eigenvalue(n);
    let inv_2lambda = 0.5 / lambda;
    let inv_2lambda_sq = inv_2lambda * inv_2lambda;

    for k in 0..n {
        // Snapshot the guide sector so we can then borrow cap mutably.
        let g_sec = *guide.sector(k);
        let o_piece = cap.sector_mut(k).piece_mut(CapPieceId::O);

        for (u_p, v_p) in [(0.0_f64, 0.0_f64), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            let u_g = u_p * inv_2lambda;
            let v_g = v_p * inv_2lambda;
            let oh = g_sec.outer_hessian(u_g, v_g);

            // Chain rule for the piece-local → scaled-guide
            // coordinate change (u_p → u_g = u_p · inv_2lambda).
            // First partials multiply by inv_2lambda; second
            // partials by inv_2lambda².
            let jet = TaylorJet2 {
                value: oh.value,
                d_u: scale_point(oh.d_x, inv_2lambda),
                d_v: scale_point(oh.d_y, inv_2lambda),
                d_uu: scale_point(oh.d_xx, inv_2lambda_sq),
                d_uv: scale_point(oh.d_xy, inv_2lambda_sq),
                d_vv: scale_point(oh.d_yy, inv_2lambda_sq),
            };

            match (u_p as i32, v_p as i32) {
                (0, 0) => set_corner_jet_00(o_piece, &jet),
                (1, 0) => set_corner_jet_10(o_piece, &jet),
                (0, 1) => set_corner_jet_01(o_piece, &jet),
                (1, 1) => set_corner_jet_11(o_piece, &jet),
                _ => unreachable!(),
            }
        }
    }
}

#[inline]
fn scale_point(p: Point3, s: f64) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

#[cfg(test)]
#[path = "kp13_bi5_cap_inner_tests.rs"]
mod tests;
