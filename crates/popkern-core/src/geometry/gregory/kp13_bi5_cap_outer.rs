//! KP13 § 5 I.14.b — tensor-border transfer to the outer CPs of
//! the auxiliary cap `a`.
//!
//! Per KP13 § 5: "The outermost layers with BB-coefficients
//! `a^{s,k}_{ij}`, `s ∈ {ℓ, m, r}` match the tensor-border
//! reparameterized by β as detailed in § 4.2." The reparameterisation
//! β is a two-piece bi-3 map on `[0, 1]²` with specific coefficients
//! (KP13 eq (15)) that ensures `b ∘ β` is `C²` internally and joins
//! adjacent sectors `G²` per eqs (11), (12).
//!
//! # Scope of this commit (β = identity, corner-jet-only)
//!
//! This module implements a **simplified** transfer: for each
//! sector `k`, the 6 CPs controlling the 2nd-order BB-jet at
//! each outer corner of the non-`o` pieces are set from the
//! tensor-border `b^k`'s jet at the corresponding corner. The full
//! paper-prescribed β reparameterisation and outer-layer (not just
//! corner) transfer remain additive work documented below.
//!
//! ## Corner mapping (sector-local → piece-local)
//!
//! | sector-local corner | piece  | piece-local corner |
//! |---------------------|--------|--------------------|
//! | `m_k   = (1, 0)`    | `R`    | `(1, 0)`           |
//! | `m_{k+1} = (0, 1)`  | `L`    | `(0, 1)`           |
//! | `v_k   = (1, 1)`    | `M`    | `(1, 1)`           |
//!
//! The origin corner `(0, 0)` = star is owned by the `O` piece
//! (piece-local `(0, 0)`) and is populated later by the I.11 / I.14.c
//! pipeline from the guide, not the tensor-border.
//!
//! # Deferred (additive refinements)
//!
//! - **Full outer-layer transfer**: the paper sets the outer *two*
//!   BB layers (depth-2 ring matching the bi-3 tensor-border of
//!   depth 2), not only the 2nd-order corner jet. Implementing this
//!   requires degree-raising `b^k` from bi-3 to bi-5 and splitting
//!   into sub-patches via de Casteljau to match the 2 × 2 piece
//!   layout.
//! - **β reparameterisation**: the paper's β^m, β^r bi-3 reparam
//!   coefs are known in closed form
//!   ([`super::kp13_bi5_reparam::Bi5Reparams`]) but not yet
//!   composed into `b ∘ β` for the outer-layer transfer.
//! - **G² connectivity**: the full outer-layer transfer is what
//!   makes `a^{r,k-1}` ↔ `a^{ℓ,k}` join `G²` after eq (18);
//!   corner-jet-only transfer is insufficient.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_cap::{Bi5Cap, Bi5CapPiece, CapPieceId};
use super::kp13_bi5_chi::ChiCorner;
use super::kp13_bi5_jet::TaylorJet2;
use crate::geometry::Point3;

/// Populate the outer corner jets of a cap's `ℓ, m, r` pieces from
/// the tensor-border's jets at `m_k`, `m_{k+1}`, `v_k` respectively.
///
/// Uses β = identity (no reparameterisation). See module header for
/// what's deferred.
///
/// Each call writes exactly 3 × 6 × `n` = `18n` CPs total.
pub fn populate_outer_corner_jets_identity(cap: &mut Bi5Cap, target: &TensorBorder) {
    assert_eq!(
        cap.n,
        target.n(),
        "populate_outer_corner_jets_identity: cap and target n mismatch"
    );
    for k in 0..cap.n {
        let j_mk = target.jet(k, ChiCorner::Mk);
        let j_mk1 = target.jet(k, ChiCorner::Mk1);
        let j_vk = target.jet(k, ChiCorner::Vk);
        let sector = cap.sector_mut(k);
        set_corner_jet_10(sector.piece_mut(CapPieceId::R), &j_mk);
        set_corner_jet_01(sector.piece_mut(CapPieceId::L), &j_mk1);
        set_corner_jet_11(sector.piece_mut(CapPieceId::M), &j_vk);
    }
}

/// Write the six bi-5 CPs that determine the 2nd-order BB-jet at
/// piece-local corner `(u, v) = (1, 0)`.
///
/// Uses the analytical bi-5 corner-jet → CP inversion (standard
/// Bernstein basis identities). The other 30 CPs of the piece are
/// untouched.
pub fn set_corner_jet_10(piece: &mut Bi5CapPiece, jet: &TaylorJet2) {
    // Bi-5 at piece corner (1, 0):
    //   value = P[5][0]
    //   ∂u   = 5 (P[5][0] − P[4][0])
    //   ∂v   = 5 (P[5][1] − P[5][0])
    //   ∂²u  = 20 (P[5][0] − 2 P[4][0] + P[3][0])
    //   ∂²v  = 20 (P[5][2] − 2 P[5][1] + P[5][0])
    //   ∂u∂v = 25 (P[5][1] − P[5][0] − P[4][1] + P[4][0])
    let p50 = jet.value;
    let p40 = sub(p50, scale(jet.d_u, 1.0 / 5.0));
    let p51 = add(p50, scale(jet.d_v, 1.0 / 5.0));
    let p30 = sub_add(p40, p40, p50, scale(jet.d_uu, 1.0 / 20.0));
    let p52 = sub_add(p51, p51, p50, scale(jet.d_vv, 1.0 / 20.0));
    // P[4][1] = P[5][1] − P[5][0] + P[4][0] − d_uv / 25
    let p41 = sub(add(sub(p51, p50), p40), scale(jet.d_uv, 1.0 / 25.0));

    piece.cps[5][0] = p50;
    piece.cps[4][0] = p40;
    piece.cps[5][1] = p51;
    piece.cps[3][0] = p30;
    piece.cps[5][2] = p52;
    piece.cps[4][1] = p41;
}

/// Write the six bi-5 CPs that determine the 2nd-order BB-jet at
/// piece-local corner `(u, v) = (0, 1)`.
pub fn set_corner_jet_01(piece: &mut Bi5CapPiece, jet: &TaylorJet2) {
    // Bi-5 at piece corner (0, 1):
    //   value = P[0][5]
    //   ∂u   = 5 (P[1][5] − P[0][5])
    //   ∂v   = 5 (P[0][5] − P[0][4])
    //   ∂²u  = 20 (P[2][5] − 2 P[1][5] + P[0][5])
    //   ∂²v  = 20 (P[0][5] − 2 P[0][4] + P[0][3])
    //   ∂u∂v = 25 (P[1][5] − P[0][5] − P[1][4] + P[0][4])
    let p05 = jet.value;
    let p15 = add(p05, scale(jet.d_u, 1.0 / 5.0));
    let p04 = sub(p05, scale(jet.d_v, 1.0 / 5.0));
    let p25 = sub_add(p15, p15, p05, scale(jet.d_uu, 1.0 / 20.0));
    let p03 = sub_add(p04, p04, p05, scale(jet.d_vv, 1.0 / 20.0));
    // P[1][4] = P[1][5] − P[0][5] + P[0][4] − d_uv / 25
    let p14 = sub(add(sub(p15, p05), p04), scale(jet.d_uv, 1.0 / 25.0));

    piece.cps[0][5] = p05;
    piece.cps[1][5] = p15;
    piece.cps[0][4] = p04;
    piece.cps[2][5] = p25;
    piece.cps[0][3] = p03;
    piece.cps[1][4] = p14;
}

/// Write the six bi-5 CPs that determine the 2nd-order BB-jet at
/// piece-local corner `(u, v) = (1, 1)`.
///
/// Mirrors [`super::kp13_bi5_guide::from_gtri5_with_eq8_bilinear`]'s
/// analytical inversion.
pub fn set_corner_jet_11(piece: &mut Bi5CapPiece, jet: &TaylorJet2) {
    // Bi-5 at piece corner (1, 1):
    //   value = P[5][5]
    //   ∂u   = 5 (P[5][5] − P[4][5])
    //   ∂v   = 5 (P[5][5] − P[5][4])
    //   ∂²u  = 20 (P[5][5] − 2 P[4][5] + P[3][5])
    //   ∂²v  = 20 (P[5][5] − 2 P[5][4] + P[5][3])
    //   ∂u∂v = 25 (P[5][5] − P[4][5] − P[5][4] + P[4][4])
    let p55 = jet.value;
    let p45 = sub(p55, scale(jet.d_u, 1.0 / 5.0));
    let p54 = sub(p55, scale(jet.d_v, 1.0 / 5.0));
    let p35 = sub_add(p45, p45, p55, scale(jet.d_uu, 1.0 / 20.0));
    let p53 = sub_add(p54, p54, p55, scale(jet.d_vv, 1.0 / 20.0));
    // P[4][4] = P[4][5] + P[5][4] − P[5][5] + d_uv / 25
    let p44 = add(sub(add(p45, p54), p55), scale(jet.d_uv, 1.0 / 25.0));

    piece.cps[5][5] = p55;
    piece.cps[4][5] = p45;
    piece.cps[5][4] = p54;
    piece.cps[3][5] = p35;
    piece.cps[5][3] = p53;
    piece.cps[4][4] = p44;
}

// ---- small vector helpers ----

#[inline]
fn add(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

#[inline]
fn sub(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

#[inline]
fn scale(a: Point3, s: f64) -> Point3 {
    Point3::new(s * a.x, s * a.y, s * a.z)
}

/// `a + b − c + delta`. Useful for expressing the 2×·adj − corner +
/// (scaled 2nd-order) pattern that shows up in all three corner
/// setters.
#[inline]
fn sub_add(a: Point3, b: Point3, c: Point3, delta: Point3) -> Point3 {
    Point3::new(
        a.x + b.x - c.x + delta.x,
        a.y + b.y - c.y + delta.y,
        a.z + b.z - c.z + delta.z,
    )
}

#[cfg(test)]
#[path = "kp13_bi5_cap_outer_tests.rs"]
mod tests;
