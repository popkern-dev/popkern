//! Phase I.11e.5 — KP13 eq (6) post-LSQ exact match at `m_k`.
//!
//! After the LSQ gives a least-squares fit, the paper's eq (6) pins
//! the jet at every `m_k` corner *exactly* to the tensor-border:
//!
//! > "Next, to improve the fit, we re-compute, locally, one by one,
//! > a part of g△^{k-1} and g△^k so that
//! >
//! >   H_{m_k}(g△^k ∘ χ) = H_{m_k}(b^k)       ... (6)
//! >
//! > Fig. 5b illustrates (6): nine BB-coefficients of two pieces of
//! > g△ (marked by filled circles) are determined by nine BB-
//! > coefficients of b (marked by black boxes) plus some coefficients
//! > of g△ (marked by circles). Not-marked coefficients are
//! > determined by enforcing C² continuity of g△."
//!
//! # Concrete derivation (bilinear χ scope)
//!
//! At corner `(u, v) = (1, 0)` (paper `m_k` in the χ-domain unit
//! square), the tri-5 BB-jet involves exactly six CPs of `g△^k`:
//!
//! | component  | linear combination                                     |
//! |------------|--------------------------------------------------------|
//! | value      | `flat 20`                                              |
//! | ∂u         | `5·(flat 20 − flat 18)`                                |
//! | ∂v         | `5·(flat 19 − flat 18)`                                |
//! | ∂²u        | `20·(flat 20 − 2·flat 18 + flat 15)`                   |
//! | ∂²uv       | `20·(flat 19 − flat 18 − flat 16 + flat 15)`           |
//! | ∂²v        | `20·(flat 17 − 2·flat 16 + flat 15)`                   |
//!
//! Flat-index legend per sector `k`:
//!   - 20 = `g^k_{050}` = `m_k`
//!   - 18 = `g^k_{140}` (one step from `m_k` toward star along the
//!     sector edge with sector `k-1`)
//!   - 19 = `g^k_{141}` (depth-1 off that edge)
//!   - 15 = `g^k_{230}` (two steps from `m_k` toward star)
//!   - 16 = `g^k_{131}`
//!   - 17 = `g^k_{032}`
//!
//! The map is lower-triangular in forward order, so we can solve
//! by substitution in six steps.
//!
//! Shared edges: `flat 20 = flat 5 of sector (k-1)` (the `m_k`
//! vertex itself), `flat 18 = flat 4 of sector (k-1)`, `flat 15 =
//! flat 3 of sector (k-1)`. After writing, we call
//! `sync_shared_edge_all_depths` to propagate.
//!
//! # Re-run of C² back-fills
//!
//! The paper's "not-marked coefs are determined by C²" hint maps
//! to re-running the circulant + eq (3) back-fills AFTER eq (6)
//! overrides. We expose both the raw "write the 6 CPs" and a
//! higher-level "run full pipeline + eq (6)" entry point. The
//! user can choose.
//!
//! # Interplay with higher-j C¹ back-fills
//!
//! `fill_g31_via_c1` would overwrite `flat 16` of the *next* sector
//! using this sector's `flat 3`/`flat 4`/`flat 9`. Since `flat 3` of
//! sector `k` equals `flat 15` of sector `k+1` (and we've just
//! written a new `flat 15`), re-running `fill_g31_via_c1` is
//! consistent: sector `k+1`'s new `flat 16` is derived from the
//! new shared-edge values. Similarly `fill_g41_via_c1`.

#![allow(dead_code)]

use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_chi::ChiCorner;
use super::kp13_bi5_circulant::{
    CirculantError, complete_derivative_layers_with_higher_j, sync_shared_edge_all_depths,
};
use super::kp13_bi5_gtri5::Gtri5;
use super::kp13_bi5_jet::TaylorJet2;
use crate::geometry::Point3;

/// Near-`m_k` tri-5 CP indices for sector `k`, in the forward-
/// substitution order used by [`solve_mk_corner_cps_per_axis`].
pub const MK_CORNER_FLATS: [usize; 6] = [20, 18, 19, 15, 16, 17];

/// Solve the 6-CP forward-substitution system **per axis** to pin
/// the jet at `(u, v) = (1, 0)` to the given target.
///
/// Returns `[flat 20, flat 18, flat 19, flat 15, flat 16, flat 17]`
/// per [`MK_CORNER_FLATS`].
pub fn solve_mk_corner_cps_per_axis(target: &TaylorJet2) -> [Point3; 6] {
    // Per-axis scalar solve; Point3 is componentwise.
    let solve = |val: f64, du: f64, dv: f64, duu: f64, duv: f64, dvv: f64| -> [f64; 6] {
        // 1. flat 20 = value
        let f20 = val;
        // 2. ∂u = 5·(flat 20 − flat 18) → flat 18 = flat 20 − ∂u/5
        let f18 = f20 - du / 5.0;
        // 3. ∂v = 5·(flat 19 − flat 18) → flat 19 = flat 18 + ∂v/5
        let f19 = f18 + dv / 5.0;
        // 4. ∂²u = 20·(flat 20 − 2·flat 18 + flat 15)
        //    → flat 15 = ∂²u/20 − flat 20 + 2·flat 18
        let f15 = duu / 20.0 - f20 + 2.0 * f18;
        // 5. ∂²uv = 20·(flat 19 − flat 18 − flat 16 + flat 15)
        //    → flat 16 = flat 19 − flat 18 + flat 15 − ∂²uv/20
        let f16 = f19 - f18 + f15 - duv / 20.0;
        // 6. ∂²v = 20·(flat 17 − 2·flat 16 + flat 15)
        //    → flat 17 = ∂²v/20 + 2·flat 16 − flat 15
        let f17 = dvv / 20.0 + 2.0 * f16 - f15;
        [f20, f18, f19, f15, f16, f17]
    };
    let x = solve(
        target.value.x,
        target.d_u.x,
        target.d_v.x,
        target.d_uu.x,
        target.d_uv.x,
        target.d_vv.x,
    );
    let y = solve(
        target.value.y,
        target.d_u.y,
        target.d_v.y,
        target.d_uu.y,
        target.d_uv.y,
        target.d_vv.y,
    );
    let z = solve(
        target.value.z,
        target.d_u.z,
        target.d_v.z,
        target.d_uu.z,
        target.d_uv.z,
        target.d_vv.z,
    );
    [
        Point3::new(x[0], y[0], z[0]),
        Point3::new(x[1], y[1], z[1]),
        Point3::new(x[2], y[2], z[2]),
        Point3::new(x[3], y[3], z[3]),
        Point3::new(x[4], y[4], z[4]),
        Point3::new(x[5], y[5], z[5]),
    ]
}

/// Write the 6 near-`m_k` CPs of every sector so that
/// `H_{m_k}(g△^k ∘ χ_bil^k)` matches the target `b^k` jet at `m_k`
/// **exactly**.
///
/// Per-sector writes happen independently; after all sectors are
/// written, `resync_mid_edges` propagates each new `flat 20` (= m_k)
/// into the adjacent sector's `flat 5`, and `sync_shared_edge_all_depths`
/// propagates `flat 15` / `flat 18` into their shared-edge partners
/// (`flat 3` / `flat 4` of the neighbour). The mid-edge direction is
/// authoritative from the higher-index sector because that's where
/// the CP naturally lives at its paper `050` index.
///
/// If the tensor-border targets are inconsistent at a shared vertex
/// (e.g., `b^k`'s value at Mk disagrees with `b^{k+1}`'s value at
/// Mk1), the last-write-wins via these sync steps resolves the
/// conflict in favour of the higher-index sector.
///
/// Preserves the Q identity near the star (untouched) but may
/// break C² at cross-sector layer-2 positions (flat 12 etc.). See
/// [`apply_eq6_with_c2_reconcile`] for the canonical pipeline.
pub fn apply_eq6_raw(gtri5: &mut Gtri5, target: &TensorBorder) {
    let n = gtri5.n;
    for k in 0..n {
        let tgt = target.jet(k, ChiCorner::Mk);
        let cps = solve_mk_corner_cps_per_axis(&tgt);
        for (slot, &flat) in MK_CORNER_FLATS.iter().enumerate() {
            gtri5.sector_mut(k).cps[flat] = cps[slot];
        }
    }
    // Restore structural invariants: mid-edge share (flat 20 ↔
    // flat 5) and depth-0 shared-edge coefs (flat 15 ↔ flat 3,
    // flat 18 ↔ flat 4).
    gtri5.resync_mid_edges();
    sync_shared_edge_all_depths(gtri5);
}

/// End-to-end eq (6) application: write the 6 near-`m_k` CPs per
/// sector, sync shared edges, then re-run the C² chain to re-
/// derive the "not-marked" coefs (circulant + eq 3 back-fills)
/// consistent with the new values.
///
/// The re-run may slightly perturb `flat 12`, `flat 8`, `flat 16`,
/// `flat 19` of each sector because those depend on the newly-
/// written near-`m_k` values. This is the paper's "not-marked
/// coefficients determined by C²" step.
pub fn apply_eq6_with_c2_reconcile(
    gtri5: &mut Gtri5,
    target: &TensorBorder,
) -> Result<(), CirculantError> {
    apply_eq6_raw(gtri5, target);
    complete_derivative_layers_with_higher_j(gtri5)?;
    // After the reconcile, the near-m_k jet may no longer match
    // exactly — we've just recomputed some of those CPs. Re-apply
    // eq (6) one more time to take precedence over the C²-derived
    // values for the 6 marked CPs. apply_eq6_raw handles its own
    // invariant resync.
    apply_eq6_raw(gtri5, target);
    Ok(())
}

#[cfg(test)]
#[path = "kp13_bi5_eq6_tests.rs"]
mod tests;
