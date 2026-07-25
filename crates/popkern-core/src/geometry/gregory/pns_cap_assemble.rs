//! Phase IV PnS — End-to-end bi-4 cap assembly.
//!
//! Composes the PnS pipeline primitives into a single function
//! [`assemble_pns_bi4_cap`] that turns a boundary ring + pole into
//! `N` bi-4 sector patches:
//!
//! ```text
//! (Bi2TensorBorderRing, β, σ, pole) → N × [[Point3; 5]; 5]
//! ```
//!
//! # Architecture
//!
//! Storage: `patches[k][u_c][v_c]` uses **canonical** indexing
//! (sector's own parameterization). Each sector's 25 BB coefs are
//! classified by source:
//!
//! | Source                          | Canonical `(u_c, v_c)`             | Count |
//! |---------------------------------|------------------------------------|-------|
//! | Pole (constant)                 | `(0, 0)`                           | 1     |
//! | Central `a∘σ^k`                 | `u_c+v_c ≤ 2` (excl. pole)         | 5     |
//! | Boundary `b∘β` (v_c=3)          | `(0..4, 3)`                        | 5     |
//! | Boundary `b∘β` (v_c=4)          | `(0..4, 4)`                        | 5     |
//! | Prev-sharing (cur-prev seam)    | `(1..4, 0)`                        | 4     |
//! | Eq (2) closed-forms             | `(2, 1), (3, 1), (4, 1)`           | 3     |
//! | Free (F_κ min)                  | `(1, 2), (2, 2), (3, 2), (4, 2)`   | 4     |
//!
//! Sector-sharing relations (in canonical coords):
//! - cur-nxt seam: `cur[0][v_c] = nxt[v_c][0]` for `v_c ∈ 0..4`
//! - cur-prev seam: `cur[u_c][0] = prev[0][u_c]` for `u_c ∈ 0..4`
//!
//! Paper vs canonical axis swap (for cur side only):
//! - paper `p_{ij}` (cur) = canonical `p^k_{j, i}` (swapped)
//! - paper `´p_{ij}` (nxt) = canonical `p^{k+1}_{i, j}` (no swap)
//!
//! Eq (2) closed-forms (derived from paper's 5 outputs):
//! - paper cur_{10}, cur_{20}: automatic via central chain rule +
//!   paper's σ_{20,x} = (1 − 4/(3c))·d_1 + (4·cos(α/2)/(3c))·d_2
//!   (transcription fix 2026-04-24: was (1 − 4c/3), should be
//!   (1 − 4/(3c)))
//! - paper cur_{40} → canonical (4, 1) for sector k:
//!   `p^k_{4, 1} = 2·p^{k-1}_{0, 4} − p^{k-1}_{1, 4}`
//! - paper cur_{30} → canonical (3, 1) for sector k:
//!   `p^k_{3, 1} = [2·(4-c)·p^{k-1}_{0, 3} + c·(p^{k-1}_{1, 4} + p^k_{4, 1}) − 4·p^{k-1}_{1, 3}] / 4`
//! - paper nxt_{21} → canonical (2, 1) for sector k:
//!   `p^k_{2, 1} = −p^{k-1}_{1, 2} + (2-c)·p^{k-1}_{0, 2} + c·p^{k-1}_{0, 3}`
//!
//! # DOFs and the F_κ solve
//!
//! The truly free variables, determined jointly by the F_κ
//! minimisation (paper eq 6), are:
//!
//! - `a`'s 5 Point3 components `(a_x, a_y, a_xx, a_xy, a_yy)` — 15
//!   scalars (shared across sectors; `a_0` pinned to pole)
//! - `4 × N` free interior CPs per sector — `12·N` scalars
//!
//! Total: `15 + 12·N` scalars. Problem is separable per coordinate
//! (x, y, z), so we solve three `(5 + 4·N)`-dimensional symmetric
//! PSD systems.
//!
//! # Valence support
//!
//! Supported: `n ∈ [3, 12] \ {4}`. The bi-4 construction is
//! numerically singular at `n=4` (σ_{20,x} has a `1/c` factor
//! with `c=cos(π/2)=0`); for that case callers should use the
//! bi-3 variant or a different cap construction.
//!
//! Even `n ∈ {6, 8, 10, 12}` is supported through the asymmetric
//! closed-form treatment: instead of symmetrically solving
//! p_{2,1} via a circulant pair-sum (rank `n−1` for even n), we
//! treat canonical `(1, 2)` as a free DOF and derive canonical
//! `(2, 1)` one-way from the previous sector's data. No rank
//! deficiency to resolve.
//!
//! # Auto-compute pole
//!
//! [`pole_from_boundary`] derives the central pole from the
//! boundary ring via paper's eq (5) (see [`super::pns_center`]).
//! Useful for callers who don't want to compute the pole
//! manually. Assumes bi-4 e_def weight.
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` §4 (eqs 2, 4, 5, 6), Fig. 8c. The
//! paper-vs-canonical axis swap was clarified by clean-room
//! feedback from the KK-Peters reference implementation author
//! (2026-04-24); see memory note `project_pns_sigma_bug.md`.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::pns_beta::BetaReparam;
use super::pns_boundary::{Bi2TensorBorderRing, compose_bi4_outer_layers};
use super::pns_center::{central_point, e_def_bi4};
use super::pns_central_quad::{QuadMapJet, rotate_sigma_jet};
use super::pns_sigma::SigmaCenterJet;
use crate::geometry::Point3;

#[path = "pns_cap_assemble_helpers.rs"]
mod pns_cap_assemble_helpers;

use pns_cap_assemble_helpers::{
    build_coef_map, reconstruct_a, reconstruct_patches, solve_dof_system,
};

/// Inputs for the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleInput {
    /// Per-sector bi-2 tensor-border strips.
    pub boundary: Bi2TensorBorderRing,
    /// β reparameterisation for `n` sectors.
    pub beta: BetaReparam,
    /// Canonical σ centre jet.
    pub sigma: SigmaCenterJet,
    /// Pole (central point, eq 5).
    pub pole: Point3,
    /// F_κ order, `κ ∈ 1..=4`.
    pub kappa: usize,
}

/// Result of the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleResult {
    /// Fitted quadratic map `a` (with `a.a_0 = pole`).
    pub a: QuadMapJet,
    /// `N` bi-4 sector patches.
    pub patches: Vec<[[Point3; 5]; 5]>,
}

/// Cap assembly error cases.
#[derive(Debug, thiserror::Error)]
pub enum CapAssembleError {
    /// Valence is out of the supported range `[3, 12]`.
    #[error("PnS cap assemble: n must be in [3, 12], got {0}")]
    InvalidValence(usize),
    /// The supplied σ center jet contains a non-finite coefficient.
    ///
    /// Assembly consumes a jet it did not build, so this is the
    /// backstop against a caller passing a divergent σ — most
    /// plausibly a hand-rolled `n = 4` jet whose `1/c` pole did not
    /// cancel. Callers building σ through
    /// [`super::pns_sigma::sigma_center_jet`] get a typed
    /// `SigmaJetError` at construction instead and never reach here.
    #[error("PnS cap assemble: σ center jet has a non-finite {coefficient} coefficient")]
    NonFiniteSigmaJet {
        /// Which jet coefficient was non-finite.
        coefficient: &'static str,
    },
    /// `κ` is out of the supported range `1..=4`.
    #[error("PnS cap assemble: kappa must be 1..=4, got {0}")]
    InvalidKappa(usize),
    /// The F_κ linear system was singular (shouldn't happen in
    /// practice for well-posed input).
    #[error("PnS cap assemble: F_κ linear system is singular")]
    SingularSystem,
}

/// Compute the central pole via paper's eq (5) from the boundary
/// ring alone. The DS-net points `d^k_1, d^k_2, d^k_4` are
/// interpreted from the bi-2 tensor-border's inner row:
///
/// ```text
///   b^k_{i, 1} = strip.pos[i] + strip.deriv[i] / 2    (standard bi-2 derivative form)
///
///   d^k_1 = b^k_{1, 1}    (inner-row middle)
///   d^k_2 = b^k_{2, 1}    (inner-row right)
///   d^k_4 = b^k_{0, 1}    (inner-row left)
/// ```
///
/// Returns `pole = e·avg_k(d^k_1) + (1−e)·avg_k((d^k_2+d^k_4)/2)`
/// with `e = e_def_bi4(n)`.
pub fn pole_from_boundary(boundary: &Bi2TensorBorderRing) -> Point3 {
    let n = boundary.n();
    let e = e_def_bi4(n);
    let mut d1 = Vec::with_capacity(n);
    let mut d24_half = Vec::with_capacity(n);
    for strip in &boundary.sectors {
        // b^k_{i, 1} = pos[i] + deriv[i] / 2 (bi-2 derivative: row 1 = row 0 + deriv/2).
        let b01 = Point3::new(
            strip.pos[0].x + 0.5 * strip.deriv[0].x,
            strip.pos[0].y + 0.5 * strip.deriv[0].y,
            strip.pos[0].z + 0.5 * strip.deriv[0].z,
        );
        let b11 = Point3::new(
            strip.pos[1].x + 0.5 * strip.deriv[1].x,
            strip.pos[1].y + 0.5 * strip.deriv[1].y,
            strip.pos[1].z + 0.5 * strip.deriv[1].z,
        );
        let b21 = Point3::new(
            strip.pos[2].x + 0.5 * strip.deriv[2].x,
            strip.pos[2].y + 0.5 * strip.deriv[2].y,
            strip.pos[2].z + 0.5 * strip.deriv[2].z,
        );
        d1.push(b11);
        d24_half.push(Point3::new(
            0.5 * (b01.x + b21.x),
            0.5 * (b01.y + b21.y),
            0.5 * (b01.z + b21.z),
        ));
    }
    central_point(&d1, &d24_half, e)
}

/// Reject a σ center jet carrying a non-finite coefficient before it
/// can poison the whole cap.
///
/// Assembly receives σ as data rather than building it, so a caller
/// that hand-rolls a divergent jet (classically an `n = 4` jet whose
/// `1/c` pole did not cancel) would otherwise propagate NaN silently
/// into every patch. Jets built via
/// [`super::pns_sigma::sigma_center_jet`] are already checked at
/// construction and always pass here.
fn check_sigma_jet_finite(sigma: &SigmaCenterJet) -> Result<(), CapAssembleError> {
    for (coefficient, value) in [
        ("σ(0,0)", sigma.s00),
        ("∂_u σ", sigma.ds_u),
        ("∂_v σ", sigma.ds_v),
        ("∂²_u σ", sigma.ds_uu),
        ("∂²_v σ", sigma.ds_vv),
        ("∂²_uv σ", sigma.ds_uv),
    ] {
        if !value[0].is_finite() || !value[1].is_finite() {
            return Err(CapAssembleError::NonFiniteSigmaJet { coefficient });
        }
    }
    Ok(())
}

/// Assemble a PnS bi-4 cap from boundary + pole + σ jet.
pub fn assemble_pns_bi4_cap(
    input: &CapAssembleInput,
) -> Result<CapAssembleResult, CapAssembleError> {
    let n = input.boundary.n();
    if !(3..=12).contains(&n) {
        return Err(CapAssembleError::InvalidValence(n));
    }
    // n=4 was long rejected here on the grounds that c = cos(2π/4) = 0
    // makes σ_{20,x}'s `1/(3c)` factor singular. That reading was
    // wrong: regrouped as `d₁ + (4/(3c))·(cos(α/2)·d₂ − d₁)`, the
    // residue vanishes bit-exactly at n=4 under the paper's Appendix B
    // scalings, so the pole is removable and the jet is finite. σ now
    // carries that limit itself, and reports an irremovable pole as a
    // typed error for non-default scalings — see `pns_sigma`. Nothing
    // else in this assembly divides by c, so n=4 needs no other guard.
    // Other even n (6, 8, 10, 12) work fine — my asymmetric
    // closed-form treatment (treating canonical (1, 2) as the
    // free DOF and (2, 1), (3, 1), (4, 1) as one-way derivations
    // from prev's data) avoids the rank-deficient symmetric
    // pair-sum issue that would otherwise plague even n.
    if !(1..=4).contains(&input.kappa) {
        return Err(CapAssembleError::InvalidKappa(input.kappa));
    }
    check_sigma_jet_finite(&input.sigma)?;

    // Precompute per-sector σ^k and boundary (row0, row1).
    //
    // σ convention (refined 2026-04-24 per follow-up clarification):
    // σ^k = R_{kα} · σ^0 coefficient-by-coefficient (image-space
    // rotation). Each sector's BB-coefs genuinely differ — σ¹_{ij}
    // = R_α·σ⁰_{ij}. The "common recipe" phrase referred to the
    // FORMULAS in d₁, d₂, d₃ being identical per sector, not the
    // BB-coef values in ambient ℝ².
    let mut sigma_k = Vec::with_capacity(n);
    let mut row0 = Vec::with_capacity(n);
    let mut row1 = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        sigma_k.push(rotate_sigma_jet(&input.sigma, theta));
        let (r0, r1) = compose_bi4_outer_layers(&input.boundary.sectors[k], &input.beta);
        row0.push(r0);
        row1.push(r1);
    }

    // Build the coefficient map: for each sector k and each CP
    // (i, j), express it as a linear combination of DOFs plus a
    // constant Point3.
    let coef_map = build_coef_map(input, &sigma_k, &row0, &row1);

    // Assemble the DOF-space Gram matrix and RHS, solve per coord,
    // reconstruct patches.
    let dofs = solve_dof_system(&coef_map, input.kappa, n)?;
    let patches = reconstruct_patches(&coef_map, &dofs, n);
    let a = reconstruct_a(&dofs, input.pole);

    Ok(CapAssembleResult { a, patches })
}

#[cfg(test)]
#[path = "pns_cap_assemble_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "pns_cap_assemble_n4_tests.rs"]
mod n4_tests;
