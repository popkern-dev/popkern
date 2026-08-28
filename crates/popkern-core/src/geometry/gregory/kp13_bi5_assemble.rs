//! End-to-end KP13 bi-5 macro-patch assembly + bi-9
//! elevation for fillet-pipeline consumption.
//!
//! Wires the full `g△ → g → a → p → bi-9` chain into a single
//! entry point, producing a commercial-grade bi-9 cap suitable for
//! integration with the existing `popkern-fillet` bi-9 dispatcher
//! (same output shape as the the bi-9 cap campaign direct bidegree-9 and the
//! KP15 bi-6 polar pipelines).
//!
//! # Pipeline
//!
//! 1. **LSQ fit `g△`** — solve the I.11e LSQ for the tri-5
//!    auxiliary map from the tensor-border + CC limit point.
//! 2. **Degree-raise to `g`** (I.12 foundation + eq (8) override)
//!    — produce the bi-5 guide surface.
//! 3. **Build auxiliary cap `a`** (I.14.a–d), sub-steps:
//!    allocate zero cap, populate `ℓ/m/r` depth-2 outer layers
//!    from target (β = identity — see I.14.b.2), populate `O`
//!    piece's corner-block CPs from the guide via eq (17)
//!    (I.14.c), propagate `C²` internally: `O → R/L`, then
//!    `R/L → M` (I.14.d).
//! 4. **Build final cap `p`** (I.14.e) — eq (18) G² strip
//!    transfer + F_1 minimisation over the 5n + 6 free DoFs.
//! 5. **Elevate to bi-9** — each of the `4n` bi-5 sub-patches is
//!    degree-raised to bi-9 via
//!    [`super::kp_degree_elevate::elevate_bi5_grid_to_bi9`].
//!
//! # Output shape
//!
//! [`Bi5CapOutput::patches_bi9`] is a row-major
//! `Vec<[[Point3; 10]; 10]>` of length `4n`, indexed by `(k,
//! piece)` as `4 * k + piece.index()` where `piece.index()` is
//! the canonical `CapPieceId` order `O, L, R, M`. This matches the
//! bi-6 polar cap's output convention, so downstream fillet code
//! can dispatch on a `Vec<[[Point3; 10]; 10]>` unconditionally.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::kp_degree_elevate::elevate_bi5_grid_to_bi9;
use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_cap::{Bi5Cap, CapPieceId};
use super::kp13_bi5_cap_connect::connect_all_sectors;
use super::kp13_bi5_cap_final::{CapMinError, solve_f1_min};
use super::kp13_bi5_cap_inner::populate_o_piece_from_guide_simplified;
use super::kp13_bi5_cap_outer::populate_outer_layers_full_identity;
use super::kp13_bi5_guide::Bi5Guide;
use super::kp13_bi5_lsq::{LsqError, solve_g_tri_lsq};
use crate::geometry::Point3;

/// Input to the full KP13 bi-5 macro-patch cap assembly.
#[derive(Debug, Clone)]
pub struct Bi5CapInput {
    /// Per-sector tensor-border `b^k` providing the boundary data
    /// the cap must interpolate. `target.n()` is the EV valence.
    pub target: TensorBorder,
    /// Catmull-Clark limit point for the star vertex of the cap-
    /// net. Pinned into `g△` at its star CP (paper `IDX_500`).
    pub cc_limit: Point3,
}

/// Error from the full KP13 bi-5 cap assembly.
#[derive(Debug, Clone, PartialEq)]
pub enum Bi5AssembleError {
    /// LSQ stage (I.11e) failed — typically singular normal
    /// equations or a shape mismatch on the tensor border.
    Lsq(LsqError),
    /// `F_1` minimisation stage (I.14.e) failed — singular normal
    /// equations on an axis.
    F1Min(CapMinError),
}

impl std::fmt::Display for Bi5AssembleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Bi5AssembleError::Lsq(e) => write!(f, "KP13 bi-5 LSQ: {e}"),
            Bi5AssembleError::F1Min(e) => write!(f, "KP13 bi-5 F_1 min: {e}"),
        }
    }
}

impl std::error::Error for Bi5AssembleError {}

impl From<LsqError> for Bi5AssembleError {
    fn from(e: LsqError) -> Self {
        Bi5AssembleError::Lsq(e)
    }
}

impl From<CapMinError> for Bi5AssembleError {
    fn from(e: CapMinError) -> Self {
        Bi5AssembleError::F1Min(e)
    }
}

/// Assembled output of the full KP13 bi-5 macro-patch pipeline.
#[derive(Debug, Clone, PartialEq)]
pub struct Bi5CapOutput {
    /// Valence of the extraordinary vertex.
    pub n: usize,
    /// Final bi-5 cap `p` with eq (18) G² strips + F_1-minimised
    /// free DoFs. Has `n` sectors × 4 bi-5 sub-patches each.
    pub cap: Bi5Cap,
    /// Degree-raised `4n` bi-9 CP grids suitable for consumption
    /// by the fillet bi-9 dispatcher. See module header for
    /// layout.
    pub patches_bi9: Vec<[[Point3; 10]; 10]>,
}

/// End-to-end KP13 bi-5 macro-patch cap assembly.
///
/// Runs the full `g△ → g → a → p → bi-9` pipeline. Returns an
/// error only if the LSQ / F_1-min stages have singular linear
/// systems — shouldn't happen for well-posed inputs at `n ≥ 3`.
pub fn assemble_bi5_macro_cap(input: &Bi5CapInput) -> Result<Bi5CapOutput, Bi5AssembleError> {
    let n = input.target.n();
    assert!(n >= 3, "assemble_bi5_macro_cap: valence n = {n} < 3");

    // Stage 1: LSQ fit g△ from tensor-border + CC limit.
    let g_tri = solve_g_tri_lsq(&input.target, input.cc_limit)?;

    // Stage 2: degree-raise g△ → g (bi-5 guide). Uses the
    // bilinear eq (8) path (paper-faithful corner-jet match at v_k
    // under bilinear χ; the CC-χ variant is also available via
    // Bi5Guide::from_gtri5_with_eq8_cc but not used here for
    // speed).
    let guide = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &input.target);

    // Stage 3: build auxiliary cap a from guide + target via
    // I.14.a–d.
    let mut a = Bi5Cap::zero(n);
    populate_outer_layers_full_identity(&mut a, &input.target);
    populate_o_piece_from_guide_simplified(&mut a, &guide);
    connect_all_sectors(&mut a);

    // Stage 4: build final cap p from a via I.14.e (eq 18 strips
    // + F_1 minimisation over 5n + 6 free DoFs).
    let p = solve_f1_min(&a)?;

    // Stage 5: elevate each bi-5 sub-patch to bi-9. Row-major:
    // for each sector k, for each piece in canonical order
    // (O, L, R, M), append the elevated grid.
    let mut patches_bi9 = Vec::with_capacity(4 * n);
    for k in 0..n {
        for id in CapPieceId::all() {
            let bi5 = p.sector(k).piece(id).cps;
            patches_bi9.push(elevate_bi5_grid_to_bi9(bi5));
        }
    }

    Ok(Bi5CapOutput {
        n,
        cap: p,
        patches_bi9,
    })
}

#[cfg(test)]
#[path = "kp13_bi5_assemble_tests.rs"]
mod tests;
