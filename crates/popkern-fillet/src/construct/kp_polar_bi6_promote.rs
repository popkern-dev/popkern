//! Phase I.10f — K-P 2015 bi-6 polar promotion for N ≥ 5 sub-
//! patch rings.
//!
//! Sibling to [`promote_decomp_to_kp_binonic`](super::three_quad_apply)
//! (Phase H, rectangular N ∈ {3, 4}). This module runs the
//! **polar** K-P 2015 bi-6 pipeline (core modules I.1–I.10d)
//! via the wrapper [`try_kp_polar_bi6_cap`](super::kp_polar_bi6_apply::try_kp_polar_bi6_cap)
//! and replaces the sub-patches in-place on success.
//!
//! # Dispatch
//!
//! Called from
//! [`apply_three_quad_decomposition_to_blend`](super::three_quad_apply::apply_three_quad_decomposition_to_blend)
//! when `n_sides >= 5 && options.use_kp_polar_gregory`. Returns
//! `Some(())` on success (sub_patches now bi-9 grids), `None` on
//! any failure (caller retains existing two-stripe / kp-binonic
//! fallback).
//!
//! # Pipeline
//!
//! 1. Reuse [`extract_direct_boundary_input`] (shared with
//!    `kp_binonic` path) to read per-spoke + per-corner data
//!    from the two-stripe sub-patches.
//! 2. Convert via
//!    [`PolarBoundaryInput::from_direct_boundary_input`].
//! 3. Run the full polar pipeline via `try_kp_polar_bi6_cap` —
//!    tensor-border extraction → bi-6 cap input → F_5 SPD min
//!    → bi-9 output.
//! 4. Wrap each `[[Point3; 10]; 10]` as
//!    [`RectangularGregoryBiNonic`] and commit.

use popkern_core::geometry::gregory::{
    GregoryPatch, PolarBoundaryInput, RectangularGregoryBiNonic, RectangularGregoryTwoStripe, Side,
    Spoke, SpokeEnd,
};

use super::kp_polar_bi6_apply::try_kp_polar_bi6_cap;
use super::three_quad_apply::extract_direct_boundary_input;

/// Attempt to promote a ring of two-stripe sub-patches to bi-6
/// polar K-P output. Overwrites `sub_patches` with N
/// `RectangularBiNonic` variants on success.
///
/// Returns `Some(())` on success; `None` on any failure (the
/// caller's existing fallback ladder then applies).
#[allow(dead_code)]
pub(super) fn promote_decomp_to_kp_polar_bi6(sub_patches: &mut [GregoryPatch]) -> Option<()> {
    let n = sub_patches.len();
    if n < 5 {
        return None;
    }

    // Every sub-patch must already be two-stripe — the biseptic-
    // promoted seed. (The ladder in `apply_three_quad_decomposition_to_blend`
    // ensures this by the time we're called.)
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for patch in sub_patches.iter() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            return None;
        };
        two_stripe.push((**ts).clone());
    }

    // Ring-topology spokes — identical convention to the binonic
    // path (spoke `k` shares sub-patch `k` UOne with sub-patch
    // `(k+1) % n` VZero reversed).
    let spokes: Vec<Spoke> = (0..n)
        .map(|k| Spoke {
            a: SpokeEnd {
                patch_idx: k,
                side: Side::UOne,
                reversed: false,
            },
            b: SpokeEnd {
                patch_idx: (k + 1) % n,
                side: Side::VZero,
                reversed: true,
            },
        })
        .collect();

    let direct = extract_direct_boundary_input(&two_stripe, &spokes)?;
    let polar = PolarBoundaryInput::from_direct_boundary_input(&direct).ok()?;
    let bi9_grids = try_kp_polar_bi6_cap(&polar).ok()?;
    if bi9_grids.len() != n {
        return None;
    }

    // Self-intersection gate — never commit a fold.
    let new_patches: Vec<GregoryPatch> = bi9_grids
        .into_iter()
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect();
    const SELF_INT_TOL: f64 = 1e-6;
    for gp in &new_patches {
        if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
            return None;
        }
    }

    for (idx, gp) in new_patches.into_iter().enumerate() {
        sub_patches[idx] = gp;
    }
    Some(())
}

#[cfg(test)]
#[path = "kp_polar_bi6_promote_tests.rs"]
mod tests;
