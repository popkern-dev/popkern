//! KP13 bi-5 macro-patch promotion — sibling of
//! [`kp_polar_bi6_promote`](super::kp_polar_bi6_promote) but
//! produces `4N` bi-9 Gregory patches (paper's 2 × 2 piece layout
//! per polar sector) instead of `N`.
//!
//! # Dispatch
//!
//! Called from
//! [`apply_three_quad_decomposition_to_blend`](super::three_quad_apply::apply_three_quad_decomposition_to_blend)
//! when `options.use_kp_polar_gregory && n_sides >= 5` AND the
//! env var `POPKERN_KP13_BI5=1` is set. Replaces the `N`-entry
//! `sub_patches` vector with `4N` entries in canonical
//! `(sector, piece)` row-major order matching the shape expected
//! by [`poke_n_gon_face_into_4n_quads`](super::poke_4n).
//!
//! # Pipeline
//!
//! Mirrors the bi-6 polar pipeline up to the
//! [`PolarBoundaryInput`] boundary extraction, then diverges:
//!
//! 1. Reuse [`extract_direct_boundary_input`] to read per-spoke
//!    + per-corner data from the two-stripe sub-patches.
//! 2. Convert via
//!    [`PolarBoundaryInput::from_direct_boundary_input`].
//! 3. Call
//!    [`try_kp_bi5_macro_cap_as_gregory`](super::kp_bi5_macro_apply::try_kp_bi5_macro_cap_as_gregory)
//!    → `4N` [`RectangularGregoryBiNonic`] wrappers.
//! 4. Screen self-intersection across all `4N` patches; reject
//!    the whole promotion on any fold.
//! 5. Wrap each as `GregoryPatch::RectangularBiNonic` and return.
//!
//! # Diagnostics
//!
//! Set env var `POPKERN_KP13_BI5_TRACE=1` to stream the reason on
//! every fallback. Silent by default.

#![allow(dead_code)]

use popkern_core::geometry::gregory::{
    GregoryPatch, PolarBoundaryInput, RectangularGregoryTwoStripe, Side, Spoke, SpokeEnd,
};

use super::kp_bi5_macro_apply::try_kp_bi5_macro_cap_as_gregory;
use super::kp_polar_bi6_promote::extract_fillet_radius;
use super::three_quad_apply::extract_direct_boundary_input;
use crate::fillet_surface::FilletSurfaceResult;

/// Self-intersection tolerance matching the bi-6 path.
const SELF_INT_TOL: f64 = 1e-6;

/// Outcome of a KP13 bi-5 promotion attempt.
pub(super) enum Bi5PromoteOutcome {
    /// Success — `patches_4n` has `4 * sub_patches.len()` entries,
    /// ready to replace the decomposition's `sub_patches` vector.
    Committed { patches_4n: Vec<GregoryPatch> },
    /// Promotion did not run; caller should fall through to the
    /// bi-6 path.
    Rejected,
}

/// Attempt to promote a ring of two-stripe sub-patches to `4N`
/// bi-9 KP13 bi-5 macro-patches.
///
/// Returns [`Committed`] with a fresh `Vec<GregoryPatch>` of length
/// `4 * sub_patches.len()` on success, or [`Rejected`] on any
/// failure (non-two-stripe input, polar conversion error, assembly
/// failure, self-intersection).
///
/// [`Committed`]: Bi5PromoteOutcome::Committed
/// [`Rejected`]: Bi5PromoteOutcome::Rejected
pub(super) fn promote_decomp_to_kp13_bi5(
    sub_patches: &[GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Bi5PromoteOutcome {
    let n = sub_patches.len();
    if n < 3 {
        trace("n < 3", n, None);
        return Bi5PromoteOutcome::Rejected;
    }

    // Extract two-stripe view of sub-patches.
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for (idx, patch) in sub_patches.iter().enumerate() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            trace(
                "non-two-stripe input",
                n,
                Some(format!("sub-patch {idx} is {}", variant_name(patch))),
            );
            return Bi5PromoteOutcome::Rejected;
        };
        two_stripe.push((**ts).clone());
    }

    // Ring-topology spokes (matches the bi-6 path).
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

    let Some(mut direct) = extract_direct_boundary_input(&two_stripe, &spokes) else {
        trace("extract_direct_boundary_input → None", n, None);
        return Bi5PromoteOutcome::Rejected;
    };

    // Thread per-sector fillet radii into the direct boundary input.
    if fillet_surfaces.len() >= n {
        for (k, spoke) in direct.spokes.iter_mut().enumerate() {
            if let Some(r) = extract_fillet_radius(&fillet_surfaces[k]) {
                spoke.radius_at_c = r;
                spoke.radius_at_v = r;
            }
        }
    }

    let polar = match PolarBoundaryInput::from_direct_boundary_input(&direct) {
        Ok(p) => p,
        Err(e) => {
            trace("from_direct_boundary_input", n, Some(e));
            return Bi5PromoteOutcome::Rejected;
        }
    };

    // Run the KP13 bi-5 macro-patch pipeline.
    let wrapped = match try_kp_bi5_macro_cap_as_gregory(&polar) {
        Ok(w) => w,
        Err(e) => {
            trace("try_kp_bi5_macro_cap_as_gregory", n, Some(e));
            return Bi5PromoteOutcome::Rejected;
        }
    };
    if wrapped.len() != 4 * n {
        trace(
            "output length mismatch",
            n,
            Some(format!("expected 4N={}, got {}", 4 * n, wrapped.len())),
        );
        return Bi5PromoteOutcome::Rejected;
    }

    // Wrap each bi-9 grid as a GregoryPatch and screen for folds.
    let candidate: Vec<GregoryPatch> = wrapped
        .into_iter()
        .map(|bi_nonic| GregoryPatch::RectangularBiNonic(Box::new(bi_nonic)))
        .collect();
    for (idx, gp) in candidate.iter().enumerate() {
        if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
            trace("self-intersection", n, Some(format!("patch index {idx}")));
            return Bi5PromoteOutcome::Rejected;
        }
    }

    trace("accepted", n, Some(format!("4N={} bi-9 patches", 4 * n)));
    Bi5PromoteOutcome::Committed {
        patches_4n: candidate,
    }
}

/// Top-level dispatch — when the KP13 bi-5 path is gated on
/// (`POPKERN_KP13_BI5=1`), attempt the promotion. Returns
/// `Some(patches_4n)` on success (caller replaces
/// `decomp.sub_patches`) or `None` to fall through to the bi-6
/// ladder.
pub(super) fn try_kp13_bi5_dispatch(
    sub_patches: &[GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<Vec<GregoryPatch>> {
    std::env::var_os("POPKERN_KP13_BI5")?;
    match promote_decomp_to_kp13_bi5(sub_patches, fillet_surfaces) {
        Bi5PromoteOutcome::Committed { patches_4n } => Some(patches_4n),
        Bi5PromoteOutcome::Rejected => None,
    }
}

fn variant_name(p: &GregoryPatch) -> &'static str {
    match p {
        GregoryPatch::Rectangular(_) => "Rectangular",
        GregoryPatch::RectangularDeg5(_) => "RectangularDeg5",
        GregoryPatch::RectangularBiquinticFull(_) => "RectangularBiquinticFull",
        GregoryPatch::RectangularBiseptic(_) => "RectangularBiseptic",
        GregoryPatch::RectangularTwoStripe(_) => "RectangularTwoStripe",
        GregoryPatch::RectangularBiNonic(_) => "RectangularBiNonic",
        GregoryPatch::Triangular(_) => "Triangular",
    }
}

fn trace(stage: &str, n: usize, detail: Option<String>) {
    if std::env::var_os("POPKERN_KP13_BI5_TRACE").is_some() {
        match detail {
            Some(d) => eprintln!("[kp13_bi5 N={n}] {stage}: {d}"),
            None => eprintln!("[kp13_bi5 N={n}] {stage}"),
        }
    }
}

#[cfg(test)]
#[path = "kp13_bi5_promote_tests.rs"]
mod tests;
