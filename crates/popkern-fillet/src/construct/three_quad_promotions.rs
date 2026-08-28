//! Degree-promotion ladder for three-quad sub-patches.
//!
//! Extracted from [`super::three_quad_apply`] to keep both files
//! under the 500-line structural limit. Will be reorganised
//! further in the bi-9 cap campaign.6 (kp_build orchestration) when the K-P
//! bidegree-9 path joins the ladder.
//!
//! The ladder promotes each sub-patch through the sequence
//!
//! ```text
//!   deg5 → biquintic-full → biseptic → two-stripe
//! ```
//!
//! with a spoke-G1 solver pass at every rung. Each step is
//! guarded: fold-detection + residual-gated revert keep the
//! pipeline from producing worse geometry than the preceding
//! rung.
//!
//! These symbols remain `pub(super)` — they are only called by
//! [`super::three_quad_apply::apply_three_quad_decomposition_to_blend`].

use popkern_core::geometry::gregory::{
    GregoryPatch, RectangularGregoryBiquinticFull, RectangularGregoryBiseptic,
    RectangularGregoryTwoStripe, Side, Spoke, SpokeEnd, SpokeG1Patch, spoke_g1_solve,
    spoke_g1_solve_with_corner_samples_damped,
};

/// Variant-independent G1 "misalignment" score: sum of
/// `1 - |n_A · n_B|` at each spoke sample, where `n_A`, `n_B`
/// are the UNIT surface normals on either side of the shared
/// spoke. `0.0` means all samples in G1; larger means worse.
/// Generic over patch type so pre- and post-promotion states
/// can be rated on the same scale — crucial for the residual-
/// gated revert logic (pure `n_target · (xtan_A − xtan_B)`
/// isn't comparable across variants because `n_target` is
/// recomputed from each state's xtan_A).
fn spoke_g1_misalignment<P: SpokeG1Patch>(
    patches: &[P],
    spokes: &[Spoke],
    samples_per_spoke: usize,
) -> f64 {
    spoke_g1_misalignment_ext(patches, spokes, samples_per_spoke, false)
}

/// Extended version with optional near-endpoint samples at
/// `t ≈ 0.02` and `t ≈ 0.98` — matches the sample set that
/// `spoke_g1_solve_with_corner_samples` optimises, so the
/// two-stripe gate can compare like-for-like.
fn spoke_g1_misalignment_ext<P: SpokeG1Patch>(
    patches: &[P],
    spokes: &[Spoke],
    samples_per_spoke: usize,
    include_corners: bool,
) -> f64 {
    let mut ts: Vec<f64> = (1..=samples_per_spoke)
        .map(|k| (k as f64) / ((samples_per_spoke + 1) as f64))
        .collect();
    if include_corners {
        ts.insert(0, 0.02);
        ts.push(0.98);
    }
    let mut total = 0.0_f64;
    for spoke in spokes {
        for &t in &ts {
            let t_b = if spoke.b.reversed { 1.0 - t } else { t };
            let pa = &patches[spoke.a.patch_idx];
            let pb = &patches[spoke.b.patch_idx];
            let (ua, va) = match spoke.a.side {
                Side::UZero => (0.0, t),
                Side::UOne => (1.0, t),
                Side::VZero => (t, 0.0),
                Side::VOne => (t, 1.0),
            };
            let (ub, vb) = match spoke.b.side {
                Side::UZero => (0.0, t_b),
                Side::UOne => (1.0, t_b),
                Side::VZero => (t_b, 0.0),
                Side::VOne => (t_b, 1.0),
            };
            let Some(na) = pa
                .partial_u(ua, va)
                .cross(&pa.partial_v(ua, va))
                .normalized()
            else {
                continue;
            };
            let Some(nb) = pb
                .partial_u(ub, vb)
                .cross(&pb.partial_v(ub, vb))
                .normalized()
            else {
                continue;
            };
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            total += 1.0 - cos;
        }
    }
    total
}

/// Promote each deg5 sub-patch in `sub_patches` to
/// [`RectangularGregoryBiquinticFull`] via `from_deg5`, then run
/// the Karciauskas-Peters spoke-G1 joint solver over the ring.
pub(super) fn promote_decomp_to_biquintic_full_with_spoke_g1(sub_patches: &mut [GregoryPatch]) {
    let n = sub_patches.len();
    if n == 0 {
        return;
    }
    let mut biquintic: Vec<RectangularGregoryBiquinticFull> = Vec::with_capacity(n);
    for patch in sub_patches.iter() {
        let GregoryPatch::RectangularDeg5(deg5) = patch else {
            return;
        };
        biquintic.push(RectangularGregoryBiquinticFull::from_deg5(deg5));
    }
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
    let pre_solve = biquintic.clone();
    let pre_solve_folded = pre_solve.iter().any(|p| {
        GregoryPatch::RectangularBiquinticFull(Box::new(p.clone()))
            .self_intersects_with_tolerance(1e-6)
    });
    if pre_solve_folded {
        return;
    }
    let _ = spoke_g1_solve(&mut biquintic, &spokes, 3);
    let any_folded = biquintic.iter().any(|p| {
        GregoryPatch::RectangularBiquinticFull(Box::new(p.clone()))
            .self_intersects_with_tolerance(1e-6)
    });
    if any_folded {
        biquintic = pre_solve;
    }
    for (idx, patch) in biquintic.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularBiquinticFull(Box::new(patch));
    }
}

/// Promote biquintic-full sub-patches to biseptic and run the
/// degree-7 spoke-G1 solver. Returns `true` on success.
pub(super) fn promote_decomp_to_biseptic_with_spoke_g1(sub_patches: &mut [GregoryPatch]) -> bool {
    let n = sub_patches.len();
    if n == 0 {
        return false;
    }
    let biquintic_pre: Vec<RectangularGregoryBiquinticFull> = sub_patches
        .iter()
        .filter_map(|p| {
            if let GregoryPatch::RectangularBiquinticFull(b) = p {
                Some((**b).clone())
            } else {
                None
            }
        })
        .collect();
    if biquintic_pre.len() != n {
        return false;
    }
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
    let biquintic_residual = spoke_g1_misalignment(&biquintic_pre, &spokes, 5);

    // Degree-elevate biquintic → biseptic (pointwise-identical,
    // same residual as biquintic). This is the unsolved baseline
    // we commit if the solver either folds or regresses — keeping
    // the ladder monotonic in quality while still making the
    // biseptic-typed grid available to downstream rungs. Without
    // this, a single rung's solver-regression short-circuits the
    // entire binonic/K-P path (observed on cube_extreme_20).
    let biseptic_pre: Vec<RectangularGregoryBiseptic> = biquintic_pre
        .iter()
        .map(RectangularGregoryBiseptic::from_biquintic_full)
        .collect();

    let mut biseptic = biseptic_pre.clone();
    let _ = spoke_g1_solve(&mut biseptic, &spokes, 5);
    let any_folded = biseptic.iter().any(|p| {
        GregoryPatch::RectangularBiseptic(Box::new(p.clone())).self_intersects_with_tolerance(1e-6)
    });
    let biseptic_post_residual = if any_folded {
        f64::INFINITY
    } else {
        spoke_g1_misalignment(&biseptic, &spokes, 5)
    };
    // Commit whichever biseptic has lower residual: pre-solve
    // (== biquintic) or post-solve. Both start from the degree-
    // elevated biquintic, so the pre-solve option is always at
    // least as good as biquintic, never worse.
    let committed = if biseptic_post_residual <= biquintic_residual {
        biseptic
    } else {
        biseptic_pre
    };
    for (idx, patch) in committed.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularBiseptic(Box::new(patch));
    }
    true
}

/// Promote each biseptic sub-patch in `sub_patches`
/// to [`RectangularGregoryTwoStripe`] via `from_biseptic` and
/// re-run the spoke-G1 solver at 6-DOF-per-spoke.
pub(super) fn promote_decomp_to_two_stripe_with_spoke_g1(sub_patches: &mut [GregoryPatch]) -> bool {
    let n = sub_patches.len();
    if n == 0 {
        return false;
    }
    let biseptic_pre: Vec<RectangularGregoryBiseptic> = sub_patches
        .iter()
        .filter_map(|p| {
            if let GregoryPatch::RectangularBiseptic(b) = p {
                Some((**b).clone())
            } else {
                None
            }
        })
        .collect();
    if biseptic_pre.len() != n {
        return false;
    }
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = biseptic_pre
        .iter()
        .map(RectangularGregoryTwoStripe::from_biseptic)
        .collect();
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

    // Sampled-damped adaptive iteration (primary
    // fallback after the abandoned polynomial-identity direct-formula
    // attempt — see git history of this file for the historical
    // comment block). 6 interior samples + corner samples at
    // t ≈ 0.02 / 0.98 is the measured sweet spot for the QP.
    let samples = 6usize;
    let mut best = two_stripe.clone();
    let mut best_residual = spoke_g1_misalignment_ext(&two_stripe, &spokes, samples, true);
    let warmup_interior = 3usize;
    for _ in 0..warmup_interior {
        let snapshot = two_stripe.clone();
        let _ = spoke_g1_solve(&mut two_stripe, &spokes, samples);
        let any_folded = two_stripe.iter().any(|p| {
            GregoryPatch::RectangularTwoStripe(Box::new(p.clone()))
                .self_intersects_with_tolerance(1e-6)
        });
        if any_folded {
            two_stripe = snapshot;
            break;
        }
        let r = spoke_g1_misalignment_ext(&two_stripe, &spokes, samples, true);
        if r < best_residual {
            best_residual = r;
            best = two_stripe.clone();
        }
    }
    let mut damping = 1.0_f64;
    let max_corner_iters = 120usize;
    let mut stagnation = 0usize;
    let mut prev_best = best_residual;
    for _ in 0..max_corner_iters {
        let snapshot = two_stripe.clone();
        let _ =
            spoke_g1_solve_with_corner_samples_damped(&mut two_stripe, &spokes, samples, damping);
        let any_folded = two_stripe.iter().any(|p| {
            GregoryPatch::RectangularTwoStripe(Box::new(p.clone()))
                .self_intersects_with_tolerance(1e-6)
        });
        if any_folded {
            two_stripe = snapshot;
            damping *= 0.5;
            if damping < 1e-4 {
                break;
            }
            continue;
        }
        let r_after = spoke_g1_misalignment_ext(&two_stripe, &spokes, samples, true);
        if r_after < best_residual {
            best_residual = r_after;
            best = two_stripe.clone();
        }
        if (prev_best - best_residual) < 1e-10 * prev_best.max(1e-12) {
            stagnation += 1;
            if stagnation >= 8 {
                break;
            }
        } else {
            stagnation = 0;
        }
        prev_best = best_residual;
        damping = (damping * 1.25).min(1.0);
    }
    two_stripe = best;
    // Previously this short-circuited with `return false` when the
    // post-solve two-stripe didn't improve on the biseptic input —
    // which blocked the downstream K-P binonic rung entirely even
    // though binonic+L-M is strong enough to recover from a
    // slightly-worse boundary. Commit unconditionally instead; the
    // binonic rung's own quality guard (`post_worst > 1.2 *
    // pre_worst`, `three_quad_apply.rs`) still catches genuine
    // regressions and reverts to the predecessor.
    for (idx, patch) in two_stripe.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularTwoStripe(Box::new(patch));
    }
    true
}
