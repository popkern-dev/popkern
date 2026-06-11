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
//! when `n_sides >= 5 && options.use_kp_polar_gregory` via the
//! [`promote_and_update_relaxed`] helper, which handles the
//! commercial-grade gate + `relaxed_coplanarity` telemetry
//! plumbing in one call.
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
//! 4. Screen self-intersection + worst-case G1 coplanarity.
//! 5. Wrap each `[[Point3; 10]; 10]` as
//!    [`RectangularGregoryBiNonic`] and commit.
//!
//! # Diagnostics
//!
//! Rejections return structured exact-boundary diagnostics so
//! downstream classification, UI, and tests can explain why the
//! polar rung declined before lower-rung fallback.

use popkern_core::geometry::gregory::{
    GregoryPatch, PolarBoundaryInput, RectangularGregoryBiNonic, RectangularGregoryTwoStripe, Side,
    Spoke, SpokeEnd, convert_polar_to_bi6_cap_input, extract_bi3_tensor_border_from_polar_input,
    f_kappa_solve_optimal_bi6_cap,
};
use popkern_core::geometry::surface::SurfaceGeometry;

use super::kp_polar_bi6_diagnostics::{self, trace_polar};
use super::kp_polar_deep_fill::{DeepFillOutcome, try_boundary_owned_deep_fill};
use super::kp_polar_deep_trace::deep_trace_bi6_output;
use super::kp_polar_lm_optimizer::try_polar_lm_with_report;
use super::kp_polar_weighted_f5::try_weighted_f5_profiles;
use super::three_quad_apply::{extract_direct_boundary_input, worst_ring_coplanarity};
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::ApexConstructionDiagnostic;

/// Commercial-grade gate on worst-case G1 coplanarity across the
/// ring's shared spokes. Output above this threshold is still
/// committed (downstream depends on the cap existing) but flagged
/// via [`CornerPatchInfo::relaxed_continuity_coplanarity`] so
/// callers can emit a `RelaxedContinuity` warning.
///
/// [`CornerPatchInfo::relaxed_continuity_coplanarity`]: crate::types::CornerPatchInfo::relaxed_continuity_coplanarity
const COMMERCIAL_GATE: f64 = 5e-2;

/// Self-intersection tolerance for the fold-revert gate.
pub(super) const SELF_INT_TOL: f64 = 1e-6;

/// Outcome of a polar bi-6 promotion attempt.
pub(super) enum PolarPromoteOutcome {
    /// Sub-patches replaced with bi-nonic grids. `worst_coplanarity`
    /// is the maximum G1 residual on any spoke of the committed
    /// ring. Values above [`COMMERCIAL_GATE`] indicate relaxed
    /// continuity.
    Committed { worst_coplanarity: f64 },
    /// Promotion did not run; sub-patches unchanged. Upstream
    /// fallback applies (existing two-stripe / kp-binonic output
    /// is preserved).
    Rejected {
        diagnostics: Vec<ApexConstructionDiagnostic>,
    },
}

/// Extract the characteristic fillet radius from an analytic
/// [`FilletSurfaceResult`]. Returns `Some(radius)` when the
/// underlying [`SurfaceGeometry`] carries a scalar radius —
/// `Cylinder { radius }`, `Torus { minor_radius }`, or
/// `Sphere { radius }`. Returns `None` for planar, conic,
/// NURBS or Gregory surfaces (those don't expose a single
/// characteristic radius; the caller falls back to the
/// placeholder average baked into the two-stripe seed).
pub(super) fn extract_fillet_radius(fs: &FilletSurfaceResult) -> Option<f64> {
    match &fs.surface {
        SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
        SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
        SurfaceGeometry::Sphere { radius, .. } => Some(*radius),
        _ => None,
    }
}

/// Attempt to promote a ring of two-stripe sub-patches to bi-6
/// polar K-P output. Overwrites `sub_patches` with N
/// `RectangularBiNonic` variants on success.
///
/// `fillet_surfaces` are the per-spoke input fillet strips used
/// to recover real radii for the K-P 2015 coefficient tables.
/// Spoke `k` connects sub-patches `k` and `(k+1) % N` and sits
/// between fillets `k` and `(k+1) % N`; its radius is the mean
/// of those two (analytic cases) or falls through to the two-
/// stripe seed's placeholder for NURBS/Gregory fillets.
///
/// Returns a [`PolarPromoteOutcome`]: [`Committed`] on success
/// (with the worst-case coplanarity for commercial-gate dispatch),
/// or [`Rejected`] on any failure. Set `POPKERN_KP_POLAR_TRACE=1`
/// to see the fallback reason.
///
/// [`Committed`]: PolarPromoteOutcome::Committed
/// [`Rejected`]: PolarPromoteOutcome::Rejected
#[allow(dead_code)]
pub(super) fn promote_decomp_to_kp_polar_bi6(
    sub_patches: &mut [GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
) -> PolarPromoteOutcome {
    let n = sub_patches.len();
    let mut diagnostics = Vec::new();
    if n < 5 {
        trace_polar("n < 5", n, None);
        record_polar_decline(
            &mut diagnostics,
            "under_valence",
            n,
            Some("n < 5".to_string()),
        );
        return PolarPromoteOutcome::Rejected { diagnostics };
    }

    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for (idx, patch) in sub_patches.iter().enumerate() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            trace_polar(
                "non-two-stripe input",
                n,
                Some(format!("sub-patch {idx} is {}", variant_name(patch))),
            );
            record_polar_decline(
                &mut diagnostics,
                "non_two_stripe_input",
                n,
                Some(format!("sub-patch {idx} is {}", variant_name(patch))),
            );
            return PolarPromoteOutcome::Rejected { diagnostics };
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

    // Phase 6c-D Stage 1d: polar bi-6 path doesn't go through
    // the bi-nonic K-P guide pipeline, so hex_boundaries
    // injection isn't needed. Pass `None` for the existing path.
    let Some(mut direct) = extract_direct_boundary_input(&two_stripe, &spokes, None) else {
        trace_polar("extract_direct_boundary_input → None", n, None);
        record_polar_decline(&mut diagnostics, "direct_boundary_unavailable", n, None);
        return PolarPromoteOutcome::Rejected { diagnostics };
    };

    // I.22 — per-sector radius threading. Math spec §6 (polar-kp-
    // math.md) treats each sector as having its own characteristic
    // radius `r_k`; the only place radii enter the K-P 2015
    // construction is through cross-derivative scaling in the
    // tensor-border of sector k, driven by `r_k`. The ρ/β/σ/τ/µ
    // tables are all radius-independent (valence-only).
    //
    // In the data structure, `polar.spokes[k]` feeds sector k's
    // tensor-border (1:1 indexing in
    // `extract_bi3_tensor_border_from_polar_input`). So spoke k's
    // radius fields semantically hold sector k's radius — which
    // for our fillet ring is `fillet[k].radius` (the analytic
    // radius of the fillet strip sector k adjoins).
    //
    // NURBS / Gregory fillets return `None` from the extractor —
    // in that case retain the upstream two-stripe-seed placeholder
    // (no information loss vs. prior behaviour).
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
            trace_polar("from_direct_boundary_input", n, Some(e.clone()));
            record_polar_decline(&mut diagnostics, "polar_boundary_input", n, Some(e));
            return PolarPromoteOutcome::Rejected { diagnostics };
        }
    };

    // Build the `Bi6CapInput` once so we can feed the same state
    // into the F_κ ladder and (on exhaustion) L-M warm-start.
    let tensor_border = match extract_bi3_tensor_border_from_polar_input(&polar) {
        Ok(tb) => tb,
        Err(e) => {
            trace_polar("tensor-border extract", n, Some(e.clone()));
            record_polar_decline(&mut diagnostics, "tensor_border_extract", n, Some(e));
            return PolarPromoteOutcome::Rejected { diagnostics };
        }
    };
    let cap_input = match convert_polar_to_bi6_cap_input(&polar, &tensor_border) {
        Ok(ci) => ci,
        Err(e) => {
            trace_polar("convert_polar_to_bi6_cap_input", n, Some(e.clone()));
            record_polar_decline(&mut diagnostics, "bi6_cap_input", n, Some(e));
            return PolarPromoteOutcome::Rejected { diagnostics };
        }
    };

    // Deep-trace diagnostic: dump bi-6 CP grids of first K=5
    // attempt for detailed fold inspection. Gated on env var.
    if std::env::var_os("POPKERN_KP_POLAR_DEEP_TRACE").is_some()
        && let Ok((_, cap_out)) = f_kappa_solve_optimal_bi6_cap(&cap_input, 5)
    {
        deep_trace_bi6_output(&cap_out.patches_bi6, n);
    }

    // F_κ fallback ladder (f5-tuning-plan.md option C). On
    // self-intersection, retry with κ ∈ {4, 3, 1}. Track the
    // last F_κ solve's free params (even if folded) as L-M warm-
    // start candidate — L-M can refine from a folded seed toward
    // a nearby fold-safe configuration if one exists.
    let kappas_to_try = [5usize, 4, 3, 1];
    let mut committed_patches: Option<Vec<GregoryPatch>> = None;
    let mut last_free_params = None;
    for (attempt_idx, &kappa) in kappas_to_try.iter().enumerate() {
        let (free_params, cap_out) = match f_kappa_solve_optimal_bi6_cap(&cap_input, kappa) {
            Ok(result) => result,
            Err(e) => {
                trace_polar(
                    &format!("F_{kappa} solve error"),
                    n,
                    Some(format!("attempt {attempt_idx}: {e}")),
                );
                continue;
            }
        };
        last_free_params = Some(free_params);
        let bi9_grids = cap_out.patches_bi9;
        if bi9_grids.len() != n {
            trace_polar(
                "grid count mismatch",
                n,
                Some(format!("κ={kappa} got {}", bi9_grids.len())),
            );
            continue;
        }
        let candidate: Vec<GregoryPatch> = bi9_grids
            .into_iter()
            .map(|grid| {
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
            })
            .collect();
        let mut folded = false;
        let mut fold_idx = 0usize;
        for (i, gp) in candidate.iter().enumerate() {
            if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
                folded = true;
                fold_idx = i;
                break;
            }
        }
        if folded {
            trace_polar(
                &format!("F_{kappa} self-intersection"),
                n,
                Some(format!("sub-patch {fold_idx}")),
            );
            record_polar_decline(
                &mut diagnostics,
                "f_kappa_self_intersection",
                n,
                Some(format!("kappa={kappa};sub_patch={fold_idx}")),
            );
            continue;
        }
        if kappa != 5 {
            trace_polar(
                &format!("F_{kappa} accepted"),
                n,
                Some(format!("κ=5 folded; ladder fell to κ={kappa}")),
            );
        }
        committed_patches = Some(candidate);
        break;
    }

    // F_κ ladder exhausted — try option E (per-valence weight
    // tuning) before falling through to L-M (option F).
    if committed_patches.is_none() {
        trace_polar(
            "ladder exhausted",
            n,
            Some("F_κ all folded; trying weighted F_5 profiles (option E)".into()),
        );
        committed_patches = try_weighted_f5_profiles(&cap_input, n, &mut last_free_params);
    }

    // Still nothing — try a boundary-owned deep-interior
    // candidate before the nonlinear L-M fallback. This owns the
    // open U=1 interface from the next sector's row-0 data and
    // declines unless the same fold/G1 gates pass.
    if committed_patches.is_none() {
        trace_polar(
            "weighted F_5 exhausted",
            n,
            Some("attempting boundary-owned deep fill before L-M".into()),
        );
        if let Some(warm_start) = last_free_params.as_ref() {
            match try_boundary_owned_deep_fill(&cap_input, warm_start, COMMERCIAL_GATE) {
                DeepFillOutcome::Accepted(patches) => {
                    trace_polar("boundary-owned deep fill accepted", n, None);
                    committed_patches = Some(patches);
                }
                DeepFillOutcome::Rejected(reason) => {
                    trace_polar("boundary-owned deep fill rejected", n, Some(reason.into()));
                    record_polar_decline(
                        &mut diagnostics,
                        "deep_fill_rejected",
                        n,
                        Some(reason.into()),
                    );
                }
            }
        }
    }

    // Still nothing — fall through to L-M (option F).
    if committed_patches.is_none() {
        trace_polar(
            "boundary-owned deep fill exhausted",
            n,
            Some("attempting L-M refinement (option F)".into()),
        );
        if let Some(warm_start) = last_free_params.as_ref() {
            let lm_report = try_polar_lm_with_report(&cap_input, warm_start);
            if let Some(lm_patches) = lm_report.patches {
                trace_polar("L-M accepted", n, Some("escape via option F".into()));
                let lm_wrapped: Vec<GregoryPatch> = lm_patches
                    .into_iter()
                    .map(|grid| {
                        GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
                            grid,
                        }))
                    })
                    .collect();
                committed_patches = Some(lm_wrapped);
            } else {
                let reason = lm_report.rejection_reason.unwrap_or("unknown");
                trace_polar("L-M rejected", n, Some(reason.to_string()));
                record_polar_decline(
                    &mut diagnostics,
                    "lm_rejected",
                    n,
                    Some(format!(
                        "F_kappa and weighted F5 candidates folded; L-M reason={reason}"
                    )),
                );
            }
        } else {
            trace_polar("L-M rejected", n, Some("missing warm start".into()));
            record_polar_decline(
                &mut diagnostics,
                "lm_rejected",
                n,
                Some("F_kappa and weighted F5 candidates folded; missing warm start".into()),
            );
        }
    }

    let Some(new_patches) = committed_patches else {
        return PolarPromoteOutcome::Rejected { diagnostics };
    };

    let worst_coplanarity = worst_ring_coplanarity(&new_patches);

    for (idx, gp) in new_patches.into_iter().enumerate() {
        sub_patches[idx] = gp;
    }
    PolarPromoteOutcome::Committed { worst_coplanarity }
}

/// Dispatch wrapper used by `apply_three_quad_decomposition_to_blend`.
/// Runs the polar promotion and, on commit, updates
/// `relaxed_coplanarity` according to the commercial-grade gate:
///
/// - commercial-grade (`worst ≤ COMMERCIAL_GATE`) → `relaxed = None`
/// - relaxed continuity (`worst > COMMERCIAL_GATE`) → `relaxed = Some(worst)`
/// - rejected → caller's `relaxed_coplanarity` preserved
///
/// Returns `true` if the polar path committed bi-9 patches
/// (sub_patches now contain `RectangularBiNonic`), `false` if
/// polar rejected (sub_patches unchanged). The caller uses this
/// to decide whether to run the binonic fallback — at N ≥ 5
/// polar is the primary rung, binonic only fires when polar
/// rejects.
#[allow(dead_code)]
pub(super) fn promote_and_update_relaxed(
    sub_patches: &mut [GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
    relaxed_coplanarity: &mut Option<f64>,
) -> bool {
    promote_and_update_relaxed_with_diagnostics(sub_patches, fillet_surfaces, relaxed_coplanarity).0
}

#[allow(dead_code)]
pub(super) fn promote_and_update_relaxed_with_diagnostics(
    sub_patches: &mut [GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
    relaxed_coplanarity: &mut Option<f64>,
) -> (bool, Vec<ApexConstructionDiagnostic>) {
    match promote_decomp_to_kp_polar_bi6(sub_patches, fillet_surfaces) {
        PolarPromoteOutcome::Committed { worst_coplanarity } => {
            *relaxed_coplanarity = if worst_coplanarity > COMMERCIAL_GATE {
                Some(worst_coplanarity)
            } else {
                None
            };
            (true, Vec::new())
        }
        PolarPromoteOutcome::Rejected { diagnostics } => (false, diagnostics),
    }
}

fn record_polar_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    reason_code: &str,
    n_sides: usize,
    detail: Option<String>,
) {
    diagnostics.push(kp_polar_bi6_diagnostics::promotion_declined(
        reason_code,
        n_sides,
        detail,
        Some(SELF_INT_TOL),
    ));
}

/// Name of a `GregoryPatch` variant for diagnostic traces.
fn variant_name(p: &GregoryPatch) -> &'static str {
    match p {
        GregoryPatch::Rectangular(_) => "Rectangular",
        GregoryPatch::RectangularDeg5(_) => "RectangularDeg5",
        GregoryPatch::RectangularBiquinticFull(_) => "RectangularBiquinticFull",
        GregoryPatch::RectangularBiseptic(_) => "RectangularBiseptic",
        GregoryPatch::RectangularTwoStripe(_) => "RectangularTwoStripe",
        GregoryPatch::RectangularBiNonic(_) => "RectangularBiNonic",
        GregoryPatch::RectangularBiNonicDualTwist(_) => "RectangularBiNonicDualTwist",
        GregoryPatch::Triangular(_) => "Triangular",
    }
}

#[cfg(test)]
#[path = "kp_polar_bi6_promote_tests.rs"]
mod tests;
