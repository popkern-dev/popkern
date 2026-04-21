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
//! Set env var `POPKERN_KP_POLAR_TRACE=1` before a fillet run to
//! stream the reason on every fallback. Silent by default.

use popkern_core::geometry::gregory::{
    GregoryPatch, PolarBoundaryInput, RectangularGregoryBiNonic, RectangularGregoryTwoStripe, Side,
    Spoke, SpokeEnd,
};
use popkern_core::geometry::surface::SurfaceGeometry;

use super::kp_polar_bi6_apply::try_kp_polar_bi6_cap;
use super::three_quad_apply::{extract_direct_boundary_input, worst_ring_coplanarity};
use crate::fillet_surface::FilletSurfaceResult;

/// Commercial-grade gate on worst-case G1 coplanarity across the
/// ring's shared spokes. Output above this threshold is still
/// committed (downstream depends on the cap existing) but flagged
/// via [`CornerPatchInfo::relaxed_continuity_coplanarity`] so
/// callers can emit a `RelaxedContinuity` warning.
///
/// [`CornerPatchInfo::relaxed_continuity_coplanarity`]: crate::types::CornerPatchInfo::relaxed_continuity_coplanarity
const COMMERCIAL_GATE: f64 = 5e-2;

/// Self-intersection tolerance for the fold-revert gate.
const SELF_INT_TOL: f64 = 1e-6;

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
    Rejected,
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
    if n < 5 {
        trace_polar("n < 5", n, None);
        return PolarPromoteOutcome::Rejected;
    }

    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for (idx, patch) in sub_patches.iter().enumerate() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            trace_polar("non-two-stripe input", n, Some(format!("sub-patch {idx}")));
            return PolarPromoteOutcome::Rejected;
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

    let Some(mut direct) = extract_direct_boundary_input(&two_stripe, &spokes) else {
        trace_polar("extract_direct_boundary_input → None", n, None);
        return PolarPromoteOutcome::Rejected;
    };

    // I.10 polish: override the two-stripe-seed placeholder radii
    // with real values pulled from the analytic fillet surfaces.
    // Spoke `k` straddles fillets `k` and `(k+1) % n`; the spoke
    // radius is their mean. NURBS / Gregory fillets return `None`
    // from the extractor — in that case retain the upstream
    // placeholder (no information loss vs. prior behaviour).
    if fillet_surfaces.len() >= n {
        for (k, spoke) in direct.spokes.iter_mut().enumerate() {
            let r_k = extract_fillet_radius(&fillet_surfaces[k]);
            let r_kp1 = extract_fillet_radius(&fillet_surfaces[(k + 1) % n]);
            if let (Some(a), Some(b)) = (r_k, r_kp1) {
                let mean = 0.5 * (a + b);
                spoke.radius_at_c = mean;
                spoke.radius_at_v = mean;
            }
        }
    }

    let polar = match PolarBoundaryInput::from_direct_boundary_input(&direct) {
        Ok(p) => p,
        Err(e) => {
            trace_polar("from_direct_boundary_input", n, Some(e));
            return PolarPromoteOutcome::Rejected;
        }
    };
    let bi9_grids = match try_kp_polar_bi6_cap(&polar) {
        Ok(g) => g,
        Err(e) => {
            trace_polar("try_kp_polar_bi6_cap", n, Some(e));
            return PolarPromoteOutcome::Rejected;
        }
    };
    if bi9_grids.len() != n {
        trace_polar(
            "grid count mismatch",
            n,
            Some(format!("got {}", bi9_grids.len())),
        );
        return PolarPromoteOutcome::Rejected;
    }

    let new_patches: Vec<GregoryPatch> = bi9_grids
        .into_iter()
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect();
    for (i, gp) in new_patches.iter().enumerate() {
        if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
            trace_polar("self-intersection", n, Some(format!("sub-patch {i}")));
            return PolarPromoteOutcome::Rejected;
        }
    }

    let worst_coplanarity = worst_ring_coplanarity(&new_patches);

    for (idx, gp) in new_patches.into_iter().enumerate() {
        sub_patches[idx] = gp;
    }
    PolarPromoteOutcome::Committed { worst_coplanarity }
}

/// Dispatch wrapper used by `apply_three_quad_decomposition_to_blend`.
/// Runs the polar promotion and, on success, updates
/// `relaxed_coplanarity` according to the commercial-grade gate:
///
/// - commercial-grade (`worst ≤ COMMERCIAL_GATE`) → `relaxed = None`
/// - relaxed continuity (`worst > COMMERCIAL_GATE`) → `relaxed = Some(worst)`
/// - rejected → caller's `relaxed_coplanarity` preserved (e.g.,
///   binonic's value remains if it set one)
#[allow(dead_code)]
pub(super) fn promote_and_update_relaxed(
    sub_patches: &mut [GregoryPatch],
    fillet_surfaces: &[FilletSurfaceResult],
    relaxed_coplanarity: &mut Option<f64>,
) {
    match promote_decomp_to_kp_polar_bi6(sub_patches, fillet_surfaces) {
        PolarPromoteOutcome::Committed { worst_coplanarity } => {
            *relaxed_coplanarity = if worst_coplanarity > COMMERCIAL_GATE {
                Some(worst_coplanarity)
            } else {
                None
            };
        }
        PolarPromoteOutcome::Rejected => {}
    }
}

/// Emit a diagnostic trace when the polar K-P promotion rejects.
/// Gated on `POPKERN_KP_POLAR_TRACE` so production builds stay
/// silent; set `=1` to see which stage rejected and why.
fn trace_polar(stage: &str, n: usize, detail: Option<String>) {
    if std::env::var_os("POPKERN_KP_POLAR_TRACE").is_some() {
        match detail {
            Some(d) => eprintln!("[kp_polar_bi6 N={n}] fallback @ {stage}: {d}"),
            None => eprintln!("[kp_polar_bi6 N={n}] fallback @ {stage}"),
        }
    }
}

#[cfg(test)]
#[path = "kp_polar_bi6_promote_tests.rs"]
mod tests;
