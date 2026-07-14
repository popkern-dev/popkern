//! PK-087 — the KP15 shared-net cap as the **primary** candidate of
//! the polar bi-6 backstop.
//!
//! The legacy `promote_decomp_to_kp_polar_bi6` path drove a single
//! per-sector-rim `Bi6CapInput` into the F_κ minimisation ladder,
//! which self-intersects (folds) on every candidate at N ≥ 5 — the
//! measured gate-decline root (PK-087 stage 3). This module runs the
//! **correct-by-construction** replacement first:
//!
//! 1. [`convert_polar_to_bi6_cap_fitted`] re-orients the polar
//!    boundary into a **shared C²-consistent** bi-6 cap net (spokes on
//!    the sector-glue edges, tensor border on the star-opposite rim
//!    edges) whose control ring is **fitted** to the real fillet
//!    tensor border. The sharing makes the net one C²-consistent
//!    B-spline surface, so `sal` vanishes and every eq (2)/(3)
//!    cross-sector constraint closes **by construction** — the cap is
//!    fold-free where the F_κ ladder folds — while the fit inverts the
//!    B-spline → Bézier shrinkage so the cap's outer boundary
//!    interpolates the surrounding surface (PK-087 stage 5).
//! 2. [`assemble_bi6_polar_cap_from_net`] populates every determined
//!    deep-interior BB-coefficient through the oracle-gated KP15 port
//!    and degree-elevates to bi-9.
//!
//! # Fit regularisation ladder
//!
//! An *exact* boundary interpolation moves interior control points far
//! enough to fold the deep interior (the B-spline boundary corner is
//! dominated by the interior symbol `a5`). The fit is therefore
//! Tikhonov-regularised: we sweep the weight `λ` from small (near-exact
//! boundary, larger interior deviation) to large (near the coarse
//! fold-free seed) and accept the **first** cap that is both fold-free
//! *and* within the cap-outer gate — the strongest boundary fidelity
//! the geometry admits without folding.
//!
//! # Acceptance is all-or-nothing
//!
//! On decline the model is untouched and the guarded F_κ ladder runs
//! as the fallback, so the model is never left holding a folded or
//! boundary-wrong cap.

use popkern_core::geometry::gregory::{
    GregoryPatch, PolarBoundaryInput, RectangularGregoryBiNonic, assemble_bi6_polar_cap_from_net,
    convert_polar_to_bi6_cap_fitted,
};

use super::SELF_INT_TOL;

/// Cap-outer G0 gate (metres): the constructed cap's outer corners
/// must lie within this distance of the real surrounding-surface
/// corners. Mirrors the user-visible `5e-3` cap-outer acceptance the
/// rational-Gregory dispatch enforces.
pub(super) const CAP_OUTER_GATE: f64 = 5e-3;

/// Deferred eq (8) ring-closure tolerance for the shared net. The
/// shared net is C²-consistent by construction, so `sal` is machine
/// zero; a violation is a genuine bridge defect, not input noise.
const SAL_TOLERANCE: f64 = 1e-6;

/// Tikhonov regularisation ladder for the boundary fit, small (near-
/// exact boundary) to large (near the coarse fold-free seed). The
/// first rung that yields a fold-free cap within the cap-outer gate
/// wins; if none does, the smallest-`cap_outer` fold-free rung is
/// reported in the decline diagnostic.
const FIT_LAMBDA_LADDER: [f64; 6] = [3e-2, 1e-1, 3e-1, 1e0, 3e0, 1e1];

/// Outcome of the shared-net primary candidate.
pub(super) enum SharedNetOutcome {
    /// Fold-free cap whose outer boundary is within the gate — commit.
    Committed(Vec<GregoryPatch>),
    /// Declined; `reason` is the structured decline code and `detail`
    /// carries the measured evidence. The caller records the decline
    /// and falls through to the guarded F_κ ladder.
    Declined {
        reason: &'static str,
        detail: String,
    },
}

/// Euclidean distance between two points.
fn dist(a: popkern_core::geometry::Point3, b: popkern_core::geometry::Point3) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

/// Worst distance of the cap's outer-boundary corners to the real
/// surrounding-surface corners. Under the KP15 script→kernel
/// re-indexing the sector's two star-opposite outer edges terminate at
/// the rim corners `V_s` (kernel `grid[0][9]`, script `[6][0]`) and
/// `V_{s-1}` (kernel `grid[9][0]`, script `[0][6]`); `grid[9][9]` is
/// the shared apex `C`, not a rim corner.
fn worst_cap_outer(
    cap: &[[[popkern_core::geometry::Point3; 10]; 10]],
    polar: &PolarBoundaryInput,
) -> f64 {
    let n = polar.n();
    let mut worst = 0.0_f64;
    for (s, grid) in cap.iter().enumerate() {
        worst = worst.max(dist(grid[0][9], polar.outer_corners[s].position));
        worst = worst.max(dist(
            grid[9][0],
            polar.outer_corners[(s + n - 1) % n].position,
        ));
    }
    worst
}

/// Attempt the KP15 shared-net cap for `polar`, sweeping the fit
/// regularisation ladder. See the module docs for the acceptance
/// contract.
pub(super) fn try_kp15_shared_net_cap(polar: &PolarBoundaryInput) -> SharedNetOutcome {
    let n = polar.n();
    // Best fold-free candidate seen so far (smallest cap_outer), for a
    // measured decline diagnostic if no rung clears the gate.
    let mut best_fold_free: Option<(f64, f64)> = None; // (cap_outer, lambda)
    let mut last_fit_err: Option<String> = None;
    let mut last_assembly_err: Option<String> = None;

    for &lambda in &FIT_LAMBDA_LADDER {
        let (params, net, free) = match convert_polar_to_bi6_cap_fitted(polar, lambda) {
            Ok(triple) => triple,
            Err(e) => {
                last_fit_err = Some(e);
                continue;
            }
        };
        let cap = match assemble_bi6_polar_cap_from_net(&params, &net, &free, SAL_TOLERANCE) {
            Ok(cap) => cap,
            Err(e) => {
                last_assembly_err = Some(e);
                continue;
            }
        };
        if cap.patches_bi9.len() != n {
            last_assembly_err = Some(format!(
                "got {} grids, expected n = {n}",
                cap.patches_bi9.len()
            ));
            continue;
        }
        let candidate: Vec<GregoryPatch> = cap
            .patches_bi9
            .iter()
            .map(|grid| {
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
                    grid: *grid,
                }))
            })
            .collect();
        let folded = candidate
            .iter()
            .any(|gp| gp.self_intersects_with_tolerance(SELF_INT_TOL));
        let cap_outer = worst_cap_outer(&cap.patches_bi9, polar);
        if folded {
            continue;
        }
        if best_fold_free.is_none_or(|(best, _)| cap_outer < best) {
            best_fold_free = Some((cap_outer, lambda));
        }
        if cap_outer <= CAP_OUTER_GATE {
            return SharedNetOutcome::Committed(candidate);
        }
    }

    // No rung cleared both gates — decline with the strongest measured
    // evidence available (mechanism honesty; the caller falls through
    // to the guarded F_κ ladder).
    if let Some((cap_outer, lambda)) = best_fold_free {
        SharedNetOutcome::Declined {
            reason: "shared_net_cap_outer",
            detail: format!(
                "best fold-free cap_outer_g0={cap_outer:.3e} (λ={lambda:.0e}) exceeds \
                 {CAP_OUTER_GATE:.3e}: boundary fit cannot clear the gate without folding \
                 on this apex"
            ),
        }
    } else if let Some(e) = last_assembly_err {
        SharedNetOutcome::Declined {
            reason: "shared_net_assembly",
            detail: e,
        }
    } else {
        SharedNetOutcome::Declined {
            reason: "shared_net_fold",
            detail: format!(
                "every fit rung folded; last fit error: {}",
                last_fit_err.unwrap_or_else(|| "none".to_string())
            ),
        }
    }
}
