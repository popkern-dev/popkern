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
//!    edges) whose control ring is **fitted** to the fillet's real
//!    outer boundary. When the wire supplies the boundary target
//!    (PK-087 stage 6) the fit aims each sector's two outer edges at
//!    the exact bi-3 rows of the real `FilletEnd`/`FaceEdge` segments
//!    (`boundary.segments`, cubic-exact) rather than the crude
//!    spoke-derived proxy the extractor produced (~2.5e-2 off). The
//!    sharing makes the net one C²-consistent B-spline surface, so
//!    `sal` vanishes and every eq (2)/(3) cross-sector constraint
//!    closes **by construction**.
//! 2. [`assemble_bi6_polar_cap_from_net`] populates every determined
//!    deep-interior BB-coefficient through the oracle-gated KP15 port
//!    and degree-elevates to bi-9. The assembled grid is then u-flipped
//!    (`grid'[i][j] = grid[9-i][j]`) into the two-stripe sub-patch slot
//!    frame: apex at `grid[9][0]`, outer boundary on `grid[0][*]`
//!    (u=0, `FilletEnd`) and `grid[*][9]` (v=1, `FaceEdge`).
//!
//! # Fit regularisation ladder and the remaining fold
//!
//! An *exact* boundary interpolation moves interior control points far
//! enough to fold the deep interior — the B-spline boundary corner
//! `b[0][0]` is a **16/36-weighted average dominated by the interior
//! symbol `a5`**, so inverting the tensor-border shrinkage displaces
//! `a5` by ~36/16× the corner error (larger than the cap scale). The
//! fit is therefore Tikhonov-regularised: we sweep `λ` small (near-
//! exact boundary) → large (near the coarse fold-free seed) and accept
//! the **first** cap that is both fold-free *and* within the whole-edge
//! boundary gate.
//!
//! On the real `hexagon6_variable` apex this is a **measured structural
//! trade-off** (PK-087 stage 6): fold-free requires `λ ≥ 1e0`, giving a
//! whole-edge residual ~1.3e-2; the residual reaches ~5e-3 only in the
//! folding regime (`λ ≤ 1e-2`); the fold onset sits between `λ` 3e-1
//! (folds) and 1e0 (fold-free). So the cap **declines** here — an
//! honest, precise decline naming the fold-free residual, a strict
//! improvement over stage 5's 2.5e-2-skewed target. Closing the gate
//! needs fold containment (relax the free set `P` — orthogonal to the
//! boundary fit — against the fitted border via the existing F_κ/LM
//! machinery), the named follow-on.
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

/// Cap-outer G0 gate (metres): the constructed cap's two outer
/// boundary edges must lie within this whole-edge distance of the real
/// fillet boundary curves. Mirrors both the user-visible `5e-3`
/// cap-outer acceptance the rational-Gregory dispatch enforces and the
/// downstream `PROMOTED_BOUNDARY_RESIDUAL_GATE` in `three_quad_apply`.
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
const FIT_LAMBDA_LADDER: [f64; 8] = [1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1e0];

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

/// Evaluate a cubic Bézier (4 BB-CPs) at parameter `t`.
fn cubic_at(cps: &[popkern_core::geometry::Point3; 4], t: f64) -> popkern_core::geometry::Point3 {
    let u = 1.0 - t;
    let b = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
    popkern_core::geometry::Point3::new(
        cps.iter().zip(b).map(|(p, w)| p.x * w).sum(),
        cps.iter().zip(b).map(|(p, w)| p.y * w).sum(),
        cps.iter().zip(b).map(|(p, w)| p.z * w).sum(),
    )
}

/// Worst whole-edge distance of the oriented cap's two outer boundary
/// edges to the real fillet boundary curves — the same metric the
/// downstream `PROMOTED_BOUNDARY_RESIDUAL_GATE` enforces (sampled at 11
/// points per edge). The cap's `u = 0` edge (`grid[0][*]`) must match
/// `FilletEnd` (`target[s][1]`), its `v = 1` edge (`grid[*][9]`) must
/// match `FaceEdge` (`target[s][0]`). Orientation of the cubic target is
/// resolved by taking the closer of the two parameterisations, so the
/// metric is a pure geometric distance independent of edge direction.
fn worst_boundary_residual(
    cap: &[[[popkern_core::geometry::Point3; 10]; 10]],
    target: &[[[popkern_core::geometry::Point3; 4]; 2]],
) -> f64 {
    let mut worst = 0.0_f64;
    for (s, grid) in cap.iter().enumerate() {
        for sample in 0..=10 {
            let t = sample as f64 / 10.0;
            let u0_pt = eval_edge(grid, 0, t);
            let v1_pt = eval_edge(grid, 1, t);
            let fe = target[s][1];
            let face = target[s][0];
            let fe_d = dist(u0_pt, cubic_at(&fe, t)).min(dist(u0_pt, cubic_at(&fe, 1.0 - t)));
            let face_d = dist(v1_pt, cubic_at(&face, t)).min(dist(v1_pt, cubic_at(&face, 1.0 - t)));
            worst = worst.max(fe_d).max(face_d);
        }
    }
    worst
}

/// Evaluate the cap's outer boundary edge at parameter `t` by
/// de Casteljau over the bi-9 row/column. `edge = 0` is the `u = 0`
/// edge (`grid[0][j]`), `edge = 1` the `v = 1` edge (`grid[i][9]`).
fn eval_edge(
    grid: &[[popkern_core::geometry::Point3; 10]; 10],
    edge: usize,
    t: f64,
) -> popkern_core::geometry::Point3 {
    let mut work: [popkern_core::geometry::Point3; 10] = if edge == 0 {
        grid[0]
    } else {
        core::array::from_fn(|i| grid[i][9])
    };
    for r in 1..10 {
        for i in 0..(10 - r) {
            let a = work[i];
            let b = work[i + 1];
            work[i] = popkern_core::geometry::Point3::new(
                (1.0 - t) * a.x + t * b.x,
                (1.0 - t) * a.y + t * b.y,
                (1.0 - t) * a.z + t * b.z,
            );
        }
    }
    work[0]
}

/// Corner-only fallback metric used when the fillet wire supplies no
/// explicit boundary target (a non-`2n`-segment boundary). Worst
/// distance of the oriented cap's rim corners `grid[0][9]` (junction)
/// and `grid[9][9]` to the nearest supplied outer corner — a coarse
/// proxy that keeps the acceptance strict on inputs where the true
/// whole-edge residual cannot be measured.
fn worst_cap_outer(
    cap: &[[[popkern_core::geometry::Point3; 10]; 10]],
    polar: &PolarBoundaryInput,
) -> f64 {
    let mut worst = 0.0_f64;
    for grid in cap.iter() {
        let mut nearest_09 = f64::INFINITY;
        let mut nearest_99 = f64::INFINITY;
        for oc in &polar.outer_corners {
            nearest_09 = nearest_09.min(dist(grid[0][9], oc.position));
            nearest_99 = nearest_99.min(dist(grid[9][9], oc.position));
        }
        worst = worst.max(nearest_09).max(nearest_99);
    }
    worst
}

/// Attempt the KP15 shared-net cap for `polar`, sweeping the fit
/// regularisation ladder. See the module docs for the acceptance
/// contract.
pub(super) fn try_kp15_shared_net_cap(
    polar: &PolarBoundaryInput,
    boundary_target: Option<&[[[popkern_core::geometry::Point3; 4]; 2]]>,
) -> SharedNetOutcome {
    let n = polar.n();
    // Best fold-free candidate seen so far (smallest cap_outer), for a
    // measured decline diagnostic if no rung clears the gate.
    let mut best_fold_free: Option<(f64, f64)> = None; // (cap_outer, lambda)
    let mut last_fit_err: Option<String> = None;
    let mut last_assembly_err: Option<String> = None;

    for &lambda in &FIT_LAMBDA_LADDER {
        let (params, net, free) =
            match convert_polar_to_bi6_cap_fitted(polar, lambda, boundary_target) {
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
        // Orientation: the assembled bi-9 grid places the apex at
        // `grid[0][0]` (script `[6][6]`) with the two tensor-border outer
        // edges on `grid[*][9]` (tad) and `grid[9][*]` (tak). The
        // two-stripe sub-patch slot it replaces expects the apex at
        // `grid[9][0]` with the outer boundary on `grid[0][*]` (u=0,
        // `FilletEnd`) and `grid[*][9]` (v=1, `FaceEdge`). The u-flip
        // `grid'[i][j] = grid[9-i][j]` re-parametrises the cap into the
        // slot's frame (PK-087 stage-6 orientation derivation).
        let oriented: Vec<[[popkern_core::geometry::Point3; 10]; 10]> = cap
            .patches_bi9
            .iter()
            .map(|grid| core::array::from_fn(|i| grid[9 - i]))
            .collect();
        let candidate: Vec<GregoryPatch> = oriented
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
        // Boundary fidelity: whole-edge residual against the real fillet
        // boundary when the wire supplied it (the honest metric that
        // matches the downstream `PROMOTED_BOUNDARY_RESIDUAL_GATE`);
        // otherwise fall back to the corner metric of the coarse frame.
        let cap_outer = match boundary_target {
            Some(target) => worst_boundary_residual(&oriented, target),
            None => worst_cap_outer(&oriented, polar),
        };
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
    // to the guarded F_κ ladder). PK-087 stage 6: with the real boundary
    // target the fold-free residual floors at ~1.3e-2 (fold-free only for
    // λ ≥ 1e0); residual ≤ 5e-3 is reachable only in the folding regime
    // (λ ≤ 1e-2). The structural fold — the B-spline corner b[0][0] is
    // a5-dominated (16/36) so inverting the tensor-border shrinkage moves
    // the interior symbol far enough to self-intersect — is the sole
    // remaining blocker; see the ADR stage-6 pin.
    if let Some((cap_outer, lambda)) = best_fold_free {
        SharedNetOutcome::Declined {
            reason: "shared_net_cap_outer",
            detail: format!(
                "best fold-free boundary residual={cap_outer:.3e} (λ={lambda:.0e}) exceeds \
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
