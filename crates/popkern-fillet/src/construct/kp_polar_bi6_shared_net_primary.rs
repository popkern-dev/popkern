//! PK-087 stage 4 — the KP15 shared-net cap as the **primary**
//! candidate of the polar bi-6 backstop.
//!
//! The legacy `promote_decomp_to_kp_polar_bi6` path drove a single
//! per-sector-rim `Bi6CapInput` into the F_κ minimisation ladder,
//! which self-intersects (folds) on every candidate at N ≥ 5 — the
//! measured gate-decline root (PK-087 stage 3). This module runs the
//! **correct-by-construction** replacement first:
//!
//! 1. [`convert_polar_to_bi6_cap_net`] re-orients the polar boundary
//!    into a **shared C²-consistent** bi-6 cap net (spokes on the
//!    sector-glue edges, tensor border on the star-opposite rim
//!    edges). Because the sharing *is* the statement that the input
//!    net is one C²-consistent B-spline surface, the eq (7)/(8) ring
//!    residual `sal` vanishes and every eq (2)/(3) cross-sector
//!    constraint closes **by construction** — the cap is fold-free
//!    where the F_κ ladder folds.
//! 2. [`assemble_bi6_polar_cap_from_net`] populates every determined
//!    deep-interior BB-coefficient through the oracle-gated KP15
//!    port and degree-elevates to bi-9.
//!
//! # Acceptance is all-or-nothing
//!
//! The cap is committed **only** when it is fold-free *and* its outer
//! boundary lies within the cap-outer gate of the real surrounding
//! surface. The coarse B-spline cap net's boundary is pulled inward
//! by the B-spline → Bézier averaging (the shrinkage measured in the
//! PK-087 ADR: ~0.5 unit on a unit fixture, ≫ 5 mm), so this
//! candidate **declines** on the cap-outer gate today and records a
//! structured diagnostic naming `cap_outer_g0`. On decline the model
//! is untouched and the guarded F_κ ladder runs as the fallback.
//!
//! The commercial-grade replacement (a C²-consistent net *fitted* to
//! the real boundary, so both continuity and boundary fidelity hold)
//! is recorded in
//! `docs/decisions/2026-06-22-kp15-bi6-deep-interior-without-maple.md`.

use popkern_core::geometry::gregory::{
    GregoryPatch, PolarBoundaryInput, RectangularGregoryBiNonic, assemble_bi6_polar_cap_from_net,
    convert_polar_to_bi6_cap_net,
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

/// Attempt the KP15 shared-net cap for `polar`. See the module docs
/// for the all-or-nothing acceptance contract.
pub(super) fn try_kp15_shared_net_cap(polar: &PolarBoundaryInput) -> SharedNetOutcome {
    let (params, net, free) = match convert_polar_to_bi6_cap_net(polar) {
        Ok(triple) => triple,
        Err(e) => {
            return SharedNetOutcome::Declined {
                reason: "shared_net_bridge",
                detail: e,
            };
        }
    };
    let cap = match assemble_bi6_polar_cap_from_net(&params, &net, &free, SAL_TOLERANCE) {
        Ok(cap) => cap,
        Err(e) => {
            return SharedNetOutcome::Declined {
                reason: "shared_net_assembly",
                detail: e,
            };
        }
    };
    let n = polar.n();
    if cap.patches_bi9.len() != n {
        return SharedNetOutcome::Declined {
            reason: "shared_net_grid_count",
            detail: format!("got {} grids, expected n = {n}", cap.patches_bi9.len()),
        };
    }

    let candidate: Vec<GregoryPatch> = cap
        .patches_bi9
        .iter()
        .map(|grid| {
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: *grid }))
        })
        .collect();

    // Fold gate — the shared net should never fold, but keep the
    // guard so a regression declines rather than ships a folded cap.
    for (i, gp) in candidate.iter().enumerate() {
        if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
            return SharedNetOutcome::Declined {
                reason: "shared_net_fold",
                detail: format!("sub_patch={i}"),
            };
        }
    }

    // Cap-outer gate — the cap's outer corner (kernel grid[9][9]) must
    // match the real surrounding-surface corner. The coarse net
    // shrinks this; a genuine C²-consistent-net fit will pass here.
    let mut worst_cap_outer = 0.0_f64;
    for (s, grid) in cap.patches_bi9.iter().enumerate() {
        let outer = grid[9][9];
        let target = polar.outer_corners[s].position;
        let d = ((outer.x - target.x).powi(2)
            + (outer.y - target.y).powi(2)
            + (outer.z - target.z).powi(2))
        .sqrt();
        worst_cap_outer = worst_cap_outer.max(d);
    }
    if worst_cap_outer > CAP_OUTER_GATE {
        return SharedNetOutcome::Declined {
            reason: "shared_net_cap_outer",
            detail: format!(
                "cap_outer_g0={worst_cap_outer:.3e} exceeds {CAP_OUTER_GATE:.3e} \
                 (B-spline cap-net boundary shrinkage; C²-consistent net fit pending)"
            ),
        };
    }

    SharedNetOutcome::Committed(candidate)
}
