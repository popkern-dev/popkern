//! `n = 4` gates for the PnS bi-4 cap, plus the division-free eq (2)
//! G¹ gate that covers the whole supported valence range.
//!
//! Split out of `pns_cap_assemble_tests.rs` to stay under the 500-line
//! structural limit. These belong together: every one of them exists
//! because `n = 4` was long believed impossible for this cap, and each
//! pins a different part of why that was wrong — the σ pole is
//! removable, the eq (2) pole is presentational, and the solved form
//! must decline rather than fabricate a row where the system is
//! underdetermined. See
//! `docs/decisions/2026-07-25-pns-n4-removable-pole.md`.

use super::super::pns_beta::BetaReparam;
use super::super::pns_sigma::sigma_center_jet_default;
use super::*;

use super::tests::{approx_pt, planar_ring, pt};

/// `n = 4` used to be rejected outright on the belief that σ's `1/c`
/// factor was structurally singular there. It is not: regrouped as
/// `d₁ + (4/(3c))·(cos(α/2)·d₂ − d₁)`, the residue vanishes exactly
/// under the Appendix B defaults, so the pole is removable and the
/// cap assembles like any other valence. This is the four-valent
/// corner — pyramid apexes and cross junctions — so the old guard cost
/// the bi-4 cap a whole class of real geometry.
#[test]
fn n4_assembles_with_paper_defaults() {
    let input = CapAssembleInput {
        boundary: planar_ring(4, 1.0),
        beta: BetaReparam::bi4(4),
        sigma: sigma_center_jet_default(4).expect("paper defaults yield a finite σ jet"),
        pole: pt(0.0, 0.0, 0.5),
        kappa: 3,
    };
    let result = assemble_pns_bi4_cap(&input).expect("n=4 must assemble; its σ pole is removable");
    assert_eq!(result.patches.len(), 4, "n=4 cap must have 4 sectors");
    for (k, patch) in result.patches.iter().enumerate() {
        assert!(
            approx_pt(patch[0][0], pt(0.0, 0.0, 0.5), 1e-12),
            "n=4 sector {k} must meet the pole"
        );
        for row in patch.iter() {
            for cp in row.iter() {
                assert!(
                    cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite(),
                    "n=4 sector {k} produced a non-finite control point {cp:?}"
                );
            }
        }
    }
}

/// The removable pole is a property of the *scalings*, not of `n = 4`
/// itself. A caller supplying `(d₁, d₂)` that leave a residue gets a
/// typed error naming the measured residue, not a NaN cap.
#[test]
fn n4_irremovable_pole_is_a_typed_error() {
    use super::super::pns_sigma::{SigmaJetError, d1_def, d2_def, sigma_center_jet};
    // Perturb d₂ off the cancelling value; d₁ is left at the default.
    let bad = sigma_center_jet(4, d1_def(4), d2_def(4) * 1.05);
    match bad {
        Err(SigmaJetError::IrremovablePole { residue, .. }) => {
            assert!(
                residue.abs() > 1e-3,
                "a 5% d₂ perturbation should leave a clearly non-zero residue, got {residue:e}"
            );
        }
        other => panic!("expected IrremovablePole for off-default n=4 scalings, got {other:?}"),
    }
}

/// Assembly consumes σ as data, so it keeps its own backstop against a
/// hand-rolled divergent jet reaching the patches.
#[test]
fn non_finite_sigma_jet_is_rejected_by_assembly() {
    let mut sigma = sigma_center_jet_default(5).expect("paper defaults yield a finite σ jet");
    sigma.ds_uu = [f64::INFINITY, 0.0];
    let input = CapAssembleInput {
        boundary: planar_ring(5, 1.0),
        beta: BetaReparam::bi4(5),
        sigma,
        pole: pt(0.0, 0.0, 0.5),
        kappa: 3,
    };
    assert!(
        matches!(
            assemble_pns_bi4_cap(&input),
            Err(CapAssembleError::NonFiniteSigmaJet { .. })
        ),
        "a non-finite σ jet must be refused, not propagated into the cap"
    );
}

/// Eq (2) G¹ gate in division-free form, across the full supported
/// valence range — **including n = 4**.
///
/// The printed eq (2) solves for the j = 0 row and divides by `2c`
/// and `6c²`, so at n = 4 (`c = 0`) it cannot be evaluated. That is a
/// property of the presentation, not the geometry: cleared of
/// denominators the same five relations hold at machine zero for
/// every valence. This gate is what proves the n = 4 cap is G¹, and
/// it is why `CapAssembleError::N4NotSupported` no longer exists.
#[test]
fn eq2_g1_holds_division_free_at_every_valence() {
    use super::super::pns_g1_constraints::bi4_g1_worst_residual;
    // Machine-zero bar. Residuals are differences of control-point
    // combinations on a radius-1.3 ring, so this is a few hundred ulp
    // of the coordinate magnitudes involved — not a fitted threshold.
    let gate = 1e-13;
    for n in [3_usize, 4, 5, 6, 7, 8, 9, 10, 11, 12] {
        let input = CapAssembleInput {
            boundary: planar_ring(n, 1.3),
            beta: BetaReparam::bi4(n),
            sigma: sigma_center_jet_default(n).expect("paper defaults yield a finite σ jet"),
            pole: pt(0.0, 0.0, 0.2),
            kappa: 3,
        };
        let result = assemble_pns_bi4_cap(&input)
            .unwrap_or_else(|e| panic!("n={n} must assemble, got {e:?}"));
        for k in 0..n {
            let nxt = (k + 1) % n;
            let p = &result.patches;
            let cur_row0 = [p[k][0][0], p[k][0][1], p[k][0][2], p[k][0][3], p[k][0][4]];
            let cur_row1 = [p[k][1][0], p[k][1][1], p[k][1][2], p[k][1][3], p[k][1][4]];
            let nxt_row1 = [
                p[nxt][0][1],
                p[nxt][1][1],
                p[nxt][2][1],
                p[nxt][3][1],
                p[nxt][4][1],
            ];
            let worst = bi4_g1_worst_residual(p[k][0][0], cur_row0, cur_row1, nxt_row1, n);
            assert!(
                worst < gate,
                "n={n} sector {k}: worst division-free eq (2) residual {worst:.3e} exceeds {gate:.0e}"
            );
        }
    }
}

/// The solved form of eq (2) must refuse n = 4 rather than divide by
/// zero and hand back a fabricated row.
#[test]
fn solved_eq2_declines_at_n4() {
    use super::super::pns_g1_constraints::{Bi4G1SolveError, bi4_g1_row_from_row1};
    let z = pt(0.0, 0.0, 0.0);
    let row = [z; 5];
    assert!(
        matches!(
            bi4_g1_row_from_row1(z, row, row, 4),
            Err(Bi4G1SolveError::UnderdeterminedAtN4)
        ),
        "n=4 must decline the solve; the G¹ system constrains rather than determines there"
    );
    assert!(
        bi4_g1_row_from_row1(z, row, row, 5).is_ok(),
        "n=5 is determined and must still solve"
    );
}
