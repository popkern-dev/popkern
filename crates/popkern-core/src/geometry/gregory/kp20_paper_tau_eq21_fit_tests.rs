//! PK-134 acceptance: the printed KP20 µ0 table is the eq (21)
//! normalization of the SOLVED paper τ.

use super::kp20_paper_tau_eq21_mu_fit;

/// With µ FIXED at the printed table, a σ-equivariant domain 3-jet ψ
/// takes the solved τ's cubic jet exactly onto the eq (22) template:
/// 20 equations vs 9 ψ parameters, so the ~1e-10 residual is 11
/// equations of genuine agreement — the printed table is reproduced.
#[test]
fn printed_mu0_table_is_the_eq21_normalized_solved_tau() {
    for n in [3_usize, 5, 6, 7, 8, 9, 10] {
        let fit = kp20_paper_tau_eq21_mu_fit(n).expect("fit");
        eprintln!(
            "n={}: residual={:.3e} mu={:?} swapped={}",
            fit.n, fit.residual_inf, fit.mu, fit.swapped,
        );
        assert!(
            fit.residual_inf < 1.0e-6,
            "n={n}: eq21 gauge fit residual {:.3e} — printed table NOT \
             reproduced by the solved tau",
            fit.residual_inf,
        );
    }
}
