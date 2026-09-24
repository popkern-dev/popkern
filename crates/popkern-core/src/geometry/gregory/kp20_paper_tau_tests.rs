use super::*;
use crate::geometry::gregory::cc_characteristic_tensor_border::cc_characteristic_tensor_border;
use crate::geometry::gregory::kp20_bi6_tau_f3_solve::kp20_planar_tau_internal_c2_constraints;

/// Measurement probe: solved system health per valence.
#[test]
fn paper_tau_measurements_are_reported() {
    for n in [3, 5, 6, 7, 8] {
        match solve_kp20_paper_tau(n) {
            Ok(report) => {
                eprintln!(
                    "n={n}: rows={rows:e} neg={neg} zero={zero} min_det={dmin:e} jet={jet:?}",
                    rows = report.max_row_residual,
                    neg = report.negative_jacobian_samples,
                    zero = report.near_zero_jacobian_samples,
                    dmin = report.min_jacobian,
                    jet = report.jet_delta,
                );
            }
            Err(err) => eprintln!("n={n}: solve failed: {err}"),
        }
    }
}

fn eval_vector_family(tau: &Kp20MacroPatch, rows: &[Kp20PlanarTauLinearConstraint]) -> f64 {
    let mut worst = 0.0_f64;
    for row in rows {
        let mut acc = [0.0_f64; 2];
        for (coeff, w) in &row.terms {
            let cp = tau.piece(coeff.piece).cps[coeff.i][coeff.j];
            acc[0] += w * cp.x;
            acc[1] += w * cp.y;
        }
        worst = worst
            .max((acc[0] - row.target[0]).abs())
            .max((acc[1] - row.target[1]).abs());
    }
    worst
}

fn eval_scalar_family(tau: &Kp20MacroPatch, rows: &[Kp20PlanarTauScalarLinearConstraint]) -> f64 {
    let mut worst = 0.0_f64;
    for row in rows {
        let mut acc = 0.0_f64;
        for (coeff, w) in &row.terms {
            let cp = tau.piece(coeff.coeff.piece).cps[coeff.coeff.i][coeff.coeff.j];
            acc += w * if coeff.axis == 0 { cp.x } else { cp.y };
        }
        worst = worst.max((acc - row.target).abs());
    }
    worst
}

/// Every row family must close at solver tolerance simultaneously —
/// the paper's τ system (hard border pins, internal C², rotational
/// eq (8)/(9), diagonal symmetry, origin/Fig-8b rows) is exactly
/// feasible on the derived characteristic border.
#[test]
fn every_row_family_closes_simultaneously() {
    for n in [3, 5, 6, 7, 8] {
        let solved = solve_kp20_paper_tau_raw(n).expect("solve");
        let tau = &solved.tau_macro_patch;
        let border = cc_characteristic_tensor_border(n).expect("border");
        let internal = kp20_planar_tau_internal_c2_constraints();
        let pins = kp20_paper_tau_border_pin_rows(n, &border.sectors[0]);
        let origin = kp20_paper_tau_origin_rows();
        let g_rows = kp20_paper_tau_rotational_g1_g2_rows(n);
        let diag = kp20_paper_tau_diagonal_symmetry_rows(n);
        for (label, residual) in [
            ("internal C2", eval_vector_family(tau, &internal)),
            ("border pins", eval_vector_family(tau, &pins)),
            ("origin", eval_vector_family(tau, &origin)),
            ("rotational G1/G2", eval_scalar_family(tau, &g_rows)),
            ("diagonal symmetry", eval_scalar_family(tau, &diag)),
        ] {
            assert!(residual < 1.0e-8, "n={n}: {label} residual {residual:e}",);
        }
    }
}

/// WP2 acceptance gates: exact feasibility, fold-free sampled
/// Jacobian, and the cubic expansion at `o` measured against the
/// printed eq (21)/(22) table.
///
/// The printed table is an analytic apex model (µ₁ = 1/3 exactly;
/// µ₂, µ₃ analytic in µ₀) produced by the paper's eq (21)
/// "reformulate and scale" normalization of the solved cubic; the
/// exact normalization/fit procedure is not printed, so the µ-level
/// comparison is bounded loosely here and recorded as a measurement.
/// Deriving that procedure (or gating through the WP4 cap-level
/// residuals instead) is the recorded open question of WP2.
#[test]
fn solved_paper_tau_is_feasible_fold_free_and_near_the_printed_table() {
    for n in [3, 5, 6, 7, 8] {
        let report = match solve_kp20_paper_tau(n) {
            Ok(report) => report,
            Err(err) => panic!("paper tau solve failed for n={n}: {err}"),
        };
        assert!(
            report.max_row_residual < 1.0e-8,
            "n={n}: row residual {:e}",
            report.max_row_residual,
        );
        assert_eq!(
            report.negative_jacobian_samples, 0,
            "n={n}: solved tau folds (min det {:e})",
            report.min_jacobian,
        );
        assert_eq!(
            report.near_zero_jacobian_samples, 0,
            "n={n}: solved tau has near-singular samples",
        );
        assert!(
            report.min_jacobian > 1.0e-3,
            "n={n}: sampled Jacobian margin too small: {:e}",
            report.min_jacobian,
        );
        let jet_delta = report
            .jet_delta
            .expect("tabulated valences carry the printed-table cross-check");
        assert!(
            jet_delta[0] < 1.0e-10,
            "n={n}: tau(o) must be the origin: {:e}",
            jet_delta[0],
        );
        // These orders carry the eq (21) reformulation GAUGE (the
        // jet_delta measure compares raw charts without the domain
        // 3-jet ψ) — the exact reproduction gate is
        // printed_mu0_table_is_the_eq21_normalized_solved_tau
        // (kp20_paper_tau_eq21_fit.rs, PK-134). Bounds here are
        // regression levels for the raw-chart deviation only.
        assert!(
            jet_delta[2] < 0.5 && jet_delta[3] < 1.5,
            "n={n}: cubic-expansion deviation from the printed table grew: {jet_delta:?}",
        );
    }
}

/// The corner-overlap symmetrization audit: the S-antisymmetric part
/// of the two strips' data is confined to P3's corner 3×3 (the
/// recorded interpretive choice), and after averaging, the assembled
/// pins are exactly diagonal-symmetric.
#[test]
fn averaged_pins_are_exactly_diagonal_symmetric() {
    for n in [3, 5, 6, 8] {
        let alpha = core::f64::consts::TAU / n as f64;
        let (wx, wy) = (alpha.cos(), alpha.sin());
        let reflect = |t: [f64; 2]| [wx * t[0] + wy * t[1], wy * t[0] - wx * t[1]];
        let diag_piece = |p: Kp20PieceId| match p {
            Kp20PieceId::NearOriginLeft => Kp20PieceId::NearOriginLeft,
            Kp20PieceId::NearMidpointRight => Kp20PieceId::NearMidpointRight,
            Kp20PieceId::NearOriginRight => Kp20PieceId::NearMidpointLeft,
            Kp20PieceId::NearMidpointLeft => Kp20PieceId::NearOriginRight,
        };
        let border = cc_characteristic_tensor_border(n).expect("border");
        let pins = kp20_paper_tau_border_pin_rows(n, &border.sectors[0]);
        let mut worst = 0.0_f64;
        for pin in &pins {
            let (coeff, _) = pin.terms[0];
            let partner = pins
                .iter()
                .find(|other| {
                    let (oc, _) = other.terms[0];
                    oc.piece == diag_piece(coeff.piece) && oc.i == coeff.j && oc.j == coeff.i
                })
                .expect("every pin has a diagonal partner");
            let expected = reflect(partner.target);
            worst = worst
                .max((pin.target[0] - expected[0]).abs())
                .max((pin.target[1] - expected[1]).abs());
        }
        assert!(
            worst < 1.0e-10,
            "n={n}: assembled pins are not diagonal-symmetric: {worst:e}",
        );
    }
}
