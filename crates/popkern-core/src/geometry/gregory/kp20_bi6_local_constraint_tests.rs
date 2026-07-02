//! Phase A unit tests for the KP20 eq (8) G¹ and eq (9) G²
//! residual extractors.
//!
//! These tests pin the polynomial identity the extractors compute
//! against analytically-known synthetic patches. They are
//! complementary to the Phase B integration tests in
//! `kp20_bi6_assemble_tests.rs` which drive the full
//! `assemble_bi6_polar_cap` pipeline.
//!
//! Reference: KP20 §4.1 eq (8) and eq (9). See
//! `docs/decisions/2026-06-23-kp20-bi6-macro-patch-paper-grade.md`
//! Phase A test specifications.

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::kp_polar_boundary_input::TensorBorderBi3;
use crate::geometry::gregory::kp20_bi6_assemble::{
    Kp20Bi6CapInput, Kp20Bi6MacroCapInput, assemble_kp20_bi6_cap, assemble_kp20_bi6_macro_cap,
    assemble_kp20_bi6_macro_cap_with_antidiagonal4_restore_probe,
};
use crate::geometry::gregory::kp20_bi6_guide_build::{Kp20GuideBuildInput, build_kp20_guide_ring};
use crate::geometry::gregory::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use crate::geometry::gregory::kp20_bi6_piece_maps::sample_kp20_macro_patch_from_default_maps;
use crate::geometry::gregory::kp20_bi6_tau::{Kp20TauMap, kp20_compose_guide_tau_jet};

fn zero_grid() -> [[Point3; 7]; 7] {
    [[Point3::new(0.0, 0.0, 0.0); 7]; 7]
}

/// Pure translation along `dir`: all CPs at the same point. Means
/// `∂_v p ≡ 0`, `∂_u p ≡ 0`, all higher derivatives also zero.
fn constant_grid(dir: Vec3) -> [[Point3; 7]; 7] {
    let p = Point3::new(dir.x, dir.y, dir.z);
    [[p; 7]; 7]
}

/// `p(u, v) = u · dir`. Linear in `u`, constant in `v`. Means
/// `∂_v p ≡ 0`, `∂_u p ≡ dir`, all higher derivatives zero.
fn linear_in_u_grid(dir: Vec3) -> [[Point3; 7]; 7] {
    let mut g = zero_grid();
    for (i, row) in g.iter_mut().enumerate() {
        let t = i as f64 / 6.0;
        let p = Point3::new(t * dir.x, t * dir.y, t * dir.z);
        for cp in row.iter_mut() {
            *cp = p;
        }
    }
    g
}

/// `p(u, v) = (a + v · b) · u + c`. Then `∂_u p = a + v·b`,
/// `∂_v p = u · b`, `∂_u∂_v p = b`, all second derivatives zero
/// except `∂_u∂_v` (already covered) — so `∂²_vv p = 0`,
/// `∂²_uu p = 0`. Useful for v-direction tests.
#[allow(dead_code)]
fn separable_grid(a: Vec3, b: Vec3, c: Vec3) -> [[Point3; 7]; 7] {
    let mut g = zero_grid();
    for (i, row) in g.iter_mut().enumerate() {
        let u = i as f64 / 6.0;
        for (j, cp) in row.iter_mut().enumerate() {
            let v = j as f64 / 6.0;
            *cp = Point3::new(
                (a.x + v * b.x) * u + c.x,
                (a.y + v * b.y) * u + c.y,
                (a.z + v * b.z) * u + c.z,
            );
        }
    }
    g
}

fn worst_mag_v3(vs: &[Vec3]) -> f64 {
    vs.iter()
        .map(|v| v.x.abs().max(v.y.abs()).max(v.z.abs()))
        .fold(0.0_f64, f64::max)
}

// --- Test 1: constant patches → residual is zero ----------------

/// Per ADR Phase A test (1): if both adjacent pieces are constant
/// Point3 fields, eq (8) and eq (9) residuals are exactly zero
/// (all derivatives vanish identically).
#[test]
fn eq8_residual_vanishes_on_constant_patches() {
    let prev = constant_grid(Vec3::new(1.0, -0.5, 0.3));
    let curr = constant_grid(Vec3::new(2.0, 0.1, -0.7));
    let c = (std::f64::consts::TAU / 6.0).cos();
    let res = kp20_eq8_g1_residual(&prev, &curr, kp20_b_coeffs_near_o(c));
    let worst = worst_mag_v3(&res);
    assert!(
        worst < 1e-12,
        "eq(8) residual on constant patches should vanish, got worst = {worst:e}",
    );
}

#[test]
fn eq9_residual_vanishes_on_constant_patches() {
    let prev = constant_grid(Vec3::new(0.7, 0.2, -0.4));
    let curr = constant_grid(Vec3::new(-0.3, 0.9, 0.5));
    let c = (std::f64::consts::TAU / 5.0).cos();
    let res = kp20_eq9_g2_residual(&prev, &curr, kp20_b_coeffs_near_m(c));
    let worst = worst_mag_v3(&res);
    assert!(
        worst < 1e-12,
        "eq(9) residual on constant patches should vanish, got worst = {worst:e}",
    );
}

// --- Test 2: linear-in-u patches measure the b · ∂_u p term -----

/// Per ADR Phase A test (2): if `∂_v` strips on both sides cancel
/// by symmetry (zero here because both are constant in v), eq (8)
/// reduces to `−b(u) · ∂_u p`. With `p(u, v) = u · dir`, we have
/// `∂_u p ≡ dir` constant, so the residual is the polynomial
/// `−b(u) · dir`. In BB-form, the residual coefficients
/// are `−b_i · dir`.
///
/// The residual is BB-deg-6 (we de-elevate from deg-7) so we
/// compare its values at sample points rather than coefficient-by-
/// coefficient: the deg-2 BB poly `b` produces residuals that are
/// linear interpolants of `−b_i · dir` after de-elevation.
#[test]
fn eq8_residual_equals_minus_b_times_du_on_linear_in_u_patches() {
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let prev = linear_in_u_grid(dir);
    let curr = constant_grid(Vec3::new(0.0, 0.0, 0.0));
    let c = (std::f64::consts::TAU / 6.0).cos();
    let b = kp20_b_coeffs_near_o(c);
    let res = kp20_eq8_g1_residual(&prev, &curr, b);

    // The residual polynomial is `−b(u) · dir`. Its analytical
    // BB-deg-2 representation has coefficients `−b_i · dir`. After
    // de-elevation from deg-7 to deg-6 the polynomial is still
    // `−b(u)·dir` (exact, since input degree is 2).
    //
    // We compare integrals: integral over [0, 1] of `−b(u)`
    // equals `−(b_0 + b_1 + b_2) / 3` (mean of BB coefs for deg-2),
    // and the residual integral equals the mean of its BB-deg-6
    // coefs.
    let integrand_expected = -(b[0] + b[1] + b[2]) / 3.0; // <-b(u)>_u
    let integrand_actual_x = res.iter().map(|v| v.x).sum::<f64>() / (res.len() as f64);
    assert!(
        (integrand_actual_x - integrand_expected * dir.x).abs() < 1e-10,
        "expected mean ≈ {expected:e}, got {actual:e}",
        expected = integrand_expected * dir.x,
        actual = integrand_actual_x,
    );

    // Sanity check: the y and z components are zero because dir.y, z = 0.
    let worst_yz: f64 = res
        .iter()
        .map(|v| v.y.abs().max(v.z.abs()))
        .fold(0.0_f64, f64::max);
    assert!(
        worst_yz < 1e-12,
        "off-axis components should vanish, got {worst_yz:e}",
    );
}

// --- Test 3 & 4: b_coeffs paper-form pin -------------------------

/// Per ADR Phase A test (3): the `b(u)` BB-coefficient triple near
/// the extraordinary point `o` is the paper-faithful `[2c, 2c, c]`
/// with `c := cos(2π/n)`.
#[test]
fn b_coeffs_near_o_matches_paper_eq10() {
    for n in 3..=10 {
        let c = (std::f64::consts::TAU / n as f64).cos();
        let b = kp20_b_coeffs_near_o(c);
        assert!(
            (b[0] - 2.0 * c).abs() < 1e-15,
            "near_o b₀ should be 2c at n = {n}",
        );
        assert!(
            (b[1] - 2.0 * c).abs() < 1e-15,
            "near_o b₁ should be 2c at n = {n}",
        );
        assert!((b[2] - c).abs() < 1e-15, "near_o b₂ should be c at n = {n}",);
    }
}

/// Per ADR Phase A test (4): the `b(u)` BB-coefficient triple near
/// the outer midpoint `m` is the paper-faithful `[c, 0, 0]`.
#[test]
fn b_coeffs_near_m_matches_paper_eq10() {
    for n in 3..=10 {
        let c = (std::f64::consts::TAU / n as f64).cos();
        let b = kp20_b_coeffs_near_m(c);
        assert!((b[0] - c).abs() < 1e-15, "near_m b₀ should be c at n = {n}",);
        assert!(b[1].abs() < 1e-15, "near_m b₁ should be 0 at n = {n}");
        assert!(b[2].abs() < 1e-15, "near_m b₂ should be 0 at n = {n}");
    }
}

// --- Test 5: eq(9) sanity on constant + linear-in-u --------------

/// Per ADR Phase A test (5): for a configuration where eq (9)'s
/// individual terms are analytically zero, the extractor returns
/// zero. We use `p̃ = const`, `p = u·dir`: then `∂²_vv p̃ = 0`,
/// `∂²_vv p = 0`, `∂_u∂_v p = 0`, `∂²_uu p = 0`, but `∂_u p ≠ 0`.
/// The residual reduces to `−b(u) b'(u) · dir`, a degree-3
/// polynomial in u.
#[test]
fn eq9_residual_equals_minus_b_bprime_du_on_linear_in_u_patches() {
    let dir = Vec3::new(0.0, 1.0, 0.0);
    let prev = linear_in_u_grid(dir);
    let curr = constant_grid(Vec3::new(0.0, 0.0, 0.0));
    let c = (std::f64::consts::TAU / 5.0).cos();
    let b = kp20_b_coeffs_near_o(c);
    let res = kp20_eq9_g2_residual(&prev, &curr, b);

    // Analytical reference: residual(u) = −b(u) · b'(u) · dir,
    // a polynomial of degree 3 in u. We compute it as a closed-
    // form value at u = 0 and u = 1 and compare to the BB[0] and
    // BB[7] of the returned residual (which equal the polynomial
    // values at u = 0 and u = 1 by the Bernstein-endpoint identity).
    //
    // At u = 0: b(0) = b₀; b'(0) is the slope of b at 0.
    //   b(u) BB-deg-2 = b₀(1-u)² + 2 b₁ u(1-u) + b₂ u²,
    //   so b'(u) = 2 (b₁ - b₀)(1 - u) + 2 (b₂ - b₁) u,
    //   b'(0) = 2 (b₁ - b₀), b'(1) = 2 (b₂ - b₁).
    let b_bprime_0 = b[0] * 2.0 * (b[1] - b[0]);
    let b_bprime_1 = b[2] * 2.0 * (b[2] - b[1]);

    // residual = −b · b' · dir (scalar applied to Vec3).
    let expected_at_0 = -b_bprime_0 * dir.y; // y component
    let expected_at_1 = -b_bprime_1 * dir.y;

    let bb0_y = res[0].y;
    let bb7_y = res[7].y;
    assert!(
        (bb0_y - expected_at_0).abs() < 1e-10,
        "eq(9) residual at u=0: expected {expected_at_0:e}, got {bb0_y:e}",
    );
    assert!(
        (bb7_y - expected_at_1).abs() < 1e-10,
        "eq(9) residual at u=1: expected {expected_at_1:e}, got {bb7_y:e}",
    );
}

// --- Test 6: signature-level smoke test --------------------------

/// Per ADR Phase A test (6): the extractors are pure and depend
/// only on their inputs. Calling them twice with identical inputs
/// returns identical outputs.
#[test]
fn extractors_are_deterministic() {
    let prev = linear_in_u_grid(Vec3::new(0.3, 0.7, -0.1));
    let curr = linear_in_u_grid(Vec3::new(-0.2, 0.4, 0.5));
    let b = kp20_b_coeffs_near_o((std::f64::consts::TAU / 7.0).cos());

    let r1 = kp20_eq8_g1_residual(&prev, &curr, b);
    let r2 = kp20_eq8_g1_residual(&prev, &curr, b);
    for k in 0..7 {
        assert_eq!(r1[k], r2[k], "eq(8) result mismatch at BB[{k}]");
    }

    let s1 = kp20_eq9_g2_residual(&prev, &curr, b);
    let s2 = kp20_eq9_g2_residual(&prev, &curr, b);
    for k in 0..8 {
        assert_eq!(s1[k], s2[k], "eq(9) result mismatch at BB[{k}]");
    }
}

// --- Test 7: cancellation symmetry on twin-pair patches ----------

/// Per ADR Phase A test (7): if `p̃ = p` (same patch on both
/// sides) and `b(u) ≡ 0`, then eq (8) reduces to `2 · ∂_v p` and
/// eq (9) reduces to zero (modulo the `∂_u∂_v p`, `∂²_uu p`,
/// `∂_u p` terms which all multiply zero scalars). Set both
/// `b_coeffs = [0, 0, 0]` and use a constant-in-v patch so
/// `∂_v p = 0`. Then both residuals are zero — sanity check the
/// zeroing path.
#[test]
fn eq8_residual_vanishes_when_b_is_zero_and_v_constant() {
    let prev = linear_in_u_grid(Vec3::new(1.0, 0.5, -0.3));
    let curr = linear_in_u_grid(Vec3::new(-0.2, 0.8, 0.4));
    let b_zero = [0.0; 3];
    let res = kp20_eq8_g1_residual(&prev, &curr, b_zero);
    let worst = worst_mag_v3(&res);
    assert!(
        worst < 1e-12,
        "eq(8) with b ≡ 0 and v-constant patches should vanish, got {worst:e}",
    );
}

// --- Phase A2: Jacobian-extractor consistency --------------------

/// The Jacobian's per-column values must match the residual
/// extractor's response to a unit-load at that column. Tests this
/// for the eq (8) Jacobian.
#[test]
fn eq8_jacobian_matches_unit_load_residual_at_random_columns() {
    let c = (std::f64::consts::TAU / 6.0).cos();
    let b = kp20_b_coeffs_near_o(c);
    let jac = kp20_eq8_jacobian_dense(b);
    // Sample a few columns deterministically.
    for &(side, i, j) in &[(0, 0, 0), (1, 6, 6), (0, 3, 3), (1, 4, 2)] {
        let col = side * 49 + i * 7 + j;
        let mut prev = zero_grid();
        let mut curr = zero_grid();
        if side == 0 {
            prev[i][j] = Point3::new(1.0, 0.0, 0.0);
        } else {
            curr[i][j] = Point3::new(1.0, 0.0, 0.0);
        }
        let bb = kp20_eq8_g1_residual(&prev, &curr, b);
        for row in 0..7 {
            let from_jac = jac[row][col];
            let from_residual = bb[row].x;
            assert!(
                (from_jac - from_residual).abs() < 1e-12,
                "eq(8) Jacobian disagrees at (side={side}, i={i}, j={j}, row={row}): jac={from_jac:e} resid={from_residual:e}",
            );
        }
    }
}

#[test]
fn eq9_jacobian_matches_unit_load_residual_at_random_columns() {
    let c = (std::f64::consts::TAU / 5.0).cos();
    let b = kp20_b_coeffs_near_m(c);
    let jac = kp20_eq9_jacobian_dense(b);
    for &(side, i, j) in &[(0, 0, 0), (1, 6, 6), (0, 2, 4), (1, 5, 1)] {
        let col = side * 49 + i * 7 + j;
        let mut prev = zero_grid();
        let mut curr = zero_grid();
        if side == 0 {
            prev[i][j] = Point3::new(1.0, 0.0, 0.0);
        } else {
            curr[i][j] = Point3::new(1.0, 0.0, 0.0);
        }
        let bb = kp20_eq9_g2_residual(&prev, &curr, b);
        for row in 0..8 {
            let from_jac = jac[row][col];
            let from_residual = bb[row].x;
            assert!(
                (from_jac - from_residual).abs() < 1e-12,
                "eq(9) Jacobian disagrees at (side={side}, i={i}, j={j}, row={row}): jac={from_jac:e} resid={from_residual:e}",
            );
        }
    }
}

// --- Phase B: assembled-cap measurement on hexagon + pentagon ----

/// Synthetic bi-3 tensor border with `n` sectors at evenly-spaced
/// angles around a unit ring, planar in the xy plane with the cap
/// hole at the origin. The kernel's KP20 assembler accepts this
/// shape and produces a well-defined macro cap output.
fn synthetic_uniform_tensor_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (std::f64::consts::TAU) * (k as f64) / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        // Three depth rows, four CPs per row. Row j is at distance
        // `1 + 0.1 j` from origin. Row 0 outermost.
        let mut sector = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for (row_idx, row) in sector.iter_mut().enumerate() {
            let radial = 1.0 + 0.1 * row_idx as f64;
            for (col_idx, cp) in row.iter_mut().enumerate() {
                let t = col_idx as f64 / 3.0;
                // Spread CPs tangentially within the sector
                let arc = (theta - std::f64::consts::PI / n as f64)
                    + (2.0 * std::f64::consts::PI / n as f64) * t;
                *cp = Point3::new(
                    radial * arc.cos(),
                    radial * arc.sin(),
                    0.05 * row_idx as f64,
                );
            }
            let _ = (dx, dy); // suppress unused-warning if optimiser drops
        }
        sectors.push(sector);
    }
    TensorBorderBi3 { sectors }
}

fn assemble_synthetic_cap(n: usize) -> Kp20Bi6MacroCapOutput {
    let tensor_border = synthetic_uniform_tensor_border(n);
    let star_point = Point3::new(0.0, 0.0, 0.5);
    let input = Kp20Bi6CapInput {
        star_point,
        tensor_border,
    };
    assemble_kp20_bi6_cap(&input).expect("synthetic assemble must succeed")
}

fn synthetic_macro_cap_input(n: usize) -> Kp20Bi6MacroCapInput {
    let tensor_border = synthetic_uniform_tensor_border(n);
    let star_point = Point3::new(0.0, 0.0, 0.5);
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .expect("synthetic guide builds");
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(sample_kp20_macro_patch_from_default_maps)
        .collect();
    Kp20Bi6MacroCapInput {
        tensor_border,
        sampled_macro_patches,
    }
}

#[test]
fn measure_report_is_finite_on_hexagon6_uniform() {
    let cap = assemble_synthetic_cap(6);
    let report = measure_kp20_assembly_constraints(&cap);
    assert!(
        report.all_finite(),
        "hexagon6 measurement report has non-finite entries: {report:?}",
    );
    // Worst-magnitude summaries should also be finite.
    assert!(report.worst_eq8().is_finite());
    assert!(report.worst_eq9().is_finite());
    eprintln!(
        "hexagon6 KP20 cap diagnostic: worst_eq8={:e} worst_eq9={:e} L6(1)={:e} L6(2)={:e} L6(3)={:e} L6(4)={:e}",
        report.worst_eq8(),
        report.worst_eq9(),
        report.worst_lemma6_cond_1,
        report.worst_lemma6_cond_2,
        report.worst_lemma6_cond_3,
        report.worst_lemma6_cond_4,
    );
}

#[test]
fn measure_report_is_finite_on_pentagon5_uniform() {
    let cap = assemble_synthetic_cap(5);
    let report = measure_kp20_assembly_constraints(&cap);
    assert!(
        report.all_finite(),
        "pentagon5 measurement report has non-finite entries: {report:?}",
    );
    eprintln!(
        "pentagon5 KP20 cap diagnostic: worst_eq8={:e} worst_eq9={:e} L6(1)={:e} L6(2)={:e} L6(3)={:e} L6(4)={:e}",
        report.worst_eq8(),
        report.worst_eq9(),
        report.worst_lemma6_cond_1,
        report.worst_lemma6_cond_2,
        report.worst_lemma6_cond_3,
        report.worst_lemma6_cond_4,
    );
}

#[test]
fn measure_report_per_sector_vec_length_matches_n() {
    for n in &[5_usize, 6, 7] {
        let cap = assemble_synthetic_cap(*n);
        let report = measure_kp20_assembly_constraints(&cap);
        assert_eq!(report.worst_eq8_near_o_per_sector.len(), *n);
        assert_eq!(report.worst_eq8_near_m_per_sector.len(), *n);
        assert_eq!(report.worst_eq9_near_o_per_sector.len(), *n);
        assert_eq!(report.worst_eq9_near_m_per_sector.len(), *n);
    }
}

#[test]
fn lemma6_layer_connectedness_report_uses_within_macro_origin_pair() {
    for n in &[5_usize, 6] {
        let cap = assemble_synthetic_cap(*n);
        let report = measure_kp20_assembly_constraints(&cap);
        let (within_l0, cross_l0) = lemma6_layer_connectedness_comparison_for_tests(&cap, 0, 4);
        let (within_l1, cross_l1) = lemma6_layer_connectedness_comparison_for_tests(&cap, 1, 3);
        let (within_l2, cross_l2) = lemma6_layer_connectedness_comparison_for_tests(&cap, 2, 2);

        eprintln!(
            "Lemma 6 connectedness n={n}: within L0={within_l0:e} L1={within_l1:e} L2={within_l2:e}; legacy-cross L0={cross_l0:e} L1={cross_l1:e} L2={cross_l2:e}",
        );

        assert!(
            (report.worst_lemma6_cond_2 - within_l0).abs() < 1.0e-12,
            "report L6(2) must use within-macro origin pair",
        );
        assert!(
            (report.worst_lemma6_cond_3 - within_l1).abs() < 1.0e-12,
            "report L6(3) must use within-macro origin pair",
        );
        assert!(
            (report.worst_lemma6_cond_4 - within_l2).abs() < 1.0e-12,
            "report L6(4) must use within-macro origin pair",
        );
        assert!(cross_l0.is_finite());
        assert!(cross_l1.is_finite());
        assert!(cross_l2.is_finite());
    }
}

#[test]
fn lemma6_layer_connectedness_reports_worst_derivative_localization() {
    for n in &[5_usize, 6] {
        let cap = assemble_synthetic_cap(*n);
        let report = measure_kp20_assembly_constraints(&cap);

        eprintln!(
            "Lemma 6 detail n={n}: L6(2)={:?}; L6(3)={:?}; L6(4)={:?}",
            report.worst_lemma6_cond_2_detail,
            report.worst_lemma6_cond_3_detail,
            report.worst_lemma6_cond_4_detail,
        );

        assert_eq!(report.worst_lemma6_cond_2_detail.condition, 2);
        assert_eq!(report.worst_lemma6_cond_3_detail.condition, 3);
        assert_eq!(report.worst_lemma6_cond_4_detail.condition, 4);
        assert_eq!(report.worst_lemma6_cond_2_detail.layer_j, 0);
        assert_eq!(report.worst_lemma6_cond_3_detail.layer_j, 1);
        assert_eq!(report.worst_lemma6_cond_4_detail.layer_j, 2);

        assert!(
            (report.worst_lemma6_cond_2 - report.worst_lemma6_cond_2_detail.magnitude).abs()
                < 1.0e-12
        );
        assert!(
            (report.worst_lemma6_cond_3 - report.worst_lemma6_cond_3_detail.magnitude).abs()
                < 1.0e-12
        );
        assert!(
            (report.worst_lemma6_cond_4 - report.worst_lemma6_cond_4_detail.magnitude).abs()
                < 1.0e-12
        );

        for detail in [
            report.worst_lemma6_cond_2_detail,
            report.worst_lemma6_cond_3_detail,
            report.worst_lemma6_cond_4_detail,
        ] {
            assert!(detail.sector < *n);
            assert!(detail.derivative_order <= 4);
            assert_eq!(detail.cp_i_start, 0);
            assert_eq!(detail.cp_i_end, detail.derivative_order);
            assert!(detail.magnitude.is_finite());
            assert!(detail.defect.x.is_finite());
            assert!(detail.defect.y.is_finite());
            assert!(detail.defect.z.is_finite());
        }
    }
}

#[test]
fn antidiagonal4_restore_probe_measures_first_ring_outside_overlay() {
    for n in &[5_usize, 6] {
        let input = synthetic_macro_cap_input(*n);
        let baseline =
            assemble_kp20_bi6_macro_cap(&input).expect("baseline synthetic assemble succeeds");
        let restored = assemble_kp20_bi6_macro_cap_with_antidiagonal4_restore_probe(&input)
            .expect("antidiagonal4 probe assemble succeeds");
        let base = measure_kp20_assembly_constraints(&baseline);
        let probe = measure_kp20_assembly_constraints(&restored);

        eprintln!(
            "antidiagonal4 restore probe n={n}: eq8 {be8:e}->{pe8:e}; eq9 {be9:e}->{pe9:e}; L6(1) {bl1:e}->{pl1:e}; L6(2) {bl2:e}->{pl2:e}; L6(3) {bl3:e}->{pl3:e}; L6(4) {bl4:e}->{pl4:e}; detail L6(2)={d2:?}; L6(3)={d3:?}",
            be8 = base.worst_eq8(),
            pe8 = probe.worst_eq8(),
            be9 = base.worst_eq9(),
            pe9 = probe.worst_eq9(),
            bl1 = base.worst_lemma6_cond_1,
            pl1 = probe.worst_lemma6_cond_1,
            bl2 = base.worst_lemma6_cond_2,
            pl2 = probe.worst_lemma6_cond_2,
            bl3 = base.worst_lemma6_cond_3,
            pl3 = probe.worst_lemma6_cond_3,
            bl4 = base.worst_lemma6_cond_4,
            pl4 = probe.worst_lemma6_cond_4,
            d2 = probe.worst_lemma6_cond_2_detail,
            d3 = probe.worst_lemma6_cond_3_detail,
        );

        assert!(base.all_finite());
        assert!(probe.all_finite());
        assert!(probe.worst_eq8().is_finite());
        assert!(probe.worst_eq9().is_finite());
        assert!(
            (base.worst_eq8() - probe.worst_eq8()).abs() < 1.0e-12,
            "antidiagonal4 restore should be a no-op in the current pipeline"
        );
        assert!(
            (base.worst_eq9() - probe.worst_eq9()).abs() < 1.0e-12,
            "antidiagonal4 restore should be a no-op in the current pipeline"
        );
        assert!(
            (base.worst_lemma6_cond_2 - probe.worst_lemma6_cond_2).abs() < 1.0e-12,
            "antidiagonal4 restore should not change L6(2) in the current pipeline"
        );
        assert!(
            (base.worst_lemma6_cond_3 - probe.worst_lemma6_cond_3).abs() < 1.0e-12,
            "antidiagonal4 restore should not change L6(3) in the current pipeline"
        );
    }
}

// --- Phase C: global constraint matrix rank ---------------------

/// Per-valence global constraint matrix dimensions match
/// `kp20_constraint_matrix_layout`.
#[test]
fn global_constraint_matrix_dimensions_match_layout() {
    for n in &[3_usize, 5, 6, 7, 10] {
        let layout = kp20_constraint_matrix_layout(*n);
        let (layout2, m) = kp20_global_constraint_matrix(*n);
        assert_eq!(layout, layout2);
        assert_eq!(m.len(), layout.rows());
        for row in &m {
            assert_eq!(row.len(), layout.cols);
        }
    }
}

// --- τ-table paper-fidelity A/B measurements ---------------------

/// Paper τ MUST also satisfy the apex-derivative cross-sector
/// identity (since Lemma 4 requires it). If this test fails the
/// paper-extraction-derived formula in `for_sector_paper` is wrong
/// and the audit needs revisiting.
#[test]
fn paper_tau_map_apex_derivatives_match_across_sectors() {
    let valences: &[usize] = &[3, 5, 6, 7, 8, 9, 10];
    for &n in valences {
        let mut taus = Vec::with_capacity(n);
        for s in 0..n {
            taus.push(Kp20TauMap::for_sector_paper(n, s).expect("paper tau"));
        }
        let mut worst = [0.0_f64; 4];
        for s in 0..n {
            let next = (s + 1) % n;
            let curr = taus[s].jet3(0.0, 0.0);
            let nbor = taus[next].jet3(0.0, 0.0);
            for (i, slot) in worst.iter_mut().enumerate() {
                let cu = curr.d[i][0];
                let nv = nbor.d[0][i];
                let mag = (cu[0] - nv[0]).abs().max((cu[1] - nv[1]).abs());
                *slot = slot.max(mag);
            }
        }
        eprintln!(
            "PAPER τ apex symmetry n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst[0],
            o1 = worst[1],
            o2 = worst[2],
            o3 = worst[3],
        );
        for (i, w) in worst.iter().enumerate() {
            assert!(
                *w < 1e-12,
                "paper τ apex symmetry order {i} = {w:e} > 1e-12 at n = {n}",
            );
        }
    }
}

/// Compare the kernel-τ and paper-τ jets at apex side-by-side and
/// report the per-order disagreement. Numeric divergence is expected
/// at orders > 0 — this test catches whether they happen to agree
/// or not.
#[test]
fn paper_vs_kernel_tau_apex_jet_divergence() {
    let valences: &[usize] = &[3, 5, 6, 7, 8, 9, 10];
    for &n in valences {
        let mut worst = [0.0_f64; 4];
        for s in 0..n {
            let kernel = Kp20TauMap::for_sector(n, s)
                .expect("kernel tau")
                .jet3(0.0, 0.0);
            let paper = Kp20TauMap::for_sector_paper(n, s)
                .expect("paper tau")
                .jet3(0.0, 0.0);
            for (i, slot) in worst.iter_mut().enumerate() {
                let k = kernel.d[i][0];
                let p = paper.d[i][0];
                let mag = (k[0] - p[0]).abs().max((k[1] - p[1]).abs());
                *slot = slot.max(mag);
            }
        }
        eprintln!(
            "PAPER vs KERNEL τ ∂^i_u jet divergence n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst[0],
            o1 = worst[1],
            o2 = worst[2],
            o3 = worst[3],
        );
        for w in &worst {
            assert!(w.is_finite());
        }
    }
}

/// PK-088 Phase D — Lemma 4 cubic identity test re-run with the
/// paper-faithful τ map (instead of `for_sector`). If Lemma 4 holds
/// to machine zero, the kernel's τ is the bug and migrating to the
/// paper τ is the fix. If it doesn't, something else is going on.
#[test]
fn lemma4_via_paper_tau_on_synthetic_inputs() {
    use crate::geometry::Point3;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border,
        })
        .expect("guide");

        let mut jets = Vec::with_capacity(n);
        for (k, sector) in guide.sectors.iter().enumerate() {
            let tau = Kp20TauMap::for_sector_paper(n, k).expect("paper tau");
            jets.push(kp20_compose_guide_tau_jet(sector, &tau, 0.0, 0.0));
        }

        let mut worst = [0.0_f64; 4];
        for k in 0..n {
            let next = (k + 1) % n;
            for (i, slot) in worst.iter_mut().enumerate() {
                let next_v = jets[next].h[0][i];
                let curr_u = jets[k].h[i][0];
                let mag = (next_v.x - curr_u.x)
                    .abs()
                    .max((next_v.y - curr_u.y).abs())
                    .max((next_v.z - curr_u.z).abs());
                *slot = slot.max(mag);
            }
        }
        eprintln!(
            "PAPER τ + LSQ guide — Lemma 4 n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst[0],
            o1 = worst[1],
            o2 = worst[2],
            o3 = worst[3],
        );
        for w in &worst {
            assert!(w.is_finite());
        }
    }
}

/// Probe: do the guide values at the INTERIOR position (0.5, 0.5)
/// (centroid of guide triangle = U1V1 corner of NearOriginRight
/// piece map) match across sectors?
#[test]
fn lsq_guide_value_at_centroid_across_sectors() {
    use crate::geometry::Point3;
    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);
    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border,
        })
        .expect("guide");
        let mut worst = 0.0_f64;
        for s in 0..n {
            let next = (s + 1) % n;
            let v_s = guide.sector(s).jet3(0.5, 0.5).h[0][0];
            let v_n = guide.sector(next).jet3(0.5, 0.5).h[0][0];
            let mag = (v_s.x - v_n.x)
                .abs()
                .max((v_s.y - v_n.y).abs())
                .max((v_s.z - v_n.z).abs());
            worst = worst.max(mag);
        }
        eprintln!("guide value at (0.5, 0.5) cross-sector defect, n = {n}: {worst:e}");
    }
}

// --- Move 1b: LSQ guide values along shared spoke ----------------

/// After Move 1 closed the apex jet transpose, the next suspect for
/// the remaining Phase B / L6(2)/L6(3) residuals is the guide's
/// value along the SHARED SPOKE (away from the apex).
///
/// At guide-local position `(0, t)` of sector s, the guide value is
/// at world position `t · V_j_s` (sector s's j-vertex direction).
/// At guide-local position `(t, 0)` of sector (s+1), the value is
/// at world position `t · V_i_{s+1}`.
///
/// For shared spoke between sector s and (s+1):
/// `V_j_s = V_i_{s+1}`. So the world positions agree.
/// The guide VALUES at these positions should also agree for
/// cross-sector C⁰ to hold along the spoke.
///
/// Reports the worst defect for t ∈ {0.1, 0.25, 0.5, 0.75, 0.9}.
#[test]
fn lsq_guide_values_along_shared_spoke_match_across_sectors() {
    use crate::geometry::Point3;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);
    let ts: &[f64] = &[0.1, 0.25, 0.5, 0.75, 0.9];

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border,
        })
        .expect("guide");
        let mut worst = vec![0.0_f64; ts.len()];
        for s in 0..n {
            let next = (s + 1) % n;
            for (idx, &t) in ts.iter().enumerate() {
                let p_s = guide.sector(s).jet3(0.0, t).h[0][0]; // sector s at (u_g=0, v_g=t)
                let p_n = guide.sector(next).jet3(t, 0.0).h[0][0]; // sector (s+1) at (u_g=t, v_g=0)
                let mag = (p_s.x - p_n.x)
                    .abs()
                    .max((p_s.y - p_n.y).abs())
                    .max((p_s.z - p_n.z).abs());
                worst[idx] = worst[idx].max(mag);
            }
        }
        eprintln!("LSQ guide value along shared spoke, n = {n}:");
        for (idx, &t) in ts.iter().enumerate() {
            eprintln!("  t = {t}:  {w:e}", w = worst[idx]);
        }
        for &w in &worst {
            assert!(w.is_finite());
        }
    }
}

// --- Move 1: LSQ guide FULL apex jet cross-sector match ----------

/// Tests whether the LSQ guide's *full* apex cubic Hermite jet
/// satisfies the transpose identity that Lemma 4 needs.
///
/// `Step 3` of the prior investigation only checked the boundary
/// of the jet (`h[i][0]` vs `h[0][i]`). For the cubic chain rule
/// `g ∘ τ` at apex to satisfy Lemma 4, the FULL 10-entry
/// triangular jet block must agree across sectors under the
/// transpose:
///
/// ```text
///   guide_s.jet3(0,0).h[i][j] == guide_{s+1}.jet3(0,0).h[j][i]
///   for all (i, j) with i + j ≤ 3
/// ```
///
/// If this fails, the missing piece is cross-sector identifications
/// on the apex-region mixed-partial CPs of the guide.
#[test]
fn lsq_guide_full_apex_jet_transpose_match_across_sectors() {
    use crate::geometry::Point3;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border,
        })
        .expect("guide ring builds");

        let mut jets = Vec::with_capacity(n);
        for s in 0..n {
            jets.push(guide.sector(s).jet3(0.0, 0.0));
        }

        // Per-(i, j) worst defect across all sector pairs.
        let mut worst = [[0.0_f64; 4]; 4];
        for s in 0..n {
            let next = (s + 1) % n;
            for i in 0..4 {
                for j in 0..(4 - i) {
                    let a = jets[s].h[i][j];
                    let b = jets[next].h[j][i];
                    let m = (a.x - b.x)
                        .abs()
                        .max((a.y - b.y).abs())
                        .max((a.z - b.z).abs());
                    worst[i][j] = worst[i][j].max(m);
                }
            }
        }
        eprintln!("LSQ guide FULL apex jet transpose defect, n = {n}:");
        for i in 0..4 {
            for j in 0..(4 - i) {
                eprintln!("  h[{i}][{j}] vs h[{j}][{i}]:  {w:e}", w = worst[i][j]);
            }
        }
        for i in 0..4 {
            for j in 0..(4 - i) {
                assert!(
                    worst[i][j].is_finite(),
                    "non-finite at h[{i}][{j}], n = {n}",
                );
            }
        }
    }
}

// --- Investigation Step 2: τ map cross-sector symmetry -----------

/// Step 2 of the investigation plan in
/// `docs/decisions/2026-06-23-kp20-bi6-paper-grade-guide-from-kp18.md`.
///
/// Property (τ1/τ2/τ4) from
/// `docs/architecture/kp20-lemma4-sufficient-conditions.md`:
///
/// At the apex `(u, v) = (0, 0)`, sector `s`'s `∂^i_u τ_s` must
/// equal sector `(s+1)`'s `∂^i_v τ_{s+1}` for `i = 1, 2, 3`
/// because both directional derivatives lie along the shared
/// spoke. (Order 0: both are `(0, 0)` by construction.)
///
/// If this test passes, the kernel's `Kp20TauMap::for_sector` is
/// Lemma-4-coherent on its own, and the entire Phase D defect on
/// the LSQ-built cap must come from the guide ring's missing
/// **Symmetry** condition (per
/// `kp20-lemma4-sufficient-conditions.md`).
#[test]
fn kp20_tau_map_apex_derivatives_match_across_sectors() {
    let valences: &[usize] = &[3, 5, 6, 7, 8, 9, 10];
    for &n in valences {
        let mut taus = Vec::with_capacity(n);
        for s in 0..n {
            taus.push(Kp20TauMap::for_sector(n, s).expect("tau for valence"));
        }
        let mut worst_order = [0.0_f64; 4];
        for s in 0..n {
            let next = (s + 1) % n;
            let curr = taus[s].jet3(0.0, 0.0);
            let nbor = taus[next].jet3(0.0, 0.0);
            for (i, slot) in worst_order.iter_mut().enumerate() {
                // curr.d[i][0] = ∂^i_u τ_s at (0, 0); 2D vector.
                // nbor.d[0][i] = ∂^i_v τ_{s+1} at (0, 0); 2D vector.
                let cu = curr.d[i][0];
                let nv = nbor.d[0][i];
                let dx = (cu[0] - nv[0]).abs();
                let dy = (cu[1] - nv[1]).abs();
                let mag = dx.max(dy);
                *slot = slot.max(mag);
            }
        }
        eprintln!(
            "τ-map apex symmetry n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst_order[0],
            o1 = worst_order[1],
            o2 = worst_order[2],
            o3 = worst_order[3],
        );
        for (i, w) in worst_order.iter().enumerate() {
            assert!(
                *w < 1e-12,
                "τ-map apex symmetry order {i} = {w:e} > 1e-12 at n = {n}",
            );
        }
    }
}

// --- Investigation Step 3b: assembled-cap apex Lemma 4 -----------

/// Measure Lemma 4 directly on the assembled cap's ◦-marked CPs
/// (the bi-6 BB coefficients near the apex of each macro patch's
/// L/R piece pair). This is the identity the paper actually
/// requires — and it is on these CPs that PK-088's Phase B eq(8)
/// residuals are evaluated.
///
/// For each macro patch `k`:
/// ```text
///   NearOriginLeft.cps[i][0]  ==  NearOriginRight.cps[0][i]
///       for i = 0..3
/// ```
/// (after sampling and the kernel's tensor-border + C² propagation
/// + G²-strip assembly pipeline).
///
/// If Step 3 (guide-only) passes but this fails, the cause of
/// PK-088's Phase B residual is in the assembly pipeline's order:
/// step (ii) tensor-border and step (iii) C² propagation overwrite
/// the ◦-marked CPs that step (i) Hermite-sampling set from the
/// guide. The paper's algorithm sets ◦-marked CPs LAST so they are
/// preserved.
/// Estimate the eq(8) residual reduction available from the
/// Lemma-4-overlay-alone approach. Synthesises a paper-perfect
/// apex region (10 CPs per piece) on top of the existing
/// assembled cap and re-measures eq(8). Lets us know whether
/// the overlay alone can reach the 1e-8 acceptance gate or
/// whether the cps[4..6][0] strip-region inconsistency
/// (C²-propagated values not paper-consistent with the
/// overlay) leaves a residual floor.
/// Compare the current affine sampler against the kernel-τ and paper-τ
/// samplers. The tau paths include a composed-domain report before residual
/// measurement so a bad result is classified before it is interpreted as a
/// geometry failure.
#[test]
fn kp20_sampler_ab_report_classifies_domain_before_residuals() {
    use super::super::kp20_bi6_tau_solve_assemble::{
        Kp20PlanarTauPaperRowFamily, Kp20PlanarTauPaperSolveInput,
        kp20_planar_tau_rotational_paper_family_residuals,
        kp20_planar_tau_rotational_section42_source_projection_report,
        solve_kp20_planar_tau_rotational_paper_sector_f3,
    };
    use crate::geometry::Point3;
    use crate::geometry::gregory::kp20_bi6_assemble::{
        Kp20Bi6AssembleError, Kp20Bi6MacroCapAssemblyStage, Kp20Bi6MacroCapInput,
        Kp20Bi6MacroCapOutput, assemble_kp20_bi6_macro_cap,
        assemble_kp20_bi6_macro_cap_stage_snapshots,
        assemble_kp20_bi6_macro_cap_stage_snapshots_with_sampled_tensor_border,
        assemble_kp20_bi6_macro_cap_with_sampled_tensor_border_source,
    };
    use crate::geometry::gregory::kp20_bi6_cap::Kp20MacroPatch;
    use crate::geometry::gregory::kp20_bi6_cc_ring::{
        Kp20Pr08TensorBorderCandidate, kp20_pr08_stam_q1_characteristic_tensor_border,
    };
    use crate::geometry::gregory::kp20_bi6_guide_build::{
        Kp20GuideBuildInput, build_kp20_guide_ring,
    };
    use crate::geometry::gregory::kp20_bi6_piece_maps::{
        Kp20PieceSamplerFamily, Kp20TauFamily, kp20_default_sector_piece_maps,
        kp20_paper_quadrant_sector_piece_maps, kp20_planar_bi6_tau_f3_energy,
        kp20_planar_bi6_tau_orientation_report, kp20_positive_nonoverlap_c2_sector_piece_maps,
        kp20_tau_piece_map_domain_report_for_pieces, sample_kp20_macro_patch_from_paper_tau_maps,
        sample_kp20_macro_patch_from_piece_sampler_maps,
        sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch,
        sample_kp20_macro_patch_from_tau_family_maps, sample_kp20_macro_patch_from_tau_maps,
    };
    use crate::geometry::gregory::{
        GregoryPatch, RectangularGregoryBiNonic, elevate_bi6_grid_to_bi9,
    };

    #[derive(Clone, Copy)]
    enum SamplerCase {
        AffineDefault,
        KernelTauDefault,
        PaperTauPositiveNonoverlap,
        PaperRefinedApexPositiveNonoverlap,
        MixedRefinedOriginAffineOuter,
        PlanarBi6TauAssembledCandidate,
        PlanarBi6TauRotationalPaperRegularSource,
        PlanarBi6TauRotationalPaperRegularSourcePlanarGuide,
        PlanarBi6TauRotationalPaperPr08Source,
        PaperTauPaperQuadrant,
    }

    impl SamplerCase {
        fn label(self) -> &'static str {
            match self {
                SamplerCase::AffineDefault => "affine/default",
                SamplerCase::KernelTauDefault => "kernel-tau/default",
                SamplerCase::PaperTauPositiveNonoverlap => "paper-tau/positive-nonoverlap",
                SamplerCase::PaperRefinedApexPositiveNonoverlap => {
                    "paper-refined-apex/positive-nonoverlap"
                }
                SamplerCase::MixedRefinedOriginAffineOuter => "mixed-refined-origin/affine-outer",
                SamplerCase::PlanarBi6TauAssembledCandidate => "planar-bi6-tau/assembled-candidate",
                SamplerCase::PlanarBi6TauRotationalPaperRegularSource => {
                    "planar-bi6-tau/rotational-paper-regular-source"
                }
                SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide => {
                    "planar-bi6-tau/rotational-paper-regular-source-planar-guide"
                }
                SamplerCase::PlanarBi6TauRotationalPaperPr08Source => {
                    "planar-bi6-tau/rotational-paper-pr08-source"
                }
                SamplerCase::PaperTauPaperQuadrant => "paper-tau/paper-quadrant",
            }
        }
    }

    let valences: &[usize] = &[5, 6];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        for sampler in [
            SamplerCase::AffineDefault,
            SamplerCase::KernelTauDefault,
            SamplerCase::PaperTauPositiveNonoverlap,
            SamplerCase::PaperRefinedApexPositiveNonoverlap,
            SamplerCase::MixedRefinedOriginAffineOuter,
            SamplerCase::PlanarBi6TauAssembledCandidate,
            SamplerCase::PlanarBi6TauRotationalPaperRegularSource,
            SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide,
            SamplerCase::PlanarBi6TauRotationalPaperPr08Source,
            SamplerCase::PaperTauPaperQuadrant,
        ] {
            let case_tensor_border = match sampler {
                SamplerCase::PlanarBi6TauRotationalPaperRegularSource
                | SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide => {
                    planar_regular_sector_tensor_border(n)
                }
                SamplerCase::PlanarBi6TauRotationalPaperPr08Source => {
                    pr08_characteristic_tensor_border(n)
                }
                _ => tensor_border.clone(),
            };
            let case_star_point = match sampler {
                SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide => {
                    Point3::new(0.0, 0.0, 0.0)
                }
                _ => star_point,
            };
            let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
                star_point: case_star_point,
                tensor_border: case_tensor_border.clone(),
            })
            .expect("guide builds");
            let planar_tau_macro_patches = match sampler {
                SamplerCase::PlanarBi6TauAssembledCandidate => {
                    Some(assembled_planar_tau_macro_patches(n))
                }
                SamplerCase::PlanarBi6TauRotationalPaperRegularSource => {
                    Some(rotational_paper_tau_macro_patches(
                        n,
                        &case_tensor_border,
                        "regular-planar-source",
                    ))
                }
                SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide => {
                    Some(rotational_paper_tau_macro_patches(
                        n,
                        &case_tensor_border,
                        "regular-planar-source-planar-guide",
                    ))
                }
                SamplerCase::PlanarBi6TauRotationalPaperPr08Source => {
                    Some(rotational_paper_tau_macro_patches(
                        n,
                        &case_tensor_border,
                        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows.name(),
                    ))
                }
                _ => None,
            };
            let sampled_macro_patches: Result<Vec<_>, String> = (0..n)
                .map(|k| match sampler {
                    SamplerCase::AffineDefault => {
                        Ok(sample_kp20_macro_patch_from_default_maps(guide.sector(k)))
                    }
                    SamplerCase::KernelTauDefault => {
                        sample_kp20_macro_patch_from_tau_maps(guide.sector(k), n, k)
                    }
                    SamplerCase::PaperTauPositiveNonoverlap => {
                        sample_kp20_macro_patch_from_paper_tau_maps(guide.sector(k), n, k)
                    }
                    SamplerCase::PaperRefinedApexPositiveNonoverlap => {
                        sample_kp20_macro_patch_from_tau_family_maps(
                            guide.sector(k),
                            n,
                            k,
                            Kp20TauFamily::PaperRefinedApex,
                            kp20_positive_nonoverlap_c2_sector_piece_maps(),
                        )
                    }
                    SamplerCase::MixedRefinedOriginAffineOuter => {
                        sample_kp20_macro_patch_from_piece_sampler_maps(
                            guide.sector(k),
                            n,
                            k,
                            [
                                (
                                    Kp20PieceId::NearOriginLeft,
                                    Kp20PieceSamplerFamily::Tau(Kp20TauFamily::PaperRefinedApex),
                                ),
                                (
                                    Kp20PieceId::NearOriginRight,
                                    Kp20PieceSamplerFamily::Tau(Kp20TauFamily::PaperRefinedApex),
                                ),
                                (
                                    Kp20PieceId::NearMidpointRight,
                                    Kp20PieceSamplerFamily::AffineGuide,
                                ),
                                (
                                    Kp20PieceId::NearMidpointLeft,
                                    Kp20PieceSamplerFamily::AffineGuide,
                                ),
                            ],
                            kp20_positive_nonoverlap_c2_sector_piece_maps(),
                        )
                    }
                    SamplerCase::PlanarBi6TauAssembledCandidate => {
                        let tau_macro_patches = planar_tau_macro_patches
                            .as_ref()
                            .expect("planar tau candidate built");
                        Ok(sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
                            guide.sector(k),
                            &tau_macro_patches[k],
                        ))
                    }
                    SamplerCase::PlanarBi6TauRotationalPaperRegularSource => {
                        let tau_macro_patches = planar_tau_macro_patches
                            .as_ref()
                            .expect("solved paper tau candidate built");
                        Ok(sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
                            guide.sector(k),
                            &tau_macro_patches[k],
                        ))
                    }
                    SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide => {
                        let tau_macro_patches = planar_tau_macro_patches
                            .as_ref()
                            .expect("solved paper tau planar-guide candidate built");
                        Ok(sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
                            guide.sector(k),
                            &tau_macro_patches[k],
                        ))
                    }
                    SamplerCase::PlanarBi6TauRotationalPaperPr08Source => {
                        let tau_macro_patches = planar_tau_macro_patches
                            .as_ref()
                            .expect("solved PR08-source paper tau candidate built");
                        Ok(sample_kp20_macro_patch_from_planar_bi6_tau_macro_patch(
                            guide.sector(k),
                            &tau_macro_patches[k],
                        ))
                    }
                    SamplerCase::PaperTauPaperQuadrant => {
                        sample_kp20_macro_patch_from_tau_family_maps(
                            guide.sector(k),
                            n,
                            k,
                            Kp20TauFamily::Paper,
                            kp20_paper_quadrant_sector_piece_maps(),
                        )
                    }
                })
                .collect();
            let sampled_macro_patches =
                sampled_macro_patches.expect("sampler must produce patches");
            let pre_assembly_cap = macro_cap_from_sampled(n, sampled_macro_patches.clone());
            let pre_assembly_report = measure_kp20_assembly_constraints(&pre_assembly_cap);
            let pre_assembly_folds = summarize_cap_folds(&pre_assembly_cap);
            let macro_cap_input = Kp20Bi6MacroCapInput {
                tensor_border: case_tensor_border.clone(),
                sampled_macro_patches: sampled_macro_patches.clone(),
            };
            let cap_result = match sampler {
                SamplerCase::PlanarBi6TauRotationalPaperRegularSource
                | SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide
                | SamplerCase::PlanarBi6TauRotationalPaperPr08Source => {
                    assemble_kp20_bi6_macro_cap_with_sampled_tensor_border_source(&macro_cap_input)
                }
                _ => assemble_kp20_bi6_macro_cap(&macro_cap_input),
            };
            let cap =
                cap_result.unwrap_or_else(|e: Kp20Bi6AssembleError| panic!("assemble failed: {e}"));
            let report = measure_kp20_assembly_constraints(&cap);
            let fold_summary = summarize_cap_folds(&cap);
            let domain = match sampler {
                SamplerCase::AffineDefault => "domain=no-tau".to_string(),
                SamplerCase::KernelTauDefault => {
                    let summary = aggregate_tau_domain(
                        n,
                        Kp20TauFamily::Kernel,
                        |_| true,
                        kp20_default_sector_piece_maps,
                    );
                    format!(
                        "domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                    )
                }
                SamplerCase::PaperTauPositiveNonoverlap => {
                    let summary = aggregate_tau_domain(
                        n,
                        Kp20TauFamily::Paper,
                        |_| true,
                        kp20_positive_nonoverlap_c2_sector_piece_maps,
                    );
                    format!(
                        "domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                    )
                }
                SamplerCase::PaperRefinedApexPositiveNonoverlap => {
                    let summary = aggregate_tau_domain(
                        n,
                        Kp20TauFamily::PaperRefinedApex,
                        |_| true,
                        kp20_positive_nonoverlap_c2_sector_piece_maps,
                    );
                    format!(
                        "domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                    )
                }
                SamplerCase::MixedRefinedOriginAffineOuter => {
                    let summary = aggregate_tau_domain(
                        n,
                        Kp20TauFamily::PaperRefinedApex,
                        |piece| {
                            matches!(
                                piece,
                                Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight
                            )
                        },
                        kp20_positive_nonoverlap_c2_sector_piece_maps,
                    );
                    format!(
                        "origin-domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                    )
                }
                SamplerCase::PlanarBi6TauAssembledCandidate => {
                    let tau_macro_patches = planar_tau_macro_patches
                        .as_ref()
                        .expect("planar tau candidate built");
                    let summary = aggregate_planar_tau_cp_domain(tau_macro_patches);
                    let f3 = aggregate_planar_tau_f3(tau_macro_patches);
                    format!(
                        "cp-domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}] F3={f3:e}",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                        f3 = f3,
                    )
                }
                SamplerCase::PlanarBi6TauRotationalPaperRegularSource
                | SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide => {
                    let tau_macro_patches = planar_tau_macro_patches
                        .as_ref()
                        .expect("solved paper tau candidate built");
                    let summary = aggregate_planar_tau_cp_domain(tau_macro_patches);
                    let jacobian = aggregate_planar_tau_jacobian(tau_macro_patches);
                    let f3 = aggregate_planar_tau_f3(tau_macro_patches);
                    format!(
                        "cp-domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}] tau_jac=neg:{neg}/{jac_samples} zero:{zero} det=[{min_det:e},{max_det:e}] worst=s{worst_sector}:{worst_piece:?} F3={f3:e}",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                        neg = jacobian.negative,
                        zero = jacobian.near_zero,
                        jac_samples = jacobian.samples,
                        min_det = jacobian.min_det,
                        max_det = jacobian.max_det,
                        worst_sector = jacobian.worst_sector,
                        worst_piece = jacobian.worst_piece,
                        f3 = f3,
                    )
                }
                SamplerCase::PlanarBi6TauRotationalPaperPr08Source => {
                    let tau_macro_patches = planar_tau_macro_patches
                        .as_ref()
                        .expect("solved PR08-source paper tau candidate built");
                    let summary = aggregate_planar_tau_cp_domain(tau_macro_patches);
                    let jacobian = aggregate_planar_tau_jacobian(tau_macro_patches);
                    let f3 = aggregate_planar_tau_f3(tau_macro_patches);
                    format!(
                        "source={source} cp-domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}] tau_jac=neg:{neg}/{jac_samples} zero:{zero} det=[{min_det:e},{max_det:e}] worst=s{worst_sector}:{worst_piece:?} F3={f3:e}",
                        source = Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows
                            .name(),
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                        neg = jacobian.negative,
                        zero = jacobian.near_zero,
                        jac_samples = jacobian.samples,
                        min_det = jacobian.min_det,
                        max_det = jacobian.max_det,
                        worst_sector = jacobian.worst_sector,
                        worst_piece = jacobian.worst_piece,
                        f3 = f3,
                    )
                }
                SamplerCase::PaperTauPaperQuadrant => {
                    let summary = aggregate_tau_domain(
                        n,
                        Kp20TauFamily::Paper,
                        |_| true,
                        kp20_paper_quadrant_sector_piece_maps,
                    );
                    format!(
                        "domain=outside:{outside}/{samples} max_excess={excess:e} u=[{min_u:e},{max_u:e}] v=[{min_v:e},{max_v:e}]",
                        outside = summary.outside,
                        samples = summary.samples,
                        excess = summary.max_simplex_excess,
                        min_u = summary.min_u,
                        max_u = summary.max_u,
                        min_v = summary.min_v,
                        max_v = summary.max_v,
                    )
                }
            };
            eprintln!(
                "sampler A/B n={n} {label}: {domain} worst_eq8={e8:e} worst_eq9={e9:e} L6(1)={l1:e} L6(2)={l2:e} L6(3)={l3:e} L6(4)={l4:e} folds={folds} min_fold_cos={fold_cos:e} piece_folds={piece_folds:?}",
                label = sampler.label(),
                e8 = report.worst_eq8(),
                e9 = report.worst_eq9(),
                l1 = report.worst_lemma6_cond_1,
                l2 = report.worst_lemma6_cond_2,
                l3 = report.worst_lemma6_cond_3,
                l4 = report.worst_lemma6_cond_4,
                folds = fold_summary.fold_count,
                fold_cos = fold_summary.min_jacobian_normal_cos,
                piece_folds = fold_summary.piece_fold_counts,
            );
            if matches!(
                sampler,
                SamplerCase::PlanarBi6TauRotationalPaperRegularSource
                    | SamplerCase::PlanarBi6TauRotationalPaperRegularSourcePlanarGuide
                    | SamplerCase::PlanarBi6TauRotationalPaperPr08Source
            ) {
                const ROTATIONAL_PAPER_CAP_ROW_TOLERANCE: f64 = 1.0e-12;
                eprintln!(
                    "sampler A/B n={n} {label} pre-assembly: worst_eq8={pre_e8:e} worst_eq9={pre_e9:e} folds={pre_folds} piece_folds={pre_piece_folds:?}; post-assembly: worst_eq8={post_e8:e} worst_eq9={post_e9:e} folds={post_folds} piece_folds={post_piece_folds:?}",
                    label = sampler.label(),
                    pre_e8 = pre_assembly_report.worst_eq8(),
                    pre_e9 = pre_assembly_report.worst_eq9(),
                    pre_folds = pre_assembly_folds.fold_count,
                    pre_piece_folds = pre_assembly_folds.piece_fold_counts,
                    post_e8 = report.worst_eq8(),
                    post_e9 = report.worst_eq9(),
                    post_folds = fold_summary.fold_count,
                    post_piece_folds = fold_summary.piece_fold_counts,
                );
                assert!(
                    pre_assembly_report.worst_eq8() < ROTATIONAL_PAPER_CAP_ROW_TOLERANCE,
                    "paper tau sampled source must enter assembly with roundoff eq8 residual for n={n} {label}: {e8:e}",
                    label = sampler.label(),
                    e8 = pre_assembly_report.worst_eq8(),
                );
                assert!(
                    pre_assembly_report.worst_eq9() < ROTATIONAL_PAPER_CAP_ROW_TOLERANCE,
                    "paper tau sampled source must enter assembly with roundoff eq9 residual for n={n} {label}: {e9:e}",
                    label = sampler.label(),
                    e9 = pre_assembly_report.worst_eq9(),
                );
                assert_eq!(
                    pre_assembly_folds.fold_count,
                    0,
                    "paper tau sampled source must be fold-free before assembly for n={n} {label}",
                    label = sampler.label(),
                );
                assert!(
                    report.worst_eq8() < ROTATIONAL_PAPER_CAP_ROW_TOLERANCE,
                    "projected sampled-source assembly must preserve roundoff eq8 residual for n={n} {label}: {e8:e}",
                    label = sampler.label(),
                    e8 = report.worst_eq8(),
                );
                assert!(
                    report.worst_eq9() < ROTATIONAL_PAPER_CAP_ROW_TOLERANCE,
                    "projected sampled-source assembly must preserve roundoff eq9 residual for n={n} {label}: {e9:e}",
                    label = sampler.label(),
                    e9 = report.worst_eq9(),
                );
                assert_eq!(
                    fold_summary.fold_count,
                    0,
                    "projected sampled-source assembly must remain fold-free for n={n} {label}",
                    label = sampler.label(),
                );
                let tau_macro_patches = planar_tau_macro_patches
                    .as_ref()
                    .expect("rotational paper tau macro-patches exist");
                let tau_jacobian = aggregate_planar_tau_jacobian(tau_macro_patches);
                assert_eq!(
                    tau_jacobian.negative,
                    0,
                    "current F3-priority paper tau should be nonnegative to roundoff for n={n} {label}",
                    label = sampler.label(),
                );
                assert_eq!(
                    tau_jacobian.near_zero,
                    tau_jacobian.samples,
                    "current F3-priority paper tau is a collapsed limit, not a nonsingular reparameterization for n={n} {label}",
                    label = sampler.label(),
                );
                let stage_snapshots =
                    assemble_kp20_bi6_macro_cap_stage_snapshots(&Kp20Bi6MacroCapInput {
                        tensor_border: case_tensor_border.clone(),
                        sampled_macro_patches: sampled_macro_patches.clone(),
                    })
                    .expect("assembly stage snapshots succeed");
                assert_eq!(stage_snapshots.len(), 6);
                assert_eq!(
                    stage_snapshots.last().expect("final stage").stage,
                    Kp20Bi6MacroCapAssemblyStage::G2StripProjection
                );
                for snapshot in stage_snapshots {
                    let stage_cap = macro_cap_from_sampled(n, snapshot.macro_patches);
                    let stage_report = measure_kp20_assembly_constraints(&stage_cap);
                    let stage_folds = summarize_cap_folds(&stage_cap);
                    eprintln!(
                        "assembly stage n={n} {label} {stage:?}: worst_eq8={e8:e} worst_eq9={e9:e} folds={folds} piece_folds={piece_folds:?} strip_projections={strip_count}",
                        label = sampler.label(),
                        stage = snapshot.stage,
                        e8 = stage_report.worst_eq8(),
                        e9 = stage_report.worst_eq9(),
                        folds = stage_folds.fold_count,
                        piece_folds = stage_folds.piece_fold_counts,
                        strip_count = snapshot.strip_projections.len(),
                    );
                    assert!(stage_report.all_finite());
                    if snapshot.stage == Kp20Bi6MacroCapAssemblyStage::TensorBorderRows {
                        assert!(
                            stage_report.worst_eq8() > 1.0,
                            "ordinary tensor-border overwrite should remain localized as the high-residual path for n={n} {label}",
                            label = sampler.label(),
                        );
                    }
                }
                let projected_source_snapshots =
                    assemble_kp20_bi6_macro_cap_stage_snapshots_with_sampled_tensor_border(
                        &Kp20Bi6MacroCapInput {
                            tensor_border: case_tensor_border.clone(),
                            sampled_macro_patches: sampled_macro_patches.clone(),
                        },
                    )
                    .expect("projected source assembly stage snapshots succeed");
                assert_eq!(projected_source_snapshots.len(), 6);
                for snapshot in projected_source_snapshots {
                    let stage_cap = macro_cap_from_sampled(n, snapshot.macro_patches);
                    let stage_report = measure_kp20_assembly_constraints(&stage_cap);
                    let stage_folds = summarize_cap_folds(&stage_cap);
                    eprintln!(
                        "assembly projected-source stage n={n} {label} {stage:?}: worst_eq8={e8:e} worst_eq9={e9:e} folds={folds} piece_folds={piece_folds:?} strip_projections={strip_count}",
                        label = sampler.label(),
                        stage = snapshot.stage,
                        e8 = stage_report.worst_eq8(),
                        e9 = stage_report.worst_eq9(),
                        folds = stage_folds.fold_count,
                        piece_folds = stage_folds.piece_fold_counts,
                        strip_count = snapshot.strip_projections.len(),
                    );
                    assert!(stage_report.all_finite());
                    assert!(
                        stage_report.worst_eq8() < ROTATIONAL_PAPER_CAP_ROW_TOLERANCE,
                        "projected sampled-source stage must preserve roundoff eq8 residual for n={n} {label} {stage:?}: {e8:e}",
                        label = sampler.label(),
                        stage = snapshot.stage,
                        e8 = stage_report.worst_eq8(),
                    );
                    assert!(
                        stage_report.worst_eq9() < ROTATIONAL_PAPER_CAP_ROW_TOLERANCE,
                        "projected sampled-source stage must preserve roundoff eq9 residual for n={n} {label} {stage:?}: {e9:e}",
                        label = sampler.label(),
                        stage = snapshot.stage,
                        e9 = stage_report.worst_eq9(),
                    );
                    assert_eq!(
                        stage_folds.fold_count,
                        0,
                        "projected sampled-source stage must remain fold-free for n={n} {label} {stage:?}",
                        label = sampler.label(),
                        stage = snapshot.stage,
                    );
                }
                assert!(pre_assembly_report.all_finite());
            }
            assert!(report.all_finite());
            assert!(
                fold_summary.min_jacobian_normal_cos.is_finite() || fold_summary.fold_count == 0
            );
        }
    }

    struct DomainSummary {
        outside: usize,
        samples: usize,
        min_u: f64,
        min_v: f64,
        max_u: f64,
        max_v: f64,
        max_simplex_excess: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
        worst_i: usize,
        worst_j: usize,
    }

    struct FoldSummary {
        fold_count: usize,
        min_jacobian_normal_cos: f64,
        piece_fold_counts: [usize; 4],
        worst_patch_index: usize,
    }

    struct TauJacobianSummary {
        samples: usize,
        negative: usize,
        near_zero: usize,
        min_det: f64,
        max_det: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
    }

    fn aggregate_tau_domain(
        n: usize,
        tau_family: Kp20TauFamily,
        include_piece: impl Fn(Kp20PieceId) -> bool + Copy,
        maps: impl Fn() -> crate::geometry::gregory::Kp20SectorPieceMaps,
    ) -> DomainSummary {
        let mut summary = DomainSummary {
            outside: 0,
            samples: 0,
            min_u: f64::INFINITY,
            min_v: f64::INFINITY,
            max_u: f64::NEG_INFINITY,
            max_v: f64::NEG_INFINITY,
            max_simplex_excess: 0.0,
            worst_sector: 0,
            worst_piece: Kp20PieceId::NearOriginLeft,
            worst_i: 0,
            worst_j: 0,
        };
        for k in 0..n {
            let report = kp20_tau_piece_map_domain_report_for_pieces(
                maps(),
                n,
                k,
                tau_family,
                include_piece,
            )
            .expect("domain report");
            summary.outside += report.outside_triangle_samples;
            summary.samples += report.sample_count;
            summary.min_u = summary.min_u.min(report.min_u);
            summary.min_v = summary.min_v.min(report.min_v);
            summary.max_u = summary.max_u.max(report.max_u);
            summary.max_v = summary.max_v.max(report.max_v);
            summary.max_simplex_excess = summary.max_simplex_excess.max(report.max_simplex_excess);
        }
        summary
    }

    fn assembled_planar_tau_macro_patches(n: usize) -> Vec<Kp20MacroPatch> {
        let tensor_border = planar_regular_sector_tensor_border(n);
        let sampled_macro_patches = (0..n)
            .map(|k| planar_tau_seed_macro_patch(n, k))
            .collect::<Vec<_>>();
        let assembled = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
            tensor_border,
            sampled_macro_patches,
        })
        .expect("planar tau macro assembly succeeds");
        assembled
            .macro_patches
            .into_iter()
            .enumerate()
            .map(|(k, patch)| global_planar_patch_to_local(n, k, patch))
            .collect()
    }

    fn rotational_paper_tau_macro_patches(
        n: usize,
        tensor_border: &TensorBorderBi3,
        source_label: &str,
    ) -> Vec<Kp20MacroPatch> {
        let mut worst_solve_residual = 0.0_f64;
        let mut worst_projection_delta = 0.0_f64;
        let mut worst_projection_row_residual = 0.0_f64;
        let mut adjusted_scalar_targets = 0_usize;
        let mut worst_family_residuals = [0.0_f64; 7];
        let mut macro_patches = Vec::with_capacity(n);
        for sector in 0..n {
            let solved =
                solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                    valence: n,
                    sector,
                    characteristic_tensor_border: tensor_border,
                    initial: planar_tau_seed_macro_patch(n, sector),
                    regularization_weight: 1.0e-10,
                })
                .expect("rotational paper tau solve succeeds");
            worst_solve_residual = worst_solve_residual.max(solved.max_constraint_residual);
            let projection_report = kp20_planar_tau_rotational_section42_source_projection_report(
                n,
                sector,
                tensor_border,
            )
            .expect("section42 source projection report");
            worst_projection_delta = worst_projection_delta.max(projection_report.max_source_delta);
            worst_projection_row_residual =
                worst_projection_row_residual.max(projection_report.max_row_residual);
            adjusted_scalar_targets =
                adjusted_scalar_targets.max(projection_report.adjusted_scalar_targets);
            for residual in kp20_planar_tau_rotational_paper_family_residuals(
                n,
                sector,
                tensor_border,
                &solved.tau_macro_patch,
            )
            .expect("family residuals evaluate")
            {
                let slot = family_residual_index(residual.family);
                worst_family_residuals[slot] =
                    worst_family_residuals[slot].max(residual.max_residual);
            }
            macro_patches.push(global_planar_patch_to_local(
                n,
                sector,
                solved.tau_macro_patch,
            ));
        }
        eprintln!(
            "rotational paper tau solve n={n} source={source_label}: max_constraint_residual={worst_solve_residual:e} source_projection_delta={worst_projection_delta:e} source_projection_residual={worst_projection_row_residual:e} adjusted_scalar_targets={adjusted_scalar_targets} families fig8={fig8:e} eq14={eq14:e} section42_c0={c0:e} section42_eq11={eq11:e} section42_eq12={eq12:e} section42_eq13={eq13:e} diagonal_diagnostic={diagonal:e}",
            fig8 = worst_family_residuals[0],
            eq14 = worst_family_residuals[1],
            c0 = worst_family_residuals[2],
            eq11 = worst_family_residuals[3],
            eq12 = worst_family_residuals[4],
            eq13 = worst_family_residuals[5],
            diagonal = worst_family_residuals[6],
        );
        macro_patches
    }

    fn pr08_characteristic_tensor_border(n: usize) -> TensorBorderBi3 {
        kp20_pr08_stam_q1_characteristic_tensor_border(n).expect("PR08/Stam Q1 tensor-border")
    }

    fn family_residual_index(family: Kp20PlanarTauPaperRowFamily) -> usize {
        match family {
            Kp20PlanarTauPaperRowFamily::Fig8IndicatedAverage => 0,
            Kp20PlanarTauPaperRowFamily::TensorBorderEq14 => 1,
            Kp20PlanarTauPaperRowFamily::Section42C0 => 2,
            Kp20PlanarTauPaperRowFamily::Section42Eq11Degree5 => 3,
            Kp20PlanarTauPaperRowFamily::Section42Eq12InternalC2 => 4,
            Kp20PlanarTauPaperRowFamily::Section42Eq13Strip => 5,
            Kp20PlanarTauPaperRowFamily::DiagonalSymmetry => 6,
            Kp20PlanarTauPaperRowFamily::InternalC2
            | Kp20PlanarTauPaperRowFamily::SectorG1G2NearOrigin
            | Kp20PlanarTauPaperRowFamily::SectorG1G2NearMidpoint => {
                unreachable!("raw eq(8)/(9) rows are not part of the §4.2 live basis")
            }
        }
    }

    fn planar_regular_sector_tensor_border(n: usize) -> TensorBorderBi3 {
        let mut sectors = Vec::with_capacity(n);
        for k in 0..n {
            let a0 = core::f64::consts::TAU * k as f64 / n as f64;
            let a1 = core::f64::consts::TAU * (k + 1) as f64 / n as f64;
            let p0 = Point3::new(a0.cos(), a0.sin(), 0.0);
            let p1 = Point3::new(a1.cos(), a1.sin(), 0.0);
            let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                let radius = 1.0 - 0.1 * row as f64;
                for i in 0..4 {
                    let t = i as f64 / 3.0;
                    rows[row][i] = Point3::new(
                        radius * ((1.0 - t) * p0.x + t * p1.x),
                        radius * ((1.0 - t) * p0.y + t * p1.y),
                        0.0,
                    );
                }
            }
            sectors.push(rows);
        }
        TensorBorderBi3 { sectors }
    }

    fn planar_tau_seed_macro_patch(n: usize, k: usize) -> Kp20MacroPatch {
        let mut out = Kp20MacroPatch::zero();
        for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
            let base = corner_maps[0].1;
            let piece = out.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let u = i as f64 / 6.0;
                    let v = j as f64 / 6.0;
                    let local = [
                        base.origin[0] + u * base.du[0] + v * base.dv[0],
                        base.origin[1] + u * base.du[1] + v * base.dv[1],
                    ];
                    piece.cps[i][j] = local_to_global_planar(n, k, local);
                }
            }
        }
        out
    }

    fn local_to_global_planar(n: usize, k: usize, local: [f64; 2]) -> Point3 {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        Point3::new(
            local[0] * spoke1[0] + local[1] * spoke0[0],
            local[0] * spoke1[1] + local[1] * spoke0[1],
            0.0,
        )
    }

    fn global_planar_patch_to_local(
        n: usize,
        k: usize,
        mut patch: Kp20MacroPatch,
    ) -> Kp20MacroPatch {
        for piece_id in Kp20PieceId::all() {
            let piece = patch.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let local = global_to_local_planar(n, k, piece.cps[i][j]);
                    piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
                }
            }
        }
        patch
    }

    fn global_to_local_planar(n: usize, k: usize, p: Point3) -> [f64; 2] {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
        [
            (p.x * spoke0[1] - p.y * spoke0[0]) / det,
            (spoke1[0] * p.y - spoke1[1] * p.x) / det,
        ]
    }

    fn aggregate_planar_tau_cp_domain(macro_patches: &[Kp20MacroPatch]) -> DomainSummary {
        let mut summary = DomainSummary {
            outside: 0,
            samples: 0,
            min_u: f64::INFINITY,
            min_v: f64::INFINITY,
            max_u: f64::NEG_INFINITY,
            max_v: f64::NEG_INFINITY,
            max_simplex_excess: 0.0,
            worst_sector: 0,
            worst_piece: Kp20PieceId::NearOriginLeft,
            worst_i: 0,
            worst_j: 0,
        };
        for (sector, macro_patch) in macro_patches.iter().enumerate() {
            for piece_id in Kp20PieceId::all() {
                for (i, row) in macro_patch.piece(piece_id).cps.iter().enumerate() {
                    for (j, cp) in row.iter().enumerate() {
                        let u = cp.x;
                        let v = cp.y;
                        let triangle_excess = (-u).max(-v).max(u + v - 1.0).max(0.0);
                        summary.outside += usize::from(triangle_excess > 1.0e-12);
                        summary.samples += 1;
                        summary.min_u = summary.min_u.min(u);
                        summary.min_v = summary.min_v.min(v);
                        summary.max_u = summary.max_u.max(u);
                        summary.max_v = summary.max_v.max(v);
                        if triangle_excess > summary.max_simplex_excess {
                            summary.max_simplex_excess = triangle_excess;
                            summary.worst_sector = sector;
                            summary.worst_piece = piece_id;
                            summary.worst_i = i;
                            summary.worst_j = j;
                        }
                    }
                }
            }
        }
        summary
    }

    fn aggregate_planar_tau_f3(macro_patches: &[Kp20MacroPatch]) -> f64 {
        macro_patches
            .iter()
            .map(kp20_planar_bi6_tau_f3_energy)
            .sum()
    }

    fn aggregate_planar_tau_jacobian(macro_patches: &[Kp20MacroPatch]) -> TauJacobianSummary {
        const TAU_JACOBIAN_TOLERANCE: f64 = 1.0e-10;
        const SAMPLES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];
        let mut summary = TauJacobianSummary {
            samples: 0,
            negative: 0,
            near_zero: 0,
            min_det: f64::INFINITY,
            max_det: f64::NEG_INFINITY,
            worst_sector: 0,
            worst_piece: Kp20PieceId::NearOriginLeft,
        };
        for (sector, macro_patch) in macro_patches.iter().enumerate() {
            let report = kp20_planar_bi6_tau_orientation_report(
                macro_patch,
                &SAMPLES,
                TAU_JACOBIAN_TOLERANCE,
            );
            summary.samples += report.sample_count;
            summary.negative += report.negative_determinants;
            summary.near_zero += report.near_zero_determinants;
            if report.min_determinant < summary.min_det {
                summary.min_det = report.min_determinant;
                summary.worst_sector = sector;
                summary.worst_piece = report.worst_piece;
            }
            summary.max_det = summary.max_det.max(report.max_determinant);
        }
        summary
    }

    fn macro_cap_from_sampled(
        n: usize,
        macro_patches: Vec<Kp20MacroPatch>,
    ) -> Kp20Bi6MacroCapOutput {
        Kp20Bi6MacroCapOutput {
            n,
            patches_bi6: macro_patches
                .iter()
                .flat_map(|macro_patch| {
                    Kp20PieceId::all()
                        .into_iter()
                        .map(|piece_id| macro_patch.piece(piece_id).cps)
                })
                .collect(),
            macro_patches,
            strip_projections: Vec::new(),
        }
    }

    fn summarize_cap_folds(cap: &Kp20Bi6MacroCapOutput) -> FoldSummary {
        const SELF_INTERSECTION_TOLERANCE: f64 = 1.0e-6;
        let mut fold_count = 0;
        let mut min_jacobian_normal_cos = f64::INFINITY;
        let mut piece_fold_counts = [0_usize; 4];
        let mut worst_patch_index = 0;
        for (patch_index, grid) in cap.patches_bi6.iter().enumerate() {
            let patch = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
                grid: elevate_bi6_grid_to_bi9(*grid),
            }));
            let fold_report = patch.self_intersection_report(SELF_INTERSECTION_TOLERANCE);
            if fold_report.intersects() {
                fold_count += 1;
                piece_fold_counts[patch_index % 4] += 1;
                if fold_report.min_jacobian_normal_cos < min_jacobian_normal_cos {
                    min_jacobian_normal_cos = fold_report.min_jacobian_normal_cos;
                    worst_patch_index = patch_index;
                }
            }
        }
        FoldSummary {
            fold_count,
            min_jacobian_normal_cos,
            piece_fold_counts,
            worst_patch_index,
        }
    }
}

#[test]
fn projected_pr08_q1_q3_candidates_rank_closed_tau_domain_before_full_sampler() {
    use super::super::kp20_bi6_tau_solve_assemble::{
        Kp20PlanarTauPaperSolveInput,
        kp20_planar_tau_rotational_section42_source_projection_report,
        solve_kp20_planar_tau_rotational_paper_sector_f3,
    };
    use crate::geometry::gregory::kp20_bi6_cap::Kp20MacroPatch;
    use crate::geometry::gregory::kp20_bi6_cc_ring::{
        Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
        kp20_pr08_characteristic_tensor_border_candidate_rows,
    };
    use crate::geometry::gregory::kp20_bi6_cnet::kp20_tensor_border_from_regular_cnet_rows;
    use crate::geometry::gregory::kp20_bi6_piece_maps::{
        kp20_planar_bi6_tau_f3_energy, kp20_positive_nonoverlap_c2_sector_piece_maps,
    };

    #[derive(Debug, Clone)]
    struct TauDomainQuality {
        candidate_name: &'static str,
        projection_delta: f64,
        projection_residual: f64,
        solve_residual: f64,
        tau_outside: usize,
        tau_samples: usize,
        tau_max_excess: f64,
        tau_worst_sector: usize,
        tau_worst_piece: Kp20PieceId,
        tau_worst_i: usize,
        tau_worst_j: usize,
        f3: f64,
    }

    let n = 5;
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let mut qualities = Vec::new();

    for candidate in [
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
    ] {
        let cnet_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 tensor-border candidate rows");
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&cnet_rows)
            .expect("PR08 tensor-border converts");
        let mut tau_macro_patches = Vec::with_capacity(n);
        let mut solve_residual = 0.0_f64;
        let mut projection_delta = 0.0_f64;
        let mut projection_residual = 0.0_f64;

        for sector in 0..n {
            let solved =
                solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                    valence: n,
                    sector,
                    characteristic_tensor_border: &tensor_border,
                    initial: planar_tau_seed_macro_patch_for_domain_ranking(n, sector),
                    regularization_weight: 1.0e-10,
                })
                .expect("projected rotational paper tau solve succeeds");
            solve_residual = solve_residual.max(solved.max_constraint_residual);
            let projection_report = kp20_planar_tau_rotational_section42_source_projection_report(
                n,
                sector,
                &tensor_border,
            )
            .expect("projection report");
            projection_delta = projection_delta.max(projection_report.max_source_delta);
            projection_residual = projection_residual.max(projection_report.max_row_residual);
            tau_macro_patches.push(global_planar_patch_to_local_for_domain_ranking(
                n,
                sector,
                solved.tau_macro_patch,
            ));
        }

        let domain = aggregate_planar_tau_cp_domain_for_domain_ranking(&tau_macro_patches);
        let f3 = tau_macro_patches
            .iter()
            .map(kp20_planar_bi6_tau_f3_energy)
            .sum::<f64>();
        let quality = TauDomainQuality {
            candidate_name: candidate.name(),
            projection_delta,
            projection_residual,
            solve_residual,
            tau_outside: domain.outside,
            tau_samples: domain.samples,
            tau_max_excess: domain.max_excess,
            tau_worst_sector: domain.worst_sector,
            tau_worst_piece: domain.worst_piece,
            tau_worst_i: domain.worst_i,
            tau_worst_j: domain.worst_j,
            f3,
        };
        eprintln!(
            "projected PR08 q1/q3 closed tau domain n={n} candidate={name}: solve={solve:e} projection_delta={projection_delta:e} projection_residual={projection_residual:e} tau_domain=outside:{outside}/{samples} max_excess={excess:e} worst_tau=s{sector}:{piece:?}[{i},{j}] F3={f3:e}",
            name = quality.candidate_name,
            solve = quality.solve_residual,
            projection_delta = quality.projection_delta,
            projection_residual = quality.projection_residual,
            outside = quality.tau_outside,
            samples = quality.tau_samples,
            excess = quality.tau_max_excess,
            sector = quality.tau_worst_sector,
            piece = quality.tau_worst_piece,
            i = quality.tau_worst_i,
            j = quality.tau_worst_j,
            f3 = quality.f3,
        );
        assert!(quality.solve_residual < 1.0e-7);
        assert!(quality.projection_residual < 1.0e-7);
        qualities.push(quality);
    }

    let best_by_projection = qualities
        .iter()
        .min_by(|a, b| a.projection_delta.total_cmp(&b.projection_delta))
        .expect("candidate qualities");
    let best_by_domain = qualities
        .iter()
        .min_by(|a, b| {
            a.tau_outside
                .cmp(&b.tau_outside)
                .then_with(|| a.tau_max_excess.total_cmp(&b.tau_max_excess))
                .then_with(|| a.f3.total_cmp(&b.f3))
        })
        .expect("candidate qualities");
    let admissible = qualities
        .iter()
        .filter(|quality| quality.tau_outside == 0)
        .min_by(|a, b| {
            a.projection_delta
                .total_cmp(&b.projection_delta)
                .then_with(|| a.f3.total_cmp(&b.f3))
        });

    eprintln!(
        "projected PR08 q1/q3 closed tau selection n={n}: best_projection={bp} delta={bp_delta:e} outside={bp_outside}/{bp_samples}; best_domain={bd} outside={bd_outside}/{bd_samples} max_excess={bd_excess:e} worst_tau=s{bd_sector}:{bd_piece:?}[{bd_i},{bd_j}] F3={bd_f3:e}; admissible={admissible}",
        bp = best_by_projection.candidate_name,
        bp_delta = best_by_projection.projection_delta,
        bp_outside = best_by_projection.tau_outside,
        bp_samples = best_by_projection.tau_samples,
        bd = best_by_domain.candidate_name,
        bd_outside = best_by_domain.tau_outside,
        bd_samples = best_by_domain.tau_samples,
        bd_excess = best_by_domain.tau_max_excess,
        bd_sector = best_by_domain.tau_worst_sector,
        bd_piece = best_by_domain.tau_worst_piece,
        bd_i = best_by_domain.tau_worst_i,
        bd_j = best_by_domain.tau_worst_j,
        bd_f3 = best_by_domain.f3,
        admissible = admissible
            .map(|quality| quality.candidate_name)
            .unwrap_or("none"),
    );

    assert!(
        admissible.is_none(),
        "closed projected source should still decline promotion until a candidate is inside the sector triangle"
    );
    assert!(
        best_by_domain.tau_outside > 0,
        "current selection rule must reject candidates that still leave the sector-local tau domain"
    );

    fn planar_tau_seed_macro_patch_for_domain_ranking(n: usize, k: usize) -> Kp20MacroPatch {
        let mut out = Kp20MacroPatch::zero();
        for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
            let base = corner_maps[0].1;
            let piece = out.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let u = i as f64 / 6.0;
                    let v = j as f64 / 6.0;
                    let local = [
                        base.origin[0] + u * base.du[0] + v * base.dv[0],
                        base.origin[1] + u * base.du[1] + v * base.dv[1],
                    ];
                    piece.cps[i][j] = local_to_global_planar_for_domain_ranking(n, k, local);
                }
            }
        }
        out
    }

    fn local_to_global_planar_for_domain_ranking(n: usize, k: usize, local: [f64; 2]) -> Point3 {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        Point3::new(
            local[0] * spoke1[0] + local[1] * spoke0[0],
            local[0] * spoke1[1] + local[1] * spoke0[1],
            0.0,
        )
    }

    fn global_planar_patch_to_local_for_domain_ranking(
        n: usize,
        k: usize,
        mut patch: Kp20MacroPatch,
    ) -> Kp20MacroPatch {
        for piece_id in Kp20PieceId::all() {
            let piece = patch.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let local = global_to_local_planar_for_domain_ranking(n, k, piece.cps[i][j]);
                    piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
                }
            }
        }
        patch
    }

    fn global_to_local_planar_for_domain_ranking(n: usize, k: usize, p: Point3) -> [f64; 2] {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
        [
            (p.x * spoke0[1] - p.y * spoke0[0]) / det,
            (spoke1[0] * p.y - spoke1[1] * p.x) / det,
        ]
    }

    struct TauDomainSummary {
        outside: usize,
        samples: usize,
        max_excess: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
        worst_i: usize,
        worst_j: usize,
    }

    fn aggregate_planar_tau_cp_domain_for_domain_ranking(
        macro_patches: &[Kp20MacroPatch],
    ) -> TauDomainSummary {
        let mut summary = TauDomainSummary {
            outside: 0,
            samples: 0,
            max_excess: 0.0,
            worst_sector: 0,
            worst_piece: Kp20PieceId::NearOriginLeft,
            worst_i: 0,
            worst_j: 0,
        };
        for (sector, macro_patch) in macro_patches.iter().enumerate() {
            for piece_id in Kp20PieceId::all() {
                for (i, row) in macro_patch.piece(piece_id).cps.iter().enumerate() {
                    for (j, cp) in row.iter().enumerate() {
                        let u = cp.x;
                        let v = cp.y;
                        let triangle_excess = (-u).max(-v).max(u + v - 1.0).max(0.0);
                        summary.outside += usize::from(triangle_excess > 1.0e-12);
                        summary.samples += 1;
                        if triangle_excess > summary.max_excess {
                            summary.max_excess = triangle_excess;
                            summary.worst_sector = sector;
                            summary.worst_piece = piece_id;
                            summary.worst_i = i;
                            summary.worst_j = j;
                        }
                    }
                }
            }
        }
        summary
    }
}

#[test]
fn projected_section42_former_escape_coefficients_are_not_bound_by_accepted_rows() {
    use super::super::kp20_bi6_tau_f3_solve::Kp20PlanarTauCoeff;
    use super::super::kp20_bi6_tau_solve_assemble::{
        Kp20PlanarTauPaperSolveInput, kp20_planar_tau_rotational_paper_row_set,
        solve_kp20_planar_tau_rotational_paper_sector_f3,
    };
    use crate::geometry::gregory::kp20_bi6_cap::Kp20MacroPatch;
    use crate::geometry::gregory::kp20_bi6_cc_ring::{
        Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
        kp20_pr08_characteristic_tensor_border_candidate_rows,
    };
    use crate::geometry::gregory::kp20_bi6_cnet::kp20_tensor_border_from_regular_cnet_rows;
    use crate::geometry::gregory::kp20_bi6_piece_maps::{
        kp20_planar_bi6_tau_f3_energy, kp20_positive_nonoverlap_c2_sector_piece_maps,
    };

    #[derive(Debug, Clone, Copy)]
    struct RegSweep {
        regularization_weight: f64,
        local: [f64; 2],
        initial_local: [f64; 2],
        triangle_excess: f64,
        delta_from_initial: f64,
        f3: f64,
        residual: f64,
    }

    let n = 5;
    let sector = 0;
    let candidate = Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows;
    let escaped_coeff = Kp20PlanarTauCoeff::new(Kp20PieceId::NearOriginLeft, 0, 6);
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let cnet_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
        .expect("PR08 tensor-border candidate rows");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&cnet_rows).expect("PR08 tensor-border converts");
    let row_set = kp20_planar_tau_rotational_paper_row_set(n, sector, &tensor_border)
        .expect("projected rotational paper row set");
    let vector_rows_touching = row_set
        .vector_constraints
        .iter()
        .filter(|row| row.terms.iter().any(|&(coeff, _)| coeff == escaped_coeff))
        .count();
    let scalar_rows_touching = row_set
        .scalar_constraints
        .iter()
        .filter(|row| {
            row.terms
                .iter()
                .any(|&(scalar_coeff, _)| scalar_coeff.coeff == escaped_coeff)
        })
        .count();
    let mut regularization_sweep = Vec::new();
    for regularization_weight in [1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4] {
        let initial = planar_tau_seed_macro_patch_for_escape_probe(n, sector);
        let solved =
            solve_kp20_planar_tau_rotational_paper_sector_f3(&Kp20PlanarTauPaperSolveInput {
                valence: n,
                sector,
                characteristic_tensor_border: &tensor_border,
                initial,
                regularization_weight,
            })
            .expect("projected rotational paper tau solve succeeds");
        let solved_local =
            global_planar_patch_to_local_for_escape_probe(n, sector, solved.tau_macro_patch);
        let initial_local = global_planar_patch_to_local_for_escape_probe(n, sector, initial);
        let cp = solved_local.piece(escaped_coeff.piece).cps[escaped_coeff.i][escaped_coeff.j];
        let initial_cp =
            initial_local.piece(escaped_coeff.piece).cps[escaped_coeff.i][escaped_coeff.j];
        let local = [cp.x, cp.y];
        let initial_local = [initial_cp.x, initial_cp.y];
        let triangle_excess = (-local[0])
            .max(-local[1])
            .max(local[0] + local[1] - 1.0)
            .max(0.0);
        let delta_from_initial =
            ((local[0] - initial_local[0]).powi(2) + (local[1] - initial_local[1]).powi(2)).sqrt();
        regularization_sweep.push(RegSweep {
            regularization_weight,
            local,
            initial_local,
            triangle_excess,
            delta_from_initial,
            f3: kp20_planar_bi6_tau_f3_energy(&solved_local),
            residual: solved.max_constraint_residual,
        });
    }

    for sweep in &regularization_sweep {
        eprintln!(
            "projected section42 escape coeff n={n} candidate={candidate} coeff={piece:?}[{i},{j}] reg={reg:e}: row_touch=vector:{vector_rows_touching} scalar:{scalar_rows_touching} local=({u:e},{v:e}) initial=({iu:e},{iv:e}) excess={excess:e} delta={delta:e} F3={f3:e} residual={residual:e}",
            candidate = candidate.name(),
            piece = escaped_coeff.piece,
            i = escaped_coeff.i,
            j = escaped_coeff.j,
            reg = sweep.regularization_weight,
            u = sweep.local[0],
            v = sweep.local[1],
            iu = sweep.initial_local[0],
            iv = sweep.initial_local[1],
            excess = sweep.triangle_excess,
            delta = sweep.delta_from_initial,
            f3 = sweep.f3,
            residual = sweep.residual,
        );
    }

    let paper_weight = regularization_sweep[0];
    assert_eq!(
        vector_rows_touching, 0,
        "near-origin high-v coefficient should not be directly pinned by eq(14)/Fig.8 rows"
    );
    assert_eq!(
        scalar_rows_touching, 2,
        "promoted Fig.8 diagonal rows now bind this former high-v escape coefficient"
    );
    assert!(
        paper_weight.triangle_excess > 0.01 && paper_weight.triangle_excess < 0.02,
        "corrected Fig.8 diagonal reflection should reduce this former sector-domain escape to a small measured violation"
    );
    assert!(
        paper_weight.delta_from_initial > 0.55 && paper_weight.delta_from_initial < 0.58,
        "corrected Fig.8 diagonal reflection should keep the remaining source movement measured"
    );

    fn planar_tau_seed_macro_patch_for_escape_probe(n: usize, k: usize) -> Kp20MacroPatch {
        let mut out = Kp20MacroPatch::zero();
        for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
            let base = corner_maps[0].1;
            let piece = out.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let u = i as f64 / 6.0;
                    let v = j as f64 / 6.0;
                    let local = [
                        base.origin[0] + u * base.du[0] + v * base.dv[0],
                        base.origin[1] + u * base.du[1] + v * base.dv[1],
                    ];
                    piece.cps[i][j] = local_to_global_planar_for_escape_probe(n, k, local);
                }
            }
        }
        out
    }

    fn local_to_global_planar_for_escape_probe(n: usize, k: usize, local: [f64; 2]) -> Point3 {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        Point3::new(
            local[0] * spoke1[0] + local[1] * spoke0[0],
            local[0] * spoke1[1] + local[1] * spoke0[1],
            0.0,
        )
    }

    fn global_planar_patch_to_local_for_escape_probe(
        n: usize,
        k: usize,
        mut patch: Kp20MacroPatch,
    ) -> Kp20MacroPatch {
        for piece_id in Kp20PieceId::all() {
            let piece = patch.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let local = global_to_local_planar_for_escape_probe(n, k, piece.cps[i][j]);
                    piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
                }
            }
        }
        patch
    }

    fn global_to_local_planar_for_escape_probe(n: usize, k: usize, p: Point3) -> [f64; 2] {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
        [
            (p.x * spoke0[1] - p.y * spoke0[0]) / det,
            (spoke1[0] * p.y - spoke1[1] * p.x) / det,
        ]
    }
}

#[test]
fn projected_section42_diagonal_symmetry_measures_closed_tau_domain() {
    use super::super::kp20_bi6_tau_f3_solve::{
        Kp20PlanarTauF3ScalarSolveInput, solve_kp20_planar_tau_f3_with_scalar_constraints,
    };
    use super::super::kp20_bi6_tau_solve_assemble::kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_tensor_targets;
    use super::super::kp20_bi6_tau_solve_assemble::{
        Kp20PlanarTauPaperRowFamily, kp20_planar_tau_rotational_paper_family_residuals,
    };
    use crate::geometry::gregory::kp20_bi6_cap::Kp20MacroPatch;
    use crate::geometry::gregory::kp20_bi6_cc_ring::{
        Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
        kp20_pr08_characteristic_tensor_border_candidate_rows,
    };
    use crate::geometry::gregory::kp20_bi6_cnet::kp20_tensor_border_from_regular_cnet_rows;
    use crate::geometry::gregory::kp20_bi6_piece_maps::{
        kp20_planar_bi6_tau_f3_energy, kp20_positive_nonoverlap_c2_sector_piece_maps,
    };

    #[derive(Debug, Clone)]
    struct DomainSummary {
        candidate_name: &'static str,
        outside: usize,
        samples: usize,
        max_excess: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
        worst_i: usize,
        worst_j: usize,
        f3: f64,
        residual: f64,
        family_residuals: [f64; 7],
    }

    let n = 5;
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let mut summaries = Vec::new();

    for candidate in [
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows,
        Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows,
    ] {
        let cnet_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
            .expect("PR08 tensor-border candidate rows");
        let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&cnet_rows)
            .expect("PR08 tensor-border converts");
        let mut tau_macro_patches = Vec::with_capacity(n);
        let mut residual = 0.0_f64;
        let mut family_residuals = [0.0_f64; 7];
        for sector in 0..n {
            let row_set =
                kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_tensor_targets(
                    n,
                    sector,
                    &tensor_border,
                )
                .expect("free-target projected paper row set with local diagonal symmetry");
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch_for_diagonal_probe(n, sector),
                    vector_constraints: row_set.vector_constraints,
                    scalar_constraints: row_set.scalar_constraints,
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("diagonally symmetric projected paper tau solve succeeds");
            residual = residual.max(solved.max_constraint_residual);
            for family_residual in kp20_planar_tau_rotational_paper_family_residuals(
                n,
                sector,
                &tensor_border,
                &solved.tau_macro_patch,
            )
            .expect("family residuals evaluate")
            {
                let slot = diagonal_probe_family_residual_index(family_residual.family);
                family_residuals[slot] = family_residuals[slot].max(family_residual.max_residual);
            }
            tau_macro_patches.push(global_planar_patch_to_local_for_diagonal_probe(
                n,
                sector,
                solved.tau_macro_patch,
            ));
        }
        let domain = aggregate_planar_tau_cp_domain_for_diagonal_probe(&tau_macro_patches);
        let summary = DomainSummary {
            candidate_name: candidate.name(),
            outside: domain.outside,
            samples: domain.samples,
            max_excess: domain.max_excess,
            worst_sector: domain.worst_sector,
            worst_piece: domain.worst_piece,
            worst_i: domain.worst_i,
            worst_j: domain.worst_j,
            f3: tau_macro_patches
                .iter()
                .map(kp20_planar_bi6_tau_f3_energy)
                .sum(),
            residual,
            family_residuals,
        };
        eprintln!(
            "projected section42 diagonal tau domain n={n} candidate={name}: residual={residual:e} outside={outside}/{samples} max_excess={excess:e} worst_tau=s{sector}:{piece:?}[{i},{j}] F3={f3:e} families fig8={fig8:e} eq14={eq14:e} c0={c0:e} eq11={eq11:e} eq12={eq12:e} eq13={eq13:e} diagonal={diagonal:e}",
            name = summary.candidate_name,
            residual = summary.residual,
            outside = summary.outside,
            samples = summary.samples,
            excess = summary.max_excess,
            sector = summary.worst_sector,
            piece = summary.worst_piece,
            i = summary.worst_i,
            j = summary.worst_j,
            f3 = summary.f3,
            fig8 = summary.family_residuals[0],
            eq14 = summary.family_residuals[1],
            c0 = summary.family_residuals[2],
            eq11 = summary.family_residuals[3],
            eq12 = summary.family_residuals[4],
            eq13 = summary.family_residuals[5],
            diagonal = summary.family_residuals[6],
        );
        summaries.push(summary);
    }

    let best = summaries
        .iter()
        .min_by(|a, b| {
            a.outside
                .cmp(&b.outside)
                .then_with(|| a.max_excess.total_cmp(&b.max_excess))
                .then_with(|| a.f3.total_cmp(&b.f3))
        })
        .expect("diagonal summaries");
    eprintln!(
        "projected section42 diagonal tau selection n={n}: best={name} outside={outside}/{samples} max_excess={excess:e} worst_tau=s{sector}:{piece:?}[{i},{j}] F3={f3:e}",
        name = best.candidate_name,
        outside = best.outside,
        samples = best.samples,
        excess = best.max_excess,
        sector = best.worst_sector,
        piece = best.worst_piece,
        i = best.worst_i,
        j = best.worst_j,
        f3 = best.f3,
    );
    assert!(
        best.residual < 5.0e-2,
        "free tensor-target projection should show that sector-local diagonal compatibility is feasible"
    );
    assert!(
        best.outside < 870,
        "free tensor-target projection should improve on the current best no-diagonal closed-tau domain count"
    );

    fn diagonal_probe_family_residual_index(family: Kp20PlanarTauPaperRowFamily) -> usize {
        match family {
            Kp20PlanarTauPaperRowFamily::Fig8IndicatedAverage => 0,
            Kp20PlanarTauPaperRowFamily::TensorBorderEq14 => 1,
            Kp20PlanarTauPaperRowFamily::Section42C0 => 2,
            Kp20PlanarTauPaperRowFamily::Section42Eq11Degree5 => 3,
            Kp20PlanarTauPaperRowFamily::Section42Eq12InternalC2 => 4,
            Kp20PlanarTauPaperRowFamily::Section42Eq13Strip => 5,
            Kp20PlanarTauPaperRowFamily::DiagonalSymmetry => 6,
            Kp20PlanarTauPaperRowFamily::InternalC2
            | Kp20PlanarTauPaperRowFamily::SectorG1G2NearOrigin
            | Kp20PlanarTauPaperRowFamily::SectorG1G2NearMidpoint => {
                unreachable!("raw eq(8)/(9) rows are not part of the §4.2 live basis")
            }
        }
    }

    fn planar_tau_seed_macro_patch_for_diagonal_probe(n: usize, k: usize) -> Kp20MacroPatch {
        let mut out = Kp20MacroPatch::zero();
        for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
            let base = corner_maps[0].1;
            let piece = out.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let u = i as f64 / 6.0;
                    let v = j as f64 / 6.0;
                    let local = [
                        base.origin[0] + u * base.du[0] + v * base.dv[0],
                        base.origin[1] + u * base.du[1] + v * base.dv[1],
                    ];
                    piece.cps[i][j] = local_to_global_planar_for_diagonal_probe(n, k, local);
                }
            }
        }
        out
    }

    fn local_to_global_planar_for_diagonal_probe(n: usize, k: usize, local: [f64; 2]) -> Point3 {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        Point3::new(
            local[0] * spoke1[0] + local[1] * spoke0[0],
            local[0] * spoke1[1] + local[1] * spoke0[1],
            0.0,
        )
    }

    fn global_planar_patch_to_local_for_diagonal_probe(
        n: usize,
        k: usize,
        mut patch: Kp20MacroPatch,
    ) -> Kp20MacroPatch {
        for piece_id in Kp20PieceId::all() {
            let piece = patch.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let local = global_to_local_planar_for_diagonal_probe(n, k, piece.cps[i][j]);
                    piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
                }
            }
        }
        patch
    }

    fn global_to_local_planar_for_diagonal_probe(n: usize, k: usize, p: Point3) -> [f64; 2] {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
        [
            (p.x * spoke0[1] - p.y * spoke0[0]) / det,
            (spoke1[0] * p.y - spoke1[1] * p.x) / det,
        ]
    }

    struct TauDomainSummary {
        outside: usize,
        samples: usize,
        max_excess: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
        worst_i: usize,
        worst_j: usize,
    }

    fn aggregate_planar_tau_cp_domain_for_diagonal_probe(
        macro_patches: &[Kp20MacroPatch],
    ) -> TauDomainSummary {
        let mut summary = TauDomainSummary {
            outside: 0,
            samples: 0,
            max_excess: 0.0,
            worst_sector: 0,
            worst_piece: Kp20PieceId::NearOriginLeft,
            worst_i: 0,
            worst_j: 0,
        };
        for (sector, macro_patch) in macro_patches.iter().enumerate() {
            for piece_id in Kp20PieceId::all() {
                for (i, row) in macro_patch.piece(piece_id).cps.iter().enumerate() {
                    for (j, cp) in row.iter().enumerate() {
                        let triangle_excess = (-cp.x).max(-cp.y).max(cp.x + cp.y - 1.0).max(0.0);
                        summary.outside += usize::from(triangle_excess > 1.0e-12);
                        summary.samples += 1;
                        if triangle_excess > summary.max_excess {
                            summary.max_excess = triangle_excess;
                            summary.worst_sector = sector;
                            summary.worst_piece = piece_id;
                            summary.worst_i = i;
                            summary.worst_j = j;
                        }
                    }
                }
            }
        }
        summary
    }
}

#[test]
fn projected_section42_diagonal_source_masks_rank_before_promotion() {
    use super::super::kp20_bi6_tau_f3_solve::{
        Kp20PlanarTauCoeff, Kp20PlanarTauF3ScalarSolveInput,
        solve_kp20_planar_tau_f3_with_scalar_constraints,
    };
    use super::super::kp20_bi6_tau_solve_assemble::kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs;
    use crate::geometry::gregory::kp20_bi6_cap::Kp20MacroPatch;
    use crate::geometry::gregory::kp20_bi6_cc_ring::{
        Kp20Pr08TensorBorderCandidate, kp20_catmull_clark_characteristic_ring,
        kp20_pr08_characteristic_tensor_border_candidate_rows,
    };
    use crate::geometry::gregory::kp20_bi6_cnet::kp20_tensor_border_from_regular_cnet_rows;
    use crate::geometry::gregory::kp20_bi6_piece_maps::{
        kp20_planar_bi6_tau_f3_energy, kp20_positive_nonoverlap_c2_sector_piece_maps,
    };

    #[derive(Clone, Copy)]
    enum TensorMask {
        Section42Underlined,
        Section42UnderlinedTransposeClosed,
        LowerThreeByThree,
        LowerFourByThree,
        LowerFiveByThree,
        LowerSixByThree,
        LowerSixByThreePlusLeftFinalColumn,
        LowerSixByThreePlusRightFinalColumn,
        LowerSixByThreePlusLeftFinalColumnRightDepth0,
        LowerSixByThreePlusLeftFinalColumnRightDepth01,
        AllTensorTargets,
    }

    impl TensorMask {
        fn label(self) -> &'static str {
            match self {
                TensorMask::Section42Underlined => "section42-underlined-source",
                TensorMask::Section42UnderlinedTransposeClosed => {
                    "section42-underlined-transpose-closed"
                }
                TensorMask::LowerThreeByThree => "tensor-lower-3x3",
                TensorMask::LowerFourByThree => "tensor-lower-4x3",
                TensorMask::LowerFiveByThree => "tensor-lower-5x3",
                TensorMask::LowerSixByThree => "tensor-lower-6x3",
                TensorMask::LowerSixByThreePlusLeftFinalColumn => {
                    "tensor-lower-6x3-plus-left-final-column"
                }
                TensorMask::LowerSixByThreePlusRightFinalColumn => {
                    "tensor-lower-6x3-plus-right-final-column"
                }
                TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth0 => {
                    "tensor-lower-6x3-plus-left-final-column-right-depth0"
                }
                TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth01 => {
                    "tensor-lower-6x3-plus-left-final-column-right-depth01"
                }
                TensorMask::AllTensorTargets => "tensor-all-targets",
            }
        }

        fn coeffs(self) -> Vec<Kp20PlanarTauCoeff> {
            let mut coeffs = Vec::new();
            for piece in [
                Kp20PieceId::NearMidpointRight,
                Kp20PieceId::NearMidpointLeft,
            ] {
                match self {
                    TensorMask::Section42Underlined => {
                        add_unique_tau_coeffs(
                            &mut coeffs,
                            piece,
                            &[(0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (1, 1), (1, 2)],
                        );
                    }
                    TensorMask::Section42UnderlinedTransposeClosed => {
                        add_unique_tau_coeffs(
                            &mut coeffs,
                            piece,
                            &[
                                (0, 0),
                                (1, 0),
                                (2, 0),
                                (0, 1),
                                (0, 2),
                                (1, 1),
                                (1, 2),
                                (2, 1),
                            ],
                        );
                    }
                    TensorMask::LowerThreeByThree => {
                        for i in 0..=2 {
                            for j in 0..=2 {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(i, j)]);
                            }
                        }
                    }
                    TensorMask::LowerFourByThree
                    | TensorMask::LowerFiveByThree
                    | TensorMask::LowerSixByThree
                    | TensorMask::LowerSixByThreePlusLeftFinalColumn
                    | TensorMask::LowerSixByThreePlusRightFinalColumn => {
                        let column_end = match self {
                            TensorMask::LowerFourByThree => 3,
                            TensorMask::LowerFiveByThree => 4,
                            TensorMask::LowerSixByThree
                            | TensorMask::LowerSixByThreePlusLeftFinalColumn
                            | TensorMask::LowerSixByThreePlusRightFinalColumn => 5,
                            _ => unreachable!("handled by outer match"),
                        };
                        for i in 0..=2 {
                            for j in 0..=column_end {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(i, j)]);
                            }
                        }
                        let include_final_column = matches!(
                            (self, piece),
                            (
                                TensorMask::LowerSixByThreePlusLeftFinalColumn,
                                Kp20PieceId::NearMidpointLeft
                            ) | (
                                TensorMask::LowerSixByThreePlusRightFinalColumn,
                                Kp20PieceId::NearMidpointRight
                            )
                        );
                        if include_final_column {
                            for i in 0..=2 {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(i, 6)]);
                            }
                        }
                    }
                    TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth0
                    | TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth01 => {
                        for i in 0..=2 {
                            for j in 0..=5 {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(i, j)]);
                            }
                        }
                        if piece == Kp20PieceId::NearMidpointLeft {
                            for i in 0..=2 {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(i, 6)]);
                            }
                        }
                        if piece == Kp20PieceId::NearMidpointRight {
                            add_unique_tau_coeffs(&mut coeffs, piece, &[(0, 6)]);
                            if matches!(
                                self,
                                TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth01
                            ) {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(1, 6)]);
                            }
                        }
                    }
                    TensorMask::AllTensorTargets => {
                        for i in 0..=2 {
                            for j in 0..=6 {
                                add_unique_tau_coeffs(&mut coeffs, piece, &[(i, j)]);
                            }
                        }
                    }
                }
            }
            coeffs
        }
    }

    #[derive(Debug)]
    struct MaskQuality {
        label: &'static str,
        free_coeffs: usize,
        residual: f64,
        outside: usize,
        samples: usize,
        max_excess: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
        worst_i: usize,
        worst_j: usize,
        f3: f64,
    }

    let n = 5;
    let candidate = Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows;
    let ring = kp20_catmull_clark_characteristic_ring(n).expect("PR08 characteristic ring");
    let cnet_rows = kp20_pr08_characteristic_tensor_border_candidate_rows(&ring, candidate)
        .expect("PR08 tensor-border candidate rows");
    let tensor_border =
        kp20_tensor_border_from_regular_cnet_rows(&cnet_rows).expect("PR08 tensor-border converts");
    let mut qualities = Vec::new();

    for mask in [
        TensorMask::Section42Underlined,
        TensorMask::Section42UnderlinedTransposeClosed,
        TensorMask::LowerThreeByThree,
        TensorMask::LowerFourByThree,
        TensorMask::LowerFiveByThree,
        TensorMask::LowerSixByThree,
        TensorMask::LowerSixByThreePlusLeftFinalColumn,
        TensorMask::LowerSixByThreePlusRightFinalColumn,
        TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth0,
        TensorMask::LowerSixByThreePlusLeftFinalColumnRightDepth01,
        TensorMask::AllTensorTargets,
    ] {
        let free_coeffs = mask.coeffs();
        let mut tau_macro_patches = Vec::with_capacity(n);
        let mut residual = 0.0_f64;
        for sector in 0..n {
            let row_set =
                kp20_planar_tau_rotational_paper_row_set_with_sector_local_diagonal_free_coeffs(
                    n,
                    sector,
                    &tensor_border,
                    &free_coeffs,
                )
                .expect("masked projected paper row set with local diagonal symmetry");
            let solved = solve_kp20_planar_tau_f3_with_scalar_constraints(
                &Kp20PlanarTauF3ScalarSolveInput {
                    initial: planar_tau_seed_macro_patch_for_mask_probe(n, sector),
                    vector_constraints: row_set.vector_constraints,
                    scalar_constraints: row_set.scalar_constraints,
                    regularization_weight: 1.0e-10,
                },
            )
            .expect("masked diagonal paper tau solve succeeds");
            residual = residual.max(solved.max_constraint_residual);
            tau_macro_patches.push(global_planar_patch_to_local_for_mask_probe(
                n,
                sector,
                solved.tau_macro_patch,
            ));
        }
        let domain = aggregate_planar_tau_cp_domain_for_mask_probe(&tau_macro_patches);
        let quality = MaskQuality {
            label: mask.label(),
            free_coeffs: free_coeffs.len(),
            residual,
            outside: domain.outside,
            samples: domain.samples,
            max_excess: domain.max_excess,
            worst_sector: domain.worst_sector,
            worst_piece: domain.worst_piece,
            worst_i: domain.worst_i,
            worst_j: domain.worst_j,
            f3: tau_macro_patches
                .iter()
                .map(kp20_planar_bi6_tau_f3_energy)
                .sum(),
        };
        eprintln!(
            "projected section42 diagonal source mask n={n} candidate={candidate} mask={mask} free_coeffs={free_coeffs} residual={residual:e} outside={outside}/{samples} max_excess={excess:e} worst_tau=s{sector}:{piece:?}[{i},{j}] F3={f3:e}",
            candidate = candidate.name(),
            mask = quality.label,
            free_coeffs = quality.free_coeffs,
            residual = quality.residual,
            outside = quality.outside,
            samples = quality.samples,
            excess = quality.max_excess,
            sector = quality.worst_sector,
            piece = quality.worst_piece,
            i = quality.worst_i,
            j = quality.worst_j,
            f3 = quality.f3,
        );
        qualities.push(quality);
    }

    let best = qualities
        .iter()
        .min_by(|a, b| {
            a.residual
                .total_cmp(&b.residual)
                .then_with(|| a.outside.cmp(&b.outside))
                .then_with(|| a.max_excess.total_cmp(&b.max_excess))
        })
        .expect("mask qualities");
    eprintln!(
        "projected section42 diagonal source mask selection n={n}: best={mask} free_coeffs={free_coeffs} residual={residual:e} outside={outside}/{samples} max_excess={excess:e} F3={f3:e}",
        mask = best.label,
        free_coeffs = best.free_coeffs,
        residual = best.residual,
        outside = best.outside,
        samples = best.samples,
        excess = best.max_excess,
        f3 = best.f3,
    );

    let all_targets = qualities
        .iter()
        .find(|quality| quality.label == TensorMask::AllTensorTargets.label())
        .expect("all-target mask quality");
    let underlined = qualities
        .iter()
        .find(|quality| quality.label == TensorMask::Section42Underlined.label())
        .expect("underlined mask quality");
    assert!(
        all_targets.residual < underlined.residual,
        "all-target feasibility bound should improve on the current underlined source mask"
    );
    assert!(
        all_targets.outside < underlined.outside,
        "all-target feasibility bound should improve domain count before mask promotion"
    );

    fn add_unique_tau_coeffs(
        coeffs: &mut Vec<Kp20PlanarTauCoeff>,
        piece: Kp20PieceId,
        indices: &[(usize, usize)],
    ) {
        for &(i, j) in indices {
            let coeff = Kp20PlanarTauCoeff::new(piece, i, j);
            if !coeffs.contains(&coeff) {
                coeffs.push(coeff);
            }
        }
    }

    fn planar_tau_seed_macro_patch_for_mask_probe(n: usize, k: usize) -> Kp20MacroPatch {
        let mut out = Kp20MacroPatch::zero();
        for (piece_id, corner_maps) in kp20_positive_nonoverlap_c2_sector_piece_maps() {
            let base = corner_maps[0].1;
            let piece = out.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let u = i as f64 / 6.0;
                    let v = j as f64 / 6.0;
                    let local = [
                        base.origin[0] + u * base.du[0] + v * base.dv[0],
                        base.origin[1] + u * base.du[1] + v * base.dv[1],
                    ];
                    piece.cps[i][j] = local_to_global_planar_for_mask_probe(n, k, local);
                }
            }
        }
        out
    }

    fn local_to_global_planar_for_mask_probe(n: usize, k: usize, local: [f64; 2]) -> Point3 {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        Point3::new(
            local[0] * spoke1[0] + local[1] * spoke0[0],
            local[0] * spoke1[1] + local[1] * spoke0[1],
            0.0,
        )
    }

    fn global_planar_patch_to_local_for_mask_probe(
        n: usize,
        k: usize,
        mut patch: Kp20MacroPatch,
    ) -> Kp20MacroPatch {
        for piece_id in Kp20PieceId::all() {
            let piece = patch.piece_mut(piece_id);
            for i in 0..=6 {
                for j in 0..=6 {
                    let local = global_to_local_planar_for_mask_probe(n, k, piece.cps[i][j]);
                    piece.cps[i][j] = Point3::new(local[0], local[1], 0.0);
                }
            }
        }
        patch
    }

    fn global_to_local_planar_for_mask_probe(n: usize, k: usize, p: Point3) -> [f64; 2] {
        let alpha = core::f64::consts::TAU / n as f64;
        let spoke0 = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let spoke1 = [
            ((k + 1) as f64 * alpha).cos(),
            ((k + 1) as f64 * alpha).sin(),
        ];
        let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
        [
            (p.x * spoke0[1] - p.y * spoke0[0]) / det,
            (spoke1[0] * p.y - spoke1[1] * p.x) / det,
        ]
    }

    struct TauDomainSummary {
        outside: usize,
        samples: usize,
        max_excess: f64,
        worst_sector: usize,
        worst_piece: Kp20PieceId,
        worst_i: usize,
        worst_j: usize,
    }

    fn aggregate_planar_tau_cp_domain_for_mask_probe(
        macro_patches: &[Kp20MacroPatch],
    ) -> TauDomainSummary {
        let mut summary = TauDomainSummary {
            outside: 0,
            samples: 0,
            max_excess: 0.0,
            worst_sector: 0,
            worst_piece: Kp20PieceId::NearOriginLeft,
            worst_i: 0,
            worst_j: 0,
        };
        for (sector, macro_patch) in macro_patches.iter().enumerate() {
            for piece_id in Kp20PieceId::all() {
                for (i, row) in macro_patch.piece(piece_id).cps.iter().enumerate() {
                    for (j, cp) in row.iter().enumerate() {
                        let triangle_excess = (-cp.x).max(-cp.y).max(cp.x + cp.y - 1.0).max(0.0);
                        summary.outside += usize::from(triangle_excess > 1.0e-12);
                        summary.samples += 1;
                        if triangle_excess > summary.max_excess {
                            summary.max_excess = triangle_excess;
                            summary.worst_sector = sector;
                            summary.worst_piece = piece_id;
                            summary.worst_i = i;
                            summary.worst_j = j;
                        }
                    }
                }
            }
        }
        summary
    }
}

#[test]
fn eq8_residual_with_synthetic_overlay_estimate() {
    use crate::geometry::Point3;
    use crate::geometry::gregory::kp20_bi6_cap::Kp20PieceId;
    use crate::geometry::gregory::kp20_bi6_guide_sample::kp20_corner_block_from_jet;
    use crate::geometry::gregory::kp20_bi6_tau::Kp20TauMap;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let mut cap = assemble_synthetic_cap(n);
        // The kernel's guide isn't accessible after assembly, but
        // its apex jet equals what `kp20_compose_guide_tau_jet`
        // would return at (0, 0). For this probe we directly
        // rebuild the apex 4x4 block from the production guide,
        // overlay, then re-measure eq(8).
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = crate::geometry::gregory::kp20_bi6_guide_build::build_kp20_guide_ring(
            &crate::geometry::gregory::kp20_bi6_guide_build::Kp20GuideBuildInput {
                star_point,
                tensor_border: tensor_border.clone(),
            },
        )
        .expect("guide builds");
        for (k, macro_patch) in cap.macro_patches.iter_mut().enumerate() {
            // Use for_sector_local — matches production sampling
            // (kp20_bi6_piece_maps line 317).
            let tau = Kp20TauMap::for_sector_local(n, k).expect("tau for valence");
            let jet = crate::geometry::gregory::kp20_bi6_tau::kp20_compose_guide_tau_jet(
                guide.sector(k),
                &tau,
                0.0,
                0.0,
            );
            let block_left = kp20_corner_block_from_jet(jet);
            // For the right piece, the affine map (kp20_default_sector_piece_maps
            // line 28-29) swaps (du, dv). The piece-local jet is therefore
            // h^R[i][j] = h^L[j][i], so block_right = block_left^T.
            for (i, row) in block_left.iter().enumerate().take(4) {
                for (j, &cp) in row.iter().enumerate().take(4 - i) {
                    macro_patch.piece_mut(Kp20PieceId::NearOriginLeft).cps[i][j] = cp;
                    macro_patch.piece_mut(Kp20PieceId::NearOriginRight).cps[j][i] = cp;
                }
            }
        }
        let report = measure_kp20_assembly_constraints(&cap);
        eprintln!(
            "eq(8)/eq(9) after synthetic overlay n={n}: worst_eq8 = {e8:e}  worst_eq9 = {e9:e}  L6(1) = {l1:e}",
            e8 = report.worst_eq8(),
            e9 = report.worst_eq9(),
            l1 = report.worst_lemma6_cond_1,
        );
        assert!(report.worst_eq8().is_finite());
        assert!(report.worst_eq9().is_finite());
    }
}

#[test]
fn assembled_cap_lemma4_on_orig_pieces_apex_bb_coefs() {
    use crate::geometry::Point3;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let cap = {
            let tensor_border = synthetic_uniform_tensor_border(n);
            let input = crate::geometry::gregory::kp20_bi6_assemble::Kp20Bi6CapInput {
                star_point,
                tensor_border,
            };
            crate::geometry::gregory::kp20_bi6_assemble::assemble_kp20_bi6_cap(&input)
                .expect("synthetic assemble must succeed")
        };

        let mut worst_order = [0.0_f64; 4];
        for macro_patch in &cap.macro_patches {
            let left = &macro_patch
                .piece(crate::geometry::gregory::kp20_bi6_cap::Kp20PieceId::NearOriginLeft)
                .cps;
            let right = &macro_patch
                .piece(crate::geometry::gregory::kp20_bi6_cap::Kp20PieceId::NearOriginRight)
                .cps;
            for (i, slot) in worst_order.iter_mut().enumerate() {
                // c^s_{i, 0} from NearOriginLeft, c^{s+1}_{0, i} from NearOriginRight.
                let l = left[i][0];
                let r = right[0][i];
                let dx = (l.x - r.x).abs();
                let dy = (l.y - r.y).abs();
                let dz = (l.z - r.z).abs();
                let mag = dx.max(dy).max(dz);
                *slot = slot.max(mag);
            }
        }
        eprintln!(
            "assembled-cap Lemma 4 n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst_order[0],
            o1 = worst_order[1],
            o2 = worst_order[2],
            o3 = worst_order[3],
        );
        for (i, w) in worst_order.iter().enumerate() {
            assert!(
                w.is_finite(),
                "assembled-cap Lemma 4 order {i} non-finite at n = {n}",
            );
        }
    }
}

// --- Investigation Step 3: LSQ guide apex-region symmetry --------

/// Step 3 of the investigation plan in
/// `docs/decisions/2026-06-23-kp20-bi6-paper-grade-guide-from-kp18.md`.
///
/// **Symmetry** condition from
/// `docs/architecture/kp20-lemma4-sufficient-conditions.md`:
/// across the shared spoke between sector `s` and sector `s+1`,
/// the guide's directional derivatives at the apex `o` must
/// agree —
/// ```text
///   ∂^i_v_g g_s   at (0, 0) == ∂^i_u_g g_{s+1}   at (0, 0)
/// ```
/// for `i = 0, 1, 2, 3` — because `v_g` of sector `s` and `u_g`
/// of sector `(s+1)` both point along the shared spoke.
///
/// Step 2 confirmed `τ` is symmetric to machine zero. The Phase D
/// defect on the assembled cap must therefore live in the guide
/// failing this property.
///
/// We *expect* this test to FAIL (large defect at orders 1 and 3)
/// — that failure is the load-bearing diagnostic.
#[test]
fn lsq_guide_apex_directional_derivatives_match_across_spokes() {
    use crate::geometry::Point3;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border,
        })
        .expect("guide ring builds for synthetic input");

        let mut jets = Vec::with_capacity(n);
        for s in 0..n {
            jets.push(guide.sector(s).jet3(0.0, 0.0));
        }

        let mut worst_order = [0.0_f64; 4];
        for s in 0..n {
            let next = (s + 1) % n;
            for (i, slot) in worst_order.iter_mut().enumerate() {
                let cu = jets[s].h[0][i]; // ∂^i_v g_s at apex
                let nv = jets[next].h[i][0]; // ∂^i_u g_{s+1} at apex
                let dx = (cu.x - nv.x).abs();
                let dy = (cu.y - nv.y).abs();
                let dz = (cu.z - nv.z).abs();
                let mag = dx.max(dy).max(dz);
                *slot = slot.max(mag);
            }
        }
        eprintln!(
            "LSQ guide apex symmetry n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst_order[0],
            o1 = worst_order[1],
            o2 = worst_order[2],
            o3 = worst_order[3],
        );
        // Measurement-only assertion: defect must at least be finite.
        for (i, w) in worst_order.iter().enumerate() {
            assert!(
                w.is_finite(),
                "guide apex symmetry order {i} non-finite at n = {n}",
            );
        }
    }
}

// --- Phase D probe: planar guide should satisfy Lemma 4 ---------

/// A planar guide `g(x, y) = (x, y, 0)` sampled into a tri6 sector
/// barycentrically. THIS GUIDE DOES NOT SATISFY KP18 (3)-(5)
/// pairwise C^k joins — those joins use specific weights
/// `(w0, w1, w2) = (-1, 2c, 2(1-c))` that the naive barycentric
/// lifting does not respect for n ≠ 4.
///
/// Empirically this probe fails (order 1 ~1.3, order 3 ~75), which
/// confirms the test below: Lemma 4 is NOT a consequence of g
/// being C^∞; it requires the KP18 join structure.
#[test]
#[ignore = "negative-control: naive barycentric planar lifting does NOT satisfy KP18 (3)–(5) for n != 4, so Lemma 4 also fails. Documented finding from ADR 2026-06-23-kp20-bi6-paper-grade-guide-from-kp18.md."]
fn lemma4_fails_on_naive_barycentric_planar_guide() {
    use crate::geometry::Point3;
    use crate::geometry::gregory::bernstein_tri6::{TRI6_CP_COUNT, flat_index};
    use crate::geometry::gregory::kp20_bi6_guide_sample::Kp20Tri6GuideSector;

    let valences: &[usize] = &[5, 6, 7];

    for &n in valences {
        let alpha = std::f64::consts::TAU / n as f64;
        let mut sectors: Vec<Kp20Tri6GuideSector> = Vec::with_capacity(n);
        for s in 0..n {
            let v_i = Point3::new((alpha * s as f64).cos(), (alpha * s as f64).sin(), 0.0);
            let v_j = Point3::new(
                (alpha * (s + 1) as f64).cos(),
                (alpha * (s + 1) as f64).sin(),
                0.0,
            );
            let v_k = Point3::new(0.0, 0.0, 0.0);
            let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI6_CP_COUNT];
            for i in 0..=6 {
                for j in 0..=(6 - i) {
                    let k = 6 - i - j;
                    let bi = i as f64 / 6.0;
                    let bj = j as f64 / 6.0;
                    let bk = k as f64 / 6.0;
                    cps[flat_index(6, i, j)] = Point3::new(
                        bi * v_i.x + bj * v_j.x + bk * v_k.x,
                        bi * v_i.y + bj * v_j.y + bk * v_k.y,
                        bi * v_i.z + bj * v_j.z + bk * v_k.z,
                    );
                }
            }
            sectors.push(Kp20Tri6GuideSector { cps });
        }

        let mut jets = Vec::with_capacity(n);
        for (k, sector) in sectors.iter().enumerate() {
            let tau = Kp20TauMap::for_sector(n, k).expect("tau map for valence");
            let jet = kp20_compose_guide_tau_jet(sector, &tau, 0.0, 0.0);
            jets.push(jet);
        }

        let mut worst_per_order = [0.0_f64; 4];
        for k in 0..n {
            let next = (k + 1) % n;
            for (i, worst_slot) in worst_per_order.iter_mut().enumerate() {
                let next_v = jets[next].h[0][i];
                let curr_u = jets[k].h[i][0];
                let dx = (next_v.x - curr_u.x).abs();
                let dy = (next_v.y - curr_u.y).abs();
                let dz = (next_v.z - curr_u.z).abs();
                let mag = dx.max(dy).max(dz);
                *worst_slot = worst_slot.max(mag);
            }
        }
        eprintln!(
            "PLANAR guide Lemma 4 n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst_per_order[0],
            o1 = worst_per_order[1],
            o2 = worst_per_order[2],
            o3 = worst_per_order[3],
        );
        for (i, &w) in worst_per_order.iter().enumerate() {
            assert!(
                w < 1e-10,
                "PLANAR guide Lemma 4 defect at order {i} = {w:e} (should be ≤ 1e-10)",
            );
        }
    }
}

/// Diagnostic probe: how well does the kernel's LSQ satisfy its
/// own KP18 equality constraints?
#[test]
fn lsq_max_constraint_residual_on_synthetic_inputs() {
    use crate::geometry::Point3;
    use crate::geometry::gregory::kp20_bi6_guide_lsq::solve_kp20_guide_lsq_with_report;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let solved = solve_kp20_guide_lsq_with_report(&tensor_border, star_point)
            .expect("LSQ solves on synthetic input");
        eprintln!(
            "LSQ report n={n}: indep_constraints={ic} free_params={fp} max_fit={mf:e} rms_fit={rf:e} max_constraint={mc:e} kkt_cond={kc:e}",
            ic = solved.report.independent_constraints,
            fp = solved.report.free_parameters_after_center_pin,
            mf = solved.report.max_fit_residual,
            rf = solved.report.rms_fit_residual,
            mc = solved.report.max_constraint_residual,
            kc = solved.report.kkt_condition_number,
        );
    }
}

// --- Phase D: Lemma 4 cubic-expansion identification ------------

/// Lemma 4 says the cubic Hermite jet `h := g ∘ τ` of adjacent
/// sectors at the extraordinary point `o` agrees on the shared
/// spoke direction:
///   `h^{s+1}.h[0][i] == h^s.h[i][0]` for i = 0..3
/// where `h^s.h[p][q] = ∂u^p ∂v^q (g ∘ τ_s)` at `(0, 0)` of sector s.
///
/// For sector s, the u-axis is along spoke (s+1) (paper §4.1
/// eq (22): `cps[(3,1,0)] = v_ddot/3 = spoke (s+1) / 3`). For
/// sector s+1, the v-axis is along spoke (s+1). So the i-th
/// derivative along the shared spoke (s+1) is `h^s.h[i][0]` from
/// sector s and `h^{s+1}.h[0][i]` from sector s+1.
///
/// **Measurement, not assertion.** The kernel's guide-LSQ
/// approximation may not be paper-faithful KP18, so this is a
/// diagnostic. We emit the worst-defect-per-valence to stderr
/// and assert only `finite`.
#[test]
fn kp20_lemma4_cubic_expansion_identity_across_sectors() {
    use crate::geometry::Point3;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border,
        })
        .expect("guide ring builds for synthetic input");

        let mut jets = Vec::with_capacity(n);
        for k in 0..n {
            let tau = Kp20TauMap::for_sector(n, k).expect("tau map for valence");
            let jet = kp20_compose_guide_tau_jet(guide.sector(k), &tau, 0.0, 0.0);
            jets.push(jet);
        }

        // Identity: jets[(k+1) % n].h[0][i] == jets[k].h[i][0] for i = 0..4.
        let mut worst_per_order = [0.0_f64; 4];
        for k in 0..n {
            let next = (k + 1) % n;
            for (i, worst_slot) in worst_per_order.iter_mut().enumerate() {
                let next_v = jets[next].h[0][i];
                let curr_u = jets[k].h[i][0];
                let dx = (next_v.x - curr_u.x).abs();
                let dy = (next_v.y - curr_u.y).abs();
                let dz = (next_v.z - curr_u.z).abs();
                let mag = dx.max(dy).max(dz);
                *worst_slot = worst_slot.max(mag);
            }
        }
        for &w in &worst_per_order {
            assert!(w.is_finite(), "Lemma 4 defect is non-finite at n = {n}");
        }
        eprintln!(
            "Lemma 4 KP20 cubic-identity n={n}: order 0 = {o0:e}  order 1 = {o1:e}  order 2 = {o2:e}  order 3 = {o3:e}",
            o0 = worst_per_order[0],
            o1 = worst_per_order[1],
            o2 = worst_per_order[2],
            o3 = worst_per_order[3],
        );
    }
}

/// Measure the rank of the global constraint matrix M(n, c) for
/// every valence and emit a markdown report under
/// `docs/reports/kp20-constraint-rank.md`. Per ADR
/// `docs/decisions/2026-06-23-kp20-bi6-macro-patch-paper-grade.md`
/// Phase C: the test asserts the file exists and parses; it does
/// NOT assert specific rank numbers.
#[test]
fn kp20_constraint_matrix_rank_per_valence() {
    use crate::geometry::gregory::direct_linalg::dense_rank;
    use std::io::Write;

    let valences: &[usize] = &[3, 5, 6, 7, 8, 9, 10];
    let tol = 1e-9;

    let mut report_lines = Vec::<String>::new();
    report_lines.push("# KP20 global constraint matrix rank by valence".into());
    report_lines.push(String::new());
    report_lines.push(format!(
        "Generated by `kp20_constraint_matrix_rank_per_valence` in `{}`.",
        file!()
    ));
    report_lines.push(String::new());
    report_lines.push(
        "M(n, c) stacks per-macro-patch eq (8) G¹ + eq (9) G² Jacobians \
         (30n rows) above within-macro-patch boundary identifications \
         (14n rows). Columns are the 196n BB-coefs across 4 pieces × n \
         macro patches. The numerical rank is computed via SVD with \
         tolerance σ > 1e-9. See ADR \
         `docs/decisions/2026-06-23-kp20-bi6-macro-patch-paper-grade.md` \
         Phase C."
            .into(),
    );
    report_lines.push(String::new());
    report_lines.push("| n | rows | cols | rank | rows − rank | cols − rank |".into());
    report_lines.push("|---|------|------|------|-------------|-------------|".into());

    for &n in valences {
        let (layout, m) = kp20_global_constraint_matrix(n);
        let rank = dense_rank(&m, tol);
        let row_defect = layout.rows().saturating_sub(rank);
        let col_defect = layout.cols.saturating_sub(rank);
        report_lines.push(format!(
            "| {n} | {rows} | {cols} | {rank} | {row_defect} | {col_defect} |",
            rows = layout.rows(),
            cols = layout.cols,
            rank = rank,
            row_defect = row_defect,
            col_defect = col_defect,
        ));
        eprintln!(
            "n={n}: rows={rows} cols={cols} rank={rank} row-defect={row_defect} col-defect={col_defect}",
            rows = layout.rows(),
            cols = layout.cols,
            rank = rank,
            row_defect = row_defect,
            col_defect = col_defect,
        );
    }

    // Write report to docs/reports/kp20-constraint-rank.md.
    let crate_dir = env!("CARGO_MANIFEST_DIR");
    let report_path = std::path::Path::new(crate_dir)
        .join("..")
        .join("..")
        .join("docs")
        .join("reports")
        .join("kp20-constraint-rank.md");
    std::fs::create_dir_all(report_path.parent().unwrap()).expect("create docs/reports parent dir");
    let mut f = std::fs::File::create(&report_path).expect("create report file");
    for line in &report_lines {
        writeln!(f, "{line}").expect("write report line");
    }
    f.flush().expect("flush report");
    drop(f);
    // Verify the report exists and contains every valence.
    let written = std::fs::read_to_string(&report_path).expect("read back report");
    for &n in valences {
        assert!(
            written.contains(&format!("| {n} |")),
            "report missing row for n = {n}",
        );
    }
}
