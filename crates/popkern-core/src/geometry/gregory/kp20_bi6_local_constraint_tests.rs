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
use crate::geometry::gregory::kp20_bi6_assemble::{Kp20Bi6CapInput, assemble_kp20_bi6_cap};
use crate::geometry::gregory::kp20_bi6_guide_build::{Kp20GuideBuildInput, build_kp20_guide_ring};
use crate::geometry::gregory::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
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
/// Probe: replace the kernel's affine-only step (iii) sampler
/// with the τ-aware variant `sample_kp20_macro_patch_from_tau_maps`
/// and re-run the rest of the assembly pipeline. If Phase B
/// residuals drop, the load-bearing fix is "switch the sampler",
/// not "overlay step (i)".
#[test]
fn eq8_residual_with_tau_based_sampler_estimate() {
    use crate::geometry::Point3;
    use crate::geometry::gregory::kp20_bi6_assemble::{
        Kp20Bi6AssembleError, Kp20Bi6MacroCapInput, assemble_kp20_bi6_macro_cap,
    };
    use crate::geometry::gregory::kp20_bi6_guide_build::{
        Kp20GuideBuildInput, build_kp20_guide_ring,
    };
    use crate::geometry::gregory::kp20_bi6_piece_maps::sample_kp20_macro_patch_from_tau_maps;

    let valences: &[usize] = &[5, 6, 7];
    let star_point = Point3::new(0.0, 0.0, 0.5);

    for &n in valences {
        let tensor_border = synthetic_uniform_tensor_border(n);
        let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
            star_point,
            tensor_border: tensor_border.clone(),
        })
        .expect("guide builds");
        let sampled_macro_patches: Result<Vec<_>, String> = (0..n)
            .map(|k| sample_kp20_macro_patch_from_tau_maps(guide.sector(k), n, k))
            .collect();
        let sampled_macro_patches = sampled_macro_patches.expect("tau sampling");
        let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
            tensor_border,
            sampled_macro_patches,
        })
        .unwrap_or_else(|e: Kp20Bi6AssembleError| panic!("assemble failed: {e}"));
        let report = measure_kp20_assembly_constraints(&cap);
        eprintln!(
            "τ-based sampler n={n}: worst_eq8 = {e8:e}  worst_eq9 = {e9:e}  L6(1) = {l1:e}",
            e8 = report.worst_eq8(),
            e9 = report.worst_eq9(),
            l1 = report.worst_lemma6_cond_1,
        );
        assert!(report.worst_eq8().is_finite());
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
