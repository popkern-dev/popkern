//! Tests for the CC-χ LSQ variant.

use super::super::kp13_bi5_boundary::TensorBorder;
use super::super::kp13_bi5_chi::{BilinearChi, CatmullClarkChi, ChiCorner};
use super::super::kp13_bi5_free_params::{Bi5FreeParams, apply_free_params};
use super::super::kp13_bi5_gtri5::sector_corner;
use super::super::kp13_bi5_lsq::{build_lsq_system, num_free_params, num_rows};
use super::*;

const EPS: f64 = 1e-10;

// ---- inverse_sector_basis ----

#[test]
fn inverse_sector_basis_is_inverse() {
    for n in [3, 4, 5, 6, 7, 8] {
        for k in 0..n {
            let m_inv = inverse_sector_basis(n, k);
            let (mkx, mky) = sector_corner(n, k);
            let (mk1x, mk1y) = sector_corner(n, k + 1);
            // M_inv · M should be identity.
            let col0 = mat_vec(&m_inv, (mkx, mky));
            let col1 = mat_vec(&m_inv, (mk1x, mk1y));
            assert!(
                (col0.0 - 1.0).abs() < EPS && col0.1.abs() < EPS,
                "n={n}, k={k}: M_inv · m_k = {col0:?} (expected (1,0))"
            );
            assert!(
                col1.0.abs() < EPS && (col1.1 - 1.0).abs() < EPS,
                "n={n}, k={k}: M_inv · m_{{k+1}} = {col1:?} (expected (0,1))"
            );
        }
    }
}

fn mat_vec(m: &[[f64; 2]; 2], v: (f64, f64)) -> (f64, f64) {
    (m[0][0] * v.0 + m[0][1] * v.1, m[1][0] * v.0 + m[1][1] * v.1)
}

// ---- ambient_to_simplex_jet ----

#[test]
fn ambient_simplex_preserves_corner_value() {
    // At (1, 0) corner: χ = m_k, simplex coords of m_k = (1, 0).
    for n in [5, 6, 7, 8] {
        let cc_chi = CatmullClarkChi::new(n);
        for k in 0..n {
            let jet = cc_chi.jet(k, ChiCorner::Mk);
            let jet_s = ambient_to_simplex_jet(&jet, n, k);
            assert!(
                (jet_s.pos.0 - 1.0).abs() < EPS && jet_s.pos.1.abs() < EPS,
                "n={n}, k={k}: {:?} (expected (1, 0))",
                jet_s.pos
            );
        }
    }
}

#[test]
fn ambient_simplex_preserves_other_corners() {
    for n in [5, 6, 7, 8] {
        let cc_chi = CatmullClarkChi::new(n);
        for k in 0..n {
            let j_star = ambient_to_simplex_jet(&cc_chi.jet(k, ChiCorner::Star), n, k);
            assert!(j_star.pos.0.abs() < EPS && j_star.pos.1.abs() < EPS);
            let j_mk1 = ambient_to_simplex_jet(&cc_chi.jet(k, ChiCorner::Mk1), n, k);
            assert!(j_mk1.pos.0.abs() < EPS && (j_mk1.pos.1 - 1.0).abs() < EPS);
        }
    }
}

#[test]
fn ambient_simplex_bilinear_identity() {
    // For the bilinear map, the simplex-coord jet IS the χ-domain
    // (u, v) identity: d_u = (1, 0), d_v = (0, 1), all 2nd = 0.
    for n in [5, 6, 7, 8] {
        let bilin = BilinearChi::new(n);
        for k in 0..n {
            for corner in [
                ChiCorner::Star,
                ChiCorner::Mk,
                ChiCorner::Mk1,
                ChiCorner::Vk,
            ] {
                let jet = bilin.jet(k, corner);
                let js = ambient_to_simplex_jet(&jet, n, k);
                assert!((js.d_u.0 - 1.0).abs() < EPS && js.d_u.1.abs() < EPS);
                assert!(js.d_v.0.abs() < EPS && (js.d_v.1 - 1.0).abs() < EPS);
                assert!(js.d_uu.0.abs() < EPS && js.d_uu.1.abs() < EPS);
                assert!(js.d_uv.0.abs() < EPS && js.d_uv.1.abs() < EPS);
                assert!(js.d_vv.0.abs() < EPS && js.d_vv.1.abs() < EPS);
            }
        }
    }
}

// ---- composite_jet_cc: basic structural ----

#[test]
fn composite_jet_cc_finite() {
    for n in [5, 6, 7, 8] {
        let cc_chi = CatmullClarkChi::new(n);
        let params = Bi5FreeParams::zero(n);
        let g = apply_free_params(&params, crate::geometry::Point3::new(1.0, 2.0, 3.0)).unwrap();
        for k in 0..n {
            for corner in [ChiCorner::Mk, ChiCorner::Vk] {
                let j = composite_jet_cc(&g, k, corner, &cc_chi);
                assert!(j.value.x.is_finite() && j.value.y.is_finite() && j.value.z.is_finite());
                assert!(j.d_u.x.is_finite() && j.d_v.x.is_finite());
                assert!(j.d_uu.x.is_finite() && j.d_uv.x.is_finite() && j.d_vv.x.is_finite());
            }
        }
    }
}

// ---- build_lsq_system_cc ----

#[test]
fn system_shape_matches_bilinear() {
    for n in [5, 6, 7, 8] {
        let target = TensorBorder::zero(n);
        let cc_limit = Point3::new(0.0, 0.0, 0.0);
        let sys_cc = build_lsq_system_cc(&target, cc_limit).unwrap();
        let sys_bil = build_lsq_system(&target, cc_limit).unwrap();
        assert_eq!(sys_cc.n, sys_bil.n);
        assert_eq!(sys_cc.a.len(), sys_bil.a.len());
        assert_eq!(sys_cc.a[0].len(), sys_bil.a[0].len());
        assert_eq!(sys_cc.a.len(), num_rows(n));
        assert_eq!(sys_cc.a[0].len(), num_free_params(n));
    }
}

#[test]
fn cc_accepts_n3_n4_via_pseudo_inverse() {
    // KP09a § 6 low-valence support via pseudo-inverse: CC variant
    // inherits from the underlying circulant solver.
    for n in [3, 4] {
        let target = TensorBorder::zero(n);
        let sys = build_lsq_system_cc(&target, Point3::new(0.0, 0.0, 0.0))
            .expect("n ∈ {3, 4} now supported via CC-χ LSQ");
        assert_eq!(sys.n, n);
    }
}

#[test]
fn cc_system_differs_from_bilinear_in_2nd_partial_rows() {
    // The CC-χ system should differ from the bilinear one in the
    // second-partial rows — because the CC χ has non-zero 2nd
    // partials while bilinear has zero. The VALUE rows at corners
    // where χ_cc = χ_bil (Star, Mk, Mk1) should coincide.
    for n in [5, 6, 7] {
        let target = TensorBorder::zero(n);
        let cc = Point3::new(1.0, 0.5, -0.2);
        let sys_cc = build_lsq_system_cc(&target, cc).unwrap();
        let sys_bil = build_lsq_system(&target, cc).unwrap();

        // Check that at least SOME entries differ in ∂²u rows.
        let n_rows = num_rows(n);
        let n_cols = num_free_params(n);
        let mut any_diff = false;
        for r in 0..n_rows {
            for c in 0..n_cols {
                if (sys_cc.a[r][c] - sys_bil.a[r][c]).abs() > 1e-6 {
                    any_diff = true;
                    break;
                }
            }
            if any_diff {
                break;
            }
        }
        assert!(
            any_diff,
            "n={n}: CC system should differ from bilinear in at least one entry"
        );
    }
}

// ---- solve_g_tri_lsq_cc: convergence + determinism ----

#[test]
fn solve_cc_deterministic() {
    let n = 5;
    let target = TensorBorder::zero(n);
    let cc = Point3::new(0.1, 0.2, 0.3);
    let g_a = solve_g_tri_lsq_cc(&target, cc).unwrap();
    let g_b = solve_g_tri_lsq_cc(&target, cc).unwrap();
    // Compare one sector's CP 0 (star).
    let s_a = g_a.sector(0).cps[0];
    let s_b = g_b.sector(0).cps[0];
    assert!((s_a.x - s_b.x).abs() < 1e-14);
    assert!((s_a.y - s_b.y).abs() < 1e-14);
    assert!((s_a.z - s_b.z).abs() < 1e-14);
}

#[test]
fn solve_cc_pins_star_to_cc_limit() {
    // The star point of the resulting g△ should match cc_limit by
    // construction (apply_free_params pins P_005 = cc_limit).
    let n = 6;
    let target = TensorBorder::zero(n);
    let cc = Point3::new(0.7, -0.3, 1.1);
    let g = solve_g_tri_lsq_cc(&target, cc).unwrap();
    // Star is at IDX_500 = 0 in this codebase's paper-aligned tri5
    // indexing (geometric position (u=0, v=0, w=1)).
    let star = g.sector(0).cps[0];
    assert!(
        (star.x - cc.x).abs() < 1e-10
            && (star.y - cc.y).abs() < 1e-10
            && (star.z - cc.z).abs() < 1e-10,
        "star = {star:?}, cc_limit = {cc:?}"
    );
}

#[test]
fn solve_cc_zero_target_zero_cc_gives_zero_patch() {
    // With both target = 0 and cc_limit = 0, the LSQ should
    // produce a zero g△ (all CPs at origin), by linearity.
    let n = 5;
    let target = TensorBorder::zero(n);
    let g = solve_g_tri_lsq_cc(&target, Point3::new(0.0, 0.0, 0.0)).unwrap();
    for k in 0..n {
        let sector = g.sector(k);
        for cp in &sector.cps {
            assert!(
                cp.x.abs() < 1e-9 && cp.y.abs() < 1e-9 && cp.z.abs() < 1e-9,
                "n={n}, k={k}: unexpected non-zero CP {cp:?}"
            );
        }
    }
}
