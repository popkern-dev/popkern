use super::tests_support::*;
use super::*;

// ---- CC-χ eq (8) variant ----

#[test]
fn eq8_cc_matches_target_value_and_partials() {
    // With full CC χ, the composite (g ∘ L⁻¹ ∘ χ) jet at the v_k
    // corner must match b^k's jet there. We verify by re-evaluating
    // the composite through the same pipeline and checking agreement.
    use super::super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::super::kp13_bi5_chi::{CatmullClarkChi, ChiCorner};
    use super::super::super::kp13_bi5_jet::compose_chain_rule;
    use super::super::super::kp13_bi5_lsq_cc::ambient_to_simplex_jet;

    for n in [5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.1, 0.2, 0.3)).unwrap();
        let mut target = TensorBorder::zero(n);
        for k in 0..n {
            let t = k as f64;
            target.patches[k].cps[3][3] = Point3::new(t + 1.0, t * 0.5, -t);
            target.patches[k].cps[2][3] = Point3::new(0.1, 0.2, t * 0.1);
            target.patches[k].cps[3][2] = Point3::new(-0.1, 0.3, 0.0);
            target.patches[k].cps[2][2] = Point3::new(0.0, 0.05, 0.0);
            target.patches[k].cps[1][3] = Point3::new(0.05, 0.0, t * 0.05);
            target.patches[k].cps[3][1] = Point3::new(0.0, t * 0.05, 0.0);
        }
        let g_full = Bi5Guide::from_gtri5_with_eq8_cc(&g_tri, &target)
            .unwrap_or_else(|e| panic!("n={n}: {e}"));

        let cc_chi = CatmullClarkChi::new(n);
        for k in 0..n {
            let jet_amb = cc_chi.jet(k, ChiCorner::Vk);
            let jet_simplex = ambient_to_simplex_jet(&jet_amb, n, k);
            let (u, v) = jet_simplex.pos;
            let outer = g_full.sector(k).outer_hessian(u, v);
            let composite = compose_chain_rule(&outer, &jet_simplex);

            let target_jet = target.jet(k, ChiCorner::Vk);
            assert!(
                approx_p(composite.value, target_jet.value, 1e-8),
                "n={n}, k={k}: value {:?} vs target {:?}",
                composite.value,
                target_jet.value
            );
            assert!(
                approx_p(composite.d_u, target_jet.d_u, 1e-8)
                    && approx_p(composite.d_v, target_jet.d_v, 1e-8),
                "n={n}, k={k}: first partials mismatch"
            );
            assert!(
                approx_p(composite.d_uu, target_jet.d_uu, 1e-7)
                    && approx_p(composite.d_uv, target_jet.d_uv, 1e-7)
                    && approx_p(composite.d_vv, target_jet.d_vv, 1e-7),
                "n={n}, k={k}: second partials mismatch"
            );
        }
    }
}

#[test]
fn eq8_cc_preserves_other_corners() {
    // Star, m_k, m_{k+1} CPs must not change between from_gtri5 and
    // from_gtri5_with_eq8_cc.
    use super::super::super::kp13_bi5_boundary::TensorBorder;

    let n = 6;
    let params = Bi5FreeParams::zero(n);
    let g_tri = apply_free_params(&params, Point3::new(0.4, -0.2, 0.7)).unwrap();
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(5.0, -3.0, 2.0);
    }

    let g_base = Bi5Guide::from_gtri5(&g_tri);
    let g_full = Bi5Guide::from_gtri5_with_eq8_cc(&g_tri, &target).unwrap();
    for k in 0..n {
        assert_eq!(g_base.sector(k).cps[0][0], g_full.sector(k).cps[0][0]);
        assert_eq!(g_base.sector(k).cps[5][0], g_full.sector(k).cps[5][0]);
        assert_eq!(g_base.sector(k).cps[0][5], g_full.sector(k).cps[0][5]);
    }
}

#[test]
fn eq8_cc_zero_target_zero_star_leaves_interior_cps_zero() {
    // With target = 0 and star = 0, apply_free_params produces the
    // zero g△, so the bi-5 guide is all-zero, and the eq (8) CC
    // override should find x = 0 as the unique solution.
    use super::super::super::kp13_bi5_boundary::TensorBorder;

    for n in [5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let target = TensorBorder::zero(n);
        let g = Bi5Guide::from_gtri5_with_eq8_cc(&g_tri, &target).unwrap();
        for k in 0..n {
            for p in 0..6 {
                for q in 0..6 {
                    let cp = g.sector(k).cps[p][q];
                    assert!(
                        cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                        "n={n}, k={k}, ({p}, {q}): non-zero CP {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn eq8_cc_deterministic() {
    use super::super::super::kp13_bi5_boundary::TensorBorder;

    let n = 6;
    let params = Bi5FreeParams::zero(n);
    let g_tri = apply_free_params(&params, Point3::new(0.5, 0.5, 0.5)).unwrap();
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(k as f64, 1.0, -1.0);
    }
    let a = Bi5Guide::from_gtri5_with_eq8_cc(&g_tri, &target).unwrap();
    let b = Bi5Guide::from_gtri5_with_eq8_cc(&g_tri, &target).unwrap();
    assert_eq!(a, b);
}

#[test]
fn outer_hessian_finite_diff_first_partials() {
    // Sanity: outer_hessian's d_x matches a finite-difference of
    // the eval at (u, v).
    let mut sector = Bi5GuideSector::zero();
    for p in 0..6 {
        for q in 0..6 {
            sector.cps[p][q] = Point3::new(
                (p as f64).sin() + (q as f64).cos(),
                (p + q) as f64 * 0.1,
                0.0,
            );
        }
    }
    let (u, v) = (0.3, 0.7);
    let h = 1e-5;
    let oh = sector.outer_hessian(u, v);
    let v_plus = sector.eval(u + h, v);
    let v_minus = sector.eval(u - h, v);
    let fd_d_u = Point3::new(
        (v_plus.x - v_minus.x) / (2.0 * h),
        (v_plus.y - v_minus.y) / (2.0 * h),
        (v_plus.z - v_minus.z) / (2.0 * h),
    );
    assert!(
        approx_p(oh.d_x, fd_d_u, 1e-4),
        "outer_hessian d_x = {:?}, FD = {:?}",
        oh.d_x,
        fd_d_u
    );
}
