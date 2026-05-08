#![cfg(test)]
use super::tests_support::*;
use super::*;

// ---- Eq (8) bilinear v_k override ----

#[test]
fn eq8_bilinear_matches_target_value() {
    use super::super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::super::kp13_bi5_chi::ChiCorner;
    let params = Bi5FreeParams::zero(5);
    let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();

    // Construct a non-zero tensor border so the v_k jet has content.
    let mut target = TensorBorder::zero(5);
    for k in 0..5 {
        // Set the outer-corner CP of each bi-3 patch to a known point.
        target.patches[k].cps[3][3] = Point3::new(k as f64, k as f64 + 1.0, k as f64 + 2.0);
    }

    let g_full = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &target);

    // After eq (8) override, each sector's g^k(1, 1) = target.jet(k, Vk).value.
    for k in 0..5 {
        let expected = target.jet(k, ChiCorner::Vk).value;
        let got = g_full.sector(k).cps[5][5];
        assert!(
            approx_p(got, expected, 1e-10),
            "sector {k}: g(1,1) = {got:?}, expected = {expected:?}"
        );
    }
}

#[test]
fn eq8_bilinear_matches_target_first_partials() {
    use super::super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::super::kp13_bi5_chi::ChiCorner;

    let params = Bi5FreeParams::zero(5);
    let g_tri = apply_free_params(&params, Point3::new(0.1, 0.1, 0.1)).unwrap();
    let mut target = TensorBorder::zero(5);
    for k in 0..5 {
        target.patches[k].cps[2][3] = Point3::new(0.5 * k as f64, 0.0, 0.0);
        target.patches[k].cps[3][2] = Point3::new(0.0, 0.5 * k as f64, 0.0);
        target.patches[k].cps[3][3] = Point3::new(k as f64, k as f64, k as f64);
    }
    let g_full = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &target);

    for k in 0..5 {
        let target_jet = target.jet(k, ChiCorner::Vk);
        let sec = g_full.sector(k);
        // ∂u g^k(1, 1) = 5 (P_{5,5} - P_{4,5})
        let d_u = Point3::new(
            5.0 * (sec.cps[5][5].x - sec.cps[4][5].x),
            5.0 * (sec.cps[5][5].y - sec.cps[4][5].y),
            5.0 * (sec.cps[5][5].z - sec.cps[4][5].z),
        );
        let d_v = Point3::new(
            5.0 * (sec.cps[5][5].x - sec.cps[5][4].x),
            5.0 * (sec.cps[5][5].y - sec.cps[5][4].y),
            5.0 * (sec.cps[5][5].z - sec.cps[5][4].z),
        );
        assert!(
            approx_p(d_u, target_jet.d_u, 1e-10),
            "sector {k} ∂u: got {d_u:?}, expected {:?}",
            target_jet.d_u
        );
        assert!(
            approx_p(d_v, target_jet.d_v, 1e-10),
            "sector {k} ∂v: got {d_v:?}, expected {:?}",
            target_jet.d_v
        );
    }
}

#[test]
fn eq8_bilinear_matches_target_second_partials() {
    use super::super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::super::kp13_bi5_chi::ChiCorner;

    let params = Bi5FreeParams::zero(6);
    let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mut target = TensorBorder::zero(6);
    for k in 0..6 {
        let t = k as f64;
        target.patches[k].cps[1][3] = Point3::new(t, t + 0.1, t - 0.2);
        target.patches[k].cps[3][1] = Point3::new(-t, t * 0.5, t * 0.3);
        target.patches[k].cps[2][3] = Point3::new(t * 0.2, 0.0, 0.0);
        target.patches[k].cps[3][2] = Point3::new(0.0, t * 0.2, 0.0);
        target.patches[k].cps[2][2] = Point3::new(0.0, 0.0, t * 0.1);
        target.patches[k].cps[3][3] = Point3::new(t, t, t);
    }
    let g_full = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &target);

    for k in 0..6 {
        let target_jet = target.jet(k, ChiCorner::Vk);
        let sec = g_full.sector(k);
        // ∂²u g^k(1, 1) = 20 (P_{5,5} - 2 P_{4,5} + P_{3,5})
        let d_uu = Point3::new(
            20.0 * (sec.cps[5][5].x - 2.0 * sec.cps[4][5].x + sec.cps[3][5].x),
            20.0 * (sec.cps[5][5].y - 2.0 * sec.cps[4][5].y + sec.cps[3][5].y),
            20.0 * (sec.cps[5][5].z - 2.0 * sec.cps[4][5].z + sec.cps[3][5].z),
        );
        let d_vv = Point3::new(
            20.0 * (sec.cps[5][5].x - 2.0 * sec.cps[5][4].x + sec.cps[5][3].x),
            20.0 * (sec.cps[5][5].y - 2.0 * sec.cps[5][4].y + sec.cps[5][3].y),
            20.0 * (sec.cps[5][5].z - 2.0 * sec.cps[5][4].z + sec.cps[5][3].z),
        );
        let d_uv = Point3::new(
            25.0 * (sec.cps[5][5].x - sec.cps[4][5].x - sec.cps[5][4].x + sec.cps[4][4].x),
            25.0 * (sec.cps[5][5].y - sec.cps[4][5].y - sec.cps[5][4].y + sec.cps[4][4].y),
            25.0 * (sec.cps[5][5].z - sec.cps[4][5].z - sec.cps[5][4].z + sec.cps[4][4].z),
        );
        assert!(approx_p(d_uu, target_jet.d_uu, 1e-10), "sector {k} ∂²u");
        assert!(approx_p(d_vv, target_jet.d_vv, 1e-10), "sector {k} ∂²v");
        assert!(approx_p(d_uv, target_jet.d_uv, 1e-10), "sector {k} ∂u∂v");
    }
}

#[test]
fn eq8_preserves_other_corners() {
    // Star, m_k, m_{k+1} corner CPs should NOT change between
    // from_gtri5 and from_gtri5_with_eq8_bilinear.
    use super::super::super::kp13_bi5_boundary::TensorBorder;

    let params = Bi5FreeParams::zero(5);
    let g_tri = apply_free_params(&params, Point3::new(0.3, 0.4, 0.5)).unwrap();
    let mut target = TensorBorder::zero(5);
    for k in 0..5 {
        target.patches[k].cps[3][3] = Point3::new(10.0, 20.0, 30.0);
    }

    let g_base = Bi5Guide::from_gtri5(&g_tri);
    let g_full = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &target);

    for k in 0..5 {
        // Star corner (0, 0):
        assert_eq!(g_base.sector(k).cps[0][0], g_full.sector(k).cps[0][0]);
        // m_k corner (5, 0):
        assert_eq!(g_base.sector(k).cps[5][0], g_full.sector(k).cps[5][0]);
        // m_{k+1} corner (0, 5):
        assert_eq!(g_base.sector(k).cps[0][5], g_full.sector(k).cps[0][5]);
    }
}

#[test]
fn eq8_bilinear_works_for_all_valences() {
    use super::super::super::kp13_bi5_boundary::TensorBorder;

    for n in [3, 4, 5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let target = TensorBorder::zero(n);
        let g = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &target);
        g.validate().unwrap_or_else(|e| panic!("n={n}: {e}"));
    }
}
