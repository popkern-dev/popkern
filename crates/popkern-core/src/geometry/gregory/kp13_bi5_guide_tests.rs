//! Tests for the bi-5 guide construction (KP13 § 3).

#![allow(clippy::needless_range_loop)]

use super::super::bernstein_tri5::{
    TRI5_CP_COUNT, bernstein_tri5_all, bernstein_tri5_all_extended,
};
use super::super::kp13_bi5_free_params::{Bi5FreeParams, apply_free_params};
use super::super::kp13_bi5_gtri5::{IDX_005, IDX_050, IDX_500};
use super::*;

const EPS: f64 = 1e-9;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// Evaluate a tri-5 polynomial at (u, v) using the codebase's
// extended basis (no simplex restriction — needed because bi-5 eval
// probes u + v > 1).
fn eval_tri5_extended(cps: &[Point3; TRI5_CP_COUNT], u: f64, v: f64) -> Point3 {
    let b = bernstein_tri5_all_extended(u, v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..TRI5_CP_COUNT {
        out.x += b[i] * cps[i].x;
        out.y += b[i] * cps[i].y;
        out.z += b[i] * cps[i].z;
    }
    out
}

// ---- tri5_to_bi5_degree_raise: polynomial equivalence ----

#[test]
fn degree_raise_preserves_polynomial_inside_simplex() {
    // Build a deterministic tri-5 CP vector and verify that the
    // degree-raised bi-5 patch evaluates to the same value as the
    // tri-5 patch, at interior points of the simplex.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        let t = i as f64 * 0.37;
        cps[i] = Point3::new(t.sin(), t.cos(), 0.5 * t.sin() * t.cos());
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    let sector = Bi5GuideSector { cps: bi5 };

    for &(u, v) in &[(0.1, 0.1), (0.25, 0.5), (0.3, 0.3), (0.5, 0.25), (0.7, 0.2)] {
        let tri_val = {
            let b = bernstein_tri5_all(u, v);
            let mut out = Point3::new(0.0, 0.0, 0.0);
            for i in 0..TRI5_CP_COUNT {
                out.x += b[i] * cps[i].x;
                out.y += b[i] * cps[i].y;
                out.z += b[i] * cps[i].z;
            }
            out
        };
        let bi_val = sector.eval(u, v);
        assert!(
            approx_p(tri_val, bi_val, 1e-10),
            "(u, v) = ({u}, {v}): tri={tri_val:?}, bi={bi_val:?}"
        );
    }
}

#[test]
fn degree_raise_preserves_polynomial_outside_simplex() {
    // Key observation: g△ is polynomial, so the identity
    //   tri5_polynomial(u, v) = bi5_polynomial(u, v)
    // holds for ALL (u, v) ∈ R², including (u, v) ∈ [0, 1]² where
    // u + v > 1 (the v_k region).
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        let t = i as f64 * 0.123;
        cps[i] = Point3::new(t.cos(), -t.sin(), t.sin() + 0.2);
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    let sector = Bi5GuideSector { cps: bi5 };

    for &(u, v) in &[(0.8, 0.7), (1.0, 1.0), (0.5, 0.9), (0.95, 0.6)] {
        let tri_val = eval_tri5_extended(&cps, u, v);
        let bi_val = sector.eval(u, v);
        assert!(
            approx_p(tri_val, bi_val, 1e-9),
            "(u, v) = ({u}, {v}) in v_k region: tri={tri_val:?}, bi={bi_val:?}"
        );
    }
}

#[test]
fn degree_raise_linearity() {
    // The degree-raise is linear in the tri-5 CPs.
    let mut cps_a = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    let mut cps_b = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps_a[i] = Point3::new(i as f64, 0.0, 0.0);
        cps_b[i] = Point3::new(0.0, i as f64, (-1.0_f64).powi(i as i32) * 0.5);
    }
    let bi_a = tri5_to_bi5_degree_raise(&cps_a);
    let bi_b = tri5_to_bi5_degree_raise(&cps_b);
    // Sum:
    let mut cps_sum = cps_a;
    for i in 0..TRI5_CP_COUNT {
        cps_sum[i].x += cps_b[i].x;
        cps_sum[i].y += cps_b[i].y;
        cps_sum[i].z += cps_b[i].z;
    }
    let bi_sum = tri5_to_bi5_degree_raise(&cps_sum);
    for p in 0..6 {
        for q in 0..6 {
            let expected_x = bi_a[p][q].x + bi_b[p][q].x;
            let expected_y = bi_a[p][q].y + bi_b[p][q].y;
            let expected_z = bi_a[p][q].z + bi_b[p][q].z;
            assert!(
                (bi_sum[p][q].x - expected_x).abs() < 1e-9
                    && (bi_sum[p][q].y - expected_y).abs() < 1e-9
                    && (bi_sum[p][q].z - expected_z).abs() < 1e-9,
                "linearity fails at ({p}, {q})"
            );
        }
    }
}

// ---- Corner CP identities ----

#[test]
fn bi5_star_corner_equals_tri5_star() {
    // At bi-5 CP (0, 0) the value is the patch value at (0, 0),
    // which in simplex coords is the star vertex (u_s = 0, v_s = 0,
    // w_s = 1). So bi5[0][0] = cps[IDX_500].
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps[i] = Point3::new(i as f64 + 1.0, i as f64 * 2.0, -(i as f64));
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    assert!(approx_p(bi5[0][0], cps[IDX_500], 1e-12));
}

#[test]
fn bi5_mk_corner_equals_tri5_mk() {
    // At bi-5 CP (5, 0): u = 1, v = 0 → simplex coords (1, 0, 0) = m_k.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps[i] = Point3::new(i as f64, -(i as f64), i as f64 / 2.0);
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    // IDX_050 = (i = 5, j = 0) is the m_k corner in this codebase's
    // convention (paper label "050" = barycentric (0, 5, 0)).
    assert!(approx_p(bi5[5][0], cps[IDX_050], 1e-12));
}

#[test]
fn bi5_mk1_corner_equals_tri5_mk1() {
    // At bi-5 CP (0, 5): u = 0, v = 1 → simplex (0, 1, 0) = m_{k+1}.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps[i] = Point3::new(-(i as f64), i as f64, i as f64 + 1.0);
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    // IDX_005 = (i = 0, j = 5) is the m_{k+1} corner.
    assert!(approx_p(bi5[0][5], cps[IDX_005], 1e-12));
}

// ---- Bi5Guide::from_gtri5: end-to-end ----

#[test]
fn from_gtri5_preserves_star_per_sector() {
    // All sectors share the same star point (I.11 invariant). After
    // degree-raise, each sector's bi-5 CP (0, 0) should equal the
    // shared star.
    let params = Bi5FreeParams::zero(6);
    let star = Point3::new(0.5, -0.2, 1.3);
    let g_tri = apply_free_params(&params, star).unwrap();
    let g_bi = Bi5Guide::from_gtri5(&g_tri);
    for k in 0..6 {
        let bi_star = g_bi.sector(k).cps[0][0];
        assert!(
            approx_p(bi_star, star, 1e-10),
            "sector {k} bi-5 star {bi_star:?} ≠ tri-5 star {star:?}"
        );
    }
}

#[test]
fn from_gtri5_c0_junction_between_sectors() {
    // Adjacent sectors share a boundary edge. Sector k's (0, 5)
    // corner = m_{k+1} = sector (k+1)'s (5, 0) corner = m_{k+1}.
    // The degree-raised guide inherits this shared edge.
    let params = Bi5FreeParams::zero(7);
    let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let g_bi = Bi5Guide::from_gtri5(&g_tri);
    for k in 0..7 {
        let k_next = (k + 1) % 7;
        let mk1_via_k = g_bi.sector(k).cps[0][5];
        let mk1_via_next = g_bi.sector(k_next).cps[5][0];
        assert!(
            approx_p(mk1_via_k, mk1_via_next, 1e-10),
            "sector {k}→{k_next} mid-edge mismatch: {mk1_via_k:?} vs {mk1_via_next:?}"
        );
    }
}

#[test]
fn from_gtri5_works_for_all_valences() {
    // Including low-valence n ∈ {3, 4} (via KP09a § 6 pseudo-inverse
    // in the circulant upstream).
    for n in [3, 4, 5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.1, 0.2, 0.3)).unwrap();
        let g_bi = Bi5Guide::from_gtri5(&g_tri);
        g_bi.validate().unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert_eq!(g_bi.n, n);
        assert_eq!(g_bi.sectors.len(), n);
    }
}

// ---- Validation ----

#[test]
fn zero_guide_validates() {
    for n in [3, 4, 5, 6, 10] {
        let g = Bi5Guide::zero(n);
        g.validate().unwrap();
    }
}

#[test]
#[should_panic]
fn zero_guide_below_three_panics() {
    let _ = Bi5Guide::zero(2);
}

#[test]
fn validate_detects_sector_count_mismatch() {
    let mut g = Bi5Guide::zero(5);
    g.sectors.pop();
    assert!(g.validate().is_err());
}

#[test]
fn validate_detects_non_finite_cp() {
    let mut g = Bi5Guide::zero(5);
    g.sectors[2].cps[3][4].x = f64::NAN;
    assert!(g.validate().is_err());
}

// ---- Evaluation ----

#[test]
fn sector_eval_at_corners_equals_corner_cps() {
    let mut sector = Bi5GuideSector::zero();
    sector.cps[0][0] = Point3::new(1.0, 2.0, 3.0);
    sector.cps[5][0] = Point3::new(4.0, 5.0, 6.0);
    sector.cps[0][5] = Point3::new(7.0, 8.0, 9.0);
    sector.cps[5][5] = Point3::new(10.0, 11.0, 12.0);
    assert!(approx_p(sector.eval(0.0, 0.0), sector.cps[0][0], 1e-12));
    assert!(approx_p(sector.eval(1.0, 0.0), sector.cps[5][0], 1e-12));
    assert!(approx_p(sector.eval(0.0, 1.0), sector.cps[0][5], 1e-12));
    assert!(approx_p(sector.eval(1.0, 1.0), sector.cps[5][5], 1e-12));
}

#[test]
fn eval_is_continuous_across_sector_boundary() {
    // Sample the shared edge at v = 0 for sector k+1 vs u = 1 for
    // sector k (both trace out the m_{k+1}-side of the boundary).
    let params = Bi5FreeParams::zero(6);
    let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let g_bi = Bi5Guide::from_gtri5(&g_tri);
    // Sector k's v = 1 edge vs sector (k + 1)'s u = 0 edge (both =
    // the shared m_{k+1} edge; note different parameter directions
    // along the edge in the two sectors, so we compare endpoints).
    for k in 0..6 {
        let k_next = (k + 1) % 6;
        let a = g_bi.sector(k).eval(0.0, 1.0); // m_{k+1} via sector k
        let b = g_bi.sector(k_next).eval(1.0, 0.0); // m_{k+1} via sector k+1
        assert!(
            approx_p(a, b, 1e-10),
            "C⁰ junction at k={k}, k_next={k_next}: {a:?} vs {b:?}"
        );
    }
}

// ---- Eq (8) bilinear v_k override ----

#[test]
fn eq8_bilinear_matches_target_value() {
    use super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::kp13_bi5_chi::ChiCorner;
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
    use super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::kp13_bi5_chi::ChiCorner;

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
    use super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::kp13_bi5_chi::ChiCorner;

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
    use super::super::kp13_bi5_boundary::TensorBorder;

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
    use super::super::kp13_bi5_boundary::TensorBorder;

    for n in [3, 4, 5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let target = TensorBorder::zero(n);
        let g = Bi5Guide::from_gtri5_with_eq8_bilinear(&g_tri, &target);
        g.validate().unwrap_or_else(|e| panic!("n={n}: {e}"));
    }
}

// ---- CC-χ eq (8) variant ----

#[test]
fn eq8_cc_matches_target_value_and_partials() {
    // With full CC χ, the composite (g ∘ L⁻¹ ∘ χ) jet at the v_k
    // corner must match b^k's jet there. We verify by re-evaluating
    // the composite through the same pipeline and checking agreement.
    use super::super::kp13_bi5_boundary::TensorBorder;
    use super::super::kp13_bi5_chi::{CatmullClarkChi, ChiCorner};
    use super::super::kp13_bi5_jet::compose_chain_rule;
    use super::super::kp13_bi5_lsq_cc::ambient_to_simplex_jet;

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
    use super::super::kp13_bi5_boundary::TensorBorder;

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
    use super::super::kp13_bi5_boundary::TensorBorder;

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
    use super::super::kp13_bi5_boundary::TensorBorder;

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

// ---- Determinism ----

#[test]
fn from_gtri5_deterministic() {
    let params = Bi5FreeParams::zero(5);
    let g_tri = apply_free_params(&params, Point3::new(1.0, 2.0, 3.0)).unwrap();
    let a = Bi5Guide::from_gtri5(&g_tri);
    let b = Bi5Guide::from_gtri5(&g_tri);
    assert_eq!(a, b);
}
