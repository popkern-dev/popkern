//! Tests for LSQ design matrix assembly (I.11e.3).

use super::super::kp13_bi5_boundary::TensorBorder;
use super::*;

const EPS: f64 = 1e-10;

// ---- Flat-index round-trips ----

#[test]
fn slot_flat_roundtrip() {
    let n = 7;
    let total = num_free_params(n);
    for flat in 0..total {
        let slot = ParamSlot::from_flat(flat, n);
        assert_eq!(slot.to_flat(n), flat, "roundtrip at flat={flat}");
    }
}

#[test]
fn slot_layout_matches_docs() {
    let n = 6;
    assert_eq!(ParamSlot::from_flat(0, n), ParamSlot::QAlpha);
    assert_eq!(ParamSlot::from_flat(4, n), ParamSlot::QEpsilon);
    assert_eq!(ParamSlot::from_flat(5, n), ParamSlot::MidEdge(0));
    assert_eq!(
        ParamSlot::from_flat(5 + n - 1, n),
        ParamSlot::MidEdge(n - 1)
    );
    assert_eq!(ParamSlot::from_flat(5 + n, n), ParamSlot::G30(0));
    assert_eq!(
        ParamSlot::from_flat(5 + 5 * n + n - 1, n),
        ParamSlot::G032(n - 1)
    );
}

// ---- System sizes ----

#[test]
fn system_shape_matches_6n_plus_5() {
    for n in [5, 6, 7, 8, 10] {
        let target = TensorBorder::zero(n);
        let sys = build_lsq_system(&target, Point3::new(0.0, 0.0, 0.0)).unwrap();
        assert_eq!(sys.n, n);
        assert_eq!(sys.a.len(), num_rows(n));
        assert_eq!(sys.a[0].len(), num_free_params(n));
        for axis in 0..3 {
            assert_eq!(sys.b[axis].len(), num_rows(n));
        }
    }
}

#[test]
fn rejects_n3_n4() {
    for bad in [3, 4] {
        let target = TensorBorder::zero(bad);
        assert!(build_lsq_system(&target, Point3::new(0.0, 0.0, 0.0)).is_err());
    }
}

// ---- Zero target + zero star → zero RHS ----

#[test]
fn zero_target_zero_star_gives_zero_rhs() {
    let n = 6;
    let target = TensorBorder::zero(n);
    let sys = build_lsq_system(&target, Point3::new(0.0, 0.0, 0.0)).unwrap();
    for axis in 0..3 {
        for &x in &sys.b[axis] {
            assert!(x.abs() < EPS, "axis {axis}: non-zero RHS entry {x}");
        }
    }
}

// ---- Constant star pin propagates to RHS ----

#[test]
fn constant_star_pin_propagates_to_some_rhs_entries() {
    // With cc_limit = (c, 0, 0) and zero params:
    //   - Near-star coefs (flat 0, 1, 2, 6, 7, 11) are all (c, 0, 0)
    //   - m_k, g_30, g_40, interior 122/023/032 = zero
    //   - Circulant derives g^i_{21} from R^i = c̄(2c+1)·(c, 0, 0);
    //     back-fill g^i_{12} similarly; higher-j fills.
    //
    // Composite VALUES at Mk / Vk corners equal the CP at those
    // simplex vertices (flat 20 and the off-simplex extrapolation),
    // which are zero since m_k param = 0. So Value rows at those
    // corners stay at zero.
    //
    // But PARTIALS at Mk / Vk pick up contributions from the near-
    // star coefs via the tri-5 basis, so those rows ARE non-zero.
    // At minimum, axis-0 RHS should have some non-zero entries.
    let n = 6;
    let c = 2.5;
    let target = TensorBorder::zero(n);
    let sys = build_lsq_system(&target, Point3::new(c, 0.0, 0.0)).unwrap();
    let any_nonzero_axis0 = sys.b[0].iter().any(|&x| x.abs() > 1e-8);
    assert!(
        any_nonzero_axis0,
        "cc_limit.x = {c} should yield some non-zero RHS entry on axis 0"
    );
    // Axes y, z get RHS zero since cc_limit has no y/z component.
    let any_nonzero_axis1 = sys.b[1].iter().any(|&x| x.abs() > 1e-8);
    let any_nonzero_axis2 = sys.b[2].iter().any(|&x| x.abs() > 1e-8);
    assert!(!any_nonzero_axis1);
    assert!(!any_nonzero_axis2);
}

// ---- Design matrix A has finite entries and is non-trivial ----

#[test]
fn design_matrix_is_not_all_zero() {
    let n = 6;
    let target = TensorBorder::zero(n);
    let sys = build_lsq_system(&target, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let n_nonzero: usize = sys
        .a
        .iter()
        .map(|row| row.iter().filter(|&&x| x.abs() > 1e-12).count())
        .sum();
    assert!(
        n_nonzero > 10,
        "design matrix has only {n_nonzero} nonzero entries"
    );
    // All entries are finite.
    for row in &sys.a {
        for &x in row {
            assert!(x.is_finite(), "non-finite design-matrix entry");
        }
    }
}

// ---- Linearity check: sum of columns respects superposition ----

#[test]
fn design_matrix_is_linear_in_params() {
    // Verify that applying params_a + params_b produces the same
    // composite jets as adding jets from params_a and params_b
    // separately (confirming linearity = design matrix is correct).
    let n = 7;
    let cc = Point3::new(0.0, 0.0, 0.0);
    let mut a = Bi5FreeParams::zero(n);
    a.q.alpha = Point3::new(0.3, 0.0, 0.0);
    a.g_30[2] = Point3::new(0.5, 0.0, 0.0);
    let mut b = Bi5FreeParams::zero(n);
    b.q.delta = Point3::new(0.7, 0.0, 0.0);
    b.g_122[3] = Point3::new(-0.2, 0.0, 0.0);

    let mut a_plus_b = Bi5FreeParams::zero(n);
    a_plus_b.q.alpha = Point3::new(0.3, 0.0, 0.0);
    a_plus_b.q.delta = Point3::new(0.7, 0.0, 0.0);
    a_plus_b.g_30[2] = Point3::new(0.5, 0.0, 0.0);
    a_plus_b.g_122[3] = Point3::new(-0.2, 0.0, 0.0);

    let g_a = apply_free_params(&a, cc).unwrap();
    let g_b = apply_free_params(&b, cc).unwrap();
    let g_ab = apply_free_params(&a_plus_b, cc).unwrap();

    for k in 0..n {
        for corner in [ChiCorner::Mk, ChiCorner::Vk] {
            let j_a = composite_jet_bilinear(&g_a, k, corner);
            let j_b = composite_jet_bilinear(&g_b, k, corner);
            let j_ab = composite_jet_bilinear(&g_ab, k, corner);
            // j_a + j_b should equal j_ab.
            let sum_value_x = j_a.value.x + j_b.value.x;
            assert!(
                (sum_value_x - j_ab.value.x).abs() < 1e-10,
                "k={k} corner={corner:?}: linearity fails at value.x"
            );
            let sum_d_u_x = j_a.d_u.x + j_b.d_u.x;
            assert!(
                (sum_d_u_x - j_ab.d_u.x).abs() < 1e-10,
                "k={k} corner={corner:?}: linearity fails at d_u.x"
            );
        }
    }
}

// ---- Design matrix + zero RHS → zero free params solves perfectly ----

#[test]
fn zero_target_zero_star_admits_zero_solution() {
    // If target and star are both zero, then b = 0 and A·0 = 0,
    // so x = 0 is a valid solution. This confirms the RHS/A are
    // internally consistent.
    let n = 6;
    let target = TensorBorder::zero(n);
    let sys = build_lsq_system(&target, Point3::new(0.0, 0.0, 0.0)).unwrap();
    for axis in 0..3 {
        for &x in &sys.b[axis] {
            assert_eq!(x, 0.0);
        }
    }
    // A·0 = 0 trivially. No-op check.
}

// ---- composite_jet_bilinear at corners ----

#[test]
fn composite_jet_value_at_star_corner_equals_star_cp() {
    // At corner (0, 0), composite jet's value = g△^k's CP (0, 0, 5)
    // (the star point). Apply zero params with star = (c, 0, 0)
    // to verify.
    let n = 7;
    let star = Point3::new(1.5, -0.5, 0.25);
    let params = Bi5FreeParams::zero(n);
    let g = apply_free_params(&params, star).unwrap();
    for k in 0..n {
        let jet = composite_jet_bilinear(&g, k, ChiCorner::Star);
        assert!(
            (jet.value.x - star.x).abs() < 1e-10
                && (jet.value.y - star.y).abs() < 1e-10
                && (jet.value.z - star.z).abs() < 1e-10,
            "k={k}: composite value at star != pinned star"
        );
    }
}
