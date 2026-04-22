//! Tests for the bi-3 tensor-border jet extractor.

use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- Zero patch ----

#[test]
fn zero_patch_has_zero_jets_at_every_corner() {
    let patch = Bi3TensorBorder::zero();
    for corner in [
        ChiCorner::Star,
        ChiCorner::Mk,
        ChiCorner::Mk1,
        ChiCorner::Vk,
    ] {
        let j = patch.jet_at_corner(corner);
        assert_eq!(j.sq_norm(), 0.0, "corner {corner:?}");
    }
}

// ---- Constant patch ----

#[test]
fn constant_patch_has_constant_value_and_zero_partials() {
    let c = Point3::new(1.5, -2.0, 3.25);
    let cps = [[c; 4]; 4];
    let patch = Bi3TensorBorder::from_cps(cps);
    for corner in [
        ChiCorner::Star,
        ChiCorner::Mk,
        ChiCorner::Mk1,
        ChiCorner::Vk,
    ] {
        let j = patch.jet_at_corner(corner);
        assert!(approx_p(j.value, c));
        assert!(approx_p(j.d_u, Point3::new(0.0, 0.0, 0.0)));
        assert!(approx_p(j.d_v, Point3::new(0.0, 0.0, 0.0)));
        assert!(approx_p(j.d_uu, Point3::new(0.0, 0.0, 0.0)));
        assert!(approx_p(j.d_vv, Point3::new(0.0, 0.0, 0.0)));
        assert!(approx_p(j.d_uv, Point3::new(0.0, 0.0, 0.0)));
    }
}

// ---- Linear patch ----

/// Build a bi-3 tensor-border that represents the affine function
/// `f(u, v) = a·u + b·v + c` (as a `Point3` value).
fn linear_cps(a: Point3, b: Point3, c: Point3) -> [[Point3; 4]; 4] {
    let mut cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in cps.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let u = i as f64 / 3.0;
            let v = j as f64 / 3.0;
            *cell = Point3::new(
                a.x * u + b.x * v + c.x,
                a.y * u + b.y * v + c.y,
                a.z * u + b.z * v + c.z,
            );
        }
    }
    cps
}

#[test]
fn linear_patch_jet_star() {
    let a = Point3::new(2.0, 0.0, 0.0);
    let b = Point3::new(0.0, 3.0, 0.0);
    let c = Point3::new(1.0, 2.0, 3.0);
    let patch = Bi3TensorBorder::from_cps(linear_cps(a, b, c));
    let j = patch.jet_at_corner(ChiCorner::Star);
    // At (0, 0): value = c, ∂u = a, ∂v = b, higher = 0.
    assert!(approx_p(j.value, c));
    assert!(approx_p(j.d_u, a));
    assert!(approx_p(j.d_v, b));
    assert!(approx_p(j.d_uu, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(j.d_vv, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(j.d_uv, Point3::new(0.0, 0.0, 0.0)));
}

#[test]
fn linear_patch_star_has_value_c() {
    let c = Point3::new(7.0, -3.0, 2.0);
    let patch = Bi3TensorBorder::from_cps(linear_cps(
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        c,
    ));
    let j = patch.jet_at_corner(ChiCorner::Star);
    assert!(approx_p(j.value, c));
}

#[test]
fn linear_patch_mk_value_is_a_plus_c() {
    // At Mk = (1, 0): f(1, 0) = a + c.
    let a = Point3::new(2.0, 0.0, 0.0);
    let b = Point3::new(0.0, 3.0, 0.0);
    let c = Point3::new(1.0, 1.0, 1.0);
    let patch = Bi3TensorBorder::from_cps(linear_cps(a, b, c));
    let j = patch.jet_at_corner(ChiCorner::Mk);
    let expected_value = Point3::new(a.x + c.x, a.y + c.y, a.z + c.z);
    assert!(approx_p(j.value, expected_value));
}

#[test]
fn linear_patch_mk1_value_is_b_plus_c() {
    // At Mk1 = (0, 1): f(0, 1) = b + c.
    let a = Point3::new(2.0, 0.0, 0.0);
    let b = Point3::new(0.0, 3.0, 0.0);
    let c = Point3::new(1.0, 1.0, 1.0);
    let patch = Bi3TensorBorder::from_cps(linear_cps(a, b, c));
    let j = patch.jet_at_corner(ChiCorner::Mk1);
    let expected_value = Point3::new(b.x + c.x, b.y + c.y, b.z + c.z);
    assert!(approx_p(j.value, expected_value));
}

#[test]
fn linear_patch_vk_value_is_a_plus_b_plus_c() {
    // At Vk = (1, 1): f(1, 1) = a + b + c.
    let a = Point3::new(2.0, 0.0, 0.0);
    let b = Point3::new(0.0, 3.0, 0.0);
    let c = Point3::new(1.0, 1.0, 1.0);
    let patch = Bi3TensorBorder::from_cps(linear_cps(a, b, c));
    let j = patch.jet_at_corner(ChiCorner::Vk);
    let expected_value = Point3::new(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z);
    assert!(approx_p(j.value, expected_value));
}

// ---- Corner-specific jet-value formulas (value-only sanity) ----

#[test]
fn corner_values_match_cps_at_respective_indices() {
    let mut cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    cps[0][0] = Point3::new(10.0, 0.0, 0.0);
    cps[3][0] = Point3::new(20.0, 0.0, 0.0);
    cps[0][3] = Point3::new(30.0, 0.0, 0.0);
    cps[3][3] = Point3::new(40.0, 0.0, 0.0);
    let patch = Bi3TensorBorder::from_cps(cps);
    assert!(approx_p(
        patch.jet_at_corner(ChiCorner::Star).value,
        Point3::new(10.0, 0.0, 0.0)
    ));
    assert!(approx_p(
        patch.jet_at_corner(ChiCorner::Mk).value,
        Point3::new(20.0, 0.0, 0.0)
    ));
    assert!(approx_p(
        patch.jet_at_corner(ChiCorner::Mk1).value,
        Point3::new(30.0, 0.0, 0.0)
    ));
    assert!(approx_p(
        patch.jet_at_corner(ChiCorner::Vk).value,
        Point3::new(40.0, 0.0, 0.0)
    ));
}

// ---- TensorBorder ----

#[test]
fn tensor_border_zero_has_right_length() {
    for n in [3, 5, 6, 7, 8, 10] {
        let tb = TensorBorder::zero(n);
        assert_eq!(tb.n(), n);
    }
}

#[test]
#[should_panic]
fn tensor_border_below_three_panics() {
    let _ = TensorBorder::zero(2);
}

#[test]
fn tensor_border_jet_per_sector() {
    let n = 6;
    let mut tb = TensorBorder::zero(n);
    // Set distinct constants per sector.
    for (k, patch) in tb.patches.iter_mut().enumerate() {
        let c = Point3::new(k as f64, 0.0, 0.0);
        for row in patch.cps.iter_mut() {
            for cp in row.iter_mut() {
                *cp = c;
            }
        }
    }
    for k in 0..n {
        let j = tb.jet(k, ChiCorner::Star);
        assert!(approx_p(j.value, Point3::new(k as f64, 0.0, 0.0)));
        assert!(approx_p(j.d_u, Point3::new(0.0, 0.0, 0.0)));
    }
}

// ---- BB corner-derivative identities numerically verified ----

#[test]
fn bi3_corner_at_star_matches_bb_formula() {
    // Specific CPs near the (0, 0) corner; check that the jet
    // matches the paper formula exactly.
    let p00 = Point3::new(1.0, 0.0, 0.0);
    let p10 = Point3::new(2.0, 0.0, 0.0);
    let p01 = Point3::new(3.0, 0.0, 0.0);
    let p20 = Point3::new(4.0, 0.0, 0.0);
    let p02 = Point3::new(5.0, 0.0, 0.0);
    let p11 = Point3::new(6.0, 0.0, 0.0);
    let mut cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    cps[0][0] = p00;
    cps[1][0] = p10;
    cps[0][1] = p01;
    cps[2][0] = p20;
    cps[0][2] = p02;
    cps[1][1] = p11;
    let patch = Bi3TensorBorder::from_cps(cps);
    let j = patch.jet_at_corner(ChiCorner::Star);
    // value = p00 = 1
    assert!(approx_p(j.value, Point3::new(1.0, 0.0, 0.0)));
    // ∂u = 3·(p10 - p00) = 3·(2 - 1) = 3
    assert!(approx_p(j.d_u, Point3::new(3.0, 0.0, 0.0)));
    // ∂v = 3·(p01 - p00) = 3·2 = 6
    assert!(approx_p(j.d_v, Point3::new(6.0, 0.0, 0.0)));
    // ∂²u = 6·(p20 - 2·p10 + p00) = 6·(4 - 4 + 1) = 6
    assert!(approx_p(j.d_uu, Point3::new(6.0, 0.0, 0.0)));
    // ∂²v = 6·(p02 - 2·p01 + p00) = 6·(5 - 6 + 1) = 0
    assert!(approx_p(j.d_vv, Point3::new(0.0, 0.0, 0.0)));
    // ∂²uv = 9·(p11 - p10 - p01 + p00) = 9·(6 - 2 - 3 + 1) = 18
    assert!(approx_p(j.d_uv, Point3::new(18.0, 0.0, 0.0)));
}
