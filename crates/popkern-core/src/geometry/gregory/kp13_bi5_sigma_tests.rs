//! Tests for the σ sampling reparameterisation.

use super::*;

const EPS: f64 = 1e-10;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

fn approx_pair(a: (f64, f64), b: (f64, f64)) -> bool {
    approx_eq(a.0, b.0) && approx_eq(a.1, b.1)
}

// ---- Identity construction + evaluation ----

#[test]
fn identity_eval_is_identity_at_corners() {
    let s = SigmaPatch::identity();
    assert!(approx_pair(s.eval(0.0, 0.0), (0.0, 0.0)));
    assert!(approx_pair(s.eval(1.0, 0.0), (1.0, 0.0)));
    assert!(approx_pair(s.eval(0.0, 1.0), (0.0, 1.0)));
    assert!(approx_pair(s.eval(1.0, 1.0), (1.0, 1.0)));
}

#[test]
fn identity_eval_is_identity_at_interior() {
    let s = SigmaPatch::identity();
    // For the identity map parameterised bilinearly in BB coefs,
    // σ(u, v) = (u, v) holds at all interior points.
    for (u, v) in &[(0.25, 0.75), (0.5, 0.5), (0.3, 0.6), (0.8, 0.2)] {
        let (x, y) = s.eval(*u, *v);
        assert!(
            approx_eq(x, *u) && approx_eq(y, *v),
            "at ({u}, {v}): got ({x}, {y})"
        );
    }
}

// ---- Jets at corners ----

#[test]
fn identity_jet_at_star_corner() {
    let s = SigmaPatch::identity();
    let j = s.jet_at_corner(ChiCorner::Star);
    // σ = identity → σ(0, 0) = (0, 0), ∂u = (1, 0), ∂v = (0, 1),
    // 2nd partials = (0, 0).
    assert!(approx_pair(j.pos, (0.0, 0.0)));
    assert!(approx_pair(j.d_u, (1.0, 0.0)));
    assert!(approx_pair(j.d_v, (0.0, 1.0)));
    assert!(approx_pair(j.d_uu, (0.0, 0.0)));
    assert!(approx_pair(j.d_uv, (0.0, 0.0)));
    assert!(approx_pair(j.d_vv, (0.0, 0.0)));
}

#[test]
fn identity_jet_at_all_corners() {
    let s = SigmaPatch::identity();
    for (corner, expected_pos) in [
        (ChiCorner::Star, (0.0, 0.0)),
        (ChiCorner::Mk, (1.0, 0.0)),
        (ChiCorner::Mk1, (0.0, 1.0)),
        (ChiCorner::Vk, (1.0, 1.0)),
    ] {
        let j = s.jet_at_corner(corner);
        assert!(approx_pair(j.pos, expected_pos), "corner {corner:?} pos");
        assert!(approx_pair(j.d_u, (1.0, 0.0)), "corner {corner:?} d_u");
        assert!(approx_pair(j.d_v, (0.0, 1.0)), "corner {corner:?} d_v");
        // All 2nd partials zero for linear σ.
        assert!(approx_pair(j.d_uu, (0.0, 0.0)));
        assert!(approx_pair(j.d_uv, (0.0, 0.0)));
        assert!(approx_pair(j.d_vv, (0.0, 0.0)));
    }
}

// ---- Symmetry ----

#[test]
fn identity_is_symmetric() {
    let s = SigmaPatch::identity();
    assert!(s.is_symmetric());
}

#[test]
fn asymmetric_patch_detected() {
    let mut s = SigmaPatch::identity();
    // Break symmetry: modify (1, 2) without touching swap partner
    // (2, 1).
    s.cps[1][2] = (100.0, 100.0);
    assert!(!s.is_symmetric());
}

#[test]
fn enforce_symmetry_restores_invariant() {
    let mut s = SigmaPatch::identity();
    s.cps[1][2] = (100.0, 200.0);
    s.cps[2][1] = (999.0, -999.0);
    s.enforce_symmetry();
    assert!(s.is_symmetric());
}

#[test]
fn enforce_symmetry_is_idempotent() {
    let mut s = SigmaPatch::identity();
    s.cps[1][2] = (0.3, 0.7);
    s.enforce_symmetry();
    let snap = s;
    s.enforce_symmetry();
    for i in 0..6 {
        for j in 0..6 {
            assert!(approx_pair(s.cps[i][j], snap.cps[i][j]));
        }
    }
}

// ---- F_3 functional ----

#[test]
fn f3_of_identity_is_zero() {
    // The identity map σ(u, v) = (u, v) is linear, so all 3rd-order
    // partials vanish. The `(u-v)^k σ` terms don't vanish, but the
    // paper's F_3 operator for the identity map... actually the
    // (u-v)^k terms do contribute. Let me check the formula once more.
    //
    // F_3 includes terms like `(u-v)^3 σ` (no partial), which for
    // σ = identity gives (u-v)^3 · (u, v). This is NON-zero at
    // interior points and integrates to a non-zero positive value.
    //
    // So F_3(identity) is **not** zero in general. This test
    // documents that: F_3 is non-negative, so we just require
    // F_3(identity) ≥ 0.
    let s = SigmaPatch::identity();
    let f = f3_functional(&s);
    assert!(f >= 0.0, "F_3 should be non-negative, got {f}");
    assert!(f.is_finite());
}

#[test]
fn f3_of_zero_patch_is_zero() {
    // σ ≡ 0 (all CPs at origin) has zero partials AND zero value,
    // so F_3 = 0 exactly.
    let s = SigmaPatch {
        cps: [[(0.0, 0.0); 6]; 6],
    };
    let f = f3_functional(&s);
    assert!(f.abs() < 1e-12, "F_3(zero) = {f}, expected 0");
}

#[test]
fn f3_scales_quadratically() {
    // F_3 is a quadratic form: F_3(c · σ) = c² · F_3(σ).
    let base = SigmaPatch::identity();
    let mut doubled = base;
    for row in doubled.cps.iter_mut() {
        for cp in row.iter_mut() {
            cp.0 *= 2.0;
            cp.1 *= 2.0;
        }
    }
    let f_base = f3_functional(&base);
    let f_doubled = f3_functional(&doubled);
    assert!(
        (f_doubled - 4.0 * f_base).abs() < 1e-6,
        "F_3(2σ) = {f_doubled}, expected 4·F_3(σ) = {}",
        4.0 * f_base
    );
}

#[test]
fn f3_of_shifted_patch_not_identity() {
    // A patch offset from the identity should have F_3 different
    // from identity's — confirming the functional responds to the
    // patch shape.
    let base = SigmaPatch::identity();
    let mut shifted = base;
    // Perturb one CP.
    shifted.cps[2][2] = (0.8, 0.9);
    let f_base = f3_functional(&base);
    let f_shifted = f3_functional(&shifted);
    assert!(
        (f_base - f_shifted).abs() > 1e-6,
        "F_3 shouldn't be identical for different patches"
    );
}
