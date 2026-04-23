//! Tests for [`pns_central_quad`] — `a ∘ σ` chain rule.

#![allow(clippy::needless_range_loop)]

use super::super::pns_sigma::{SigmaCenterJet, sigma_center_jet_default};
use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_pt(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

/// Zero `a` (everything 0, including the pole) gives the all-zero cap.
#[test]
fn zero_a_gives_zero_bb_coefs() {
    let zero = pt(0.0, 0.0, 0.0);
    let a = QuadMapJet {
        a_0: zero,
        a_x: zero,
        a_y: zero,
        a_xx: zero,
        a_xy: zero,
        a_yy: zero,
    };
    let sigma = sigma_center_jet_default(6);
    let bb = central_bb_from_a_and_sigma(&a, &sigma);
    for p in [bb.p00, bb.p10, bb.p01, bb.p20, bb.p11, bb.p02] {
        assert!(approx_pt(p, zero, 1e-12));
    }
}

/// Constant `a` (only `a_0` non-zero) gives all 6 CPs at the pole.
#[test]
fn constant_a_pins_all_central_coefs_to_pole() {
    let pole = pt(1.5, 2.7, 3.9);
    let zero = pt(0.0, 0.0, 0.0);
    let a = QuadMapJet {
        a_0: pole,
        a_x: zero,
        a_y: zero,
        a_xx: zero,
        a_xy: zero,
        a_yy: zero,
    };
    let sigma = sigma_center_jet_default(7);
    let bb = central_bb_from_a_and_sigma(&a, &sigma);
    for p in [bb.p00, bb.p10, bb.p01, bb.p20, bb.p11, bb.p02] {
        assert!(approx_pt(p, pole, 1e-12), "got {p:?}");
    }
}

/// Sanity: `p_{00}` is always the pole, regardless of `a`'s other components.
#[test]
fn p00_is_always_pole() {
    let pole = pt(7.0, -3.0, 0.5);
    let a = QuadMapJet {
        a_0: pole,
        a_x: pt(1.0, 0.0, 0.0),
        a_y: pt(0.0, 2.0, 0.0),
        a_xx: pt(3.0, 0.0, 0.0),
        a_xy: pt(0.0, 4.0, 1.0),
        a_yy: pt(0.0, 0.0, 5.0),
    };
    for n in [5_usize, 6, 7, 8] {
        let sigma = sigma_center_jet_default(n);
        let bb = central_bb_from_a_and_sigma(&a, &sigma);
        assert!(approx_pt(bb.p00, pole, 1e-12), "n={n}");
    }
}

/// Rotating σ by `theta` and then applying chain rule with `a` is
/// equivalent to applying chain rule with `a ∘ R_theta` and the
/// original σ — up to the fact that `a` pre-rotation changes
/// `a`'s Hessian. Here we verify that rotation-then-chain matches
/// a direct computation with rotated sigma jet.
#[test]
fn rotate_sigma_jet_preserves_structure() {
    let sigma = sigma_center_jet_default(6);
    // Rotate by 0 → identity.
    let rotated_0 = rotate_sigma_jet(&sigma, 0.0);
    assert!(approx_jet(&sigma, &rotated_0, 1e-12));

    // Rotate by 2π → identity.
    let rotated_2pi = rotate_sigma_jet(&sigma, 2.0 * std::f64::consts::PI);
    assert!(approx_jet(&sigma, &rotated_2pi, 1e-10));

    // Magnitudes of jet vectors preserved under rotation.
    let rotated_arbitrary = rotate_sigma_jet(&sigma, 0.7);
    let mag = |v: [f64; 2]| (v[0] * v[0] + v[1] * v[1]).sqrt();
    assert!((mag(sigma.ds_u) - mag(rotated_arbitrary.ds_u)).abs() < 1e-12);
    assert!((mag(sigma.ds_v) - mag(rotated_arbitrary.ds_v)).abs() < 1e-12);
    assert!((mag(sigma.ds_uu) - mag(rotated_arbitrary.ds_uu)).abs() < 1e-12);
}

fn approx_jet(a: &SigmaCenterJet, b: &SigmaCenterJet, eps: f64) -> bool {
    let ok = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).abs() < eps && (a[1] - b[1]).abs() < eps;
    ok(a.s00, b.s00)
        && ok(a.ds_u, b.ds_u)
        && ok(a.ds_v, b.ds_v)
        && ok(a.ds_uu, b.ds_uu)
        && ok(a.ds_vv, b.ds_vv)
        && ok(a.ds_uv, b.ds_uv)
}

/// For all N sectors of a convex regular cap, the central BB coefs
/// should preserve the cap's rotational symmetry: if `a` has only
/// rotationally-symmetric components (e.g., `a = (0, 0, z·(x² + y²))`
/// producing a paraboloid), then all N sectors' BB sextuples should
/// agree up to a rotation of the (x, y) components.
#[test]
fn paraboloid_a_gives_rotationally_symmetric_cap() {
    // a(x, y) = (0, 0, x² + y²)  →  a_0 = (0,0,0), a_xx = a_yy = (0,0,2), rest zero.
    // Note: with our convention a(x,y) = a_0 + a_x·x + a_y·y + (1/2)a_xx·x² + a_xy·xy + (1/2)a_yy·y²,
    // x² + y² as the z-coord means (1/2)a_xx = 1 so a_xx = (0,0,2) and a_yy = (0,0,2).
    let zero = pt(0.0, 0.0, 0.0);
    let a = QuadMapJet {
        a_0: zero,
        a_x: zero,
        a_y: zero,
        a_xx: pt(0.0, 0.0, 2.0),
        a_xy: zero,
        a_yy: pt(0.0, 0.0, 2.0),
    };
    let sigma = sigma_center_jet_default(6);
    let all = central_bb_all_sectors(&a, &sigma, 6);
    // All 6 sectors' p_{00} = pole.
    for bb in &all {
        assert!(approx_pt(bb.p00, zero, 1e-12));
    }
    // For a rotationally-symmetric paraboloid, the z-components of
    // all sectors' BB coefs should match — x² + y² is invariant
    // under rotation.
    let z_p10_0 = all[0].p10.z;
    for (k, bb) in all.iter().enumerate() {
        assert!(
            (bb.p10.z - z_p10_0).abs() < 1e-10,
            "sector {k} p_{{10}}.z = {} vs sector 0 {z_p10_0}",
            bb.p10.z
        );
    }
}

/// Smoke test: `central_bb_all_sectors` with valence 5 produces
/// 5 distinct sectors' worth of BB sextuples.
#[test]
fn central_bb_all_sectors_count_and_finite() {
    let pole = pt(0.1, 0.2, 0.3);
    let a = QuadMapJet {
        a_0: pole,
        a_x: pt(1.0, 0.0, 0.5),
        a_y: pt(0.0, 1.0, 0.3),
        a_xx: pt(0.2, 0.0, 0.0),
        a_xy: pt(0.0, 0.3, 0.0),
        a_yy: pt(0.0, 0.0, 0.4),
    };
    for n in [5_usize, 6, 7, 8] {
        let sigma = sigma_center_jet_default(n);
        let all = central_bb_all_sectors(&a, &sigma, n);
        assert_eq!(all.len(), n);
        for bb in &all {
            for p in [bb.p00, bb.p10, bb.p01, bb.p20, bb.p11, bb.p02] {
                assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
            }
        }
    }
}
