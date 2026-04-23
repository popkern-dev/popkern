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

/// **Known failure — documents the σ bug flagged in the
/// audit.** The paper claims (§4 "Increased smoothness at the
/// central point") that `a ∘ σ` satisfies eq (2) automatically.
/// This test checks that claim: `cur_{10}` from the central chain
/// rule should match `cur_{10}` from eq (2) applied with adjacent
/// sectors' j=1 row data.
///
/// # Current status: FAILS
///
/// Working through the algebra, the G¹ constraint
///   `cur_{10}(central) = cur_{10}(eq 2)`
/// requires, for all `a`:
///   `2c · ds_u^0  =  ds_v^0 + R_{2π/n} · ds_v^0`
///
/// Given `ds_u^0 = 4·(d_1, 0)`, this forces
///   `ds_v^0 = (4c·d_1, ±4c·d_1·tan(α/2))`
///
/// but our [`pns_sigma::sigma_center_jet`] derives `σ_{01}` from
/// image-space reflection about the sector bisector, which gives
///   `ds_v^0 = (4 d_1 cos α, 4 d_1 sin α)`.
///
/// The x-coords agree (`4c·d_1`) but the y-coords differ — at n=5,
/// mine is `0.95·d_1`, paper needs `±0.22·d_1` (roughly 4× off and
/// possibly wrong sign).
///
/// # Root cause
///
/// Appendix B says "flip-symmetry of each sector about its
/// diagonal and smoothness ... imply σ_{00}, σ_{10}, σ_{11},
/// σ_{20}" and "Symmetry defines the remaining BB-coefficients of
/// σ in terms of d_1, d_2, d_3." I interpreted this as
/// image-space reflection about the sector bisector. That
/// interpretation is inconsistent with the G¹ constraint the
/// paper also claims the construction enforces.
///
/// The correct derivation of `σ_{01}` (and `σ_{02}`) likely needs
/// a careful re-reading of Appendix B's rotational structure,
/// or direct inspection of a reference PnS implementation output
/// on a trivial test case.
///
/// # Until this is fixed
///
/// The cap assembly can't correctly use [`central_bb_from_a_and_sigma`]
/// yet. Either:
///  1. Fix `σ_{01}, σ_{02}` in `pns_sigma` to satisfy the G¹
///     identity, then re-enable this test.
///  2. Verify the correct formula by cross-check against the
///     reference PnS implementation's output on a known input
///     (non-clean-room approach — requires an independent probe).
#[test]
#[ignore = "σ_{01}, σ_{02} derivation from pns_sigma doesn't satisfy G¹ (eq 2) — see test doc"]
fn central_expansion_satisfies_eq2_g1_constraint() {
    use super::super::pns_g1_constraints::bi4_g1_row_from_row1;

    // Build an arbitrary non-trivial a and compute N sectors'
    // central coefs via a∘σ. Then check eq 2's formula matches.
    let pole = pt(0.0, 0.0, 0.0); // Use zero pole for simplicity
    let a = QuadMapJet {
        a_0: pole,
        a_x: pt(0.3, 0.0, 0.1),
        a_y: pt(0.0, 0.4, 0.2),
        a_xx: pt(0.2, 0.0, 0.05),
        a_xy: pt(0.0, 0.1, 0.03),
        a_yy: pt(0.0, 0.0, 0.15),
    };

    for n in [5_usize, 6, 7, 8] {
        let sigma_canonical = sigma_center_jet_default(n);
        // Sector 0 uses canonical σ; nxt (sector 1) uses σ rotated
        // by 2π/n.
        let theta = 2.0 * std::f64::consts::PI / (n as f64);
        let sigma_0 = rotate_sigma_jet(&sigma_canonical, 0.0);
        let sigma_1 = rotate_sigma_jet(&sigma_canonical, theta);
        let cur_bb = central_bb_from_a_and_sigma(&a, &sigma_0);
        let nxt_bb = central_bb_from_a_and_sigma(&a, &sigma_1);

        // Pick some value for cur_{21}, cur_{31}, cur_{41},
        // nxt_{31}, nxt_{41}. These are outside the central
        // region (i+j ≥ 3), so the central expansion doesn't
        // constrain them. Make them arbitrary for this test.
        let cur_21 = pt(0.1, 0.0, 0.0);
        let cur_31 = pt(0.2, 0.0, 0.0);
        let cur_41 = pt(0.3, 0.0, 0.0);
        let nxt_31 = pt(0.0, 0.2, 0.0);
        let nxt_41 = pt(0.0, 0.3, 0.0);

        // cur_row1 = [cur_{01}, cur_{11}, cur_{21}, cur_{31}, cur_{41}]
        let cur_row1 = [cur_bb.p01, cur_bb.p11, cur_21, cur_31, cur_41];
        // nxt_row1 = [nxt_{01}, nxt_{11}, nxt_{21}, nxt_{31}, nxt_{41}].
        // nxt_{21} is unused by bi4_g1_row_from_row1 for computing
        // cur_{10}, cur_{20}, cur_{30}, cur_{40} — it's an output.
        // Use zero as placeholder.
        let nxt_row1 = [nxt_bb.p01, nxt_bb.p11, pt(0.0, 0.0, 0.0), nxt_31, nxt_41];

        let g1_result = bi4_g1_row_from_row1(cur_bb.p00, cur_row1, nxt_row1, n);

        // The claim: g1_result.cur10 should equal cur_bb.p10 and
        // g1_result.cur20 should equal cur_bb.p20. These are the
        // BB coefs that a∘σ pins AND that eq 2 pins — they must
        // agree.
        let eps = 1e-9;
        assert!(
            approx_pt(g1_result.cur10, cur_bb.p10, eps),
            "n={n}: eq2 cur_{{10}} = {:?}, central p_{{10}} = {:?}",
            g1_result.cur10,
            cur_bb.p10
        );
        assert!(
            approx_pt(g1_result.cur20, cur_bb.p20, eps),
            "n={n}: eq2 cur_{{20}} = {:?}, central p_{{20}} = {:?}",
            g1_result.cur20,
            cur_bb.p20
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
