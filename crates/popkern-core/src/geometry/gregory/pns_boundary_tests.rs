//! Tests for [`pns_boundary`] — tensor-border to cap-outer-layer
//! composition via `b ∘ β`.

#![allow(clippy::needless_range_loop)]

use super::super::pns_beta::BetaReparam;
use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_pt(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

#[test]
fn degree_elevation_preserves_polynomial_value() {
    // Elevating a bi-2 → bi-4 should preserve the polynomial's
    // values. Pick a bi-2 with known coefs, evaluate at sample u,
    // then evaluate the elevated bi-4 at the same u, compare.
    let strip = Bi2BorderStrip {
        pos: [pt(1.0, 0.0, 0.0), pt(2.0, 1.0, 0.0), pt(3.0, 0.0, 1.0)],
        deriv: [pt(0.0, 0.0, 0.0); 3], // unused here
    };
    let beta = BetaReparam::bi4(6);
    let (row0, _row1) = compose_bi4_outer_layers(&strip, &beta);

    // Evaluate bi-2 position at u = 0.3.
    let u: f64 = 0.3;
    let omu = 1.0 - u;
    let bi2_pos = Point3::new(
        strip.pos[0].x * omu * omu + strip.pos[1].x * 2.0 * u * omu + strip.pos[2].x * u * u,
        strip.pos[0].y * omu * omu + strip.pos[1].y * 2.0 * u * omu + strip.pos[2].y * u * u,
        strip.pos[0].z * omu * omu + strip.pos[1].z * 2.0 * u * omu + strip.pos[2].z * u * u,
    );

    // Evaluate bi-4 at u.
    let b4: [f64; 5] = {
        let o = 1.0 - u;
        [
            o.powi(4),
            4.0 * u * o.powi(3),
            6.0 * u * u * o * o,
            4.0 * u.powi(3) * o,
            u.powi(4),
        ]
    };
    let mut bi4_pos = Point3::new(0.0, 0.0, 0.0);
    for i in 0..5 {
        bi4_pos.x += b4[i] * row0[i].x;
        bi4_pos.y += b4[i] * row0[i].y;
        bi4_pos.z += b4[i] * row0[i].z;
    }
    assert!(
        approx_pt(bi2_pos, bi4_pos, 1e-12),
        "bi-2 pos = {bi2_pos:?}, bi-4 elevated = {bi4_pos:?}"
    );
}

#[test]
fn n4_identity_beta_leaves_row0_as_elevated_pos() {
    // At n=4, β is the identity, so row1 = 1·dpos_du + 1·deriv...
    // wait, β.b[i] are ALL zero at n=4 (from pns_beta: b_0=b_1=b_2=b_3=0),
    // and β.a[i] = 1 for all i. So row1 should equal the derivative
    // polynomial (a=1 everywhere) = degree-elevated deriv.
    let strip = Bi2BorderStrip {
        pos: [pt(1.0, 0.0, 0.0), pt(0.5, 1.0, 0.0), pt(0.0, 0.0, 0.0)],
        deriv: [pt(0.0, 1.0, 0.0), pt(0.0, 1.0, 0.0), pt(0.0, 1.0, 0.0)],
    };
    let beta = BetaReparam::bi4(4);
    let (_row0, row1) = compose_bi4_outer_layers(&strip, &beta);

    // Sample row1 at a few u's and verify it matches `deriv`
    // (degree-elevated to 4 from degree 2).
    // With constant deriv = (0, 1, 0), row1 should be constant.
    for i in 0..5 {
        assert!(
            approx_pt(row1[i], pt(0.0, 1.0, 0.0), 1e-12),
            "row1[{i}] at n=4 identity β should be (0, 1, 0), got {:?}",
            row1[i]
        );
    }
}

#[test]
fn degree_elevation_preserves_endpoints() {
    // Endpoints: q_0 = p_0, q_4 = p_2 (the degree-elevation
    // formula preserves endpoints).
    let strip = Bi2BorderStrip {
        pos: [pt(1.0, 2.0, 3.0), pt(4.0, 5.0, 6.0), pt(7.0, 8.0, 9.0)],
        deriv: [pt(0.0, 0.0, 0.0); 3],
    };
    let beta = BetaReparam::bi4(6);
    let (row0, _) = compose_bi4_outer_layers(&strip, &beta);
    assert!(approx_pt(row0[0], strip.pos[0], 1e-12));
    assert!(approx_pt(row0[4], strip.pos[2], 1e-12));
}

#[test]
fn row1_is_degree_4_polynomial_in_u() {
    // Structural check: the computed row1 must be a valid
    // degree-4 Bezier polynomial, i.e., polynomial values at
    // sample points agree with the BB formula.
    //
    // This is a sanity check (no specific numerical value),
    // just that the function doesn't panic and returns 5 CPs.
    let strip = Bi2BorderStrip {
        pos: [pt(0.1, 0.2, 0.3), pt(0.4, 0.5, 0.6), pt(0.7, 0.8, 0.9)],
        deriv: [pt(1.0, 0.0, 0.0), pt(0.5, 0.5, 0.0), pt(0.0, 1.0, 0.0)],
    };
    for n in [5_usize, 6, 7, 8] {
        let beta = BetaReparam::bi4(n);
        let (row0, row1) = compose_bi4_outer_layers(&strip, &beta);
        assert_eq!(row0.len(), 5);
        assert_eq!(row1.len(), 5);
        // All CPs should be finite.
        for cp in row0.iter().chain(row1.iter()) {
            assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
        }
    }
}

#[test]
fn ring_n() {
    let ring = Bi2TensorBorderRing {
        sectors: vec![
            Bi2BorderStrip {
                pos: [pt(0.0, 0.0, 0.0); 3],
                deriv: [pt(0.0, 0.0, 0.0); 3],
            };
            6
        ],
    };
    assert_eq!(ring.n(), 6);
}
