//! Tests for the KP15 bi-6 tensor-border reparameterisation port.
//!
//! The authoritative numeric gate against the Maple-script oracle
//! is the full-cap fixture test in `kp15_bi6_deep_interior_tests`
//! (which exercises this border path end to end). These tests
//! pin the border stages' structural invariants that hold
//! independent of the oracle: affine invariance (constant nets map
//! to constant borders) and the row-0 position-row endpoints.

use super::*;
use crate::geometry::Point3;

const EPS: f64 = 1e-12;

fn const_net(v: f64) -> [[Point3; 4]; 4] {
    [[Point3::new(v, 0.0, 0.0); 4]; 4]
}

#[test]
fn btobb_preserves_constants() {
    // BtoBB is a convex reparameterisation of the net; a constant
    // net must map to a constant bi-3 border (weights sum to 1).
    let b = bspline_net_to_bi3_bb(&const_net(2.5));
    for row in b.iter() {
        for cp in row.iter() {
            assert!(
                (cp.x - 2.5).abs() < EPS,
                "BtoBB not affine-invariant: {}",
                cp.x
            );
        }
    }
}

#[test]
fn reparam_row_preserves_constants() {
    // The β tensor-border reparameterisation preserves position,
    // so a constant field maps to constant CPs on all three depth
    // rows — including row 2, whose d2/c/τ terms must cancel.
    let val = -1.75;
    let rows = reparam_tensor_border_row(|_, _| Point3::new(val, 0.0, 0.0), 0.5, 0.0891, 0.85);
    for (j, row) in rows.iter().enumerate() {
        for (i, cp) in row.iter().enumerate() {
            assert!(
                (cp.x - val).abs() < EPS,
                "reparam row {j}[{i}] = {} not affine-invariant",
                cp.x
            );
        }
    }
}

#[test]
fn reparam_row0_endpoints_interpolate_edge() {
    // Row 0 (position row) endpoints are the bi-3 edge endpoints:
    // rows[0][0] = p(0,0), rows[0][6] = p(3,0).
    let rows = reparam_tensor_border_row(
        |i, j| Point3::new((10 * i + j) as f64, 0.0, 0.0),
        0.5,
        0.09,
        0.85,
    );
    assert!((rows[0][0].x - 0.0).abs() < EPS);
    assert!((rows[0][6].x - 30.0).abs() < EPS);
}
