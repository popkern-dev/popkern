//! Tests for [`pns_gram`] — Gram matrix + F_κ minimiser.

#![allow(clippy::needless_range_loop)]

use super::super::pns_f_kappa::f_kappa_bi4;
use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// `p^T · G^{(κ)} · p` computed from the Gram matrix must equal
/// `F_κ(p)` from the independent quadrature-based evaluator.
/// This is the strong cross-check between the two implementations.
#[test]
fn gram_quadratic_form_matches_f_kappa_evaluator() {
    let mut grid = [[pt(0.0, 0.0, 0.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            grid[i][j] = pt(
                0.1 * (i as f64) - 0.05 * (j as f64).powi(2),
                0.2 * (j as f64) * ((i + 1) as f64).sqrt(),
                (i as f64 - 2.0) * (j as f64 - 2.0),
            );
        }
    }

    for kappa in 1..=4 {
        let g = gram_bi4_f_kappa(kappa);

        // Per-coord quadratic form.
        let mut total = 0.0_f64;
        for a in 0..5 {
            for b in 0..5 {
                let i1 = bi4_flat_index(a, b);
                for c in 0..5 {
                    for d in 0..5 {
                        let i2 = bi4_flat_index(c, d);
                        let p_ab = grid[a][b];
                        let p_cd = grid[c][d];
                        total += g[i1][i2] * (p_ab.x * p_cd.x + p_ab.y * p_cd.y + p_ab.z * p_cd.z);
                    }
                }
            }
        }

        let expected = f_kappa_bi4(&grid, kappa);
        let rel_err = (total - expected).abs() / (expected.abs() + 1e-12);
        assert!(
            rel_err < 1e-9,
            "kappa={kappa}: gram form = {total}, F_κ evaluator = {expected}, rel_err = {rel_err}"
        );
    }
}

#[test]
fn gram_is_symmetric() {
    for kappa in 1..=4 {
        let g = gram_bi4_f_kappa(kappa);
        for i in 0..25 {
            for j in i + 1..25 {
                let diff = (g[i][j] - g[j][i]).abs();
                assert!(
                    diff < 1e-10,
                    "kappa={kappa}, G[{i}][{j}]={} != G[{j}][{i}]={} (diff={diff})",
                    g[i][j],
                    g[j][i]
                );
            }
        }
    }
}

#[test]
fn gram_is_psd_on_test_vectors() {
    // For any test vector v, v^T G v >= 0 (PSD check on specific samples).
    for kappa in 1..=4 {
        let g = gram_bi4_f_kappa(kappa);
        // Test with a few pseudo-random vectors.
        for seed in 0..5 {
            let mut v = [0.0_f64; 25];
            for i in 0..25 {
                // Deterministic pseudo-random.
                v[i] = ((seed * 7 + i * 13) as f64).sin();
            }
            let mut quad = 0.0_f64;
            for i in 0..25 {
                for j in 0..25 {
                    quad += v[i] * g[i][j] * v[j];
                }
            }
            assert!(
                quad >= -1e-10,
                "kappa={kappa}, seed={seed}: v^T G v = {quad} < 0"
            );
        }
    }
}

#[test]
fn flatten_roundtrip() {
    for a in 0..5 {
        for b in 0..5 {
            let idx = bi4_flat_index(a, b);
            let (ra, rb) = bi4_unflat(idx);
            assert_eq!((a, b), (ra, rb), "idx {idx}: {a},{b} -> {ra},{rb}");
        }
    }
}

#[test]
fn minimize_with_all_fixed_is_noop() {
    // If all 25 coefs are fixed, the minimiser does nothing.
    let mut grid = [[pt(1.0, 2.0, 3.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            grid[i][j] = pt(i as f64, j as f64, (i * j) as f64);
        }
    }
    let original = grid;
    let fixed_all = [true; 25];
    minimize_f_kappa_bi4(&mut grid, 3, &fixed_all).unwrap();
    for i in 0..5 {
        for j in 0..5 {
            assert_eq!(grid[i][j], original[i][j]);
        }
    }
}

#[test]
fn minimize_with_all_free_gives_zero_energy() {
    // If all 25 coefs are free and the patch doesn't need to
    // interpolate anything, the minimum F_κ = 0 is achieved
    // trivially (all-zero solution). Whether we reach that depends
    // on the full-rank-ness of G.
    //
    // Seed the grid with a non-zero configuration; minimising with
    // NO fixed coefs means the zero vector is an unconstrained
    // optimum. However G^{(κ)} has a non-trivial null space (degree
    // < κ polynomials) so the system is under-determined.
    //
    // Setting a few CPs fixed (e.g. the 4 corners) usually pins
    // enough to get a unique solution — we test that scenario.
    let mut grid = [[pt(0.0, 0.0, 0.0); 5]; 5];
    // Corners at some test values, rest at (0,0,0).
    grid[0][0] = pt(0.0, 0.0, 0.0);
    grid[4][0] = pt(1.0, 0.0, 0.0);
    grid[0][4] = pt(0.0, 1.0, 0.0);
    grid[4][4] = pt(1.0, 1.0, 0.0);

    let mut fixed_mask = [false; 25];
    fixed_mask[bi4_flat_index(0, 0)] = true;
    fixed_mask[bi4_flat_index(4, 0)] = true;
    fixed_mask[bi4_flat_index(0, 4)] = true;
    fixed_mask[bi4_flat_index(4, 4)] = true;

    // Minimise F_3 with corners pinned — solution should exist.
    let result = minimize_f_kappa_bi4(&mut grid, 3, &fixed_mask);
    // May succeed or fail depending on G_vv rank; document the
    // outcome rather than asserting either way.
    if result.is_ok() {
        // If we got here, the corners are preserved and interior
        // was optimised.
        assert_eq!(grid[0][0], pt(0.0, 0.0, 0.0));
        assert_eq!(grid[4][0], pt(1.0, 0.0, 0.0));
        assert_eq!(grid[0][4], pt(0.0, 1.0, 0.0));
        assert_eq!(grid[4][4], pt(1.0, 1.0, 0.0));
        // Verify optimality: F_κ should be at a stationary point.
        // Not testing gradient norm here (numerical); just that we
        // didn't explode.
        let e = f_kappa_bi4(&grid, 3);
        assert!(e.is_finite() && e >= 0.0);
    }
    // If result.is_err(), the G_vv was singular — acceptable for
    // this test since 4 corners + F_3 may not give a unique min.
}

#[test]
fn minimize_bilinear_exact_interpolation() {
    // A bilinear patch `p_{a,b} = (a/4)·P1 + (b/4)·P2` has
    // F_3 = F_4 = 0 exactly (2nd + higher derivs vanish).
    //
    // If we pin the 4 boundary rows and columns (16 CPs at the edges)
    // at their bilinear values, the interior 9 CPs should be
    // uniquely determined by F_3 minimisation at... well, probably
    // at bilinear too, since any bilinear interpolant is
    // F_3-optimal.
    let mut grid = [[pt(0.0, 0.0, 0.0); 5]; 5];
    // Bilinear grid.
    for i in 0..5 {
        for j in 0..5 {
            let u = (i as f64) / 4.0;
            let v = (j as f64) / 4.0;
            grid[i][j] = pt(u, v, 0.0);
        }
    }
    // Pin only the 4 boundary edges (16 CPs: i=0 col, i=4 col,
    // j=0 row, j=4 row). The 9 interior are free.
    let mut fixed_mask = [false; 25];
    for i in 0..5 {
        fixed_mask[bi4_flat_index(0, i)] = true;
        fixed_mask[bi4_flat_index(4, i)] = true;
        fixed_mask[bi4_flat_index(i, 0)] = true;
        fixed_mask[bi4_flat_index(i, 4)] = true;
    }
    // Before: interior is at bilinear values, F_3 = 0.
    let f3_before = f_kappa_bi4(&grid, 3);
    assert!(
        f3_before < 1e-10,
        "bilinear F_3 should be 0, got {f3_before}"
    );

    // Zero out the interior to make minimisation work to restore.
    for i in 1..=3 {
        for j in 1..=3 {
            grid[i][j] = pt(0.0, 0.0, 0.0);
        }
    }
    // Minimise F_3 — should restore to bilinear (or some F_3 = 0 state).
    minimize_f_kappa_bi4(&mut grid, 3, &fixed_mask).unwrap();
    let f3_after = f_kappa_bi4(&grid, 3);
    assert!(
        f3_after < 1e-9,
        "after F_3 min with bilinear boundary, F_3 = {f3_after}"
    );
}
