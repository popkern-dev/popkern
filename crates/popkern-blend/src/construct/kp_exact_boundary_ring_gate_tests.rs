//! PK-283 wiring tests: the exact-boundary optimizer's entry/accept
//! feasibility must decline rings whose patches cross EACH OTHER even
//! when every patch is individually fold-free — the class the
//! per-patch gates cannot see.

use super::*;
use popkern_core::geometry::Point3;

/// Synthetic pyramid cap ring obeying the promotion ladder's spoke
/// convention (spoke `k`: patch `k` UOne ↔ patch `k+1` VZero reversed;
/// apex at every patch's `(1, 0)` corner). Bilinear sheets are exactly
/// represented by degree-9 Bernstein grids sampled at `i/9`.
fn pyramid_ring(n: usize) -> Vec<RectangularGregoryBiNonic> {
    let apex = Point3::new(0.0, 0.0, 1.0);
    let outer: Vec<Point3> = (0..n)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    (0..n)
        .map(|k| {
            let q_prev = outer[(k + n - 1) % n];
            let q_next = outer[k];
            let far = Point3::new(q_prev.x + q_next.x, q_prev.y + q_next.y, 0.0);
            let corner = |u: f64, v: f64| -> Point3 {
                let w00 = (1.0 - u) * (1.0 - v);
                let w01 = (1.0 - u) * v;
                let w10 = u * (1.0 - v);
                let w11 = u * v;
                Point3::new(
                    w00 * q_prev.x + w01 * far.x + w10 * apex.x + w11 * q_next.x,
                    w00 * q_prev.y + w01 * far.y + w10 * apex.y + w11 * q_next.y,
                    w00 * q_prev.z + w01 * far.z + w10 * apex.z + w11 * q_next.z,
                )
            };
            let mut patch = RectangularGregoryBiNonic::zeros();
            for i in 0..10 {
                for j in 0..10 {
                    patch.grid[i][j] = corner(i as f64 / 9.0, j as f64 / 9.0);
                }
            }
            patch
        })
        .collect()
}

/// A per-patch-clean ring whose patch 0 has been translated into patch
/// 1's sector — the deterministic stand-in for the PK-275 pyramid4
/// escape (the reverted analytic-Jacobian trajectory moved exactly this
/// DOF class, the spoke rows, into a neighbouring patch).
fn crossing_ring() -> Vec<RectangularGregoryBiNonic> {
    let mut ring = pyramid_ring(4);
    for row in ring[0].grid.iter_mut() {
        for p in row.iter_mut() {
            p.y += 0.4;
        }
    }
    ring
}

#[test]
fn ring_gate_declines_what_the_per_patch_gate_cannot_see() {
    let ring = crossing_ring();
    // The per-patch fold gate is blind to this ring…
    assert!(
        !any_self_intersects(&ring),
        "fixture must be per-patch clean to exercise the ring gate"
    );
    // …the ring gates are not: the crossing is macro-scale, so both
    // the coarse feasible-start screen and the full accept semantics
    // decline it.
    assert!(
        !feasible_start(&ring),
        "feasible-start screen must decline a macro-scale crossing"
    );
    assert!(
        !acceptable_ring(&ring),
        "accept gate must decline the pairwise crossing"
    );
    // Clean control: the same fixture without the translation passes.
    assert!(feasible_start(&pyramid_ring(4)));
    assert!(acceptable_ring(&pyramid_ring(4)));
}

#[test]
fn optimizer_entry_declines_a_crossing_ring() {
    // End-to-end through the optimizer entry: with per-patch
    // self-intersection proven absent (above), a `None` here is the
    // ring-crossing gate declining the infeasible start.
    assert!(optimize_exact_boundary_binonic(&crossing_ring()).is_none());
}
