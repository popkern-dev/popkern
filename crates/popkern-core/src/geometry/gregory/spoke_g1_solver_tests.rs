//! Tests for [`super::spoke_g1_solver::spoke_g1_solve`].

use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
use super::spoke_g1_solver::{Side, Spoke, SpokeEnd, spoke_g1_solve};
use crate::geometry::{Point3, Vec3};

/// Build a flat biquintic-full patch over the rectangle
/// `[x0, x1] × [0, 1]` in the `z = 0` plane.
fn flat_patch(x0: f64, x1: f64) -> RectangularGregoryBiquinticFull {
    let pt =
        |i: usize, j: usize| Point3::new(x0 + (x1 - x0) * (i as f64) / 5.0, (j as f64) / 5.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    let corner_inner_cps = [
        (pt(1, 1), pt(1, 1)),
        (pt(1, 4), pt(1, 4)),
        (pt(4, 4), pt(4, 4)),
        (pt(4, 1), pt(4, 1)),
    ];
    RectangularGregoryBiquinticFull {
        grid,
        corner_inner_cps,
    }
}

/// Residual of the linearised G1 condition at a sample `t` on a
/// shared spoke, using the SAME normal-target formulation as the
/// solver (average along-spoke tangent × patch-A cross-tangent).
fn g1_residual(
    patch_a: &RectangularGregoryBiquinticFull,
    side_a: Side,
    patch_b: &RectangularGregoryBiquinticFull,
    side_b: Side,
    reversed: bool,
    t: f64,
) -> f64 {
    let t_b = if reversed { 1.0 - t } else { t };
    let xtan_a = cross_tangent(patch_a, side_a, t);
    let xtan_b = cross_tangent(patch_b, side_b, t_b);
    let along_a = spoke_tangent(patch_a, side_a, t);
    let along_b = spoke_tangent(patch_b, side_b, t_b);
    let along = Vec3::new(
        0.5 * (along_a.x + along_b.x),
        0.5 * (along_a.y + along_b.y),
        0.5 * (along_a.z + along_b.z),
    );
    let n_unnorm = along.cross(&xtan_a);
    let n = n_unnorm.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    (xtan_a.x - xtan_b.x) * n.x + (xtan_a.y - xtan_b.y) * n.y + (xtan_a.z - xtan_b.z) * n.z
}

fn cross_tangent(patch: &RectangularGregoryBiquinticFull, side: Side, t: f64) -> Vec3 {
    match side {
        Side::UZero => patch.partial_u(0.0, t),
        Side::UOne => patch.partial_u(1.0, t),
        Side::VZero => patch.partial_v(t, 0.0),
        Side::VOne => patch.partial_v(t, 1.0),
    }
}

fn spoke_tangent(patch: &RectangularGregoryBiquinticFull, side: Side, t: f64) -> Vec3 {
    match side {
        Side::UZero => patch.partial_v(0.0, t),
        Side::UOne => patch.partial_v(1.0, t),
        Side::VZero => patch.partial_u(t, 0.0),
        Side::VOne => patch.partial_u(t, 1.0),
    }
}

/// Baseline: two already-G1 flat patches share a spoke. The
/// solver should leave them (essentially) unchanged — no spurious
/// perturbation of the free interior CPs.
#[test]
fn no_op_on_already_g1_pair() {
    let mut patches = vec![flat_patch(0.0, 1.0), flat_patch(1.0, 2.0)];
    let before = patches[0].grid[4][2];
    let spokes = [Spoke {
        a: SpokeEnd {
            patch_idx: 0,
            side: Side::UOne,
            reversed: false,
        },
        b: SpokeEnd {
            patch_idx: 1,
            side: Side::UZero,
            reversed: false,
        },
    }];
    let improved = spoke_g1_solve(&mut patches, &spokes, 3).unwrap();
    let after = patches[0].grid[4][2];
    let drift = ((after.x - before.x).powi(2)
        + (after.y - before.y).powi(2)
        + (after.z - before.z).powi(2))
    .sqrt();
    assert!(
        drift < 1e-10,
        "solver perturbed already-G1 patch by {drift:.4e}"
    );
    assert!(improved >= 1, "already-G1 samples must register as met");
}

/// Core test: patch B's free interior CP `(1, 2)` is lifted in
/// `z` so the initial cross-tangent at the shared `u = 0`
/// boundary has a pure `B_2^5(v)` z-component that lies within
/// the solver's reachable span (the free interior CPs at `(1, 2)
/// / (1, 3) / (4, 2) / (4, 3)` span exactly the `B_2^5 / B_3^5`
/// quintic basis). After one pass, the G1 residual must drop
/// to machine precision.
#[test]
fn in_span_mismatch_drops_to_machine_precision() {
    let mut patches = vec![flat_patch(0.0, 1.0), flat_patch(1.0, 2.0)];
    patches[1].grid[1][2].z += 0.2;
    let spokes = [Spoke {
        a: SpokeEnd {
            patch_idx: 0,
            side: Side::UOne,
            reversed: false,
        },
        b: SpokeEnd {
            patch_idx: 1,
            side: Side::UZero,
            reversed: false,
        },
    }];
    let ts = [0.25_f64, 0.5, 0.75];
    let max_before = ts
        .iter()
        .map(|&t| g1_residual(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_before > 0.1,
        "test scaffolding: expected nonzero initial G1 residual, got {max_before:.4e}"
    );
    let _ = spoke_g1_solve(&mut patches, &spokes, 3).expect("solve should not be singular");
    let max_after = ts
        .iter()
        .map(|&t| g1_residual(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_after < 1e-6,
        "post-solve G1 residual {max_after:.4e} exceeds 1e-6 (was {max_before:.4e})"
    );
}

/// Out-of-span test: the mismatch comes from lifting patch B's
/// corner inner-CP pairs, which injects a `B_1^5 + B_4^5` z-
/// pattern that is NOT in the `span(B_2^5, B_3^5)` that the
/// free interior CPs can reach. The solver minimises the
/// residual in least-squares sense; we only assert that it
/// reduces the residual (does not diverge), documenting the
/// structural DOF limitation that motivates the degree-7
/// or two-stripe extension.
#[test]
fn out_of_span_mismatch_reduces_but_does_not_fully_cancel() {
    let mut patches = vec![flat_patch(0.0, 1.0), flat_patch(1.0, 2.0)];
    for pair in patches[1].corner_inner_cps.iter_mut() {
        pair.0.z += 0.15;
        pair.1.z += 0.15;
    }
    let spokes = [Spoke {
        a: SpokeEnd {
            patch_idx: 0,
            side: Side::UOne,
            reversed: false,
        },
        b: SpokeEnd {
            patch_idx: 1,
            side: Side::UZero,
            reversed: false,
        },
    }];
    let ts = [0.25_f64, 0.5, 0.75];
    let max_before = ts
        .iter()
        .map(|&t| g1_residual(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t).abs())
        .fold(0.0_f64, f64::max);
    let _ = spoke_g1_solve(&mut patches, &spokes, 3).expect("solve should not be singular");
    let max_after = ts
        .iter()
        .map(|&t| g1_residual(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_after < max_before,
        "solver must reduce residual even on out-of-span mismatch: before {max_before:.4e} after {max_after:.4e}"
    );
}

/// Free interior CPs that are NOT adjacent to any shared spoke
/// (the 4 "central" positions `(2,2), (2,3), (3,2), (3,3)`) must
/// not be modified by the solver.
#[test]
fn central_interior_cps_are_preserved() {
    let mut patches = vec![flat_patch(0.0, 1.0), flat_patch(1.0, 2.0)];
    for pair in patches[1].corner_inner_cps.iter_mut() {
        pair.0.z += 0.1;
        pair.1.z += 0.1;
    }
    let central_before_b = [
        patches[1].grid[2][2],
        patches[1].grid[2][3],
        patches[1].grid[3][2],
        patches[1].grid[3][3],
    ];
    let spokes = [Spoke {
        a: SpokeEnd {
            patch_idx: 0,
            side: Side::UOne,
            reversed: false,
        },
        b: SpokeEnd {
            patch_idx: 1,
            side: Side::UZero,
            reversed: false,
        },
    }];
    let _ = spoke_g1_solve(&mut patches, &spokes, 3);
    let central_after_b = [
        patches[1].grid[2][2],
        patches[1].grid[2][3],
        patches[1].grid[3][2],
        patches[1].grid[3][3],
    ];
    for (b, a) in central_before_b.iter().zip(central_after_b.iter()) {
        let d = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2) + (b.z - a.z).powi(2)).sqrt();
        assert!(d < 1e-14, "central CP drifted by {d:.4e}");
    }
}
