//! Biseptic (degree-7) and Peters two-stripe spoke-G1 solver
//! tests. Extracted from `spoke_g1_solver_tests.rs` to keep
//! each test file under the 500-line structural cap.

use super::spoke_g1_solver::{Side, Spoke, SpokeEnd, spoke_g1_solve};
use crate::geometry::{Point3, Vec3};

/// Biseptic analogue of `flat_patch` — degree-7 patches
/// spanning `[x0, x1] × [0, 1]` in the `z = 0` plane.
fn flat_biseptic_patch(
    x0: f64,
    x1: f64,
) -> super::rectangular_biseptic::RectangularGregoryBiseptic {
    use super::rectangular_biseptic::RectangularGregoryBiseptic;
    let pt =
        |i: usize, j: usize| Point3::new(x0 + (x1 - x0) * (i as f64) / 7.0, (j as f64) / 7.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    let corner_inner_cps = [
        (pt(1, 1), pt(1, 1)),
        (pt(1, 6), pt(1, 6)),
        (pt(6, 6), pt(6, 6)),
        (pt(6, 1), pt(6, 1)),
    ];
    RectangularGregoryBiseptic {
        grid,
        corner_inner_cps,
    }
}

fn g1_residual_biseptic(
    patch_a: &super::rectangular_biseptic::RectangularGregoryBiseptic,
    side_a: Side,
    patch_b: &super::rectangular_biseptic::RectangularGregoryBiseptic,
    side_b: Side,
    reversed: bool,
    t: f64,
) -> f64 {
    let t_b = if reversed { 1.0 - t } else { t };
    let xtan_a = match side_a {
        Side::UZero => patch_a.partial_u(0.0, t),
        Side::UOne => patch_a.partial_u(1.0, t),
        Side::VZero => patch_a.partial_v(t, 0.0),
        Side::VOne => patch_a.partial_v(t, 1.0),
    };
    let xtan_b = match side_b {
        Side::UZero => patch_b.partial_u(0.0, t_b),
        Side::UOne => patch_b.partial_u(1.0, t_b),
        Side::VZero => patch_b.partial_v(t_b, 0.0),
        Side::VOne => patch_b.partial_v(t_b, 1.0),
    };
    let along_a = match side_a {
        Side::UZero => patch_a.partial_v(0.0, t),
        Side::UOne => patch_a.partial_v(1.0, t),
        Side::VZero => patch_a.partial_u(t, 0.0),
        Side::VOne => patch_a.partial_u(t, 1.0),
    };
    let along_b = match side_b {
        Side::UZero => patch_b.partial_v(0.0, t_b),
        Side::UOne => patch_b.partial_v(1.0, t_b),
        Side::VZero => patch_b.partial_u(t_b, 0.0),
        Side::VOne => patch_b.partial_u(t_b, 1.0),
    };
    let along = Vec3::new(
        0.5 * (along_a.x + along_b.x),
        0.5 * (along_a.y + along_b.y),
        0.5 * (along_a.z + along_b.z),
    );
    let n_unnorm = along.cross(&xtan_a);
    let n = n_unnorm.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    (xtan_a.x - xtan_b.x) * n.x + (xtan_a.y - xtan_b.y) * n.y + (xtan_a.z - xtan_b.z) * n.z
}

/// **Biseptic (degree-7) solver on out-of-span mismatch.**
///
/// Lifts all 4 corner inner-CP pairs of patch B in z. Corner-
/// adjacent positions `(1,1), (1,6)` contribute to the spoke
/// cross-tangent via the `{B_1^7, B_6^7}` Bernstein modes —
/// **outside** the free-interior reachable span
/// `{B_2^7, B_3^7, B_4^7, B_5^7}`. The same structural gap as
/// biquintic-full (the `{B_1^5, B_{n-1}^5}` vs `{B_2^5, B_3^5}`
/// separation), now at degree 7.
///
/// Biseptic **does** measurably outperform biquintic-full on
/// this mismatch — 4 free DOF instead of 2 gives tighter best-
/// fit — but the residual does not reach machine precision.
/// This test records the observed reduction ratio (≈5.8× on
/// this fixture) and asserts at least 4× reduction to catch
/// solver regressions. Full cancellation requires the Peters
/// two-stripe or Loop-Schaefer construction (later solver work).
#[test]
fn biseptic_reduces_out_of_span_corner_mismatch_by_4x() {
    let mut patches = vec![flat_biseptic_patch(0.0, 1.0), flat_biseptic_patch(1.0, 2.0)];
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
    let ts = [0.15_f64, 0.3, 0.5, 0.7, 0.85];
    let residual_max = |patches: &[super::rectangular_biseptic::RectangularGregoryBiseptic; 2]| {
        ts.iter()
            .map(|&t| {
                g1_residual_biseptic(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t)
                    .abs()
            })
            .fold(0.0_f64, f64::max)
    };
    let mut pair_before: [super::rectangular_biseptic::RectangularGregoryBiseptic; 2] =
        [patches[0].clone(), patches[1].clone()];
    let max_before = residual_max(&pair_before);
    // Silence the `mut unused` warning — we only use pair_before for metrics.
    let _ = &mut pair_before;
    assert!(
        max_before > 0.1,
        "test scaffolding: expected nonzero initial residual, got {max_before:.4e}"
    );
    let _ = spoke_g1_solve(&mut patches, &spokes, 5).expect("solve should not be singular");
    let pair_after = [patches[0].clone(), patches[1].clone()];
    let max_after = residual_max(&pair_after);
    let ratio = max_before / max_after.max(1e-18);
    assert!(
        ratio > 4.0,
        "biseptic solver improvement ratio {ratio:.2}× below 4× floor: before {max_before:.4e} after {max_after:.4e}"
    );
}

/// Build a flat two-stripe patch spanning `[x0, x1] × [0, 1]`
/// in the `z = 0` plane.
fn flat_two_stripe_patch(
    x0: f64,
    x1: f64,
) -> super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe {
    use super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;
    let pt =
        |i: usize, j: usize| Point3::new(x0 + (x1 - x0) * (i as f64) / 7.0, (j as f64) / 7.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    RectangularGregoryTwoStripe { grid }
}

fn g1_residual_two_stripe(
    patch_a: &super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe,
    side_a: Side,
    patch_b: &super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe,
    side_b: Side,
    reversed: bool,
    t: f64,
) -> f64 {
    let t_b = if reversed { 1.0 - t } else { t };
    let xtan_a = match side_a {
        Side::UZero => patch_a.partial_u(0.0, t),
        Side::UOne => patch_a.partial_u(1.0, t),
        Side::VZero => patch_a.partial_v(t, 0.0),
        Side::VOne => patch_a.partial_v(t, 1.0),
    };
    let xtan_b = match side_b {
        Side::UZero => patch_b.partial_u(0.0, t_b),
        Side::UOne => patch_b.partial_u(1.0, t_b),
        Side::VZero => patch_b.partial_v(t_b, 0.0),
        Side::VOne => patch_b.partial_v(t_b, 1.0),
    };
    let along_a = match side_a {
        Side::UZero => patch_a.partial_v(0.0, t),
        Side::UOne => patch_a.partial_v(1.0, t),
        Side::VZero => patch_a.partial_u(t, 0.0),
        Side::VOne => patch_a.partial_u(t, 1.0),
    };
    let along_b = match side_b {
        Side::UZero => patch_b.partial_v(0.0, t_b),
        Side::UOne => patch_b.partial_v(1.0, t_b),
        Side::VZero => patch_b.partial_u(t_b, 0.0),
        Side::VOne => patch_b.partial_u(t_b, 1.0),
    };
    let along = Vec3::new(
        0.5 * (along_a.x + along_b.x),
        0.5 * (along_a.y + along_b.y),
        0.5 * (along_a.z + along_b.z),
    );
    let n_unnorm = along.cross(&xtan_a);
    let n = n_unnorm.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    (xtan_a.x - xtan_b.x) * n.x + (xtan_a.y - xtan_b.y) * n.y + (xtan_a.z - xtan_b.z) * n.z
}

/// **Two-stripe critical test.** On the scenario
/// where biseptic plateaus at residual 0.071 (83% reduction of
/// corner-inner-CP mismatch), the two-stripe solver has 6 free
/// DOF per spoke side giving full `{B_1^7..B_6^7}` Bernstein-
/// basis coverage. A pure-z perturbation of the corner-adjacent
/// positions `(1,1), (1,6)` on patch B must drop to machine
/// precision — this is the structural fix we claimed.
#[test]
fn two_stripe_cancels_corner_adjacent_perturbation_to_machine_precision() {
    let mut patches = vec![
        flat_two_stripe_patch(0.0, 1.0),
        flat_two_stripe_patch(1.0, 2.0),
    ];
    patches[1].grid[1][1].z += 0.15;
    patches[1].grid[1][6].z += 0.15;
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
    let ts = [0.15_f64, 0.3, 0.5, 0.7, 0.85];
    let max_before = ts
        .iter()
        .map(|&t| {
            g1_residual_two_stripe(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t)
                .abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_before > 0.05,
        "test scaffolding: expected nonzero initial residual; got {max_before:.4e}"
    );
    let _ = spoke_g1_solve(&mut patches, &spokes, 6).expect("solve should not be singular");
    let max_after = ts
        .iter()
        .map(|&t| {
            g1_residual_two_stripe(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t)
                .abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_after < 1e-6,
        "two-stripe must fully cancel B_1^7 + B_6^7 perturbation: before {max_before:.4e} after {max_after:.4e}"
    );
}

/// Biseptic **in-span** test: perturbation inside the solver's
/// reachable span (`{B_2^7, B_3^7, B_4^7, B_5^7}`) must drop to
/// machine precision. Lifting `grid[1][3].z` gives a pure
/// `B_3^7(v)` cross-tangent perturbation.
#[test]
fn biseptic_in_span_mismatch_drops_to_machine_precision() {
    let mut patches = vec![flat_biseptic_patch(0.0, 1.0), flat_biseptic_patch(1.0, 2.0)];
    patches[1].grid[1][3].z += 0.2;
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
    let ts = [0.2_f64, 0.4, 0.5, 0.6, 0.8];
    let max_before = ts
        .iter()
        .map(|&t| {
            g1_residual_biseptic(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t).abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_before > 0.05,
        "test scaffolding: expected nonzero initial residual, got {max_before:.4e}"
    );
    let _ = spoke_g1_solve(&mut patches, &spokes, 5).expect("solve should not be singular");
    let max_after = ts
        .iter()
        .map(|&t| {
            g1_residual_biseptic(&patches[0], Side::UOne, &patches[1], Side::UZero, false, t).abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_after < 1e-6,
        "biseptic in-span residual did not drop to precision: before {max_before:.4e} after {max_after:.4e}"
    );
}
