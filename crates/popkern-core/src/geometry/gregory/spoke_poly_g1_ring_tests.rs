//! Tests for the N-spoke ring polynomial-identity solver.

use super::*;
use crate::geometry::gregory::{RectangularGregoryTwoStripe, Side, Spoke, SpokeEnd};
use crate::geometry::{Point3, Vec3};

fn flat_patch(
    origin: (f64, f64),
    u_dir: (f64, f64),
    v_dir: (f64, f64),
) -> RectangularGregoryTwoStripe {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            let s = (i as f64) / 7.0;
            let t = (j as f64) / 7.0;
            grid[i][j] = Point3::new(
                origin.0 + s * u_dir.0 + t * v_dir.0,
                origin.1 + s * u_dir.1 + t * v_dir.1,
                0.0,
            );
        }
    }
    RectangularGregoryTwoStripe { grid }
}

fn unit_normal(p: &RectangularGregoryTwoStripe, u: f64, v: f64) -> Option<Vec3> {
    p.partial_u(u, v).cross(&p.partial_v(u, v)).normalized()
}

fn spoke_min_cos_dense(
    patches: &[RectangularGregoryTwoStripe],
    spokes: &[Spoke],
    n_samples: usize,
) -> f64 {
    let mut min_cos = 1.0_f64;
    for spoke in spokes {
        for k in 1..=n_samples {
            let t = (k as f64) / ((n_samples + 1) as f64);
            let (ua, va) = uv_for_side(spoke.a.side, t);
            let t_b = if spoke.b.reversed { 1.0 - t } else { t };
            let (ub, vb) = uv_for_side(spoke.b.side, t_b);
            let Some(na) = unit_normal(&patches[spoke.a.patch_idx], ua, va) else {
                return 0.0;
            };
            let Some(nb) = unit_normal(&patches[spoke.b.patch_idx], ub, vb) else {
                return 0.0;
            };
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            if cos < min_cos {
                min_cos = cos;
            }
        }
    }
    min_cos
}

fn uv_for_side(side: Side, t: f64) -> (f64, f64) {
    match side {
        Side::UZero => (0.0, t),
        Side::UOne => (1.0, t),
        Side::VZero => (t, 0.0),
        Side::VOne => (t, 1.0),
    }
}

#[test]
fn two_patch_one_spoke_flat_residual_zero() {
    // Two coplanar flat patches sharing x=0. Already G1.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut patches = vec![a, b];
    let spokes = vec![Spoke {
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
    let result = solve_ring_poly_g1(&mut patches, &spokes).expect("solved");
    assert_eq!(result.per_spoke_residual_max.len(), 1);
    assert!(
        result.per_spoke_residual_max[0] < 1e-10,
        "flat residual should be zero, got {}",
        result.per_spoke_residual_max[0]
    );
}

#[test]
fn two_patch_one_spoke_perturbed_recovers_g1() {
    // Perturb patch B's row 1 so the initial state is non-G1.
    // Only perturb j ∈ 1..7 — j = 0 / j = 7 are endpoint CPs that
    // also sit on adjacent-spoke boundaries, and the ring solver
    // intentionally excludes them from its free set to preserve
    // G0 continuity. Ring solver should recover G1 to machine
    // precision on the rest of the spoke.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    for j in 1..7 {
        b.grid[1][j].z += 0.12;
    }
    let mut patches = vec![a, b];
    let spokes = vec![Spoke {
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

    let min_cos_before = spoke_min_cos_dense(&patches, &spokes, 50);
    assert!(min_cos_before < 0.99, "expect non-G1 before solve");

    let result = solve_ring_poly_g1(&mut patches, &spokes).expect("solved");
    assert!(
        result.per_spoke_residual_max[0] < 1e-9,
        "ring residual should vanish, got {}",
        result.per_spoke_residual_max[0]
    );

    let min_cos_after = spoke_min_cos_dense(&patches, &spokes, 100);
    assert!(
        (min_cos_after - 1.0).abs() < 1e-9,
        "post-solve dense min_cos = {min_cos_after}, want 1.0"
    );
}

#[test]
fn alpha_sign_flips_when_cross_tangents_oppose() {
    // Build a setup where xtan_A and xtan_B point in OPPOSITE
    // directions at the mid-spoke. The ring solver should detect
    // α = +1 (not −1).
    //
    //   A: u-dir = (1, 0, 0), v-dir = (0, 1, 0), origin (-1, 0, 0).
    //      UOne at x=0.
    //   B: u-dir = (1, 0, 0), v-dir = (0, 1, 0), origin (0, 0, 0).
    //      But flip the sense so B's UZero's partial_u at u=0
    //      points in −x.  This needs u-direction reversed.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((1.0, 0.0), (-1.0, 0.0), (0.0, 1.0));
    // B: u-dir = −x, so partial_u points in −x. A's partial_u at
    // UOne points in +x. Opposite direction ⇒ α = +1 expected.
    let mut patches = vec![a, b];
    let spokes = vec![Spoke {
        a: SpokeEnd {
            patch_idx: 0,
            side: Side::UOne,
            reversed: false,
        },
        b: SpokeEnd {
            patch_idx: 1,
            side: Side::UOne,
            reversed: false,
        },
    }];
    let result = solve_ring_poly_g1(&mut patches, &spokes).expect("solved");
    assert_eq!(result.per_spoke_alpha[0], 1.0, "expected α = +1");
    assert!(
        result.per_spoke_residual_max[0] < 1e-9,
        "residual should still vanish for the correct α sign"
    );
}

#[test]
fn shared_cp_is_deduplicated() {
    // Two spokes both incident on patch 1, one on each of its two
    // adjacent sides (UZero and VZero). They share grid[1][1] on
    // patch 1. The ring solver must assign a single variable index
    // for that grid position — writing back after the solve must
    // produce ONE consistent value, not two.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let c = flat_patch((0.0, -1.0), (1.0, 0.0), (0.0, 1.0));
    let mut patches = vec![a, b, c];
    let spokes = vec![
        // spoke 0: A UOne ↔ B UZero (x=0 boundary between patch 0 and 1)
        Spoke {
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
        },
        // spoke 1: B VZero ↔ C VOne (y=0 boundary between patch 1 and 2)
        Spoke {
            a: SpokeEnd {
                patch_idx: 1,
                side: Side::VZero,
                reversed: false,
            },
            b: SpokeEnd {
                patch_idx: 2,
                side: Side::VOne,
                reversed: false,
            },
        },
    ];
    // grid[1][1] on patch 1 is shared between spoke 0 (UZero's
    // row-1 at j=1) and spoke 1 (VZero's row-1 at i=1).

    let result = solve_ring_poly_g1(&mut patches, &spokes).expect("solved");
    for (k, r) in result.per_spoke_residual_max.iter().enumerate() {
        assert!(
            *r < 1e-10,
            "spoke {k} residual {r} — flat coplanar should be zero",
        );
    }
    // Verify the value of grid[1][1] on patch 1 is CONSISTENT —
    // equal to what the flat geometry prescribed (it shouldn't
    // have moved, as the system was already G1).
    let expected = Point3::new(1.0 / 7.0, 1.0 / 7.0, 0.0);
    let actual = patches[1].grid[1][1];
    assert!(
        (actual.x - expected.x).abs() < 1e-10
            && (actual.y - expected.y).abs() < 1e-10
            && (actual.z - expected.z).abs() < 1e-10,
        "grid[1][1] moved: expected {expected:?}, got {actual:?}"
    );
}
