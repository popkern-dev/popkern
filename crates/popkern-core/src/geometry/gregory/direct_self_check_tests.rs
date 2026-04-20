//! Tests for [`super::measure_ring_quality`].

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::spoke_g1_solver::SpokeEnd;

/// A flat unit square two-stripe patch on the z=0 plane with
/// configurable origin and basis. Used to build simple ring
/// arrangements whose G1 / G0 can be reasoned about directly.
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

#[test]
fn two_coplanar_flat_patches_report_exact_g1_and_g0() {
    // Two adjacent unit squares sharing x=0. Exactly G1
    // everywhere on the shared boundary; G0 gap should be 0.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let patches = vec![a, b];
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

    let result = measure_ring_quality(&patches, &spokes, &SelfCheckOptions::default());
    assert!(
        (result.dense_min_cos_g1 - 1.0).abs() < 1e-12,
        "expected exact G1, got min_cos = {}",
        result.dense_min_cos_g1
    );
    assert!(
        result.max_position_gap_g0 < 1e-12,
        "expected exact G0, got max gap = {}",
        result.max_position_gap_g0
    );
    assert!(result.dense_min_cos_g2.is_none());
}

#[test]
fn non_g1_patch_ring_reports_low_min_cos() {
    // Second patch is tilted up at z = 0.05 on its row-1 inner CPs —
    // breaks G1 at the shared boundary.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    for j in 1..7 {
        b.grid[1][j].z += 0.3;
    }
    let patches = vec![a, b];
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

    let result = measure_ring_quality(&patches, &spokes, &SelfCheckOptions::default());
    assert!(
        result.dense_min_cos_g1 < 0.99,
        "expected non-G1, got min_cos = {}",
        result.dense_min_cos_g1
    );
    // G0 should still be exact (boundary CPs unchanged).
    assert!(result.max_position_gap_g0 < 1e-12);
}

#[test]
fn per_spoke_quality_reports_individual_residuals() {
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let patches = vec![a, b];
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
    let per_spoke = measure_ring_quality_per_spoke(&patches, &spokes, &SelfCheckOptions::default());
    assert_eq!(per_spoke.len(), 1);
    assert!((per_spoke[0].min_cos_g1 - 1.0).abs() < 1e-12);
    assert!(per_spoke[0].max_position_gap < 1e-12);
}

#[test]
fn g2_measurement_on_flat_coplanar_returns_agreement() {
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let patches = vec![a, b];
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
    let options = SelfCheckOptions {
        measure_g2: true,
        ..SelfCheckOptions::default()
    };
    let result = measure_ring_quality(&patches, &spokes, &options);
    // Flat coplanar patches have zero curvature — G2 agreement
    // trivially 1.0.
    let g2 = result.dense_min_cos_g2.expect("G2 measured");
    assert!(
        (g2 - 1.0).abs() < 1e-9,
        "expected G2 agreement ≈ 1.0 on flat, got {g2}"
    );
}

#[test]
fn empty_ring_reports_trivial_perfect_quality() {
    let patches: Vec<RectangularGregoryTwoStripe> = Vec::new();
    let spokes: Vec<Spoke> = Vec::new();
    let result = measure_ring_quality(&patches, &spokes, &SelfCheckOptions::default());
    // No spokes means nothing to measure — fold-identity returns 1.0 / 0.0.
    assert!((result.dense_min_cos_g1 - 1.0).abs() < 1e-12);
    assert!(result.max_position_gap_g0 < 1e-12);
}
