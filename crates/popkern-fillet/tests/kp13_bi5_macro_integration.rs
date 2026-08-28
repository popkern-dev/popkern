//! Integration test: exercise the KP13 bi-5 macro-patch pipeline
//! from within the fillet crate.

#![allow(clippy::needless_range_loop)]
//!
//! Confirms the bi-5 macro wiring:
//! - popkern-fillet can see and call
//!   [`popkern_core::geometry::gregory::assemble_bi5_macro_cap`]
//!   and friends.
//! - The end-to-end pipeline produces a finite `4n`-patch bi-9
//!   output for a well-posed input.
//!
//! This is the **public-API reach-through** demo; production
//! dispatch through the full fillet pipeline (selecting KP13 bi-5
//! vs. bi-9 direct vs. KP15 bi-6 polar based on valence
//! or input heuristics) is a separate concern that builds on top.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Bi5CapInput, Bi5CapOutput, Bi5CapTensorBorder, assemble_bi5_macro_cap,
};

/// Build a minimal synthetic tensor-border with a single non-zero
/// CP per sector.
fn build_test_input(n: usize) -> Bi5CapInput {
    let mut target = Bi5CapTensorBorder::zero(n);
    for k in 0..n {
        // Put a distinctive value in each sector's outer-corner CP.
        target.patches[k].cps[3][3] = Point3::new(k as f64, k as f64 * 0.5, k as f64 * 0.25);
    }
    Bi5CapInput {
        target,
        cc_limit: Point3::new(0.0, 0.0, 0.0),
    }
}

#[test]
fn kp13_bi5_pipeline_reaches_through_fillet_crate() {
    for n in [5, 6, 7] {
        let input = build_test_input(n);
        let out: Bi5CapOutput =
            assemble_bi5_macro_cap(&input).expect("KP13 bi-5 pipeline should succeed");
        assert_eq!(out.n, n);
        assert_eq!(out.patches_bi9.len(), 4 * n);
        // All CPs finite.
        for (idx, patch) in out.patches_bi9.iter().enumerate() {
            for i in 0..10 {
                for j in 0..10 {
                    let cp = patch[i][j];
                    assert!(
                        cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite(),
                        "n={n}, patch {idx} ({i}, {j}): non-finite {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn kp13_bi5_pipeline_supports_low_valence() {
    // KP09a § 6 / Deferred #4 coverage: n ∈ {3, 4} via pseudo-
    // inverse in the upstream circulant.
    for n in [3, 4] {
        let input = build_test_input(n);
        let out = assemble_bi5_macro_cap(&input).expect("n ∈ {3, 4} should succeed");
        assert_eq!(out.n, n);
        assert_eq!(out.patches_bi9.len(), 4 * n);
    }
}

#[test]
fn kp13_bi5_pipeline_deterministic_across_fillet_crate_calls() {
    let input = build_test_input(5);
    let a = assemble_bi5_macro_cap(&input).unwrap();
    let b = assemble_bi5_macro_cap(&input).unwrap();
    assert_eq!(a, b);
}
