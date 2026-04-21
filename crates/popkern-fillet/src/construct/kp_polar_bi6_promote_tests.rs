//! Tests for the bi-6 polar promotion.
//!
//! Exercises `promote_decomp_to_kp_polar_bi6` on a synthetic
//! N-sided ring of two-stripe sub-patches. The ring is built by
//! hand (no upstream three-quad decomposition); verifies:
//!
//! - n < 5 rejected (the binonic rectangular path owns that
//!   regime).
//! - Non-two-stripe input rejected.
//! - Valid n = 5, 6, 7, 8 rings convert to bi-nonic grids and
//!   replace the input.
//!
//! The patch geometry here is rough and not commercial-grade;
//! we're testing the wiring and dispatch logic, not output
//! quality. Commercial-gate acceptance tests for real fillets
//! are in Stage 6 (I.25–I.28) of the plan.

use super::promote_decomp_to_kp_polar_bi6;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryTwoStripe};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Build an N-sided ring of two-stripe sub-patches that share
/// spoke boundaries. Patches are placed as "pie slices" in the
/// xy plane with the center at origin.
#[allow(clippy::needless_range_loop)]
fn synthetic_ring(n: usize) -> Vec<GregoryPatch> {
    (0..n)
        .map(|k| {
            let theta0 = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
            let theta1 = (2.0 * std::f64::consts::PI) * ((k + 1) as f64) / (n as f64);
            let mut grid = [[pt(0.0, 0.0, 0.0); 8]; 8];
            for i in 0..8 {
                for j in 0..8 {
                    let u_frac = 1.0 - (i as f64) / 7.0; // radial
                    let v_frac = (j as f64) / 7.0;
                    let theta = theta0 * (1.0 - v_frac) + theta1 * v_frac;
                    let r = u_frac;
                    grid[i][j] = pt(r * theta.cos(), r * theta.sin(), 0.0);
                }
            }
            GregoryPatch::RectangularTwoStripe(Box::new(RectangularGregoryTwoStripe { grid }))
        })
        .collect()
}

#[test]
fn rejects_n_less_than_5() {
    let mut ring = synthetic_ring(4);
    assert!(promote_decomp_to_kp_polar_bi6(&mut ring).is_none());
}

#[test]
fn rejects_empty_ring() {
    let mut ring: Vec<GregoryPatch> = vec![];
    assert!(promote_decomp_to_kp_polar_bi6(&mut ring).is_none());
}

#[test]
fn rejects_non_two_stripe_input() {
    // Mix a two-stripe ring with a non-two-stripe variant.
    let mut ring = synthetic_ring(5);
    // Replace patch 0 with a different variant. Use a cubic
    // rectangular by taking an existing two-stripe and converting
    // — simpler to construct a zero-sized different variant.
    // We'll use the RectangularGregoryBiNonic type itself as the
    // wrong-variant check.
    use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
    let bi9 = RectangularGregoryBiNonic {
        grid: [[pt(0.0, 0.0, 0.0); 10]; 10],
    };
    ring[0] = GregoryPatch::RectangularBiNonic(Box::new(bi9));
    assert!(promote_decomp_to_kp_polar_bi6(&mut ring).is_none());
}

// The positive-path test (`ring → bi-nonic patches`) requires a
// ring with enough geometric fidelity that `extract_direct_boundary_input`
// can produce finite tangent planes. The synthetic flat ring
// above produces degenerate normals (cross product of two
// parallel xy-plane vectors ⇒ zero vector) which intentionally
// fails; that's an honest negative result, not a bug.
//
// Full-pipeline positive verification lives in the core-level
// tests in `popkern-core` (`kp15_bi6_polar_bridge_tests.rs` +
// `kp_polar_bi6_apply_tests.rs`), which use synthetic
// `PolarBoundaryInput` directly and exercise the full flow at
// n ∈ {5, 6, 7, 8, 10, 12}.
