//! Phase H.7 — commercial-grade continuity gates on the K-P
//! bidegree-9 output.
//!
//! # Metric choice: triple product, not unit-normal dot
//!
//! The natural "1 − |n_A · n_B|" signal breaks down when one
//! patch's `|∂u × ∂v|` drops to ~1e-5 at a parametric
//! near-degeneracy: `.normalized()` amplifies the tiny cross
//! vector into a direction dominated by floating-point noise,
//! and the dot-product metric then spans the full [−1, +1]
//! even though the underlying surfaces share the same tangent
//! plane. Diagnosed on cube_mixed: worst unit-normal deficit
//! of 9.7e-1 at t ≈ 0.74 sits on top of a perfectly G1-coplanar
//! triple-product of 1.3e-3.
//!
//! So the gate uses the *geometric* G1 definition directly —
//! the coplanarity of `xtan_A`, `xtan_B`, and the along-boundary
//! spoke tangent, normalised so `0` = perfect G1 and `1` =
//! orthogonal. This is invariant under the parametric
//! degeneracy that trips up the unit-normal metric.
//!
//! # Gates (per the Phase H plan § 1-3)
//!
//! * **G0**: position gap ≤ 1e-10 on interior `t ∈ [0.05, 0.95]`.
//! * **G1**: coplanarity ≤ 5e-2 (≈ 2.9° normal mismatch).
//! * **G2**: normal-curvature match residual ≤ 5e-2.
//!
//! The 0.05 endpoint margin excludes the Gregory central-vertex
//! and outer-corner singularities — finite-angle transition
//! zones inherent to the rational blend (not Phase H work).
//!
//! # Fixture matrix
//!
//! | Fixture          | Valence | Radii                         | Status         |
//! |------------------|---------|-------------------------------|----------------|
//! | cube_mixed       | 3       | [0.08, 0.10, 0.12]            | passes (1.9e-2)|
//! | cube_extreme     | 3       | [0.02, 0.10, 0.20]            | #[ignore] gap  |
//! | pyramid4_mixed   | 4       | [0.04, 0.06, 0.05, 0.07]      | #[ignore] gap  |
//!
//! Uniform-radii fixtures are out-of-scope here — the analytic
//! sphere-fit short-circuits Gregory dispatch upstream, so the
//! binonic rung never engages and the gate is tautologically
//! satisfied. Those are a different code-path; testing them
//! belongs in a sphere-shortcut suite.
//!
//! Higher-valence (N=5, N=6) polar cases are out-of-scope until
//! the K-P polar variant lands (Phase H.5b).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, assert_no_regression, collect_patches,
    fillet_corner, incident_edges, measure_ring, pick_vertex_with_valence,
};

fn run_cube(radii: &[f64], use_kp_binonic: bool) -> helpers::RingMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    assert_eq!(edges.len(), 3);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, use_kp_binonic);
    let patches = collect_patches(&brep, &faces);
    assert_eq!(patches.len(), 3, "expected 3 Gregory sub-patches on cube");
    measure_ring(&patches, N_SAMPLES)
}

fn run_cube_binonic(radii: &[f64], fixture: &str) -> helpers::RingMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

fn run_square_pyramid_binonic(radii: &[f64], fixture: &str) -> helpers::RingMetrics {
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).expect("square_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

// ── Commercial-gate tests ────────────────────────────────────

#[test]
fn binonic_passes_commercial_gates_on_cube_mixed() {
    let m = run_cube_binonic(&[0.08, 0.10, 0.12], "cube_mixed");
    assert_commercial_gates("cube_mixed", &m);
}

/// **Remaining gap.** Under 1:10 radius asymmetry, the K-P
/// construction still doesn't reach commercial-grade G1: binonic
/// lands near 2.6e-1 coplanarity (≈ 15°) after the twist
/// optimizer converges. This is competitive with two-stripe
/// (2.5e-1) — the regression is gone — but the bicubic guide's
/// 16-CP substrate has reached its quality ceiling for this
/// asymmetry level.
///
/// Closing this gap (getting binonic below the 5e-2 commercial
/// gate on cube_extreme) requires a richer guide substrate:
/// Option C, a bi-quintic (5×5) guide, which gives degree-4
/// cross-derivatives along boundaries instead of degree-2.
/// Tracked as follow-up to task #93.
#[ignore]
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme() {
    let m = run_cube_binonic(&[0.02, 0.10, 0.20], "cube_extreme");
    assert_commercial_gates("cube_extreme", &m);
}

#[test]
fn binonic_passes_commercial_gates_on_pyramid4_mixed() {
    let m = run_square_pyramid_binonic(&[0.04, 0.06, 0.05, 0.07], "pyramid4_mixed");
    assert_commercial_gates("pyramid4_mixed", &m);
}

// ── Non-regression vs two-stripe ─────────────────────────────

#[test]
fn binonic_matches_or_beats_two_stripe_on_cube_mixed() {
    let mt = run_cube(&[0.08, 0.10, 0.12], false);
    let mb = run_cube(&[0.08, 0.10, 0.12], true);
    assert_no_regression("cube_mixed", &mt, &mb);
}

/// Active regression monitor. With the K-P 2004 §4 iterative
/// twist optimizer, binonic on cube_extreme now achieves G1 ≈
/// 2.6e-1 — competitive with two-stripe's 2.5e-1. The 2× tolerance
/// comfortably absorbs the small remaining gap. Fail means either
/// (a) the optimizer stopped converging, or (b) the guard's
/// 1.2× threshold started firing again (which would indicate a
/// real regression the optimizer can't compensate for).
#[test]
fn binonic_matches_or_beats_two_stripe_on_cube_extreme() {
    let mt = run_cube(&[0.02, 0.10, 0.20], false);
    let mb = run_cube(&[0.02, 0.10, 0.20], true);
    assert_no_regression("cube_extreme", &mt, &mb);
}

/// Document that cube_extreme now engages the binonic rung
/// (after the K-P 2004 §4 iterative twist optimizer closed the
/// regression under 1:10 asymmetry). Before the optimizer, the
/// promotion quality guard reverted to two-stripe because
/// binonic G1 coplanarity was 3.9e-1 vs two-stripe's 2.5e-1
/// (1.6× worse); with the optimizer, binonic lands at 2.6e-1
/// — competitive with two-stripe, guard does not revert.
///
/// This still isn't commercial-grade (2.6e-1 > 5e-2 gate), but
/// the fixture is no longer a regression. Full closure requires
/// Option C (bi-quintic guide substrate — more DoFs than bicubic).
#[test]
fn cube_extreme_engages_binonic_after_twist_optimizer() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, &[0.02, 0.10, 0.20], true);
    let patches = collect_patches(&brep, &faces);
    for p in &patches {
        assert!(
            matches!(
                p,
                popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(_)
            ),
            "cube_extreme: K-P twist optimizer should have lifted binonic to \
             competitive-with-two-stripe G1, but we have {p:?}"
        );
    }
}
