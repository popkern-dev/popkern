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

/// **Gap closed** (Phase H.8 Option C step 4 — L-M phase 3).
/// Under 1:10 radius asymmetry, the three-phase K-P stack now
/// drives G1 coplanarity to ~1.1e-2 (≈ 0.63°), well below the
/// 5e-2 commercial gate. Earlier progression:
///
///   3.85e-1 one-shot
///   → 2.56e-1 bicubic C-twist (Option B)
///   → 2.19e-1 bi-quintic FD gradient descent (Option C phase 2)
///   → 1.06e-1 bi-nonic FD gradient descent (Option C phase 3 GD floor)
///   → **1.12e-2 bi-nonic Levenberg-Marquardt (Option C phase 3 LM)**
///
/// The L-M gain came from switching the phase-3 optimizer from
/// naive gradient descent to Gauss-Newton damping on the per-
/// sample G1 coplanarity Jacobian: far better conditioning in
/// the over-parameterised 192N-dim search space.
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

/// cube_extreme now hits the commercial G1 gate (≤ 5e-2) after
/// the Phase H.8 Option C L-M phase-3 optimizer landed, so no
/// relaxed-continuity warning is expected. Kept as a regression
/// monitor: if the L-M stack ever fails to converge on this
/// fixture the binonic variant must still ship (no fallback to
/// two-stripe).
#[test]
fn cube_extreme_ships_binonic_no_warning() {
    use popkern_core::tolerance::Tolerance;
    use popkern_fillet::error::BlendWarning;
    use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.02_f64, 0.10, 0.20];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();

    // (1) Binonic shipped, not two-stripe.
    let faces: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
    let patches = collect_patches(&result.brep, &faces);
    for p in &patches {
        assert!(
            matches!(
                p,
                popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(_)
            ),
            "cube_extreme should ship binonic variant, got {p:?}"
        );
    }

    // (2) No CornerContinuityRelaxed warning — L-M brought G1 under
    // the commercial gate, so graceful-degradation shouldn't fire.
    let relaxed_warnings: Vec<_> = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CornerContinuityRelaxed { .. }))
        .collect();
    assert!(
        relaxed_warnings.is_empty(),
        "expected no CornerContinuityRelaxed warnings, got {}",
        relaxed_warnings.len()
    );
}
