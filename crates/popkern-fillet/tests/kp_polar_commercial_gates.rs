//! Phase I — commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (Phase H, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.
//!
//! # Fixture matrix
//!
//! | Fixture             | Valence | Radii                                      | Status |
//! |---------------------|---------|--------------------------------------------|--------|
//! | pentagon5_uniform   | 5       | all 0.05                                       | `#[ignore]` |
//! | pentagon5_mixed     | 5       | [0.04, 0.06, 0.05, 0.07, 0.05]                 | `#[ignore]` |
//! | pentagon5_variable  | 5       | [0.02, 0.04, 0.06, 0.04, 0.02] (1:3 ratio)     | `#[ignore]` I.22 |
//! | hexagon6_uniform    | 6       | all 0.05                                       | `#[ignore]` |
//! | hexagon6_mixed      | 6       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06]           | `#[ignore]` |
//! | hexagon6_variable   | 6       | [0.02, 0.04, 0.06, 0.06, 0.04, 0.02] (1:3)     | `#[ignore]` I.22 |
//! | octagon8_mixed      | 8       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06] | `#[ignore]` |
//!
//! All fixtures are `#[ignore]` at commit time for **two separate
//! reasons**:
//!
//! 1. **Upstream N ≥ 5 gap (blocking).** The fillet pipeline's
//!    three-quad vertex-blend path does not currently reach any
//!    polar rung for N ≥ 5 apex vertices — the `compute_blend_
//!    boundary` / kef / three-quad dispatch bails out before
//!    `apply_three_quad_decomposition_to_blend` is invoked.
//!    Confirmed empirically in `kp_polar_pipeline_probe.rs` — on
//!    every tested pyramid geometry the pipeline returns a single
//!    polar-NURBS corner patch, not an N-sector Gregory ring.
//!    Fixing this is Phase II scope.
//! 2. **Bi-6 output stabilisation (gate threshold).** Even once
//!    the upstream gap closes, the bi-6 output may need tuning
//!    to land under the 5e-2 commercial gates on adversarial
//!    fixtures. That's the normal progression we expect to track
//!    by flipping `#[ignore]` off per fixture as it passes.
//!
//! When (1) closes, `probe_pentagon5_geometry_sweep` in the probe
//! file will report `gregory_surface=5` on at least one row —
//! that's the green-light to start flipping these gates on.
//!
//! The uniform-radii fixtures are deliberately included to catch
//! regressions where a sphere-fit shortcut or upstream dispatch
//! diverts the ring away from the polar rung; `assert_all_binonic`
//! will fail loudly if a fixture stops producing `RectangularBiNonic`
//! patches.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, collect_patches, fillet_corner_full,
    incident_edges, measure_ring, pick_vertex_with_valence,
};

fn run_n_gonal_pyramid_polar(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    fixture: &str,
) -> helpers::RingMetrics {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner_full(&pyr, &edges, radii, true, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_uniform() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 0.8, &[0.05; 5], "pentagon5_uniform");
    assert_commercial_gates("pentagon5_uniform", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_mixed() {
    let m = run_n_gonal_pyramid_polar(
        5,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05],
        "pentagon5_mixed",
    );
    assert_commercial_gates("pentagon5_mixed", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_uniform() {
    let m = run_n_gonal_pyramid_polar(6, 1.0, 0.8, &[0.05; 6], "hexagon6_uniform");
    assert_commercial_gates("hexagon6_uniform", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_mixed() {
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        "hexagon6_mixed",
    );
    assert_commercial_gates("hexagon6_mixed", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_variable() {
    // Phase I.22 variable-radius fixture. Wider 1:3 ratio spread
    // (0.02..0.06) than `pentagon5_mixed` (0.04..0.07, ratio ~1.75)
    // stresses the per-sector radius threading. Core-level tests
    // (`variable_radius_*` in `kp15_bi6_polar_bridge_tests`)
    // already verify the construction is homogeneous + finite on
    // asymmetric radii; this fixture verifies end-to-end through
    // the fillet pipeline + commercial G1/G2 gates.
    let m = run_n_gonal_pyramid_polar(
        5,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.04, 0.02],
        "pentagon5_variable",
    );
    assert_commercial_gates("pentagon5_variable", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_variable() {
    // Phase I.22 variable-radius — hexagonal apex with 1:3 spread.
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02],
        "hexagon6_variable",
    );
    assert_commercial_gates("hexagon6_variable", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_octagon8_mixed() {
    let m = run_n_gonal_pyramid_polar(
        8,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
        "octagon8_mixed",
    );
    assert_commercial_gates("octagon8_mixed", &m);
}

// ── Dispatch smoke test (flag flip must not crash the pipeline) ──
//
// Flipping `use_kp_polar_gregory = true` must not panic `fillet_edges`
// end-to-end. This proves the flag's plumbing is wired, NOT that the
// polar engine actually runs on the fixture — with the upstream
// N ≥ 5 gap (see module header) the polar promote function is never
// invoked, so this test is compatible with that gap. When the gap
// closes, this test will still pass and the commercial-gate fixtures
// above will start being flip-able.

#[test]
fn polar_flag_flip_does_not_crash_pipeline_at_n5() {
    let pyr = BRep::n_gonal_pyramid(5, 1.0, 0.8, Tolerance::default()).unwrap();
    let Some(vertex) = pick_vertex_with_valence(&pyr, 5) else {
        return;
    };
    let edges = incident_edges(&pyr, vertex);
    let (_brep, _faces) =
        fillet_corner_full(&pyr, &edges, &[0.04, 0.06, 0.05, 0.07, 0.05], true, true);
}
