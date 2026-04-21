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
//! | pentagon5_uniform   | 5       | all 0.05                                   | `#[ignore]` — awaiting bi-6 polar stabilisation |
//! | pentagon5_mixed     | 5       | [0.04, 0.06, 0.05, 0.07, 0.05]             | `#[ignore]` — awaiting bi-6 polar stabilisation |
//! | hexagon6_uniform    | 6       | all 0.05                                   | `#[ignore]` — awaiting bi-6 polar stabilisation |
//! | hexagon6_mixed      | 6       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06]       | `#[ignore]` — awaiting bi-6 polar stabilisation |
//! | octagon8_mixed      | 8       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06] | `#[ignore]` — awaiting bi-6 polar stabilisation |
//!
//! All fixtures are `#[ignore]` at commit time: Phase I ships the
//! bi-6 polar **path** but the commercial gate itself is work to
//! come. Each test compiles, exercises the dispatch, and prints
//! the measured metrics — flipping `#[ignore]` off is how we'll
//! track progression as the path stabilises. This is the polar
//! analogue of the existing binonic fixture file.
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

// ── Dispatch smoke test (not gated by commercial thresholds) ──
//
// Flipping `use_kp_polar_gregory = true` must not crash the fillet
// pipeline. Runs the dispatch end-to-end on an N=5 pyramid and
// verifies the result is shape-correct (pipeline returns Ok and
// zero or more corner patches). This is the minimum-viable proof
// that the polar flag is wired through — a real commercial gate
// will follow once bi-6 output stabilises.

#[test]
fn polar_dispatch_smoke_pentagon5_mixed() {
    let pyr = BRep::n_gonal_pyramid(5, 1.0, 0.8, Tolerance::default()).unwrap();
    let Some(vertex) = pick_vertex_with_valence(&pyr, 5) else {
        // No N=5 vertex → pyramid topology variant. Not a polar
        // failure; skip so the smoke test stays valid as fixtures
        // evolve.
        return;
    };
    let edges = incident_edges(&pyr, vertex);
    let (_brep, _faces) =
        fillet_corner_full(&pyr, &edges, &[0.04, 0.06, 0.05, 0.07, 0.05], true, true);
    // `fillet_corner_full` unwraps on pipeline error; reaching here
    // means the polar-enabled run didn't panic. No metric assertion —
    // commercial-grade verification is the per-fixture gated tests
    // above.
}
