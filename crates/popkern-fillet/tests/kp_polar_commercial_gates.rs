//! Phase I — commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (Phase H, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "kp_polar_support_cell_diagnostics.rs"]
mod support_cell_diagnostics;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, collect_patches, fillet_corner_full,
    incident_edges, measure_ring, pick_vertex_with_valence,
};
use support_cell_diagnostics::assert_support_cell_quality_gate;

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

fn fillet_n_gonal_pyramid_result(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
) -> popkern_fillet::BlendResult {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let specs = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None).unwrap()
}

/// Blunt-apex positive control (apex half-angle 18°, typical CAD).
/// Expected to take sphere-fit path and bypass polar entirely.
/// The `assert_all_binonic` assertion in `run_n_gonal_pyramid_polar`
/// will fail at the moment because the sphere-fit path produces a
/// single non-Gregory face, not N Gregory sub-patches. This is the
/// intended commercial behaviour for gentle corners — the test is
/// `#[ignore]` until the probe is refactored to assert "valid
/// commercial output" (sphere-fit OR polar-bi9) rather than
/// specifically "N Gregory bi-nonic patches."
///
/// Documented here as a reminder: pentagon5 at h=3.0 (18° apex) is
/// the REAL commercial-CAD test case; the `_uniform` / `_mixed` /
/// `_variable` sharp-apex tests below stress the polar path
/// specifically.
#[ignore]
#[test]
fn pentagon5_blunt_takes_sphere_fit_path() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 3.0, &[0.05; 5], "pentagon5_blunt");
    assert_commercial_gates("pentagon5_blunt", &m);
}
#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_uniform() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 0.8, &[0.05; 5], "pentagon5_uniform");
    assert_commercial_gates("pentagon5_uniform", &m);
}

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

#[test]
fn pentagon5_mixed_no_longer_records_polar_dispatch_decline() {
    let result = fillet_n_gonal_pyramid_result(5, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    result
        .brep
        .validate()
        .expect("pentagon5_mixed polar route must leave a valid BRep");
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .all(|diagnostic| !matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "rational_gregory.dispatch.declined"
            )),
        "pentagon5_mixed should reach N>=5 Gregory/polar dispatch without fallback diagnostics"
    );
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_uniform() {
    let m = run_n_gonal_pyramid_polar(6, 1.0, 0.8, &[0.05; 6], "hexagon6_uniform");
    assert_commercial_gates("hexagon6_uniform", &m);
}

#[ignore]
#[test]
fn polar_support_cell_quality_gate_on_hexagon6_mixed() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    assert_support_cell_quality_gate("hexagon6_mixed", &result, 6);
}

#[test]
fn hexagon6_mixed_accepts_support_cell_topology() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    assert_support_cell_quality_gate("hexagon6_mixed", &result, 6);
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
