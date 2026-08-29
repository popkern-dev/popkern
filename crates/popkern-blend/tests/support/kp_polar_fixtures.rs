//! Shared N-gonal-pyramid fixture builders for the KP-polar gate
//! binaries (`kp_polar_commercial_gates.rs`,
//! `kp_polar_rim_tangency_gate.rs`).
//!
//! Consumed via `#[path = "support/kp_polar_fixtures.rs"] mod fixtures;`;
//! the including crate must also declare `mod kp_gate_helpers;` at
//! its root (same convention as `support/pns_apex_outcome.rs`).
#![allow(dead_code)]

use crate::kp_gate_helpers as helpers;
use helpers::{incident_edges, pick_vertex_with_valence};
use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

pub fn fillet_n_gonal_pyramid_result(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
) -> popkern_blend::BlendResult {
    fillet_n_gonal_pyramid_result_with_continuity(
        n,
        base_radius,
        height,
        radii,
        popkern_blend::ContinuityLevel::G1,
    )
}

pub fn fillet_n_gonal_pyramid_result_with_continuity(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    continuity: popkern_blend::ContinuityLevel,
) -> popkern_blend::BlendResult {
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
        continuity,
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
