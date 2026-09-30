//! Shared unit-cube corner fixture for the binonic ladder test
//! binaries (`kp_binonic_commercial_gates`, `kp_binonic_extreme_ladder`).
#![allow(dead_code)]

use popkern_core::topology::BRep;

use crate::helpers::incident_edges;

/// Fillet the unit cube's first corner with the full binonic-ladder
/// flag set (shared by the commercial-gate and extreme-ladder test
/// binaries).
pub fn fillet_cube_result(radii: &[f64], use_kp_binonic: bool) -> popkern_blend::BlendResult {
    try_fillet_cube_result(radii, use_kp_binonic).unwrap()
}

/// Fallible variant of [`fillet_cube_result`].
pub fn try_fillet_cube_result(
    radii: &[f64],
    use_kp_binonic: bool,
) -> Result<popkern_blend::BlendResult, popkern_blend::BlendError> {
    use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
    let cube = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    assert_eq!(edges.len(), radii.len());
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
        use_kp_binonic_gregory: use_kp_binonic,
        // Candidate-quality instrument: the constructed candidate, not
        // the serving re-route (PK-250, and PK-597's crease re-run).
        serving_self_intersection_retry: false,
        ..FilletOptions::default()
    };
    fillet_edges(
        &cube,
        &specs,
        &opts,
        popkern_core::tolerance::Tolerance::default(),
        None,
    )
}
