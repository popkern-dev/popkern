//! N >= 5 production-route all-or-nothing gate helper (PK-075/PK-089).

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::helpers::{incident_edges, pick_vertex_with_valence};

pub(crate) fn assert_n_gonal_production_route_all_or_nothing(
    n: usize,
    radii: &[f64],
    fixture: &str,
) {
    let pyr = BRep::n_gonal_pyramid(n, 1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, apex);
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
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None)
        .unwrap_or_else(|e| panic!("[{fixture}] fillet_edges failed: {e}"));
    let apex_patches: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|patch| patch.n_sides == n)
        .collect();
    // Full tilings: n patches (one per sector, pre-setback macro caps) or
    // 4n patches (the PK-087 setback vertex fan commits 4 bi-9 pieces
    // per sector as individual cap faces).
    assert!(
        apex_patches.len() == n || apex_patches.len() == 4 * n || apex_patches.is_empty(),
        "[{fixture}] N={n} production dispatch must commit a full tiling ({n} macro or \
         {} setback-fan apex patches) or roll back to zero; got {}",
        4 * n,
        apex_patches.len()
    );
    if apex_patches.is_empty() {
        let has_decline = result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|d| format!("{d:?}").contains(".dispatch."));
        assert!(
            has_decline,
            "[{fixture}] all-or-nothing rollback must record structured dispatch-decline \
             diagnostics"
        );
    }
    result
        .brep
        .validate()
        .unwrap_or_else(|e| panic!("[{fixture}] production-route model must validate: {e:?}"));
}
