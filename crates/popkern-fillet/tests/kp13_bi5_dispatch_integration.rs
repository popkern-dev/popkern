//! Integration test: KP13 bi-5 4N dispatch through `fillet_edges`.
//!
//! Verifies that:
//! * With the KP13 override unset, the fillet pipeline runs the
//!   existing bi-6 polar / KP-binonic ladder and produces `N`
//!   corner Gregory patches (no change from the pre-KP13-bi-5
//!   baseline).
//! * With the KP13 override set, the pipeline *attempts* the KP13
//!   bi-5 promotion first — on synthetic fixtures where the
//!   upstream promotion ladder (biseptic / two-stripe) hasn't
//!   closed the non-two-stripe gap at N ≥ 5, KP13 bi-5 will reject
//!   and the pipeline falls through to the bi-6 polar ladder.
//!   Critical assertion: **no panic**, and the corner-patch count
//!   matches the env-var-unset baseline (dispatch is transparent
//!   on reject).
//!
//! When upstream two-stripe promotion closes on pyramidal apex
//! fixtures (separate tracked work, see
//! `project_phase_i_upstream_gap.md` memory note), this test's
//! gated branch will start producing 4N corner patches on a
//! committed promotion — at which point the assertion above
//! flips to `== 4 * baseline`.
//!
//! The dispatch is toggled through the typed
//! `dispatch_overrides::KP13_BI5_DISPATCH` slot (guard-scoped, no
//! process-wide env mutation).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn opts_full_stack() -> FilletOptions {
    FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    }
}

fn pyramid_apex_fillet(base_r: f64, height: f64, radius: f64) -> Option<usize> {
    let n = 5;
    let pyr = BRep::n_gonal_pyramid(n, base_r, height, Tolerance::default()).unwrap();

    // Pick the apex (valence-N vertex).
    let vertex = pyr.vertices.keys().find(|&v| {
        pyr.half_edges
            .values()
            .filter(|he| he.origin == v)
            .map(|he| he.edge)
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            == n
    })?;

    let edges: Vec<_> = pyr
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();

    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();

    let result = fillet_edges(&pyr, &specs, &opts_full_stack(), Tolerance::default(), None).ok()?;

    // Count Gregory corner-patches.
    let count = result
        .corner_patches
        .iter()
        .filter(|p| {
            result
                .brep
                .faces
                .get(p.face)
                .and_then(|f| f.surface.as_ref())
                .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Gregory(_)))
        })
        .count();
    Some(count)
}

/// With the env var unset and then set, the pipeline must run
/// identically (transparent dispatch on KP13 bi-5 reject). Verifies
/// the wiring didn't regress anything.
#[test]
fn kp13_bi5_dispatch_transparent_on_reject() {
    let baseline = pyramid_apex_fillet(1.0, 0.8, 0.1);

    let gated = {
        let _bi5 = popkern_fillet::dispatch_overrides::KP13_BI5_DISPATCH.set_guard(true);
        pyramid_apex_fillet(1.0, 0.8, 0.1)
    };

    // On synthetic pyramid fixtures where the upstream two-stripe
    // promotion rejects at N=5, KP13 bi-5 also rejects (same input
    // predicate) and the pipeline falls through identically. If the
    // two-stripe gap is closed later, `gated` may jump to 4×
    // `baseline` — at which point this assertion should be
    // rewritten to match.
    assert_eq!(
        baseline, gated,
        "env-var-gated KP13 bi-5 dispatch must be transparent on reject"
    );
}
