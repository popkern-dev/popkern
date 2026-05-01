//! Integration test for the §6c-A2.3 apex-construction dispatch
//! flag (`FilletOptions::use_a2_apex_construction`).
//!
//! When the flag is `true`, the legacy
//! `apply_blend_geometry`'s apex-vertex-blend loop attempts the
//! A2 planner + builder before falling back to the legacy
//! `insert_vertex_blend` path.  Failure surfaces as
//! `BlendWarning::A2FellBackToLegacy`.
//!
//! In the current A2.2.b/c → A2.2.d/e/f transition window the
//! A2 builder errors on Sphere/Gregory paths by design, so
//! every cube/pyramid apex hits the fallback.  The test
//! asserts:
//!   * fillet_edges still succeeds (per-apex fallback never
//!     causes a hard failure);
//!   * the warning fires;
//!   * the result is geometrically equivalent to the
//!     flag=false legacy run (face / vertex / edge counts).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use popkern_fillet::error::BlendWarning;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{BlendResult, fillet_edges};

fn cube_with_three_corner_fillets(use_a2: bool) -> BlendResult {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let apex = cube.vertices.keys().next().unwrap();
    let edges: Vec<EdgeHandle> = {
        let mut edge_set: std::collections::BTreeSet<EdgeHandle> =
            std::collections::BTreeSet::new();
        for he in cube.half_edges.values() {
            if he.origin == apex {
                edge_set.insert(he.edge);
            }
        }
        edge_set.into_iter().collect()
    };
    assert_eq!(edges.len(), 3, "cube apex must have 3 incident edges");
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.10),
        })
        .collect();
    let options = FilletOptions {
        use_a2_apex_construction: use_a2,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &options, tol, None).expect("fillet_edges succeeds")
}

#[test]
fn a2_dispatch_falls_back_to_legacy_with_typed_warning() {
    // Run the A2 path.  Today A2.2.b/c are landed but A2.2.d/e/f
    // are not, so construct_apex errors on Sphere/Gregory and
    // the dispatch falls back per-apex.
    let result_a2 = cube_with_three_corner_fillets(true);

    // The fallback warning must fire at least once (the cube's
    // single trihedral apex).
    let fallback_warnings: Vec<_> = result_a2
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }))
        .collect();
    assert!(
        !fallback_warnings.is_empty(),
        "A2 dispatch with use_a2_apex_construction=true should fire \
         A2FellBackToLegacy at least once on a cube apex; got \
         {} warning(s) total",
        result_a2.edge_warnings.len(),
    );
    // The reason should mention construct_apex (since plan succeeds
    // for cube apex; the failure is at the build stage).
    let reasons: Vec<&str> = fallback_warnings
        .iter()
        .filter_map(|(_, w)| match w {
            BlendWarning::A2FellBackToLegacy { reason, .. } => Some(reason.as_str()),
            _ => None,
        })
        .collect();
    // Reason must be non-empty.  The exact stage that fails
    // (descendant resolution / plan_apex_resolution /
    // construct_apex) depends on transient state and is not
    // pinned by this test — the contract is "the fallback
    // fires AND a reason is recorded."
    assert!(
        !reasons.is_empty() && reasons.iter().all(|r| !r.is_empty()),
        "fallback warnings must carry non-empty reason strings; got: {reasons:?}",
    );
    // Surface the actual reason for diagnostic visibility — the
    // string lets future maintainers see at a glance which A2
    // stage is failing today (plan_apex_resolution vs construct_apex)
    // without re-running the dispatch under a debugger.
    eprintln!("[a2_dispatch] fallback reasons observed: {reasons:?}");
}

#[test]
fn a2_dispatch_geometry_matches_legacy_after_fallback() {
    // After the per-apex fallback runs, the BRep should match
    // the pure-legacy run.  This proves the dispatch flag adds
    // no behavioural change beyond the warning.
    let result_legacy = cube_with_three_corner_fillets(false);
    let result_a2 = cube_with_three_corner_fillets(true);

    assert_eq!(
        result_legacy.brep.vertices.len(),
        result_a2.brep.vertices.len(),
        "vertex counts must match between legacy and A2-fallback runs",
    );
    assert_eq!(
        result_legacy.brep.edges.len(),
        result_a2.brep.edges.len(),
        "edge counts must match between legacy and A2-fallback runs",
    );
    assert_eq!(
        result_legacy.brep.faces.len(),
        result_a2.brep.faces.len(),
        "face counts must match between legacy and A2-fallback runs",
    );
    assert_eq!(
        result_legacy.brep.half_edges.len(),
        result_a2.brep.half_edges.len(),
        "half-edge counts must match between legacy and A2-fallback runs",
    );
}

#[test]
fn a2_dispatch_disabled_emits_no_a2_warning() {
    // Sanity: when the flag is off, the dispatch is dormant.
    let result = cube_with_three_corner_fillets(false);
    let any_a2 = result
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }));
    assert!(
        !any_a2,
        "with use_a2_apex_construction=false, A2FellBackToLegacy must not fire",
    );
}
