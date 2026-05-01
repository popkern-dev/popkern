//! Integration test for the §6c-A2.3 apex-construction dispatch
//! flag (`FilletOptions::use_a2_apex_construction`).
//!
//! When the flag is `true`, the legacy
//! `apply_blend_geometry`'s apex-vertex-blend loop attempts the
//! A2 planner + builder before falling back to the legacy
//! `insert_vertex_blend` path.  Failure surfaces as
//! `BlendWarning::A2FellBackToLegacy`.
//!
//! After §6c-A2.2.e/f landed Sphere-cap construction:
//!   * **Sphere caps (uniform-radius cube apex)** — A2 succeeds
//!     end-to-end, no fallback fires.  Geometry matches legacy
//!     (face / vertex / edge counts).
//!   * **Gregory caps (variable-radius cube apex)** — A2 errors
//!     at the cap-surface assignment stage (A2.4 territory) and
//!     the dispatch falls back per-apex with
//!     `BlendWarning::A2FellBackToLegacy`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use popkern_fillet::error::BlendWarning;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{BlendResult, fillet_edges};

fn cube_with_corner_fillets(use_a2: bool, radii: [f64; 3]) -> BlendResult {
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
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        use_a2_apex_construction: use_a2,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &options, tol, None).expect("fillet_edges succeeds")
}

fn cube_with_three_corner_fillets(use_a2: bool) -> BlendResult {
    cube_with_corner_fillets(use_a2, [0.10, 0.10, 0.10])
}

fn cube_with_variable_corner_fillets(use_a2: bool) -> BlendResult {
    cube_with_corner_fillets(use_a2, [0.08, 0.10, 0.12])
}

#[test]
fn a2_dispatch_succeeds_on_uniform_cube_without_fallback() {
    // After A2.2.e/f, Sphere caps build natively under A2.  The
    // uniform-radius cube apex dispatches to Sphere — A2 must
    // complete without firing `A2FellBackToLegacy` AND legacy
    // `insert_vertex_blend` must NOT run after A2 (else it would
    // emit `VertexModeFallback` because V is already gone).
    let result_a2 = cube_with_three_corner_fillets(true);

    let fallback_warnings: Vec<_> = result_a2
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }))
        .collect();
    assert!(
        fallback_warnings.is_empty(),
        "A2 should succeed end-to-end on uniform-radius cube (Sphere cap); \
         got {} A2FellBackToLegacy warning(s)",
        fallback_warnings.len(),
    );
    // §6c-A2.3.d skip-legacy regression guard: if A2 success
    // didn't actually skip the legacy `insert_vertex_blend` call,
    // legacy would re-enter for the now-removed apex vertex,
    // fail at `repair_vertex_outgoing`, and emit a
    // `VertexModeFallback` warning.
    let vertex_mode_fallbacks: Vec<_> = result_a2
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::VertexModeFallback { .. }))
        .collect();
    assert!(
        vertex_mode_fallbacks.is_empty(),
        "A2 success path must skip legacy; got {} VertexModeFallback warning(s) — \
         indicates legacy `insert_vertex_blend` ran after A2 already removed V",
        vertex_mode_fallbacks.len(),
    );
}

#[test]
fn a2_dispatch_falls_back_on_variable_radius_with_typed_warning() {
    // Variable radii cause A2 to fall back somewhere in the
    // pipeline.  The exact failure stage depends on plan vs.
    // post-Phase-C state coherence:
    //   * Today's variable-radius plan emits per-fillet sphere-fit
    //     positions for `apex_vertices`, but the post-Phase-C BRep
    //     has perp-edge-split vertices at each fillet's own radius
    //     — for non-uniform radii these don't coincide, so
    //     `resolve_apex_vertices` errors first.
    //   * Once A2.2.b's create-via-edge-subdivision branch lands,
    //     this should advance to A2.2.f's Gregory stub.
    // Either way, the dispatch must fall back per-apex with a
    // typed `A2FellBackToLegacy` warning carrying a non-empty
    // reason string.
    let result_a2 = cube_with_variable_corner_fillets(true);

    let fallback_warnings: Vec<_> = result_a2
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }))
        .collect();
    assert!(
        !fallback_warnings.is_empty(),
        "Variable-radius cube → A2 must fall back; got {} warning(s) total, \
         none of which are A2FellBackToLegacy",
        result_a2.edge_warnings.len(),
    );
    let reasons: Vec<&str> = fallback_warnings
        .iter()
        .filter_map(|(_, w)| match w {
            BlendWarning::A2FellBackToLegacy { reason, .. } => Some(reason.as_str()),
            _ => None,
        })
        .collect();
    assert!(
        !reasons.is_empty() && reasons.iter().all(|r| !r.is_empty()),
        "fallback warnings must carry non-empty reason strings; got: {reasons:?}",
    );
    eprintln!("[a2_dispatch] variable-radius fallback reasons: {reasons:?}");
}

#[test]
fn a2_dispatch_uniform_geometry_matches_legacy() {
    // A2 Sphere build on uniform cube must produce the same
    // topology as the pure-legacy run — vertex/edge/face/half-edge
    // counts identical.
    let result_legacy = cube_with_three_corner_fillets(false);
    let result_a2 = cube_with_three_corner_fillets(true);

    assert_eq!(
        result_legacy.brep.vertices.len(),
        result_a2.brep.vertices.len(),
        "vertex counts must match between legacy and A2 Sphere runs",
    );
    assert_eq!(
        result_legacy.brep.edges.len(),
        result_a2.brep.edges.len(),
        "edge counts must match between legacy and A2 Sphere runs",
    );
    assert_eq!(
        result_legacy.brep.faces.len(),
        result_a2.brep.faces.len(),
        "face counts must match between legacy and A2 Sphere runs",
    );
    assert_eq!(
        result_legacy.brep.half_edges.len(),
        result_a2.brep.half_edges.len(),
        "half-edge counts must match between legacy and A2 Sphere runs",
    );
}

#[test]
fn a2_dispatch_variable_geometry_matches_legacy_after_fallback() {
    // Variable-radius cube falls back to legacy per-apex.  After
    // the fallback runs, the BRep must match the pure-legacy run.
    let result_legacy = cube_with_variable_corner_fillets(false);
    let result_a2 = cube_with_variable_corner_fillets(true);

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
