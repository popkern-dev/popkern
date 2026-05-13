//! Integration test for the §6c-A2.3 apex-construction dispatch
//! flag (`FilletOptions::use_a2_apex_construction`).
//!
//! When the flag is `true`, the legacy
//! `apply_blend_geometry`'s apex-vertex-blend loop attempts the
//! A2 planner + builder before falling back to the legacy
//! `insert_vertex_blend` path.  Failure surfaces as
//! `BlendWarning::A2FellBackToLegacy`.
//!
//! After §6c-A2.2.e/f and the constructive Gregory repeated-apex
//! normalization path:
//!   * **Sphere caps (uniform-radius cube apex)** — A2 succeeds
//!     end-to-end, no fallback fires.  Geometry matches legacy
//!     (face / vertex / edge counts).
//!   * **Gregory caps (variable-radius cube apex)** — current A2
//!     construction may decline with a structured
//!     `A2FellBackToLegacy` diagnostic while unsupported
//!     repeated-apex trim topology is still behind the fallback
//!     boundary.  The containment gate requires the fallback to be
//!     explicit and the resulting BRep to remain structurally valid.

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
fn a2_dispatch_variable_radius_records_containment_fallback() {
    let result_a2 = cube_with_variable_corner_fillets(true);

    let fallback_warnings: Vec<_> = result_a2
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }))
        .collect();
    assert_eq!(
        fallback_warnings.len(),
        1,
        "variable-radius cube should decline exactly once through the recorded A2 containment boundary; got {fallback_warnings:#?}",
    );
    let (_, BlendWarning::A2FellBackToLegacy { reason, .. }) = fallback_warnings[0] else {
        unreachable!("filtered to A2FellBackToLegacy above");
    };
    assert!(
        !reason.trim().is_empty() && reason.contains("construct_apex"),
        "A2 containment fallback must carry a useful constructive-failure reason; got {reason:?}",
    );
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
fn a2_dispatch_variable_geometry_is_structurally_valid_after_containment() {
    let result_a2 = cube_with_variable_corner_fillets(true);
    let errors = result_a2.brep.validate_all();
    let structural_errors = errors
        .iter()
        .filter(|err| {
            matches!(
                err,
                popkern_core::error::ValidationError::EulerPoincare { .. }
                    | popkern_core::error::ValidationError::BrokenTwin { .. }
                    | popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::InconsistentVertex { .. }
                    | popkern_core::error::ValidationError::NonManifoldEdge { .. }
                    | popkern_core::error::ValidationError::LowValence { .. }
                    | popkern_core::error::ValidationError::OpenLoop { .. }
                    | popkern_core::error::ValidationError::ShellGap { .. }
            )
        })
        .collect::<Vec<_>>();

    assert!(
        structural_errors.is_empty(),
        "variable-radius A2 containment result must have no structural/open-loop/shell-gap errors; got {structural_errors:#?}",
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
