//! Integration test for the §6c-direct-apex dispatch apex-construction dispatch
//! flag (`FilletOptions::use_direct_apex_construction`).
//!
//! When the flag is `true`, the legacy
//! `apply_blend_geometry`'s apex-vertex-blend loop attempts the
//! direct-apex planner + builder before falling back to the legacy
//! `insert_vertex_blend` path.  Failure surfaces as
//! `BlendWarning::DirectApexFellBackToVertexBlend`.
//!
//! After §6c-direct-apex build.e/f and the constructive Gregory repeated-apex
//! normalization path:
//!   * **Sphere caps (uniform-radius cube apex)** — direct-apex succeeds
//!     end-to-end, no fallback fires.  Geometry matches legacy
//!     (face / vertex / edge counts).
//!   * **Gregory caps (mixed-radius cube apex)** — direct-apex now completes
//!     natively for the baseline mixed cube. Broader variable-radius
//!     fixtures still have containment coverage in the direct-apex matrix.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use popkern_fillet::error::BlendWarning;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{BlendResult, fillet_edges};

fn cube_with_corner_fillets(use_direct_apex: bool, radii: [f64; 3]) -> BlendResult {
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
        use_direct_apex_construction: use_direct_apex,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &options, tol, None).expect("fillet_edges succeeds")
}

fn cube_with_three_corner_fillets(use_direct_apex: bool) -> BlendResult {
    cube_with_corner_fillets(use_direct_apex, [0.10, 0.10, 0.10])
}

fn cube_with_variable_corner_fillets(use_direct_apex: bool) -> BlendResult {
    cube_with_corner_fillets(use_direct_apex, [0.08, 0.10, 0.12])
}

#[test]
fn direct_apex_dispatch_succeeds_on_uniform_cube_without_fallback() {
    // After direct-apex build.e/f, Sphere caps build natively under direct-apex.  The
    // uniform-radius cube apex dispatches to Sphere — direct-apex must
    // complete without firing `DirectApexFellBackToVertexBlend` AND legacy
    // `insert_vertex_blend` must NOT run after direct-apex (else it would
    // emit `VertexModeFallback` because V is already gone).
    let result_direct_apex = cube_with_three_corner_fillets(true);

    let fallback_warnings: Vec<_> = result_direct_apex
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::DirectApexFellBackToVertexBlend { .. }))
        .collect();
    assert!(
        fallback_warnings.is_empty(),
        "direct-apex should succeed end-to-end on uniform-radius cube (Sphere cap); \
         got {} DirectApexFellBackToVertexBlend warning(s)",
        fallback_warnings.len(),
    );
    // §6c-direct-apex dispatch.d skip-legacy regression guard: if direct-apex success
    // didn't actually skip the legacy `insert_vertex_blend` call,
    // legacy would re-enter for the now-removed apex vertex,
    // fail at `repair_vertex_outgoing`, and emit a
    // `VertexModeFallback` warning.
    let vertex_mode_fallbacks: Vec<_> = result_direct_apex
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::VertexModeFallback { .. }))
        .collect();
    assert!(
        vertex_mode_fallbacks.is_empty(),
        "direct-apex success path must skip legacy; got {} VertexModeFallback warning(s) — \
         indicates legacy `insert_vertex_blend` ran after direct-apex already removed V",
        vertex_mode_fallbacks.len(),
    );
}

#[test]
fn direct_apex_dispatch_mixed_radius_cube_completes_without_fallback() {
    let result_direct_apex = cube_with_variable_corner_fillets(true);

    let fallback_warnings: Vec<_> = result_direct_apex
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::DirectApexFellBackToVertexBlend { .. }))
        .collect();
    assert!(
        fallback_warnings.is_empty(),
        "mixed-radius cube should complete through native direct-apex; got {fallback_warnings:#?}",
    );
}

#[test]
fn direct_apex_dispatch_uniform_geometry_matches_legacy() {
    // direct-apex Sphere build on uniform cube must produce the same
    // topology as the pure-legacy run — vertex/edge/face/half-edge
    // counts identical.
    let result_legacy = cube_with_three_corner_fillets(false);
    let result_direct_apex = cube_with_three_corner_fillets(true);

    assert_eq!(
        result_legacy.brep.vertices.len(),
        result_direct_apex.brep.vertices.len(),
        "vertex counts must match between legacy and direct-apex Sphere runs",
    );
    assert_eq!(
        result_legacy.brep.edges.len(),
        result_direct_apex.brep.edges.len(),
        "edge counts must match between legacy and direct-apex Sphere runs",
    );
    assert_eq!(
        result_legacy.brep.faces.len(),
        result_direct_apex.brep.faces.len(),
        "face counts must match between legacy and direct-apex Sphere runs",
    );
    assert_eq!(
        result_legacy.brep.half_edges.len(),
        result_direct_apex.brep.half_edges.len(),
        "half-edge counts must match between legacy and direct-apex Sphere runs",
    );
}

#[test]
fn direct_apex_dispatch_mixed_radius_geometry_is_structurally_valid() {
    let result_direct_apex = cube_with_variable_corner_fillets(true);
    let errors = result_direct_apex.brep.validate_all();
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
        "mixed-radius native direct-apex result must have no structural/open-loop/shell-gap errors; got {structural_errors:#?}",
    );
}

#[test]
fn direct_apex_dispatch_disabled_emits_no_direct_apex_warning() {
    // Sanity: when the flag is off, the dispatch is dormant.
    let result = cube_with_three_corner_fillets(false);
    let any_direct_apex = result
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::DirectApexFellBackToVertexBlend { .. }));
    assert!(
        !any_direct_apex,
        "with use_direct_apex_construction=false, DirectApexFellBackToVertexBlend must not fire",
    );
}
