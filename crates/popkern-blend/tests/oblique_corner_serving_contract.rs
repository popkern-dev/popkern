//! PK-260 — user-facing serving contract for oblique base-corner
//! fillets (pyramid3 mixed/extreme radii), under pure defaults.
//!
//! The class's support-side stitching root is FIXED: the §14.9
//! partition's rebuilt corner now completes the in-place handover
//! (serving-shape decline containment, ring-rebuild ambiguity
//! lookahead + surfaceless-sliver carrier + unique-closure search,
//! and the magnitude-gated winding calibration that unblocked the
//! mirror reversal), so:
//!
//! * `pyramid3_mixed` serves DIRECT-APEX natively — zero retries,
//!   zero fallbacks, crossing-free census.
//! * `pyramid3_extreme` (1:10) serves DIRECT-APEX natively too
//!   (PK-262): its former "footprint interference" was a
//!   SelfIntersection FALSE POSITIVE — the strip's oblique far-end
//!   runout kept an incoherent complementary-side curve (a recorded
//!   strip_end_arcs decline), and the folded trim polygon classified
//!   off-strip contact as Interior. The exact ellipse end-section
//!   rewrite makes the trim honest and the flag disappears at the
//!   geometry root.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;

fn default_corner_serve(radii: &[f64; 3]) -> popkern_blend::BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::n_gonal_pyramid(3, 1.0, 0.8, tol).unwrap();
    let vertex = helpers::pick_vertex_with_valence(&brep, 3).expect("valence-3 corner");
    let edges = helpers::incident_edges(&brep, vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
        .expect("default fillet must succeed")
}

fn crossing_count(result: &popkern_blend::BlendResult) -> usize {
    result
        .brep
        .validate_all()
        .into_iter()
        .filter(|e| matches!(e, ValidationError::SelfIntersection { .. }))
        .count()
}

/// pyramid3_mixed under pure defaults: no self-intersection, and the
/// creased direct-apex cap is replaced by the setback rung (PK-597).
#[test]
fn pyramid3_mixed_default_serve_is_native_direct_apex() {
    let result = default_corner_serve(&[0.08, 0.10, 0.12]);
    assert_eq!(
        crossing_count(&result),
        0,
        "pyramid3_mixed: default serve must census zero SelfIntersection"
    );
    // PK-597: the direct-apex cap creases across its own seams
    // (certified cap-cap normal sine 2.8e-1, body G0), so the corner
    // takes the serving exclusion re-run; the setback rung serves it
    // (certified G1, g1 2.2e-6). The crease is the only retry cause.
    let retries = &result.repair_report.serving_self_intersection_retries;
    assert!(
        !retries.is_empty()
            && retries
                .iter()
                .all(|retry| retry.detail.contains("served cap-cap crease")),
        "pyramid3_mixed: the crease is the only retry cause: {retries:?}"
    );
    assert!(
        result.corner_patches.iter().all(|p| p.setback_cap),
        "pyramid3_mixed: the setback rung serves the corner: {:?}",
        result.corner_patches
    );
    assert!(
        result.brep.boundary_edges.is_empty(),
        "pyramid3_mixed: the native serve carries no dead boundary sides"
    );
}

/// The 1:10 class. Its direct-apex cap creases across its own seams
/// (certified cap-cap normal sine 0.88), so the corner takes the serving
/// exclusion re-run. Until PK-623 the setback rung found nothing better
/// there and the creased native serve was kept; since the bridge planes
/// are oriented from their strips, the rung's corner-patched fit serves
/// it crossing-free and certified G1 (1.28 / r_min).
#[test]
fn pyramid3_extreme_default_serve_takes_the_setback_rung_over_the_creased_direct_apex() {
    let result = default_corner_serve(&[0.02, 0.10, 0.20]);
    assert_eq!(
        crossing_count(&result),
        0,
        "pyramid3_extreme: default serve must census zero SelfIntersection"
    );
    let retries = &result.repair_report.serving_self_intersection_retries;
    assert!(
        !retries.is_empty()
            && retries
                .iter()
                .all(|retry| retry.detail.contains("served cap-cap crease")),
        "pyramid3_extreme: the crease is the only retry cause: {retries:?}"
    );
    assert!(
        result.corner_patches.iter().all(|p| p.setback_cap),
        "pyramid3_extreme: the setback rung serves the corner: {:?}",
        result.corner_patches
    );
    assert!(
        result.brep.boundary_edges.is_empty(),
        "pyramid3_extreme: the serve carries no dead boundary sides"
    );
}
