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

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::types::ApexConstructionRoute;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;

fn default_corner_serve(radii: &[f64; 3]) -> popkern_fillet::BlendResult {
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

fn crossing_count(result: &popkern_fillet::BlendResult) -> usize {
    result
        .brep
        .validate_all()
        .into_iter()
        .filter(|e| matches!(e, ValidationError::SelfIntersection { .. }))
        .count()
}

/// The fixed class: pyramid3_mixed serves direct-apex natively under
/// pure defaults.
#[test]
fn pyramid3_mixed_default_serve_is_native_direct_apex() {
    let result = default_corner_serve(&[0.08, 0.10, 0.12]);
    assert_eq!(
        crossing_count(&result),
        0,
        "pyramid3_mixed: default serve must census zero SelfIntersection"
    );
    assert!(
        result
            .repair_report
            .serving_self_intersection_retries
            .is_empty(),
        "pyramid3_mixed: the native serve needs no retry"
    );
    assert!(
        result
            .repair_report
            .apex_routes
            .iter()
            .all(|route| route.route == ApexConstructionRoute::DirectApex),
        "pyramid3_mixed: the corner serves direct-apex natively"
    );
    assert!(
        result.brep.boundary_edges.is_empty(),
        "pyramid3_mixed: the native serve carries no dead boundary sides"
    );
}

/// The 1:10 class, native since PK-262: the honest far-end trim
/// retires the false-positive crossing, so the steep corner serves
/// direct-apex with no serving machinery at all.
#[test]
fn pyramid3_extreme_default_serve_is_native_direct_apex() {
    let result = default_corner_serve(&[0.02, 0.10, 0.20]);
    assert_eq!(
        crossing_count(&result),
        0,
        "pyramid3_extreme: default serve must census zero SelfIntersection"
    );
    assert!(
        result
            .repair_report
            .serving_self_intersection_retries
            .is_empty(),
        "pyramid3_extreme: the native serve needs no retry"
    );
    assert!(
        result
            .repair_report
            .apex_routes
            .iter()
            .all(|route| route.route == ApexConstructionRoute::DirectApex),
        "pyramid3_extreme: the corner serves direct-apex natively"
    );
    assert!(
        result.brep.boundary_edges.is_empty(),
        "pyramid3_extreme: the native serve carries no dead boundary sides"
    );
}
