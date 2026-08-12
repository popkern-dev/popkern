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
//! * `pyramid3_extreme` (1:10) still carries the residual
//!   cap-into-strip footprint interference (the PK-250/251 family):
//!   the serving retry detects the crossing and the exclusion re-run
//!   serves the legacy dispatch, recorded. Its half of this pin flips
//!   to the native contract when the footprint interference clears.

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

/// The residual class: pyramid3_extreme's footprint interference is
/// detected by the retry and the corner serves the legacy dispatch,
/// recorded and crossing-free.
#[test]
fn pyramid3_extreme_default_serve_is_clean_legacy_with_recorded_retry() {
    let result = default_corner_serve(&[0.02, 0.10, 0.20]);
    assert_eq!(
        crossing_count(&result),
        0,
        "pyramid3_extreme: default serve must census zero SelfIntersection"
    );
    let retries = &result.repair_report.serving_self_intersection_retries;
    assert!(
        retries
            .iter()
            .any(|retry| retry.detail.contains("Self-intersection")),
        "pyramid3_extreme: the footprint-crossing retry must be recorded; got {retries:?}"
    );
    assert!(
        result
            .repair_report
            .apex_routes
            .iter()
            .all(|route| route.route == ApexConstructionRoute::Legacy),
        "pyramid3_extreme: the corner serves the legacy dispatch while the footprint \
         interference is open"
    );
    assert!(
        result.brep.boundary_edges.is_empty(),
        "pyramid3_extreme: the legacy serve carries no dead boundary sides"
    );
}
