//! PK-230/PK-260 — user-facing serving contract for oblique
//! base-corner fillets (pyramid3 mixed/extreme radii).
//!
//! The direct-apex partition on this class currently leaves the big
//! support faces unstitched from the corner complex (PK-260): the
//! direct body carried dead boundary sides, a coplanar-overlapped
//! pocket complex, and an OCC-unreadable STEP export. The
//! serving-shape gate now DECLINES the route when a defect survives
//! every rung, and the serving retry treats the recorded decline as
//! a validity defect: the pipeline re-runs ONCE with the corner
//! excluded from direct-apex promotion, so the legacy dispatch
//! serves from clean legacy anchors (the mid-pipeline fallback body
//! is built on direct-apex-convention anchors and measures 8.5e-2
//! edge/surface deviations — never served).
//!
//! This pin locks that honest contract. When PK-260 fixes the
//! partition's support-side reconnection, the route flips to
//! DirectApex with zero retries — update this pin then.

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;

#[test]
fn pyramid3_corner_default_serve_is_clean_with_recorded_retry() {
    for (name, radii) in [
        ("pyramid3_mixed", [0.08, 0.10, 0.12]),
        ("pyramid3_extreme", [0.02, 0.10, 0.20]),
    ] {
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
        let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
            .unwrap_or_else(|e| panic!("{name} default fillet must succeed: {e:?}"));

        // The served body is crossing-free — the hard user contract.
        let crossings = result
            .brep
            .validate_all()
            .into_iter()
            .filter(|e| matches!(e, ValidationError::SelfIntersection { .. }))
            .count();
        assert_eq!(
            crossings, 0,
            "{name}: default serve must census zero SelfIntersection"
        );

        // The decline is visible, not silent: the serving retry
        // records it (with the serving-shape reason) and the adopted
        // body is the exclusion re-run, whose corner serves the
        // legacy dispatch natively. PK-260 flips this to DirectApex
        // with zero retries — update these asserts then.
        let retries = &result.repair_report.serving_self_intersection_retries;
        assert!(
            retries.iter().any(|retry| retry
                .detail
                .contains(popkern_fillet::types::SERVING_SHAPE_DECLINE_PREFIX)),
            "{name}: the serving-shape decline retry must be recorded while PK-260 is \
             open; got {retries:?}"
        );
        assert!(
            result
                .repair_report
                .apex_routes
                .iter()
                .filter(|route| route.vertex == vertex)
                .all(|route| route.route == popkern_fillet::types::ApexConstructionRoute::Legacy),
            "{name}: the corner serves the legacy dispatch while PK-260 is open"
        );
        assert!(
            result.brep.boundary_edges.is_empty(),
            "{name}: the legacy serve must carry no dead boundary sides"
        );
    }
}
