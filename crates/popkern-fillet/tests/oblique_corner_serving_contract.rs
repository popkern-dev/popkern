//! PK-230/PK-260 — user-facing serving contract for oblique
//! base-corner fillets (pyramid3 mixed/extreme radii).
//!
//! The direct-apex partition on this class currently leaves the big
//! support faces unstitched from the corner complex (PK-260): the
//! flag-matrix direct serve censuses SelfIntersection and OCC cannot
//! read its STEP export. Under DEFAULT options the PK-250 serving
//! retry detects the crossing and re-routes the corner to the legacy
//! dispatch, so what USERS receive is a validate_all-clean legacy
//! serve with the retry recorded.
//!
//! This pin locks that honest contract. When PK-260 fixes the
//! partition's support-side reconnection, the route flips to
//! DirectApex with zero retries — update this pin then (and re-enable
//! the serving retry in `build_under_direct_apex`).

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

        // The re-route is visible, not silent: the retry is recorded
        // on the repair report and the corner serves Legacy. PK-260
        // flips this to DirectApex with zero retries — update both
        // asserts then.
        assert!(
            !result
                .repair_report
                .serving_self_intersection_retries
                .is_empty(),
            "{name}: the serving retry must be recorded while PK-260 is open"
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
    }
}
