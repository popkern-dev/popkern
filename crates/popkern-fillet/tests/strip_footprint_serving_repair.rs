//! PK-251 serving contract for mixed-radius direct-apex corners under
//! the order-honest self-intersection validator (analytic-first pair
//! canonicalization + the analytic × free-form implicit-distance fast
//! path).
//!
//! Fixture map (measured 2026-08-08, final configuration):
//! - pyramid4 [0.06, 0.08, 0.10, 0.08] and cube [0.02, 0.03, 0.04]:
//!   the served direct-apex caps measure CLEAN — the crossings the
//!   PK-245 triage attributed to them were artifacts of the marching
//!   solver's fitted-curve classification; the exact grid verdicts
//!   clear them. They serve DirectApex with no repair and no retry.
//! - pyramid4 [0.04, 0.06, 0.05, 0.07], pyramid4
//!   [0.04, 0.06, 0.05, 0.055], cube [0.15, 0.1, 0.1]: the served cap
//!   crosses a neighbour support along LOCKED boundary structure
//!   (bridge-boundary region or diving spoke) — boundary
//!   participation, which interior control-point pushes structurally
//!   cannot move. The serving repair attempts, stalls with the
//!   evidence recorded, and the PK-250 exclusion retry serves legacy
//!   validate-clean. The constructive partition fix is PK-253.
//! - cube [0.05, 0.1, 0.2]: the construct itself falls back to legacy
//!   for this set (pre-existing route); clean, no repair, no retry.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionRoute, BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

fn pyramid_default_result(radii: [f64; 4]) -> BlendResult {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(4, 1.0, 0.8, tol).expect("pyramid fixture");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.z > 0.1)
        .map(|(h, _)| h)
        .expect("pyramid apex");
    let slant_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        pyramid
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .filter(|e| seen.insert(*e))
            .collect()
    };
    let specs: Vec<_> = slant_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(&pyramid, &specs, &FilletOptions::default(), tol, None).expect("pyramid fillet")
}

fn cube_default_result(radii: [f64; 3]) -> BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let corner = brep
        .vertices
        .iter()
        .find(|(_, v)| {
            (v.position.x - 1.0).abs() < 1e-9
                && (v.position.y - 1.0).abs() < 1e-9
                && (v.position.z - 1.0).abs() < 1e-9
        })
        .map(|(vh, _)| vh)
        .expect("corner");
    let mut edges = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = he.origin;
        let b = brep.half_edges[he.twin].origin;
        if a == corner || b == corner {
            edges.push(e);
        }
    }
    assert_eq!(edges.len(), 3);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii)
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None).expect("cube fillet")
}

fn assert_direct_clean_serve(label: &str, result: &BlendResult) {
    assert!(
        result.repair_report.serving_footprint_repairs.is_empty(),
        "{label}: honest-validator-clean corner must not need a repair: {:?}",
        result.repair_report.serving_footprint_repairs
    );
    assert!(
        result
            .repair_report
            .serving_self_intersection_retries
            .is_empty(),
        "{label}: clean corner must not retry"
    );
    let direct = result
        .repair_report
        .apex_routes
        .iter()
        .any(|route| matches!(route.route, ApexConstructionRoute::DirectApex));
    assert!(
        direct,
        "{label}: the corner must keep its direct-apex serve: {:?}",
        result
            .repair_report
            .apex_routes
            .iter()
            .map(|route| route.route)
            .collect::<Vec<_>>()
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "{label}: direct serve must validate clean: {errors:?}"
    );
}

fn assert_boundary_participation_serve(label: &str, result: &BlendResult) {
    let repairs = &result.repair_report.serving_footprint_repairs;
    assert_eq!(
        repairs.len(),
        1,
        "{label}: repair attempt must be recorded: {repairs:?}"
    );
    assert!(
        !repairs[0].cleared && repairs[0].detail.contains("boundary participation"),
        "{label}: the boundary-participation class must stall with evidence: {repairs:?}"
    );
    let retries = &result.repair_report.serving_self_intersection_retries;
    assert_eq!(
        retries.len(),
        1,
        "{label}: the exclusion retry serves: {retries:?}"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "{label}: retried legacy serve must validate clean: {errors:?}"
    );
}

/// Corners the order-honest validator measures clean serve DIRECT.
#[test]
fn honest_clean_corners_serve_direct() {
    let result = pyramid_default_result([0.06, 0.08, 0.10, 0.08]);
    assert_direct_clean_serve("pyramid4 [0.06,0.08,0.10,0.08]", &result);
    let result = cube_default_result([0.02, 0.03, 0.04]);
    assert_direct_clean_serve("cube [0.02,0.03,0.04]", &result);
}

/// Boundary-participation corners stall the repair (recorded), retry,
/// and serve legacy validate-clean (PK-253 = the constructive fix).
#[test]
fn boundary_participation_corners_record_repair_and_retry() {
    let result = pyramid_default_result([0.04, 0.06, 0.05, 0.07]);
    assert_boundary_participation_serve("pyramid4 [0.04,0.06,0.05,0.07]", &result);
    let result = pyramid_default_result([0.04, 0.06, 0.05, 0.055]);
    assert_boundary_participation_serve("pyramid4 [0.04,0.06,0.05,0.055]", &result);
    let result = cube_default_result([0.15, 0.1, 0.1]);
    assert_boundary_participation_serve("cube [0.15,0.1,0.1]", &result);
}

/// The construct-level fallback set stays clean with no serving
/// machinery engaged.
#[test]
fn construct_fallback_set_serves_clean() {
    let result = cube_default_result([0.05, 0.1, 0.2]);
    assert!(
        result.repair_report.serving_footprint_repairs.is_empty()
            && result
                .repair_report
                .serving_self_intersection_retries
                .is_empty(),
        "no serving machinery should engage"
    );
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "must validate clean: {errors:?}");
}
