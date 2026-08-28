//! PK-251 serving contract for mixed-radius direct-apex corners under
//! the order-honest self-intersection validator (analytic-first pair
//! canonicalization + the analytic × free-form implicit-distance fast
//! path).
//!
//! Fixture map (measured 2026-08-09, PK-253 closed: material-sense +
//! scan-depth repair, sign-change-honest fast-path verdicts):
//! - pyramid4 [0.06, 0.08, 0.10, 0.08] and cube [0.02, 0.03, 0.04]:
//!   the served direct-apex caps measure CLEAN — the crossings the
//!   PK-245 triage attributed to them were artifacts of the marching
//!   solver's fitted-curve classification; the exact grid verdicts
//!   clear them. They serve DirectApex with no repair and no retry.
//! - pyramid4 [0.04, 0.06, 0.05, 0.07] and [0.04, 0.06, 0.05, 0.055]:
//!   the historical "boundary participation" was one-sided tolerance
//!   HOVER, not a crossing — the cap's bridge boundary rides its own
//!   facet 6.5e-6 over the neighbour facet's plane near the shared
//!   raw-edge corner, with the whole corner window's minimum at
//!   machine zero. The sign-change-honest fast path no longer calls
//!   hover an intersection; these corners serve DirectApex with no
//!   repair and no retry.
//! - cube [0.15, 0.1, 0.1]: the one genuinely crossing fixture — the
//!   cap sector crosses the neighbour r=0.1 strip's retained footprint
//!   (sign changes through the cylinder). The serving repair CLEARS it
//!   (material sense orients pushes out of material; the dip is
//!   measured by scanning the crossing band's UV window, not inferred
//!   at Greville abscissae). Serves DirectApex with one cleared repair
//!   record and no retry.
//! - cube [0.05, 0.1, 0.2]: serves DirectApex clean with no repair
//!   and no retry (PK-255: terminal-sliver knot hygiene + period-
//!   honest acceptance sampling; quality bars pinned in-crate by
//!   `steep_asymmetry_acceptance_tests`).

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

/// Corners the sign-change-honest validator measures clean serve
/// DIRECT — including the mixed-radius pyramids whose bridge-boundary
/// hover was historically misread as a crossing.
#[test]
fn honest_clean_corners_serve_direct() {
    let result = pyramid_default_result([0.06, 0.08, 0.10, 0.08]);
    assert_direct_clean_serve("pyramid4 [0.06,0.08,0.10,0.08]", &result);
    let result = cube_default_result([0.02, 0.03, 0.04]);
    assert_direct_clean_serve("cube [0.02,0.03,0.04]", &result);
}

/// The tolerance-hover pyramids serve DIRECT with no repair machinery
/// engaged (PK-253).
#[test]
fn hover_pyramid_corners_serve_direct() {
    let result = pyramid_default_result([0.04, 0.06, 0.05, 0.07]);
    assert_direct_clean_serve("pyramid4 [0.04,0.06,0.05,0.07]", &result);
    let result = pyramid_default_result([0.04, 0.06, 0.05, 0.055]);
    assert_direct_clean_serve("pyramid4 [0.04,0.06,0.05,0.055]", &result);
}

/// The cube's interior crossing of a neighbour strip is repairable:
/// the serving repair clears it (material-sense-oriented, scan-depth
/// pushes) and the corner keeps its direct serve with no retry.
#[test]
fn cube_interior_crossing_repairs_to_direct_serve() {
    let result = cube_default_result([0.15, 0.1, 0.1]);
    let repairs = &result.repair_report.serving_footprint_repairs;
    assert_eq!(
        repairs.len(),
        1,
        "cube [0.15,0.1,0.1]: repair must be recorded: {repairs:?}"
    );
    assert!(
        repairs[0].cleared && repairs[0].detail == "cleared",
        "cube [0.15,0.1,0.1]: the repair must clear the crossing: {repairs:?}"
    );
    assert!(
        result
            .repair_report
            .serving_self_intersection_retries
            .is_empty(),
        "cube [0.15,0.1,0.1]: a cleared repair must not retry"
    );
    let direct = result
        .repair_report
        .apex_routes
        .iter()
        .any(|route| matches!(route.route, ApexConstructionRoute::DirectApex));
    assert!(direct, "cube [0.15,0.1,0.1]: direct serve kept");
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "cube [0.15,0.1,0.1]: repaired direct serve must validate clean: {errors:?}"
    );
}

/// PK-255 CLOSED (2026-08-09): the steep-asymmetry (1:4) fixture
/// serves DirectApex with no serving machinery, validate-clean. The
/// former vertex-blend-fan fallback was a terminal-sliver knot in one
/// fillet-end source curve (propagated by the knot union, collapsing
/// the derivative-row lever arm) plus period-clamp measurement
/// artifacts in the acceptance measures. The commercial quality bars
/// (cap-outer g0 5e-3 / g1 5e-2; measured 1.7e-8 / 1.4e-5) are pinned
/// in-crate by `steep_asymmetry_acceptance_tests` at both mirrored
/// corners.
#[test]
fn steep_asymmetry_corner_serves_direct_clean() {
    let result = cube_default_result([0.05, 0.1, 0.2]);
    assert_direct_clean_serve("cube [0.05,0.1,0.2]", &result);
}
