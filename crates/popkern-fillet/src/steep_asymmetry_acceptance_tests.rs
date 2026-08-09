//! PK-255 acceptance pins: steep-asymmetry (1:4) cube corners serve
//! DirectApex at commercial cap-outer quality, at BOTH mirrored
//! corners.
//!
//! Root chain (2026-08-09): a fillet-end source curve carried a
//! terminal sliver knot (full-multiplicity break 2.98e-8 inside its
//! end); the exact-boundary knot union propagated it into every
//! prepared curve, collapsing the derivative-row lever arm to 3.3e-9
//! and killing every sector's seam G1 at the (1,1,1) corner. On top,
//! the acceptance and strip-band measures clamped adjacent-period-
//! branch support parameters to the chart edge (PK-125 class,
//! `chart_sampling`), turning the healthy pole-searched serve into a
//! false 0.93 G1 decline. Fixes: terminal-sliver hygiene in boundary
//! preparation, the degree-elevation multi-span knot bookkeeping
//! (popkern-core), and period-honest sampling in every cap/support
//! normal measure.

use crate::construct::cap_outer_acceptance::cap_support_acceptance_error;
use crate::types::{ApexConstructionRoute, BlendResult, FilletOptions, FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const MAX_CAP_OUTER_G0: f64 = 5.0e-3;
const MAX_CAP_OUTER_G1: f64 = 5.0e-2;

fn cube_corner_result(radii: [f64; 3], corner_position: [f64; 3]) -> BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let corner = brep
        .vertices
        .iter()
        .find(|(_, v)| {
            (v.position.x - corner_position[0]).abs() < 1e-9
                && (v.position.y - corner_position[1]).abs() < 1e-9
                && (v.position.z - corner_position[2]).abs() < 1e-9
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
    crate::fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None).expect("cube fillet")
}

fn assert_direct_commercial_serve(label: &str, result: &BlendResult) {
    let direct = result
        .repair_report
        .apex_routes
        .iter()
        .any(|route| matches!(route.route, ApexConstructionRoute::DirectApex));
    assert!(
        direct,
        "{label}: must serve DirectApex: {:?}",
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
        "{label}: must validate clean: {errors:?}"
    );
    let cap_faces: Vec<_> = result
        .repair_report
        .apex_routes
        .iter()
        .flat_map(|route| route.built_faces.iter().copied())
        .collect();
    assert!(
        !cap_faces.is_empty(),
        "{label}: direct serve must record built faces"
    );
    let quality =
        cap_support_acceptance_error(&result.brep, &cap_faces, MAX_CAP_OUTER_G0, MAX_CAP_OUTER_G1);
    assert!(
        quality.is_none(),
        "{label}: served cap must hold cap-outer g0 <= {MAX_CAP_OUTER_G0:.0e} and g1 <= {MAX_CAP_OUTER_G1:.0e}: {quality:?}"
    );
}

/// The PK-255 ledger fixture serves DirectApex at commercial quality
/// with no serving machinery, at both mirrored corners (measured
/// g0 = 1.7e-8, g1 = 1.4e-5 at each).
#[test]
fn steep_asymmetry_corner_serves_direct_at_commercial_quality() {
    for corner in [[1.0, 1.0, 1.0], [0.0, 0.0, 0.0]] {
        let result = cube_corner_result([0.05, 0.1, 0.2], corner);
        let label = format!("cube [0.05,0.1,0.2] @ {corner:?}");
        assert!(
            result.repair_report.serving_footprint_repairs.is_empty()
                && result
                    .repair_report
                    .serving_self_intersection_retries
                    .is_empty(),
            "{label}: must serve without serving machinery"
        );
        assert_direct_commercial_serve(&label, &result);
    }
}

/// The [0.15, 0.1, 0.1] neighbor keeps its cleared-repair DirectApex
/// serve AND now holds the commercial bars (its pre-PK-255 blessed
/// serve measured worst_g1 0.959 — the family-worst crease).
#[test]
fn interior_crossing_neighbor_serves_repaired_at_commercial_quality() {
    let result = cube_corner_result([0.15, 0.1, 0.1], [1.0, 1.0, 1.0]);
    let label = "cube [0.15,0.1,0.1] @ (1,1,1)";
    let repairs = &result.repair_report.serving_footprint_repairs;
    assert_eq!(
        repairs.len(),
        1,
        "{label}: repair must be recorded: {repairs:?}"
    );
    assert!(
        repairs[0].cleared,
        "{label}: the repair must clear: {repairs:?}"
    );
    assert!(
        result
            .repair_report
            .serving_self_intersection_retries
            .is_empty(),
        "{label}: a cleared repair must not retry"
    );
    assert_direct_commercial_serve(label, &result);
}
