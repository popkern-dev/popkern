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
    cube_corner_result_with_retry(radii, corner_position, true)
}

/// The corner fillet with the serving retry on or off (off observes
/// the construction).
fn cube_corner_result_with_retry(
    radii: [f64; 3],
    corner_position: [f64; 3],
    retry: bool,
) -> BlendResult {
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
    let options = FilletOptions {
        serving_self_intersection_retry: retry,
        ..FilletOptions::default()
    };
    crate::fillet_edges(&brep, &specs, &options, tol, None).expect("cube fillet")
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

/// PK-597: a direct-apex cap that creases across its own seams (the
/// certified cap-cap normal jump over the commercial G1 gate) takes the
/// serving exclusion re-run, and the setback rung serves the corner.
/// The cap-outer bars these fixtures were pinned on never measured the
/// cap's own seams.
fn assert_crease_rerun_serve(label: &str, result: &crate::types::BlendResult) {
    assert!(
        result
            .repair_report
            .serving_self_intersection_retries
            .iter()
            .any(|retry| retry.detail.contains("served cap-cap crease")),
        "{label}: the creased cap takes the exclusion re-run: {:?}",
        result.repair_report.serving_self_intersection_retries
    );
    assert!(
        result
            .corner_patches
            .iter()
            .all(|p| p.setback_cap && p.achieved_continuity >= crate::types::ContinuityLevel::G1),
        "{label}: the setback rung serves at G1 or better: {:?}",
        result.corner_patches
    );
    assert!(
        result.brep.validate().is_ok(),
        "{label}: the re-run serve validates"
    );
}

/// The PK-255 ledger fixture at (1,1,1) served DirectApex with no
/// serving machinery and cap-outer quality g0 1.7e-8 / g1 1.4e-5, but
/// its sector seams crease toward the cap centre (certified normal-sine
/// lower bounds 0.48-0.95; PK-597's 9.3e-3 read only the seams the
/// 30° intent bound gated, PK-609). Both it and the mirrored (0,0,0)
/// corner (certified 7.3e-2) take the exclusion re-run; the construction
/// keeps the commercial cap-outer quality it was pinned on.
#[test]
fn steep_asymmetry_corner_serves_direct_at_commercial_quality() {
    let label = "cube [0.05,0.1,0.2] @ (1,1,1)";
    let construction = cube_corner_result_with_retry([0.05, 0.1, 0.2], [1.0, 1.0, 1.0], false);
    assert_direct_commercial_serve(label, &construction);
    let result = cube_corner_result([0.05, 0.1, 0.2], [1.0, 1.0, 1.0]);
    assert_crease_rerun_serve(label, &result);
    let mirrored = cube_corner_result([0.05, 0.1, 0.2], [0.0, 0.0, 0.0]);
    assert_crease_rerun_serve("cube [0.05,0.1,0.2] @ (0,0,0)", &mirrored);
}

/// PK-256 made the mirrored (0,0,0) corner of [0.15, 0.1, 0.1] serve
/// DirectApex with no serving machinery (winding-honest bridge sides,
/// the five-side coincident-old-tangent shape, reversed-order seam
/// twins). Its cap creases across its own seams (certified cap-cap
/// normal sine 2.4e-1, PK-597), so it takes the exclusion re-run.
#[test]
fn mirrored_corner_bridge_winding_serves_direct_clean() {
    let result = cube_corner_result([0.15, 0.1, 0.1], [0.0, 0.0, 0.0]);
    assert_crease_rerun_serve("cube [0.15,0.1,0.1] @ (0,0,0)", &result);
}

/// The [0.15, 0.1, 0.1] neighbour's cap creases across its own seams,
/// so the corner takes the PK-597 crease re-run with no footprint
/// repair: its historical "cleared crossing" was a cap rim wrapped
/// round its strip by a reversed-range end-section ellipse converted to
/// a whole circle (PK-822/824).
#[test]
fn mixed_neighbor_serves_the_crease_rerun_at_commercial_quality() {
    let result = cube_corner_result([0.15, 0.1, 0.1], [1.0, 1.0, 1.0]);
    let label = "cube [0.15,0.1,0.1] @ (1,1,1)";
    assert_crease_rerun_serve(label, &result);
}
