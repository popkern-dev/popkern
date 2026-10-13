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
//! - cube [0.15, 0.1, 0.1]: the historical "crossing the footprint
//!   repair clears" was a cap rim wrapping three quarters round the
//!   r=0.1 strip, built from reversed-range end-section ellipses that
//!   converted to whole circles (PK-822/824). With exact conversion
//!   the first pass bounds its faces and only creases: the PK-597
//!   crease re-run serves the setback rung, no repair.
//! - cube [0.05, 0.1, 0.2]: serves DirectApex clean with no repair
//!   and no retry (PK-255: terminal-sliver knot hygiene + period-
//!   honest acceptance sampling; quality bars pinned in-crate by
//!   `steep_asymmetry_acceptance_tests`).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use popkern_blend::{
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
    cube_result_with(radii, &FilletOptions::default())
}

fn cube_result_with(radii: [f64; 3], options: &FilletOptions) -> BlendResult {
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
    fillet_edges(&brep, &specs, options, tol, None).expect("cube fillet")
}

/// Whether the serve took PK-597's crease re-run (a direct-apex cap
/// creasing across its own seams; not a self-intersection decision).
fn crease_rerun(result: &BlendResult) -> bool {
    result
        .repair_report
        .serving_self_intersection_retries
        .iter()
        .any(|retry| retry.detail.contains("served cap-cap crease"))
}

/// The self-intersection subject of this file: no footprint repair and no
/// self-intersection retry, and the served body validates clean. The
/// corner keeps its direct-apex serve unless its cap creases across its
/// own seams (PK-597), in which case the setback rung serves it.
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
            .iter()
            .all(|retry| retry.detail.contains("served cap-cap crease")),
        "{label}: clean corner must not take a self-intersection retry: {:?}",
        result.repair_report.serving_self_intersection_retries
    );
    let direct = result
        .repair_report
        .apex_routes
        .iter()
        .any(|route| matches!(route.route, ApexConstructionRoute::DirectApex));
    let setback =
        !result.corner_patches.is_empty() && result.corner_patches.iter().all(|p| p.setback_cap);
    assert!(
        direct || (crease_rerun(result) && setback),
        "{label}: direct serve, or the crease re-run on the setback rung: {:?}",
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
        "{label}: served body must validate clean: {errors:?}"
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

/// The cube [0.15, 0.1, 0.1] first pass bounds every face: its strip
/// end sections were reversed-range ellipse arcs that converted to the
/// WHOLE r = 0.1 circle, so the cap rim wrapped three quarters round the
/// strip behind its axis (PK-822; the conversion is PK-824). That wrap
/// was what three pins called a crossing the PK-251 footprint repair
/// cleared. Now the first pass has no winding defect and no crossing:
/// its cap only creases across its own seams, so PK-597's crease
/// re-run serves the setback rung with no footprint repair.
#[test]
fn cube_mixed_corner_first_pass_bounds_its_faces_and_crease_reruns() {
    let label = "cube [0.15,0.1,0.1]";
    let first = cube_result_with(
        [0.15, 0.1, 0.1],
        &FilletOptions {
            serving_self_intersection_retry: false,
            ..FilletOptions::default()
        },
    );
    let winding: Vec<_> = first
        .brep
        .validate_all()
        .into_iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::LoopPeriodWinding { .. }
            )
        })
        .collect();
    assert!(
        winding.is_empty(),
        "{label}: first pass trims bound regions: {winding:?}"
    );
    let result = cube_default_result([0.15, 0.1, 0.1]);
    let repairs = &result.repair_report.serving_footprint_repairs;
    assert!(
        repairs.is_empty(),
        "{label}: no crossing to repair: {repairs:?}"
    );
    assert!(
        crease_rerun(&result),
        "{label}: the crease re-run serves: {:?}",
        result.repair_report.serving_self_intersection_retries
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "{label}: the re-run serve must validate clean: {errors:?}"
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
