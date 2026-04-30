//! Unit tests for [`super::construct::construct_apex`].
//!
//! A2.2.a covers the Open rung (transactional no-op) and the
//! typed "not yet implemented" error for Sphere / Gregory.
//! A2.2.b through A2.2.f extend this with end-to-end tests that
//! actually mutate the BRep.

use popkern_core::topology::BRep;

use crate::types::{FilletRadius, FilletSpec};

use super::construct::construct_apex;
use super::error::ApexBuildError;
use super::plan::plan_apex_resolution;
use super::plan_e2e_test_helpers::{
    build_fillet_surfaces, force_unsupported_kind, pick_apex_and_incident_edges, test_tolerance,
};

#[test]
fn open_rung_construct_is_a_clean_noop() {
    // Build the same Open-rung scenario as
    // plan_e2e_tests::open_rung_when_too_few_fillets_succeed_on_cube
    // (cube apex with 2 of 3 fillets forced to plane-substituted
    // surfaces → both Direct and Setback fail → Open rung).
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    let radius = 0.1;
    let mut surfaces = build_fillet_surfaces(&cube, &edges, radius);
    force_unsupported_kind(&mut surfaces, 2);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let plan = plan_apex_resolution(&cube, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Open rung plan");
    assert!(plan.cap_kind.is_open(), "test setup expects Open plan");

    // Snapshot the BRep before construction; Open should be a
    // no-op so the post-build BRep must match exactly.
    let mut brep = cube.clone();
    let pre_v = brep.vertices.len();
    let pre_e = brep.edges.len();
    let pre_f = brep.faces.len();
    let pre_h = brep.half_edges.len();

    let result = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect("Open rung construct succeeds");

    // No new entities, no new handles, no telemetry mutation.
    assert!(result.created_faces.is_empty());
    assert!(result.apex_vertices.is_empty());
    assert!(result.apex_edges.is_empty());
    assert_eq!(brep.vertices.len(), pre_v, "Open must not add vertices");
    assert_eq!(brep.edges.len(), pre_e, "Open must not add edges");
    assert_eq!(brep.faces.len(), pre_f, "Open must not add faces");
    assert_eq!(brep.half_edges.len(), pre_h, "Open must not add half-edges");
    // Telemetry passes through.
    assert!(result.telemetry.fallback_reason.is_some());
}

#[test]
fn sphere_or_gregory_construct_returns_typed_unimplemented() {
    // Same fixture as direct_rung_succeeds_on_uniform_cube.
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    let radius = 0.1;
    let surfaces = build_fillet_surfaces(&cube, &edges, radius);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let plan = plan_apex_resolution(&cube, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Direct plan");
    assert!(
        plan.cap_kind.is_sphere() || plan.cap_kind.is_gregory(),
        "Direct rung produces Sphere or Gregory"
    );

    // BRep snapshot — on the not-yet-implemented error,
    // construct_apex must restore.
    let mut brep = cube.clone();
    let pre_v = brep.vertices.len();
    let pre_e = brep.edges.len();
    let pre_f = brep.faces.len();
    let pre_h = brep.half_edges.len();

    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("Sphere/Gregory paths not yet implemented in A2.2.a");
    assert!(
        matches!(err, ApexBuildError::PlanMismatch { .. }),
        "expected PlanMismatch, got {err:?}",
    );
    // Snapshot rollback: BRep must be byte-identical post-error.
    assert_eq!(brep.vertices.len(), pre_v);
    assert_eq!(brep.edges.len(), pre_e);
    assert_eq!(brep.faces.len(), pre_f);
    assert_eq!(brep.half_edges.len(), pre_h);
}

#[test]
fn cancellation_token_short_circuits_construct() {
    use crate::cancellation::CancellationToken;

    // Build any valid plan; cancellation should fire before any
    // build work runs.
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    let radius = 0.1;
    let surfaces = build_fillet_surfaces(&cube, &edges, radius);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let plan = plan_apex_resolution(&cube, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("plan succeeds");

    let mut brep = cube.clone();
    let token = CancellationToken::new();
    token.cancel();
    let err = construct_apex(&mut brep, &plan, &test_tolerance(), Some(&token))
        .expect_err("pre-cancelled token");
    assert!(
        matches!(err, ApexBuildError::Cancelled),
        "expected Cancelled, got {err:?}",
    );
}
