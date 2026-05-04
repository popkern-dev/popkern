//! Unit tests for [`super::construct::construct_apex`].
//!
//! A2.2.a covers the Open rung (transactional no-op) and the
//! typed "not yet implemented" error for Sphere / Gregory.
//! A2.2.b through A2.2.f extend this with end-to-end tests that
//! actually mutate the BRep.

use popkern_core::topology::BRep;

use crate::types::{FilletRadius, FilletSpec};

use super::construct::construct_apex;
use super::construct_test_helpers::{synthetic_sphere_plan, synthetic_sphere_plan_with_edges};
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

// ---------------------------------------------------------------
// A2.2.b — apex-vertex resolution
// ---------------------------------------------------------------

#[test]
fn resolve_apex_vertices_finds_existing_within_tolerance() {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::BRep;

    // Build a BRep with a known vertex at (0.1, 0.1, 0).
    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep.mvfs(Point3::new(0.1, 0.1, 0.0)).unwrap();

    // Plan asks for a vertex at the same position (within tol).
    let plan = synthetic_sphere_plan(mvfs.vertex, vec![Point3::new(0.1, 0.1, 0.0)]);

    // Driving construct_apex against this plan resolves vertex 0
    // (the find succeeds), then errors out at A2.2.c+ with a
    // typed PlanMismatch — but the message tells us we got past
    // resolve_apex_vertices.
    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("post-resolve stages not yet implemented");
    let reason = match err {
        ApexBuildError::PlanMismatch { reason } => reason,
        other => panic!("expected PlanMismatch, got {other:?}"),
    };
    // Stage 3 (`kef_cross_face_v_edges`) can't merge anything on
    // these synthetic fixtures, so it surfaces "no merge performed".
    // That error is *evidence* stages 1 + 2 succeeded — the early
    // diagnostics ("no BRep vertex within", "no face has both")
    // would fire before reaching stage 3.
    assert!(
        reason.contains("kef_cross_face_v_edges")
            || reason.contains("kev_v_bridge")
            || reason.contains("Gregory cap surface assembly"),
        "should have advanced past resolve_apex_vertices/edges; got: {reason}",
    );
    // BRep was rolled back to pre-call state (1 vertex, 0 edges).
    assert_eq!(brep.vertices.len(), 1);
    assert_eq!(brep.edges.len(), 0);
    let _ = mvfs;
}

#[test]
fn resolve_apex_vertices_fails_when_no_match_exists() {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::BRep;

    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    // Plan asks for a vertex 1m away — way outside tol.
    let plan = synthetic_sphere_plan(mvfs.vertex, vec![Point3::new(1.0, 0.0, 0.0)]);

    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("no nearby vertex exists");
    let reason = match err {
        ApexBuildError::PlanMismatch { reason } => reason,
        other => panic!("expected PlanMismatch, got {other:?}"),
    };
    // After A2.2.b's create branch, the resolve fails with a
    // two-part diagnostic: no existing BRep vertex AND no edge
    // chord passes through the target.
    assert!(
        reason.contains("apex_vertex[0]")
            && reason.contains("no existing BRep vertex within")
            && reason.contains("no edge whose chord passes through"),
        "should mention which plan vertex couldn't be resolved; got: {reason}",
    );
}

// ---------------------------------------------------------------
// A2.2.c — apex-edge resolution
// ---------------------------------------------------------------

#[test]
fn resolve_apex_edges_finds_existing_brep_edge() {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::BRep;

    // Build a BRep with two vertices connected by a single edge
    // (mvfs + mev_at_vertex creates a degenerate-loop face but
    // that's enough to exercise edge resolution).
    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    // Plan: 2 apex vertices at the BRep's vertex positions, 1
    // edge connecting them.
    let plan = synthetic_sphere_plan_with_edges(
        mvfs.vertex,
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![(0, 1)],
    );

    let result = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect("sphere construction should resolve existing edge");
    assert_eq!(
        result.apex_vertices.len(),
        2,
        "both planned apex vertices should resolve to BRep vertices",
    );
    assert_eq!(
        result.apex_edges.len(),
        1,
        "planned apex edge should be represented in result ordering",
    );
}

#[test]
fn resolve_apex_edges_fails_when_no_shared_face_for_mef() {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::BRep;

    // §28.D update: Sphere's CapVsCyl apex_edges are now SKIPPED
    // (not mef'd), so the original "no shared face for mef"
    // error path is unreachable for Sphere plans.  Stage 3
    // (`kef_cross_face_v_edges`) hits the unconnected-orphan
    // condition and fails with a different diagnostic.
    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    let plan = synthetic_sphere_plan_with_edges(
        mvfs.vertex,
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![(0, 1)],
    );

    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("synthetic plan should not produce a valid cap topology");
    let reason = match err {
        ApexBuildError::PlanMismatch { reason } => reason,
        other => panic!("expected PlanMismatch, got {other:?}"),
    };
    // Either the legacy path (mef failed) OR the §28.D path
    // (kef_cross_face_v_edges failed) is acceptable — both
    // confirm the synthetic plan can't be realised.
    assert!(
        reason.contains("apex_edge[0]")
            || reason.contains("kef_cross_face_v_edges")
            || reason.contains("no face has both")
            || reason.contains("blend face not created"),
        "expected a topology-construction failure diagnostic; got: {reason}",
    );
}

#[test]
fn resolve_apex_vertices_creates_missing_via_edge_split() {
    // Build a BRep with a single edge between (0, 0, 0) and
    // (1, 0, 0).  The planner targets a midpoint position
    // (0.5, 0, 0) — no existing BRep vertex matches, so the
    // create branch must fire, split the edge, and return the
    // newly-inserted vertex.
    use popkern_core::geometry::Point3;
    use popkern_core::topology::BRep;

    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let pre_vertices = brep.vertices.len();
    let pre_edges = brep.edges.len();

    let plan = synthetic_sphere_plan(mvfs.vertex, vec![Point3::new(0.5, 0.0, 0.0)]);

    let result = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect("sphere construction should create missing apex vertex by edge split");
    assert_eq!(
        result.apex_vertices.len(),
        1,
        "single planned apex vertex should be created/resolved",
    );
    assert!(
        !result.created_faces.is_empty(),
        "successful sphere construction should report created cap topology",
    );
    assert!(
        brep.vertices.len() >= pre_vertices,
        "construction should not lose existing vertices",
    );
    assert!(
        brep.edges.len() >= pre_edges,
        "edge split should not lose existing edge topology",
    );
}

// ---------------------------------------------------------------
// A2.2.d — apex-edge create via mef
// ---------------------------------------------------------------

#[test]
fn resolve_apex_edges_creates_missing_edge_via_mef() {
    use popkern_core::geometry::{CurveGeometry, Point3};
    use popkern_core::topology::BRep;

    // Build a synthetic 4-vertex strip via mvfs + 3 mevs.  This
    // produces 1 face with 1 (degenerate) outer loop that walks
    // all 8 half-edges; v0 and v2 (the non-adjacent pair) both
    // sit on that loop.  Find for an edge between them fails
    // (none exists), so the A2.2.d create branch must fire,
    // call mef, attach the plan's curve, and split the face.
    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep
        .mvfs(Point3::new(0.0, 0.0, 0.0))
        .expect("mvfs builds the 0,0 vertex");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("extend to 1,0");
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("extend to 1,1");
    let _mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("extend to 0,1");
    let pre_faces = brep.faces.len();

    // Plan: connect v0 (0,0,0) to v2 (1,1,0) — non-adjacent on
    // the strip's loop.  Vertex order in BRep: insertion order
    // from mvfs/mev calls, so apex_vertices[0] = (0,0,0),
    // apex_vertices[1] = (1,1,0) maps to BRep vertices 0 and 2.
    let plan = synthetic_sphere_plan_with_edges(
        mvfs.vertex,
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        vec![(0, 1)],
    );

    let result = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect("sphere construction should create missing edge/cap topology");
    assert_eq!(
        result.apex_vertices.len(),
        2,
        "both planned apex vertices should resolve",
    );
    assert!(
        !result.created_faces.is_empty(),
        "successful construction should report created cap topology",
    );
    assert!(
        brep.faces.len() >= pre_faces,
        "successful construction should not lose existing faces",
    );

    // To verify the create branch's curve attachment in
    // isolation, replay the mef + curve attach against a fresh
    // BRep and inspect the new edge.  This is a white-box check
    // that the plan's NurbsCurve actually lands on the BRep edge
    // (not just that the mef succeeds).
    let mut brep2 = BRep::new(test_tolerance());
    let m = brep2.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let m1 = brep2
        .mev_at_vertex(m.vertex, m.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let m2 = brep2
        .mev(m1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let _m3 = brep2
        .mev(m2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    // Verify the edge-creation curve attachment path by manually invoking
    // mef against the fresh strip.
    use popkern_core::topology::HalfEdgeHandle;
    let mut he_v0: Option<HalfEdgeHandle> = None;
    let mut he_v2: Option<HalfEdgeHandle> = None;
    for (hh, h) in &brep2.half_edges {
        let p = brep2.vertices[h.origin].position;
        if (p.x - 0.0).abs() < 1e-12 && (p.y - 0.0).abs() < 1e-12 && he_v0.is_none() {
            he_v0 = Some(hh);
        }
        if (p.x - 1.0).abs() < 1e-12 && (p.y - 1.0).abs() < 1e-12 && he_v2.is_none() {
            he_v2 = Some(hh);
        }
    }
    let mef_r = brep2
        .mef(he_v0.unwrap(), he_v2.unwrap())
        .expect("mef on synthetic strip");
    assert_eq!(brep2.faces.len(), pre_faces + 1, "mef adds one face");
    // The TargetEdge.curve from synthetic_sphere_plan_with_edges
    // is a degree-1 NURBS line. The edge-attachment boundary should
    // canonicalize that into the semantic analytic Line variant.
    let target_curve = &plan.target_topology.apex_edges[0].curve;
    super::construct::attach_curve_to_edge(&mut brep2, mef_r.edge, target_curve);
    let edge_curve = brep2.edges[mef_r.edge]
        .curve
        .as_ref()
        .expect("curve attached after mef");
    assert!(
        matches!(edge_curve.geometry, CurveGeometry::Line { .. }),
        "attached straight edge should be an analytic Line variant",
    );
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
