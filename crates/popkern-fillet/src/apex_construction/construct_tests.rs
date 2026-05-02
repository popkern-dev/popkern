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

// ---------------------------------------------------------------
// A2.2.b — apex-vertex resolution
// ---------------------------------------------------------------

/// Build a synthetic ApexPlan whose target_topology has the
/// supplied apex vertex positions and a Sphere cap kind.
/// Used by the A2.2.b find-only tests where we don't need the
/// geometry of the real sphere-fit pipeline — just a plan whose
/// apex_vertices match positions we control.
fn synthetic_sphere_plan(
    apex_vertices: Vec<popkern_core::geometry::Point3>,
) -> super::plan::ApexPlan {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::VertexHandle;

    use super::cap_kind::CapKind;
    use super::plan::ApexPlan;
    use super::target_topology::{TargetCapFace, TargetTopology};
    use super::telemetry::{ConstructionRung, PlanTelemetry};

    let n = apex_vertices.len();
    let cap_kind = CapKind::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        tangent_arc_endpoints: apex_vertices.clone(),
    };
    let cap_face = TargetCapFace::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        outer_loop_edge_indices: (0..n).collect(),
    };
    let target_topology = TargetTopology {
        apex_vertices,
        apex_edges: Vec::new(),
        cap_face,
        face_boundary_updates: Vec::new(),
    };
    let mut telemetry = PlanTelemetry::new();
    telemetry.rung = ConstructionRung::Direct;
    ApexPlan {
        apex_vertex: VertexHandle::default(),
        cap_kind,
        target_topology,
        telemetry,
    }
}

#[test]
fn resolve_apex_vertices_finds_existing_within_tolerance() {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::BRep;

    // Build a BRep with a known vertex at (0.1, 0.1, 0).
    let mut brep = BRep::new(test_tolerance());
    let mvfs = brep.mvfs(Point3::new(0.1, 0.1, 0.0)).unwrap();

    // Plan asks for a vertex at the same position (within tol).
    let plan = synthetic_sphere_plan(vec![Point3::new(0.1, 0.1, 0.0)]);

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
    // these synthetic fixtures (the plan's `apex_vertex` is a
    // default handle), so it surfaces "no merge performed".  That
    // error is *evidence* stages 1 + 2 succeeded — the early
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
    brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

    // Plan asks for a vertex 1m away — way outside tol.
    let plan = synthetic_sphere_plan(vec![Point3::new(1.0, 0.0, 0.0)]);

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

/// Build a synthetic ApexPlan with explicit apex_edges connecting
/// each pair of consecutive apex_vertices via a degree-1 NURBS.
/// Used by A2.2.c find-only tests.
fn synthetic_sphere_plan_with_edges(
    apex_vertices: Vec<popkern_core::geometry::Point3>,
    edge_pairs: Vec<(usize, usize)>,
) -> super::plan::ApexPlan {
    use popkern_core::geometry::{NurbsCurve, Point3};
    use popkern_core::topology::VertexHandle;

    use super::cap_kind::CapKind;
    use super::plan::ApexPlan;
    use super::target_topology::{TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology};
    use super::telemetry::{ConstructionRung, PlanTelemetry};

    let mut apex_edges: Vec<TargetEdge> = Vec::new();
    for &(origin_idx, dest_idx) in &edge_pairs {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![apex_vertices[origin_idx], apex_vertices[dest_idx]],
            weights: vec![1.0, 1.0],
        };
        apex_edges.push(TargetEdge {
            origin_vertex_idx: origin_idx,
            dest_vertex_idx: dest_idx,
            curve,
            role: TargetEdgeRole::CapVsCyl { fillet_idx: 0 },
        });
    }
    let n_v = apex_vertices.len();
    let cap_kind = CapKind::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        tangent_arc_endpoints: apex_vertices.clone(),
    };
    let cap_face = TargetCapFace::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.1,
        outer_loop_edge_indices: (0..apex_edges.len()).collect(),
    };
    let target_topology = TargetTopology {
        apex_vertices,
        apex_edges,
        cap_face,
        face_boundary_updates: Vec::new(),
    };
    let mut telemetry = PlanTelemetry::new();
    telemetry.rung = ConstructionRung::Direct;
    let _ = n_v;
    ApexPlan {
        apex_vertex: VertexHandle::default(),
        cap_kind,
        target_topology,
        telemetry,
    }
}

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
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        vec![(0, 1)],
    );

    // construct_apex resolves both vertices (stage 1) and the
    // edge (stage 2), then errors at stage 3+ — but the message
    // confirms we got past edge resolution.
    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("post-edge-resolution stages not yet implemented");
    let reason = match err {
        ApexBuildError::PlanMismatch { reason } => reason,
        other => panic!("expected PlanMismatch, got {other:?}"),
    };
    assert!(
        reason.contains("kef_cross_face_v_edges")
            || reason.contains("kev_v_bridge")
            || reason.contains("Gregory cap surface assembly"),
        "should have advanced past resolve_apex_edges; got: {reason}",
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
    brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    brep.vertices.insert(popkern_core::topology::Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    let plan = synthetic_sphere_plan_with_edges(
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

    let plan = synthetic_sphere_plan(vec![Point3::new(0.5, 0.0, 0.0)]);

    // construct_apex's stage 1 must succeed (via the create
    // branch); stage 3 then errors against the synthetic plan's
    // default apex_vertex handle.  Snapshot rolls back, so the
    // post-call BRep is byte-identical to pre-call — vertex /
    // edge counts confirm.
    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("post-resolve stages fail on synthetic plan");
    let reason = match err {
        ApexBuildError::PlanMismatch { reason } => reason,
        other => panic!("expected PlanMismatch, got {other:?}"),
    };
    assert!(
        reason.contains("kef_cross_face_v_edges")
            || reason.contains("kev_v_bridge")
            || reason.contains("Gregory cap surface assembly"),
        "should advance past stage 1 (create branch); got: {reason}",
    );
    // Snapshot rollback: vertex / edge counts unchanged.
    assert_eq!(
        brep.vertices.len(),
        pre_vertices,
        "snapshot rollback restores vertex count after create branch",
    );
    assert_eq!(
        brep.edges.len(),
        pre_edges,
        "snapshot rollback restores edge count after create branch",
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
    let pre_edges = brep.edges.len();

    // Plan: connect v0 (0,0,0) to v2 (1,1,0) — non-adjacent on
    // the strip's loop.  Vertex order in BRep: insertion order
    // from mvfs/mev calls, so apex_vertices[0] = (0,0,0),
    // apex_vertices[1] = (1,1,0) maps to BRep vertices 0 and 2.
    let plan = synthetic_sphere_plan_with_edges(
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        vec![(0, 1)],
    );

    // construct_apex fires the find/create stages, then errors
    // at the cap-face stub — the snapshot rolls everything back,
    // so post-call BRep equals pre-call BRep.  We verify the
    // create branch ran by checking the error message advanced
    // past edge resolution.
    let err = construct_apex(&mut brep, &plan, &test_tolerance(), None)
        .expect_err("cap face assembly stub fires after edge create");
    let reason = match err {
        ApexBuildError::PlanMismatch { reason } => reason,
        other => panic!("expected PlanMismatch, got {other:?}"),
    };
    assert!(
        reason.contains("kef_cross_face_v_edges")
            || reason.contains("kev_v_bridge")
            || reason.contains("Gregory cap surface assembly"),
        "should advance past edge create; got: {reason}",
    );
    // Snapshot rollback: face/edge counts unchanged.
    assert_eq!(
        brep.faces.len(),
        pre_faces,
        "snapshot rollback restores face count",
    );
    assert_eq!(
        brep.edges.len(),
        pre_edges,
        "snapshot rollback restores edge count",
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
    super::construct::construct_apex(&mut brep2, &plan, &test_tolerance(), None).ok();
    // After the failed call, snapshot has rolled back.  Re-run
    // just the edge-creation by manually invoking find/mef.
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
    // is a degree-1 NURBS line; when the create branch attaches
    // it correctly, the new edge's curve should be Some(Nurbs(_)).
    // Replicate the attach manually to verify the wrapping logic.
    let target_curve = &plan.target_topology.apex_edges[0].curve;
    super::construct::attach_curve_to_edge(&mut brep2, mef_r.edge, target_curve);
    let edge_curve = brep2.edges[mef_r.edge]
        .curve
        .as_ref()
        .expect("curve attached after mef");
    assert!(
        matches!(edge_curve.geometry, CurveGeometry::Nurbs(_)),
        "attached curve should be a NURBS variant",
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
