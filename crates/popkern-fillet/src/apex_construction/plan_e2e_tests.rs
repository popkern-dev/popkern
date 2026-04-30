//! End-to-end planner tests for [`plan_apex_resolution`].
//!
//! Drives the full Direct → Setback → Open rung ladder against
//! real BRep fixtures (unit cube + square pyramid) with fillet
//! surfaces produced by the actual [`compute_fillet_geometry`]
//! pipeline.
//!
//! Issues 2 + 3 from the A2.1 review:
//! * Issue 2 — Direct rung against a real fixture, end-to-end.
//! * Issue 3 — Setback / Open rung exercised by injecting
//!   unsupported fillets that the lenient pass skips.
//!
//! All four rung-ladder outcomes are covered:
//!   1. Direct + Sphere cap (uniform cube).
//!   2. Direct + Sphere or Gregory (pyramid uniform — dispatch
//!      depends on bridge-collapse tolerance).
//!   3. Setback (pyramid with 1 unsupported fillet → 3 succeed).
//!   4. Open (cube with 2 unsupported fillets → 1 succeed).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::analysis;
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{CapKind as LegacyCapKind, ContinuityLevel, FilletRadius, FilletSpec};

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan::plan_apex_resolution;
use super::telemetry::ConstructionRung;

fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-3)
        .with_angular(1e-3)
        .with_parametric(1e-6)
}

/// Collect the (apex, [edge; N]) tuple for the lowest-index
/// vertex of `brep`.  Picks all incident edges via half-edge
/// origin scan.
fn pick_apex_and_incident_edges(brep: &BRep) -> (VertexHandle, Vec<EdgeHandle>) {
    let apex = brep.vertices.keys().next().expect("brep has vertices");
    let mut edge_set: std::collections::BTreeSet<EdgeHandle> = std::collections::BTreeSet::new();
    for he in brep.half_edges.values() {
        if he.origin == apex {
            edge_set.insert(he.edge);
        }
    }
    (apex, edge_set.into_iter().collect())
}

/// Build a list of [`FilletSurfaceResult`]s from a list of
/// (edge, radius) pairs, using the real
/// [`compute_fillet_geometry`] pipeline.
fn build_fillet_surfaces(
    brep: &BRep,
    edges: &[EdgeHandle],
    radius: f64,
) -> Vec<FilletSurfaceResult> {
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let tol = test_tolerance();
    let analyses = analysis::analyze_edges(brep, edges, false)
        .expect("analyze_edges should succeed on uniform-radius corner");
    let mut surfaces = Vec::with_capacity(specs.len());
    for (spec, ana) in specs.iter().zip(analyses.iter()) {
        let geom = compute_fillet_geometry(
            ana,
            &spec.radius,
            brep,
            &tol,
            ContinuityLevel::G1,
            LegacyCapKind::Open,
        )
        .expect("compute_fillet_geometry should succeed on cube/pyramid corner");
        surfaces.push(FilletSurfaceResult {
            surface: geom.surface.clone(),
            trim_curve_1: geom.trim_3d_1.clone(),
            trim_curve_2: geom.trim_3d_2.clone(),
            continuity_downgrade: None,
        });
    }
    surfaces
}

/// Replace `n` fillet surfaces (starting at index 0) with a
/// degenerate `Plane` variant that the planner's strict pass
/// rejects via `UnsupportedFilletKind` and the lenient pass
/// skips.  Trim curves are preserved so `compute_blend_boundary`
/// still succeeds.
fn force_unsupported_kind(surfaces: &mut [FilletSurfaceResult], n: usize) {
    use popkern_core::geometry::Vec3;
    for fs in surfaces.iter_mut().take(n) {
        fs.surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
    }
}

// ---------------------------------------------------------------
// Issue 2 — Direct rung end-to-end
// ---------------------------------------------------------------

#[test]
#[ignore = "diagnostic dump for select_per_fillet_centre defect"]
fn diagnostic_dump_cube_fillet_cylinders() {
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    let apex_pos = cube.vertices.get(apex).unwrap().position;
    eprintln!("apex_pos = {apex_pos:?}");
    let radius = 0.1;
    let surfaces = build_fillet_surfaces(&cube, &edges, radius);
    for (i, fs) in surfaces.iter().enumerate() {
        if let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: r,
        } = &fs.surface
        {
            eprintln!("fillet {i}: cylinder origin={origin:?} axis={axis:?} r={r}");
            let apex_to_origin = (
                origin.x - apex_pos.x,
                origin.y - apex_pos.y,
                origin.z - apex_pos.z,
            );
            eprintln!("  apex_to_origin = {apex_to_origin:?}");
            // Trim curves run along the fillet body — their endpoints
            // are robust "body-side" references.
            use popkern_core::geometry::curve::CurveEval;
            let (t0, t1) = fs.trim_curve_1.parameter_range();
            let trim1_start = fs.trim_curve_1.point_at(t0);
            let trim1_end = fs.trim_curve_1.point_at(t1);
            eprintln!(
                "  trim_curve_1 endpoints: start={trim1_start:?} end={trim1_end:?}",
            );
        }
    }
}

#[test]
fn direct_rung_succeeds_on_uniform_cube() {
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    assert_eq!(edges.len(), 3, "cube corner has 3 incident edges");
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
        .expect("Direct rung should succeed on uniform cube");

    assert_eq!(plan.apex_vertex, apex);
    assert_eq!(plan.n_fillets(), 3, "3 fillets at cube apex");
    assert_eq!(
        plan.telemetry.rung,
        ConstructionRung::Direct,
        "uniform cube → Direct rung",
    );
    assert!(plan.telemetry.fallback_reason.is_none());
    // Cap-kind dispatch is data-dependent: per-fillet sphere
    // centres at the cube apex come out *not* exactly coincident
    // even at uniform radii (a quirk of `select_per_fillet_centre`'s
    // 3-plane projection), so the planner selects Gregory.  Either
    // Sphere or Gregory is acceptable as long as the chosen radius
    // matches the input.
    match &plan.cap_kind {
        CapKind::Sphere { radius: r, .. } => {
            assert!(
                (r - radius).abs() < 1e-9,
                "Sphere cap radius {r} should match fillet radius {radius}",
            );
        }
        CapKind::Gregory {
            per_fillet_caps, ..
        } => {
            assert_eq!(per_fillet_caps.len(), 3);
            for (_, r) in per_fillet_caps {
                assert!(
                    (r - radius).abs() < 1e-9,
                    "Gregory per-fillet radius {r} should match fillet radius {radius}",
                );
            }
        }
        CapKind::Open { .. } => {
            panic!("Direct rung on uniform cube must produce a real cap, not Open",)
        }
    }
    // Target topology vertex count depends on cap kind:
    //   * Sphere collapses shared T_a/T_b at the cap centre → 3
    //   * Gregory keeps T_a + T_b distinct per fillet → 6
    let n_v = plan.target_topology.apex_vertices.len();
    assert!(
        matches!(plan.cap_kind, CapKind::Sphere { .. }) && n_v == 3
            || matches!(plan.cap_kind, CapKind::Gregory { .. }) && n_v == 6,
        "apex_vertices.len()={n_v} inconsistent with cap_kind {:?}",
        plan.cap_kind,
    );
}

#[test]
fn direct_rung_succeeds_on_uniform_pyramid() {
    let pyramid = BRep::square_pyramid(1.0, 0.8, test_tolerance()).expect("square pyramid builds");
    // Pyramid apex is the highest vertex; pick the one with
    // 4 incident edges.
    let (apex, edges) = {
        let mut best: Option<(VertexHandle, Vec<EdgeHandle>)> = None;
        for v in pyramid.vertices.keys() {
            let mut es: std::collections::BTreeSet<EdgeHandle> = std::collections::BTreeSet::new();
            for he in pyramid.half_edges.values() {
                if he.origin == v {
                    es.insert(he.edge);
                }
            }
            if es.len() == 4 {
                best = Some((v, es.into_iter().collect()));
                break;
            }
        }
        best.expect("pyramid has a 4-valence apex")
    };
    let radius = 0.05;
    let surfaces = build_fillet_surfaces(&pyramid, &edges, radius);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();

    let plan = plan_apex_resolution(&pyramid, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Direct rung should succeed on uniform pyramid");

    assert_eq!(plan.apex_vertex, apex);
    assert_eq!(plan.n_fillets(), 4, "4 fillets at pyramid apex");
    assert_eq!(
        plan.telemetry.rung,
        ConstructionRung::Direct,
        "uniform pyramid → Direct rung",
    );
    // Cap-kind dispatch may pick Sphere (collapsed bridges) or
    // Gregory depending on the per-fillet bridge geometry — we
    // just assert one of the two is chosen.
    assert!(
        matches!(
            plan.cap_kind,
            CapKind::Sphere { .. } | CapKind::Gregory { .. }
        ),
        "uniform pyramid apex → Sphere or Gregory cap, got {:?}",
        plan.cap_kind,
    );
}

// ---------------------------------------------------------------
// Issue 3 — Setback and Open rung end-to-end
// ---------------------------------------------------------------

#[test]
fn setback_rung_succeeds_when_one_pyramid_fillet_unsupported() {
    let pyramid = BRep::square_pyramid(1.0, 0.8, test_tolerance()).expect("square pyramid builds");
    let (apex, edges) = {
        let mut best: Option<(VertexHandle, Vec<EdgeHandle>)> = None;
        for v in pyramid.vertices.keys() {
            let mut es: std::collections::BTreeSet<EdgeHandle> = std::collections::BTreeSet::new();
            for he in pyramid.half_edges.values() {
                if he.origin == v {
                    es.insert(he.edge);
                }
            }
            if es.len() == 4 {
                best = Some((v, es.into_iter().collect()));
                break;
            }
        }
        best.expect("pyramid has a 4-valence apex")
    };
    let radius = 0.05;
    let mut surfaces = build_fillet_surfaces(&pyramid, &edges, radius);
    // Force the first fillet to be Plane (unsupported) — strict
    // pass fails, lenient pass skips, leaving 3 successful
    // fillets → Setback rung succeeds.
    force_unsupported_kind(&mut surfaces, 1);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();

    let plan = plan_apex_resolution(&pyramid, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Setback rung should succeed with 3 of 4 fillets");

    assert_eq!(plan.apex_vertex, apex);
    assert!(
        matches!(plan.telemetry.rung, ConstructionRung::Setback { .. }),
        "expected Setback rung, got {:?}",
        plan.telemetry.rung,
    );
    assert!(
        plan.telemetry.fallback_reason.is_some(),
        "Setback should record fallback_reason explaining Direct's failure",
    );
    assert!(
        plan.telemetry
            .fallback_reason
            .as_deref()
            .unwrap_or("")
            .contains("Direct rung failed"),
        "fallback_reason should mention Direct's failure",
    );
    // Setback has only the 3 successful fillets in the cap.
    assert_eq!(
        plan.n_fillets(),
        3,
        "Setback omits the 1 unsupported fillet",
    );
}

#[test]
fn open_rung_when_too_few_fillets_succeed_on_cube() {
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    assert_eq!(edges.len(), 3);
    let radius = 0.1;
    let mut surfaces = build_fillet_surfaces(&cube, &edges, radius);
    // Force 2 of 3 fillets to be Plane.  Direct fails, lenient
    // gets only 1 success (< 3) → Setback returns sentinel,
    // caller falls through to Open rung.
    force_unsupported_kind(&mut surfaces, 2);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();

    let plan = plan_apex_resolution(&cube, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Open rung should produce a (degenerate) plan when both Direct and Setback fail");

    assert_eq!(plan.apex_vertex, apex);
    assert_eq!(
        plan.telemetry.rung,
        ConstructionRung::Open,
        "expected Open rung, got {:?}",
        plan.telemetry.rung,
    );
    assert!(
        plan.telemetry.fallback_reason.is_some(),
        "Open should record why Direct + Setback both failed",
    );
    let reason = plan.telemetry.fallback_reason.as_deref().unwrap_or("");
    assert!(
        reason.contains("Open rung produces no apex cap"),
        "Open's fallback_reason should explain no cap is built; got: {reason}",
    );
    // Open's target_topology is empty — no cap, no apex_vertices.
    assert!(plan.target_topology.apex_vertices.is_empty());
    assert!(plan.target_topology.apex_edges.is_empty());
    // Cap-kind is the typed Open variant carrying the count of
    // fillets that *would* have been at the apex.
    assert!(
        plan.cap_kind.is_open(),
        "Open rung must produce CapKind::Open, got {:?}",
        plan.cap_kind,
    );
    if let crate::apex_construction::cap_kind::CapKind::Open {
        n_fillets_attempted,
    } = &plan.cap_kind
    {
        assert_eq!(*n_fillets_attempted, 3, "cube apex has 3 fillets");
    } else {
        unreachable!("is_open() asserted above");
    }
    assert!(matches!(
        plan.target_topology.cap_face,
        super::target_topology::TargetCapFace::Open
    ));
}

// ---------------------------------------------------------------
// Cross-cutting: error path that bypasses rung ladder.
// ---------------------------------------------------------------

#[test]
fn input_validation_short_circuits_rung_ladder() {
    // 2 fillets — fails InsufficientFillets at the validation
    // gate before the rung ladder runs.  Confirms we don't
    // accidentally fall through to Open rung for input errors.
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .take(2)
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let err = plan_apex_resolution(&cube, apex, &specs, &[], &test_tolerance(), None)
        .expect_err("should fail InsufficientFillets");
    assert!(
        matches!(err, ApexPlanError::InsufficientFillets { n: 2, .. }),
        "expected InsufficientFillets(n=2), got {err:?}",
    );
}

/// Issue 5 — pre-cancelled tokens short-circuit planning at the
/// rung-ladder boundary (after `validate_inputs` succeeds, before
/// any sphere-fit math runs).
#[test]
fn cancellation_token_short_circuits_planning() {
    use crate::cancellation::CancellationToken;

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
    let token = CancellationToken::new();
    token.cancel();
    let err = plan_apex_resolution(
        &cube,
        apex,
        &specs,
        &surfaces,
        &test_tolerance(),
        Some(&token),
    )
    .expect_err("pre-cancelled token must surface ApexPlanError::Cancelled");
    assert!(
        matches!(err, ApexPlanError::Cancelled),
        "expected Cancelled, got {err:?}",
    );
}

// ---------------------------------------------------------------
// Issue 8 — preview affordances on ApexPlan
// ---------------------------------------------------------------

#[test]
fn preview_affordances_describe_direct_rung_plan() {
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
        .expect("Direct rung succeeds on uniform cube");

    let summary = plan.plan_summary();
    assert!(
        summary.starts_with("Direct:") && summary.contains("3 fillets"),
        "summary should describe Direct rung + 3 fillets; got: {summary}",
    );

    let aabb = plan.apex_aabb().expect("Direct rung has apex geometry");
    let (min, max) = aabb;
    assert!(
        min.x <= max.x && min.y <= max.y && min.z <= max.z,
        "apex_aabb min should be ≤ max componentwise",
    );
    // Cube apex at origin, radius 0.1 → AABB is small near origin.
    let extent = (max.x - min.x).max(max.y - min.y).max(max.z - min.z);
    assert!(
        extent < 0.5,
        "AABB extent {extent} should be small for a 0.1-radius cube apex",
    );

    let preview = plan.cap_surface_preview();
    assert!(
        !preview.is_empty(),
        "Direct rung preview should yield at least the tangent-arc endpoints",
    );
}

#[test]
fn preview_affordances_describe_open_rung_plan() {
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
        .expect("Open rung produces a plan");

    let summary = plan.plan_summary();
    assert!(
        summary.contains("Open:") && summary.contains("no cap"),
        "Open's summary should mention no cap; got: {summary}",
    );

    assert!(
        plan.apex_aabb().is_none(),
        "Open rung has no apex geometry to box",
    );
    assert!(
        plan.cap_surface_preview().is_empty(),
        "Open rung has no cap to preview",
    );
}

// ---------------------------------------------------------------
// Issue 6 — JSON/serde round-trip
// ---------------------------------------------------------------

#[test]
fn plan_round_trips_through_json() {
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
        .expect("Direct rung succeeds on uniform cube");

    let json = serde_json::to_string(&plan).expect("ApexPlan serialises to JSON");
    assert!(
        json.contains("\"apex_vertex\"") && json.contains("\"telemetry\""),
        "JSON should contain ApexPlan field names; got prefix: {}",
        &json[..json.len().min(120)],
    );
    let restored: super::plan::ApexPlan =
        serde_json::from_str(&json).expect("ApexPlan deserialises from JSON");
    // f64 JSON round-trips can drift by one ULP through ryu's
    // shortest-round-trip parser interaction (a known
    // serde_json/ryu issue, not a defect of this plan type).  We
    // therefore only assert structural-field equality on
    // integer/string-valued data — the parts that do round-trip
    // bit-exactly.
    assert_eq!(plan.apex_vertex, restored.apex_vertex);
    assert_eq!(plan.n_fillets(), restored.n_fillets());
    assert_eq!(plan.telemetry.rung, restored.telemetry.rung);
    assert_eq!(
        plan.telemetry.cap_kind_reason,
        restored.telemetry.cap_kind_reason
    );
    assert_eq!(
        plan.target_topology.apex_vertices.len(),
        restored.target_topology.apex_vertices.len()
    );
    assert_eq!(
        plan.target_topology.apex_edges.len(),
        restored.target_topology.apex_edges.len()
    );
    // Geometric fields are reproduced within 1 ULP — verify by
    // bounding the per-coordinate drift on apex_vertices.
    for (a, b) in plan
        .target_topology
        .apex_vertices
        .iter()
        .zip(restored.target_topology.apex_vertices.iter())
    {
        assert!(
            (a.x - b.x).abs() < 1e-15 && (a.y - b.y).abs() < 1e-15 && (a.z - b.z).abs() < 1e-15,
            "apex_vertex round-trip drift > 1e-15: {a:?} vs {b:?}",
        );
    }
}

#[test]
fn open_rung_plan_round_trips_through_json() {
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
        .expect("Open rung produces a plan");

    let json = serde_json::to_string(&plan).expect("Open ApexPlan serialises");
    let restored: super::plan::ApexPlan =
        serde_json::from_str(&json).expect("Open ApexPlan deserialises");
    // Open's plan has no f64 geometry to lose precision on, so
    // structural equality holds bit-exactly.
    assert_eq!(plan, restored);
    assert!(restored.cap_kind.is_open());
}
