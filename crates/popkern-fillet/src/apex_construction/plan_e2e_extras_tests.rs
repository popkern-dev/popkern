//! End-to-end planner tests for the affordances added in the
//! A2.1 review (Issues 5, 6, 8) plus variable-radius coverage.
//!
//! Split out of [`super::plan_e2e_tests`] so each test file stays
//! under the 500-line structural cap.  Helpers live in
//! [`super::plan_e2e_test_helpers`].

use popkern_core::topology::BRep;

use crate::types::{FilletRadius, FilletSpec};

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan::plan_apex_resolution;
use super::plan_e2e_test_helpers::{
    build_fillet_surfaces, build_fillet_surfaces_per_edge, force_unsupported_kind,
    pick_apex_and_incident_edges, test_tolerance,
};
use super::telemetry::ConstructionRung;

// ---------------------------------------------------------------
// Issue 5 — cancellation
// ---------------------------------------------------------------

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

// ---------------------------------------------------------------
// Variable-radius coverage — Gregory dispatch on real fixture
// ---------------------------------------------------------------

#[test]
fn variable_radii_on_cube_dispatches_to_gregory() {
    let cube = BRep::unit_cube(test_tolerance()).expect("cube builds");
    let (apex, edges) = pick_apex_and_incident_edges(&cube);
    assert_eq!(edges.len(), 3);
    // Three distinct radii spanning ~10× — well outside any
    // tolerance-based "uniform" classification.
    let radii = [0.05_f64, 0.10, 0.20];
    let surfaces = build_fillet_surfaces_per_edge(&cube, &edges, &radii);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();

    let plan = plan_apex_resolution(&cube, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Direct rung succeeds on variable-radius cube apex");

    assert_eq!(
        plan.telemetry.rung,
        ConstructionRung::Direct,
        "variable radii still resolve via Direct rung",
    );
    // Variable radii must dispatch to Gregory (Sphere requires
    // uniform radii within tolerance.linear()).
    assert!(
        plan.cap_kind.is_gregory(),
        "variable radii (0.05, 0.10, 0.20) → Gregory cap, got {:?}",
        plan.cap_kind,
    );
    if let CapKind::Gregory {
        per_fillet_caps, ..
    } = &plan.cap_kind
    {
        assert_eq!(per_fillet_caps.len(), 3);
        // Each per-fillet cap radius matches its input radius.
        // Caps are in cap-loop walk order, which may not match
        // the input edge order — check by radius set.
        let mut got: Vec<f64> = per_fillet_caps.iter().map(|(_, r)| *r).collect();
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut want = radii.to_vec();
        want.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (g, w) in got.iter().zip(want.iter()) {
            assert!(
                (g - w).abs() < 1e-9,
                "per-fillet radius {g} should match an input radius (sorted: got {got:?}, want {want:?})",
            );
        }
    }
    // Telemetry should attribute the Gregory dispatch to variable
    // radii, not bridge-collapse failure.
    let reason = &plan.telemetry.cap_kind_reason;
    assert!(
        reason.contains("variable radii"),
        "Gregory reason should mention variable radii; got: {reason}",
    );
    // Diagnostic numerics: radii_spread should be ≥ 0.15 (max 0.20 - min 0.05).
    let spread = plan
        .telemetry
        .diagnostic_values
        .get("radii_spread")
        .copied();
    assert!(
        spread.is_some_and(|s| (s - 0.15).abs() < 1e-9),
        "radii_spread should be 0.15 (= 0.20 - 0.05); got {spread:?}",
    );
}
