//! End-to-end planner tests covering the rung-ladder dispatch
//! against real BRep fixtures (unit cube + square pyramid) with
//! fillet surfaces produced by the actual
//! [`compute_fillet_geometry`] pipeline.
//!
//! Issues 2 + 3 from the DirectApex.1 review:
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
//!
//! Cancellation / preview / serde / variable-radius tests live in
//! [`super::plan_e2e_extras_tests`].  Shared helpers live in
//! [`super::plan_e2e_test_helpers`].

use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::types::{FilletRadius, FilletSpec};

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan::plan_apex_resolution;
use super::plan_e2e_test_helpers::{
    build_fillet_surfaces, force_unsupported_kind, pick_apex_and_incident_edges, test_tolerance,
};
use super::telemetry::ConstructionRung;

// ---------------------------------------------------------------
// Issue 2 — Direct rung end-to-end
// ---------------------------------------------------------------

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
    // Uniform-radius cube apex: per-fillet sphere centres all
    // coincide at (r, r, r), bridges collapse, dispatch lands on
    // Sphere cap.  This is the commercial-grade outcome — the
    // earlier origin-based side check used to spuriously reject
    // one of the three inscribed centres and force Gregory.
    assert!(
        plan.cap_kind.is_sphere(),
        "uniform cube → Sphere cap, got {:?}",
        plan.cap_kind,
    );
    if let CapKind::Sphere {
        radius: r, center, ..
    } = &plan.cap_kind
    {
        assert!(
            (r - radius).abs() < 1e-9,
            "Sphere cap radius {r} should match fillet radius {radius}",
        );
        // Cap centre at (r, r, r) — the diagonal offset from the apex.
        assert!(
            (center.x - radius).abs() < 1e-9
                && (center.y - radius).abs() < 1e-9
                && (center.z - radius).abs() < 1e-9,
            "Sphere cap centre should be at (r,r,r) for cube corner; got {center:?}",
        );
    }
    // Sphere cap has 3 apex vertices (one per fillet's tangent
    // arc endpoint, shared with the next fillet's start).
    assert_eq!(
        plan.target_topology.apex_vertices.len(),
        3,
        "Sphere cap on cube apex has 3 apex vertices",
    );
}

/// Re-pointed with PK-181: the "mirror-era baseline" ignore was
/// stale — the strict Direct pass succeeds on real inscribed
/// strips at N=4 (measured 2026-07-28).
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

/// Re-pointed with PK-181 to a hexagonal pyramid — the only
/// end-to-end Setback-rung coverage in the repo. Measured on real
/// inscribed strips (2026-07-28): forcing one fillet to an
/// unsupported Plane also fails sphere-fit for its two RING
/// NEIGHBOURS (the forced surface corrupts the trim curves
/// bracketing their boundary segments → NoIntersection), so a
/// forced apex loses 3 terminations, not 1. N=4/5 therefore
/// bottom out at Open (< 3 survivors); N=6 leaves exactly the 3
/// survivors `plan_setback_rung` requires.
#[test]
fn setback_rung_succeeds_when_one_pyramid_fillet_unsupported() {
    let pyramid =
        BRep::n_gonal_pyramid(6, 1.0, 0.8, test_tolerance()).expect("hexagon pyramid builds");
    let (apex, edges) = {
        let mut best: Option<(VertexHandle, Vec<EdgeHandle>)> = None;
        for v in pyramid.vertices.keys() {
            let mut es: std::collections::BTreeSet<EdgeHandle> = std::collections::BTreeSet::new();
            for he in pyramid.half_edges.values() {
                if he.origin == v {
                    es.insert(he.edge);
                }
            }
            if es.len() == 6 {
                best = Some((v, es.into_iter().collect()));
                break;
            }
        }
        best.expect("pyramid has a 6-valence apex")
    };
    let radius = 0.05;
    let mut surfaces = build_fillet_surfaces(&pyramid, &edges, radius);
    // Force one fillet to be Plane (unsupported) — strict pass
    // fails; the lenient pass skips it AND loses its two ring
    // neighbours to the corrupted boundary brackets, leaving 3
    // successful fillets → Setback rung succeeds at its minimum.
    force_unsupported_kind(&mut surfaces, 1);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();

    let plan = plan_apex_resolution(&pyramid, apex, &specs, &surfaces, &test_tolerance(), None)
        .expect("Setback rung should succeed with 3 surviving fillets of 6");

    assert_eq!(plan.apex_vertex, apex);
    assert!(
        matches!(plan.telemetry.rung, ConstructionRung::Setback { .. }),
        "expected Setback rung, got {:?}; telemetry: {:?}",
        plan.telemetry.rung,
        plan.telemetry,
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
    // Setback keeps only the 3 surviving fillets in the cap (the
    // forced fillet plus its two boundary-corrupted neighbours are
    // omitted — see the doc comment).
    assert_eq!(
        plan.n_fillets(),
        3,
        "Setback must omit the unsupported fillet and its two affected neighbours",
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
    if let CapKind::Open {
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
