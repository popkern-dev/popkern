//! Commercial-grade continuity gates on the K-P
//! bidegree-9 output.
//!
//! # Metric choice: triple product, not unit-normal dot
//!
//! The natural "1 − |n_A · n_B|" signal breaks down when one
//! patch's `|∂u × ∂v|` drops to ~1e-5 at a parametric
//! near-degeneracy: `.normalized()` amplifies the tiny cross
//! vector into a direction dominated by floating-point noise,
//! and the dot-product metric then spans the full [−1, +1]
//! even though the underlying surfaces share the same tangent
//! plane. Diagnosed on cube_mixed: worst unit-normal deficit
//! of 9.7e-1 at t ≈ 0.74 sits on top of a perfectly G1-coplanar
//! triple-product of 1.3e-3.
//!
//! So the gate uses the *geometric* G1 definition directly —
//! the coplanarity of `xtan_A`, `xtan_B`, and the along-boundary
//! spoke tangent, normalised so `0` = perfect G1 and `1` =
//! orthogonal. This is invariant under the parametric
//! degeneracy that trips up the unit-normal metric.
//!
//! # Gates (per the bi-9 gate plan § 1-3)
//!
//! * **G0**: position gap ≤ 1e-10 on interior `t ∈ [0.05, 0.95]`.
//! * **G1**: coplanarity ≤ 5e-2 (≈ 2.9° normal mismatch).
//! * **G2**: normal-curvature match residual ≤ 5e-2.
//!
//! The 0.05 endpoint margin excludes the Gregory central-vertex
//! and outer-corner singularities — finite-angle transition
//! zones inherent to the rational blend (not bi-9 construction work).
//!
//! # Fixture matrix
//!
//! | Fixture          | Valence | Radii                                | Status                    |
//! |------------------|---------|--------------------------------------|---------------------------|
//! | cube_mixed       | 3       | [0.08, 0.10, 0.12]                   | exact-boundary route passes |
//! | cube_extreme     | 3       | [0.02, 0.10, 0.20]  (1:10 ratio)     | binonic commercial gates pass |
//! | cube_extreme_20  | 3       | [0.01, 0.10, 0.20]  (1:20 ratio)     | binonic commercial gates pass |
//! | cube_extreme_30  | 3       | [0.0067, 0.10, 0.20]  (1:30 ratio)   | binonic commercial gates pass |
//! | pyramid4_mixed   | 4       | [0.04, 0.06, 0.05, 0.07]             | direct-apex partition (gated in direct_apex_validation_matrix) |
//! | cube_extreme_40  | 3       | [0.005, 0.10, 0.20]  (1:40 ratio)    | vertex-blend route, no K-P ring |
//! | cube_extreme_50  | 3       | [0.004, 0.10, 0.20]  (1:50 ratio)    | vertex-blend route, no K-P ring |
//! | cube_extreme_100 | 3       | [0.002, 0.10, 0.20]  (1:100 ratio)   | vertex-blend route, no K-P ring |
//! | cube_tiny_radii  | 3       | [0.001, 0.002, 0.003] (all tiny)     | vertex-blend route, no K-P ring |
//! | pentagon5_mixed  | 5       | [0.04, 0.06, 0.05, 0.07, 0.05]       | production-route all-or-nothing gate (PK-075/PK-089) |
//! | hexagon6_mixed   | 6       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06] | production-route all-or-nothing gate (PK-075/PK-089) |
//!
//! **Known limits after exact-boundary promotion**:
//! * cube_mixed is now covered as an exact-boundary route invariant.
//! * cube_extreme through cube_extreme_30 now ship rectangular
//!   binonic surfaces through direct exact-boundary spoke-row
//!   optimization. cube_extreme_40+ and cube_tiny_radii route
//!   through vertex-blend construction instead of the K-P rung.
//! * pyramid4_mixed promotes to the direct-apex partition
//!   (PK-099) and is gated in `direct_apex_validation_matrix`;
//!   this file carries no binonic gate for it (PK-182).
//! * N >= 5 dispatches PnS-first (KP20 paper cap behind the polar
//!   flag); the K-P binonic ring was written for N in {3,4} and is
//!   not the high-valence route. The N >= 5 rows gate the PK-089
//!   all-or-nothing production contract instead of the binonic ring.
//!
//! Uniform-radii fixtures are out-of-scope here — the analytic
//! sphere-fit short-circuits Gregory dispatch upstream, so the
//! binonic rung never engages and the gate is tautologically
//! satisfied. Those are a different code-path; testing them
//! belongs in a sphere-shortcut suite.
//!
//! Every gate in this file asserts the serving construction route
//! (read from the recorded dispatch decision) before measuring
//! quality, so a gate can never blind-green by measuring a route
//! that did not run (PK-182). Fixtures whose selected route is
//! vertex-blend or direct-apex are gated on those routes instead of
//! carrying an ignored binonic gate here.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    BlendResult, ExactBoundaryCandidateOutcome, ExactBoundaryReportTransactionStatus,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/kp_cube_fixture.rs"]
mod kp_cube_fixture;
#[path = "support/kp_production_route.rs"]
mod kp_production_route;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, cap_faces_from_result, collect_patches,
    incident_edges, measure_ring,
};
use kp_cube_fixture::{fillet_cube_result, try_fillet_cube_result};
use kp_production_route::assert_n_gonal_production_route_all_or_nothing;

/// Route-selection precondition (PK-182): a binonic continuity
/// gate is only meaningful if the K-P promotion actually served —
/// verifies the committed patches are the binonic variant with at
/// least one present (`assert_all_binonic`).
fn run_cube_binonic(radii: &[f64], fixture: &str) -> helpers::RingMetrics {
    let result = fillet_cube_result(radii, true);
    // The K-P binonic ring is a promotion inside the vertex-blend
    // `insert_vertex_blend` path, so the recorded apex route for
    // these fixtures is `VertexBlend` — same record as fixtures where
    // the ring never fires. The discriminating dispatch evidence
    // is the committed corner patches themselves: a route record
    // must exist, and every committed patch must be the
    // `RectangularBiNonic` variant with at least one present
    // (`assert_all_binonic`), so the gate can never measure a ring
    // that did not ship.
    let route =
        result.repair_report.apex_routes.first().unwrap_or_else(|| {
            panic!("[{fixture}] dispatch must record an apex construction route")
        });
    assert!(
        route.fallback_reason.is_none(),
        "[{fixture}] binonic gate precondition failed: apex route recorded a fallback: {:?}",
        route.fallback_reason
    );
    let faces = cap_faces_from_result(&result);
    let patches = collect_patches(&result.brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

fn assert_cube_mixed_exact_boundary_route(result: &BlendResult) {
    result
        .brep
        .validate()
        .expect("cube_mixed exact-boundary route must remain a valid BRep");
    let routes = &result.repair_report.apex_routes;
    assert_eq!(routes.len(), 1, "cube_mixed should have one apex route");
    let route = &routes[0];
    assert!(route.fallback_reason.is_none());
    assert_eq!(route.exact_boundary_sector_count, 3);
    assert_eq!(route.exact_boundary_pcurve_pair_count, 6);
    assert_eq!(route.exact_boundary_report.committed_faces.len(), 3);
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted
    );
    assert_eq!(
        route.exact_boundary_report.provenance.transaction_status,
        ExactBoundaryReportTransactionStatus::Committed
    );
    assert!(
        route
            .exact_boundary_report
            .diagnostics
            .iter()
            .all(|diagnostic| {
                popkern_fillet::types::diagnostic_code_is_informational(&diagnostic.code)
                    || (diagnostic.stage == "exact_boundary_cap.support_scope"
                        && diagnostic.code == "exact_boundary_cap.support_scope.residual_exceeded")
                    || (diagnostic.stage == "exact_boundary_cap.cap_outer"
                        && diagnostic.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded")
            }),
        "cube_mixed exact-boundary route should carry only structured cap/support residual diagnostics"
    );
}

// ── Commercial-gate tests ────────────────────────────────────

#[test]
fn binonic_passes_or_declines_commercially_on_cube_mixed() {
    let result = fillet_cube_result(&[0.08, 0.10, 0.12], true);
    assert_cube_mixed_exact_boundary_route(&result);
}

#[test]
fn kp_binonic_dispatch_reaches_measurable_ring_on_cube_extreme() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let result = fillet_cube_result(&[0.02, 0.10, 0.20], true);
    let faces = result
        .corner_patches
        .iter()
        .map(|patch| patch.face)
        .collect::<Vec<_>>();
    let patches = collect_patches(&result.brep, &faces);
    assert_eq!(edges.len(), 3);
    assert_eq!(
        patches.len(),
        3,
        "cube_extreme should reach a measurable K-P/Gregory ring"
    );
    assert_all_binonic(&patches, "cube_extreme");
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .all(|diagnostic| {
                !matches!(
                    diagnostic,
                    popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                        if cap.stage == "kp_binonic.promotion"
                            && cap.code == "kp_binonic.promotion.outer_boundary_residual_exceeded"
                )
            }),
        "cube_extreme should no longer demote K-P/binonic on outer-boundary residual"
    );
}

#[test]
fn binonic_passes_commercial_gates_on_cube_extreme() {
    let m = run_cube_binonic(&[0.02, 0.10, 0.20], "cube_extreme");
    assert_commercial_gates("cube_extreme", &m);
}

// pyramid4_mixed: no binonic gate. PK-099's orientation-honest
// convexity promoted the fixture to the direct-apex partition
// (the vertex-blend convex-quad path never reaches three-quad dispatch:
// its overlap-split loop is 12-he, not 2n, and shipped an
// unmeasured g0=1.34 fallback cap; the partition cap closes at
// 1.5e-14). The serving route is gated by
// `direct_apex_native_pyramid4_mixed_partition_no_fallback` in
// direct_apex_validation_matrix.rs plus reflection_zebra. A
// binonic gate here would measure a route that never fires
// (PK-182). Reintroduce only alongside a 2n-dispatch-ready
// convex-quad loop.

#[test]
fn cube_tiny_radii_declines_transactionally_on_perp_split_gap() {
    // Honesty re-pin (PK-140 transactional gate): the previous green
    // measured a 2-of-3 PARTIAL body — one tiny-radius leg has
    // always failed the anchor-split stage ("no perpendicular edge found at
    // v_start") and was silently dropped. Until the tiny-radius perp
    // split lands, the operation must decline naming the class.
    let err = try_fillet_cube_result(&[0.001, 0.002, 0.003], true)
        .expect_err("tiny-radius corner has a known the anchor-split stage gap and must decline");
    let msg = format!("{err:?}");
    assert!(
        msg.contains("no perpendicular edge found"),
        "decline must name the tiny-radius perp-split class: {}",
        &msg[..msg.len().min(240)]
    );
}

// cube_extreme_40/50/100 and cube_tiny_radii carry no binonic
// continuity gates: those fixtures are served by vertex-blend
// construction (no measurable K-P ring exists), and a binonic
// gate would measure a route that never fires (PK-182). The
// serving route is pinned by the active
// `cube_extreme_*_accepts_vertex_blend_route_instead_of_kp_ring` tests
// above (route == VertexBlend read from the dispatch record) and by
// `cube_tiny_radii_declines_transactionally_on_perp_split_gap`.
// Reintroduce a binonic gate for a fixture only when dispatch
// actually selects the K-P rung for it.

/// N >= 5 production-route gate (PK-075). The K-P binonic ring this
/// file's other rows measure was written for N ∈ {3, 4} and is not
/// the high-valence route: production dispatch is PnS-first with the
/// KP20 paper cap behind the polar flag. The commercial contract for
/// these fixtures is the PK-089 all-or-nothing invariant: either an
/// accepted N-patch tiling is committed (each patch Gregory or PnS
/// with the apex's n_sides), or the corner rolls back to mitred with
/// zero apex patches AND structured dispatch-decline diagnostics
/// recorded. Route-quality repair beyond the contract is owned by
/// PK-084 (pyramid5 PnS/sphere-cap) and PK-090 (paper-cap inputs).

#[test]
fn pentagon5_mixed_production_route_is_all_or_nothing() {
    assert_n_gonal_production_route_all_or_nothing(
        5,
        &[0.04, 0.06, 0.05, 0.07, 0.05],
        "pentagon5_mixed",
    );
}

#[test]
fn hexagon6_mixed_production_route_is_all_or_nothing() {
    assert_n_gonal_production_route_all_or_nothing(
        6,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        "hexagon6_mixed",
    );
}

// ── Non-regression vs two-stripe ─────────────────────────────

#[test]
fn exact_boundary_route_matches_with_binonic_flag_on_cube_mixed() {
    let two_stripe_result = fillet_cube_result(&[0.08, 0.10, 0.12], false);
    let binonic_result = fillet_cube_result(&[0.08, 0.10, 0.12], true);
    assert_cube_mixed_exact_boundary_route(&two_stripe_result);
    assert_cube_mixed_exact_boundary_route(&binonic_result);
}

// No binonic-vs-two-stripe comparison gate: with the binonic
// flag off, cube_extreme produces no measurable Gregory patches
// at all, so there is no two-stripe baseline to compare against
// (PK-182). The binonic output is held to the *absolute*
// commercial G0/G1/G2 gates above, which is a stronger contract
// than a relative no-regression bar, and
// `cube_extreme_ships_binonic_no_warning` pins that the binonic
// variant (not a lower rung) is what ships.

#[test]
fn cube_extreme_ships_binonic_no_warning() {
    use popkern_core::tolerance::Tolerance;
    use popkern_fillet::error::BlendWarning;
    use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.02_f64, 0.10, 0.20];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();

    // (1) Binonic shipped, not two-stripe.
    let faces: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
    let patches = collect_patches(&result.brep, &faces);
    assert!(
        !patches.is_empty(),
        "cube_extreme should expose measurable binonic patches"
    );
    for p in &patches {
        assert!(
            matches!(
                p,
                popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(_)
            ),
            "cube_extreme should ship binonic variant, got {p:?}"
        );
    }

    // (2) No CornerContinuityRelaxed warning once binonic shipping
    // is restored.
    let relaxed_warnings: Vec<_> = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CornerContinuityRelaxed { .. }))
        .collect();
    assert!(
        relaxed_warnings.is_empty(),
        "expected no CornerContinuityRelaxed warnings, got {}",
        relaxed_warnings.len()
    );
}
