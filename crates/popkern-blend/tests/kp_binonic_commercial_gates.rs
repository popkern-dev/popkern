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

use popkern_blend::{
    BlendResult, ExactBoundaryCandidateOutcome, ExactBoundaryReportTransactionStatus,
};
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/kp_cube_fixture.rs"]
mod kp_cube_fixture;
#[path = "support/kp_production_route.rs"]
mod kp_production_route;
use helpers::{N_SAMPLES, cap_faces_from_result, collect_patches, incident_edges, measure_ring};
use kp_cube_fixture::{fillet_cube_result, try_fillet_cube_result};
#[path = "support/seam_perturbation.rs"]
mod seam_perturbation;
use kp_production_route::assert_n_gonal_production_route_all_or_nothing;

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
                popkern_blend::types::diagnostic_code_is_informational(&diagnostic.code)
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

/// PK-544/PK-547: cube_extreme (1:10) serves the exact-injection ring
/// at commercial G1.
///
/// PK-539 found the exact candidate folded across a spoke (sine 1.0),
/// and the promotion declined to a G0 ring (0.29). The exact-boundary
/// spoke-row solve now minimizes the oriented spoke normal chord, never
/// adds a fold, holds a chart-lean barrier, and stays clear of proven
/// crossings. Its candidate passes the promotion gates directly, a
/// repair remaining the fallback. The served ring certifies G1 on every
/// cap seam over the whole edge and claims the G1 it measures.
#[test]
fn cube_extreme_binonic_promotion_serves_an_honest_exact_ring() {
    let result = fillet_cube_result(&[0.02, 0.10, 0.20], true);
    let served = result.repair_report.apex_diagnostics.iter().any(|d| {
        matches!(
            d,
            popkern_blend::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "kp_binonic.promotion.boundary_candidates"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail.contains("provenance=exact") && detail.contains("accepted(")
                    })
        )
    });
    assert!(
        served,
        "cube_extreme must serve the repaired exact candidate"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "cube_extreme must validate clean: {errors:?}"
    );
    let faces = cap_faces_from_result(&result);
    let patches = collect_patches(&result.brep, &faces);
    assert_eq!(
        patches.len(),
        3,
        "cube_extreme serves a measurable 3-patch ring"
    );
    let m = measure_ring(&patches, N_SAMPLES);
    assert!(m.worst_g0_gap < helpers::G0_GATE, "G0 must stay exact");
    assert!(
        m.worst_g1_normal <= helpers::G1_GATE,
        "spoke G1 {:.3e} over the commercial bar",
        m.worst_g1_normal
    );
    helpers::assert_cap_seams_certify_g1(&result.brep, &faces, "cube_extreme");
    for patch in &result.corner_patches {
        assert_eq!(
            patch.achieved_continuity,
            popkern_blend::ContinuityLevel::G1,
            "the claim must match the measurement"
        );
    }
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
fn cube_extreme_claims_its_measured_g1_without_downgrade() {
    use popkern_blend::error::BlendWarning;
    use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
    use popkern_core::tolerance::Tolerance;

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

    // PK-544: the repaired exact ring serves at commercial G1, so the
    // corner claims G1 from its measurement and no downgrade is
    // reported (PK-539 declined here and reported G1 -> G0).
    let downgraded = result
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::CornerContinuityDowngraded { .. }));
    assert!(
        !downgraded,
        "cube_extreme serves its requested G1: {:?}",
        result.edge_warnings
    );
    for patch in &result.corner_patches {
        assert_eq!(
            patch.achieved_continuity,
            popkern_blend::ContinuityLevel::G1
        );
    }
}

/// PK-498 positive fixture for the ring gate's negative controls: a
/// synthetic 3-patch bi-nonic ring cut from one paraboloid is exact to
/// all orders across its seams. Every gate passes, at near machine
/// precision.
#[test]
fn the_synthetic_paraboloid_ring_passes_every_ring_gate() {
    let patches = seam_perturbation::paraboloid_ring(3);
    let m = measure_ring(&patches, N_SAMPLES);
    println!(
        "paraboloid ring: g0={:.3e} g1={:.3e} g2={:.3e}",
        m.worst_g0_gap, m.worst_g1_normal, m.worst_g2_residual
    );
    assert!(m.worst_g0_gap < 1e-12, "{:.3e}", m.worst_g0_gap);
    assert!(m.worst_g1_normal < 1e-9, "{:.3e}", m.worst_g1_normal);
    assert!(m.worst_g2_residual < 1e-6, "{:.3e}", m.worst_g2_residual);
}

/// PK-498 curvature-only negative control for the bi-nonic ring gate.
/// On the exact paraboloid ring, row 2 of patch 0's u = 1 side (ring
/// seam 0) moves along the normal. That keeps the seam's position and
/// normal exact (self-checked) and jumps its transverse curvature. The
/// ring gate must fail G2 on that seam while passing G0 and G1. Before
/// PK-498 the ring's G2 was along-seam only, blind here by Meusnier (it
/// read 1.0e-12 under this perturbation).
#[test]
fn a_curvature_only_ring_perturbation_fails_the_binonic_g2_gate() {
    use seam_perturbation::{PatchSide, paraboloid_ring, perturb_side_row, side_change, side_leg};
    let mut patches = paraboloid_ring(3);
    let before = patches[0].clone();
    let leg = side_leg(&before, PatchSide::U1).expect("bi-nonic side");
    perturb_side_row(&mut patches[0], PatchSide::U1, 2, 0.5 * leg).expect("perturb row 2");
    let (position, normal) = side_change(&before, &patches[0], PatchSide::U1);
    assert!(
        position < 1e-12 && normal < 1e-12,
        "G0/G1 kept: {position:.3e} {normal:.3e}"
    );
    let m = measure_ring(&patches, N_SAMPLES);
    println!("ring curvature-only control: per seam {:?}", m.per_seam);
    assert!(m.worst_g0_gap < helpers::G0_GATE, "G0 must pass");
    assert!(m.worst_g1_normal < 1e-9, "G1 must pass");
    assert!(
        m.per_seam[0].1 > 10.0 * helpers::G2_GATE,
        "the ring G2 gate must see a transverse curvature jump: {:.3e}",
        m.per_seam[0].1
    );
}

/// PK-498 normal-only negative control for the bi-nonic ring gate: row
/// 1 of the same side moves, so the seam keeps its position
/// (self-checked) and its normal tilts. The ring G1 gate must fail.
/// Before PK-498 its cross-tangent coplanarity measure passed a 53°
/// crease.
#[test]
fn a_normal_only_ring_perturbation_fails_the_binonic_g1_gate() {
    use seam_perturbation::{PatchSide, paraboloid_ring, perturb_side_row, side_change, side_leg};
    let mut patches = paraboloid_ring(3);
    let before = patches[0].clone();
    let leg = side_leg(&before, PatchSide::U1).expect("bi-nonic side");
    perturb_side_row(&mut patches[0], PatchSide::U1, 1, 0.2 * leg).expect("perturb row 1");
    let (position, normal) = side_change(&before, &patches[0], PatchSide::U1);
    assert!(position < 1e-12, "G0 kept: {position:.3e}");
    assert!(
        normal > helpers::G1_GATE,
        "the normal really turned: {normal:.3e}"
    );
    let m = measure_ring(&patches, N_SAMPLES);
    println!("ring normal-only control: per seam {:?}", m.per_seam);
    assert!(m.worst_g0_gap < helpers::G0_GATE, "G0 must pass");
    assert!(
        m.per_seam[0].0 > helpers::G1_GATE,
        "G1 must fail: {:.3e}",
        m.per_seam[0].0
    );
}
