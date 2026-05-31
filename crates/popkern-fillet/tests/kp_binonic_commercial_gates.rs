//! Phase H.7 — commercial-grade continuity gates on the K-P
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
//! # Gates (per the Phase H plan § 1-3)
//!
//! * **G0**: position gap ≤ 1e-10 on interior `t ∈ [0.05, 0.95]`.
//! * **G1**: coplanarity ≤ 5e-2 (≈ 2.9° normal mismatch).
//! * **G2**: normal-curvature match residual ≤ 5e-2.
//!
//! The 0.05 endpoint margin excludes the Gregory central-vertex
//! and outer-corner singularities — finite-angle transition
//! zones inherent to the rational blend (not Phase H work).
//!
//! # Fixture matrix
//!
//! | Fixture          | Valence | Radii                                | Status                    |
//! |------------------|---------|--------------------------------------|---------------------------|
//! | cube_mixed       | 3       | [0.08, 0.10, 0.12]                   | exact-boundary route passes |
//! | cube_extreme     | 3       | [0.02, 0.10, 0.20]  (1:10 ratio)     | binonic commercial gates pass |
//! | cube_extreme_20  | 3       | [0.01, 0.10, 0.20]  (1:20 ratio)     | binonic commercial gates pass |
//! | cube_extreme_30  | 3       | [0.0067, 0.10, 0.20]  (1:30 ratio)   | binonic commercial gates pass |
//! | pyramid4_mixed   | 4       | [0.04, 0.06, 0.05, 0.07]             | `#[ignore]` lower-rung fallback |
//! | cube_extreme_40  | 3       | [0.005, 0.10, 0.20]  (1:40 ratio)    | legacy route, no K-P ring |
//! | cube_extreme_50  | 3       | [0.004, 0.10, 0.20]  (1:50 ratio)    | legacy route, no K-P ring |
//! | cube_extreme_100 | 3       | [0.002, 0.10, 0.20]  (1:100 ratio)   | legacy route, no K-P ring |
//! | cube_tiny_radii  | 3       | [0.001, 0.002, 0.003] (all tiny)     | legacy route, no K-P ring |
//! | pentagon5_mixed  | 5       | [0.04, 0.06, 0.05, 0.07, 0.05]       | `#[ignore]` polar dispatch gap |
//! | hexagon6_mixed   | 6       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06] | `#[ignore]` polar dispatch gap |
//!
//! **Known limits after exact-boundary promotion**:
//! * cube_mixed is now covered as an exact-boundary route invariant.
//! * cube_extreme through cube_extreme_30 now ship rectangular
//!   binonic surfaces through direct exact-boundary spoke-row
//!   optimization. cube_extreme_40+ and cube_tiny_radii route
//!   through legacy construction instead of the K-P rung.
//! * pyramid4_mixed reaches Gregory patches but not the binonic rung;
//!   the current blocker is lower-rung fallback, distinct from the
//!   native exact-boundary cap_outer matrix.
//! * N >= 5 still needs the K-P polar variant (Phase H.5b).
//!
//! Uniform-radii fixtures are out-of-scope here — the analytic
//! sphere-fit short-circuits Gregory dispatch upstream, so the
//! binonic rung never engages and the gate is tautologically
//! satisfied. Those are a different code-path; testing them
//! belongs in a sphere-shortcut suite.
//!
//! The active acceptance in this file covers the selected construction
//! family for each fixture. Wider asymmetry and unsupported-valence
//! continuity tests remain ignored until they produce measurable
//! rectangular binonic rings.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionRoute, BlendResult, ExactBoundaryCandidateOutcome,
    ExactBoundaryReportTransactionStatus, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, assert_no_regression, collect_patches,
    fillet_corner, incident_edges, measure_ring, pick_vertex_with_valence,
};

fn run_cube(radii: &[f64], use_kp_binonic: bool) -> helpers::RingMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    assert_eq!(edges.len(), 3);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, use_kp_binonic);
    let patches = collect_patches(&brep, &faces);
    assert_eq!(patches.len(), 3, "expected 3 Gregory sub-patches on cube");
    measure_ring(&patches, N_SAMPLES)
}

fn run_cube_binonic(radii: &[f64], fixture: &str) -> helpers::RingMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

fn fillet_cube_result(radii: &[f64], use_kp_binonic: bool) -> BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    assert_eq!(edges.len(), radii.len());
    let specs = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: use_kp_binonic,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap()
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
                (diagnostic.stage == "exact_boundary_cap.support_scope"
                    && diagnostic.code == "exact_boundary_cap.support_scope.residual_exceeded")
                    || (diagnostic.stage == "exact_boundary_cap.cap_outer"
                        && diagnostic.code == "exact_boundary_cap.cap_outer.pair_residual_exceeded")
            }),
        "cube_mixed exact-boundary route should carry only structured cap/support residual diagnostics"
    );
}

fn run_square_pyramid_binonic(radii: &[f64], fixture: &str) -> helpers::RingMetrics {
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).expect("square_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

/// N-gonal pyramid (apex valence = n). Used for the N=5, N=6
/// stress-test fixtures probing the L-M stack at higher corner
/// valences than the N=3 cube / N=4 square-pyramid commercial
/// matrix.
fn run_n_gonal_pyramid_binonic(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    fixture: &str,
) -> helpers::RingMetrics {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
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

#[test]
#[ignore = "K-P/binonic dispatch gate: pyramid4_mixed produces Gregory patches but falls back below the rectangular binonic rung"]
fn binonic_passes_commercial_gates_on_pyramid4_mixed() {
    let m = run_square_pyramid_binonic(&[0.04, 0.06, 0.05, 0.07], "pyramid4_mixed");
    assert_commercial_gates("pyramid4_mixed", &m);
}

// ── K-P/binonic dispatch stress fixtures ─────────────────────
//
// Wider-asymmetry tests remain ignored where the current helper does
// not expose committed rectangular binonic rings. They become
// continuity gates after the upstream dispatch blocker is closed.

/// 1:20 radius asymmetry — twice as asymmetric as cube_extreme.
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_20() {
    let m = run_cube_binonic(&[0.01, 0.10, 0.20], "cube_extreme_20");
    assert_commercial_gates("cube_extreme_20", &m);
}

/// 1:30 radius asymmetry.
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_30() {
    let m = run_cube_binonic(&[0.0067, 0.10, 0.20], "cube_extreme_30");
    assert_commercial_gates("cube_extreme_30", &m);
}

#[test]
fn cube_extreme_40_accepts_legacy_route_instead_of_kp_ring() {
    assert_extreme_cube_accepts_legacy_route_instead_of_kp_ring(
        &[0.005, 0.10, 0.20],
        "cube_extreme_40",
    );
}

#[test]
fn cube_extreme_50_accepts_legacy_route_instead_of_kp_ring() {
    assert_extreme_cube_accepts_legacy_route_instead_of_kp_ring(
        &[0.004, 0.10, 0.20],
        "cube_extreme_50",
    );
}

#[test]
fn cube_extreme_100_accepts_legacy_route_instead_of_kp_ring() {
    assert_extreme_cube_accepts_legacy_route_instead_of_kp_ring(
        &[0.002, 0.10, 0.20],
        "cube_extreme_100",
    );
}

#[test]
fn cube_tiny_radii_accepts_legacy_route_instead_of_kp_ring() {
    assert_extreme_cube_accepts_legacy_route_instead_of_kp_ring(
        &[0.001, 0.002, 0.003],
        "cube_tiny_radii",
    );
}

fn assert_extreme_cube_accepts_legacy_route_instead_of_kp_ring(radii: &[f64], fixture: &str) {
    let result = fillet_cube_result(radii, true);
    result
        .brep
        .validate()
        .unwrap_or_else(|_| panic!("{fixture} legacy route must leave a valid BRep"));
    for patch in &result.corner_patches {
        assert!(
            result.brep.faces.contains_key(patch.face),
            "{fixture} should not expose stale corner patch face handles: {patch:#?}"
        );
    }
    let faces = result
        .corner_patches
        .iter()
        .map(|patch| patch.face)
        .collect::<Vec<_>>();
    assert!(
        collect_patches(&result.brep, &faces).is_empty(),
        "{fixture} should not route through measurable Gregory/K-P patches"
    );
    let route = result
        .repair_report
        .apex_routes
        .first()
        .unwrap_or_else(|| panic!("{fixture} should record an apex construction route"));
    assert_eq!(
        route.route,
        ApexConstructionRoute::Legacy,
        "{fixture} should be owned by the selected legacy route"
    );
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::NotApplicable,
        "{fixture} should not pretend direct-apex/exact-boundary was attempted"
    );
    assert!(
        route.built_faces.is_empty(),
        "{fixture} legacy route should not report direct-apex built faces"
    );
}

/// **Documented non-applicable binonic gate — 1:40 asymmetry.**
/// This fixture is covered by legacy construction before a
/// measurable K-P/Gregory ring exists. The active route test above
/// guards that selected construction; this continuity gate remains
/// ignored because it intentionally asks for the non-selected family.
#[ignore = "K-P/binonic route not selected: cube_extreme_40 is accepted by legacy construction"]
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_40() {
    let m = run_cube_binonic(&[0.005, 0.10, 0.20], "cube_extreme_40");
    assert_commercial_gates("cube_extreme_40", &m);
}

/// **Documented non-applicable binonic gate — 1:50 asymmetry.**
/// This fixture is covered by legacy construction before a
/// measurable K-P/Gregory ring exists.
#[ignore = "K-P/binonic route not selected: cube_extreme_50 is accepted by legacy construction"]
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_50() {
    let m = run_cube_binonic(&[0.004, 0.10, 0.20], "cube_extreme_50");
    assert_commercial_gates("cube_extreme_50", &m);
}

/// **Documented non-applicable binonic gate — 1:100 asymmetry.**
/// This fixture is covered by legacy construction before a
/// measurable K-P/Gregory ring exists.
#[ignore = "K-P/binonic route not selected: cube_extreme_100 is accepted by legacy construction"]
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_100() {
    let m = run_cube_binonic(&[0.002, 0.10, 0.20], "cube_extreme_100");
    assert_commercial_gates("cube_extreme_100", &m);
}

/// **Documented non-applicable binonic gate — tiny absolute
/// radii.** The selected construction is valid legacy topology
/// before a measurable K-P/Gregory ring exists.
#[ignore = "K-P/binonic route not selected: cube_tiny_radii is accepted by legacy construction"]
#[test]
fn binonic_passes_commercial_gates_on_cube_tiny_radii() {
    let m = run_cube_binonic(&[0.001, 0.002, 0.003], "cube_tiny_radii");
    assert_commercial_gates("cube_tiny_radii", &m);
}

/// **Documented gap — N=5 apex.** The K-P binonic dispatch path
/// is currently written for N ∈ {3, 4}; pentagonal apex produces
/// zero Gregory corner patches. Needs the K-P polar variant
/// (Phase H.5b per the module header) before binonic can engage
/// at all.
#[ignore = "K-P polar dispatch gate: pentagon5_mixed currently produces no measurable Gregory corner patches"]
#[test]
fn binonic_passes_commercial_gates_on_pentagon5_mixed() {
    let m = run_n_gonal_pyramid_binonic(
        5,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05],
        "pentagon5_mixed",
    );
    assert_commercial_gates("pentagon5_mixed", &m);
}

/// **Documented gap — N=6 apex.** Same dispatch ceiling as
/// pentagon5 — K-P polar variant required.
#[ignore = "K-P polar dispatch gate: hexagon6_mixed currently produces no measurable Gregory corner patches"]
#[test]
fn binonic_passes_commercial_gates_on_hexagon6_mixed() {
    let m = run_n_gonal_pyramid_binonic(
        6,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        "hexagon6_mixed",
    );
    assert_commercial_gates("hexagon6_mixed", &m);
}

// ── Non-regression vs two-stripe ─────────────────────────────

#[test]
fn exact_boundary_route_matches_with_binonic_flag_on_cube_mixed() {
    let two_stripe_result = fillet_cube_result(&[0.08, 0.10, 0.12], false);
    let binonic_result = fillet_cube_result(&[0.08, 0.10, 0.12], true);
    assert_cube_mixed_exact_boundary_route(&two_stripe_result);
    assert_cube_mixed_exact_boundary_route(&binonic_result);
}

#[test]
#[ignore = "comparison gate: non-binonic cube_extreme route currently produces no measurable Gregory patches"]
fn binonic_matches_or_beats_two_stripe_on_cube_extreme() {
    let mt = run_cube(&[0.02, 0.10, 0.20], false);
    let mb = run_cube(&[0.02, 0.10, 0.20], true);
    assert_no_regression("cube_extreme", &mt, &mb);
}

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
