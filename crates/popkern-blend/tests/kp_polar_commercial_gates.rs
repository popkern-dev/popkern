//! Commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (bi-9, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

mod kp_gate_helpers;
use kp_gate_helpers as helpers;

#[path = "support/kp_polar_fixtures.rs"]
mod fixtures;

use fixtures::{fillet_n_gonal_pyramid_result, fillet_n_gonal_pyramid_result_with_continuity};
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates_live, collect_patches,
    fillet_corner_full, incident_edges, measure_live_cap_seams, pick_vertex_with_valence,
};

/// Fillet the N-gonal pyramid apex and measure the committed cap
/// across its LIVE seam topology (the setback fan ships 4 bi-9
/// pieces per sector; the one-patch-per-sector `measure_ring` frame does not
/// describe that layout — PK-087).
fn run_n_gonal_pyramid_polar(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    fixture: &str,
) -> helpers::LiveSeamMetrics {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner_full(&pyr, &edges, radii, true, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_live_cap_seams(&brep, &faces, N_SAMPLES)
}

/// Blunt-apex positive control (apex half-angle 18°, typical CAD):
/// pentagon5 at h=3.0 is the REAL commercial-CAD case; the sharp-apex
/// fixtures below stress the polar/setback path specifically. A
/// gentle corner may legitimately ship a sphere-fit cap (single
/// non-Gregory face) instead of a fan, so the assertion is "valid
/// output, and IF a multi-patch cap shipped its live seams pass the
/// commercial gates" rather than a specific patch census.
#[test]
fn pentagon5_blunt_ships_valid_commercial_cap() {
    let result = fillet_n_gonal_pyramid_result(5, 1.0, 3.0, &[0.05; 5]);
    result
        .brep
        .validate()
        .expect("pentagon5_blunt must leave a valid BRep");
    let faces: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.n_sides == 5)
        .map(|p| p.face)
        .collect();
    let m = measure_live_cap_seams(&result.brep, &faces, N_SAMPLES);
    if m.seam_count > 0 {
        // G0/G1 commercial gates only: this control requests the
        // default G1 fillet class, so the cap's seam G2 is recorded
        // but not gated (measured at re-pin 2026-07-15: g0=1.6e-13,
        // g1=2.5e-2, g2=7.0e-2 — the G2 jump across seams of a
        // G1-request blend on the wide flat blunt fan).
        println!(
            "live-seam pentagon5_blunt: seams={} g0={:.3e} g1={:.3e} g2={:.3e}",
            m.seam_count, m.worst_g0_gap, m.worst_g1_normal, m.worst_g2_residual
        );
        assert!(
            m.worst_g0_gap < helpers::G0_GATE,
            "[pentagon5_blunt] live-seam G0 gap {:.3e} exceeds {:.1e} gate",
            m.worst_g0_gap,
            helpers::G0_GATE
        );
        assert!(
            m.worst_g1_normal < helpers::G1_GATE,
            "[pentagon5_blunt] live-seam G1 misalignment {:.3e} exceeds {:.1e} gate",
            m.worst_g1_normal,
            helpers::G1_GATE
        );
    } else {
        println!(
            "pentagon5_blunt: single-face cap ({} apex patches), no live seams",
            faces.len()
        );
    }
}
/// Blunt-apex G2-request acceptance (PK-122 → PK-187): an
/// 18°-half-angle N-gonal pyramid requesting G2 must ACCEPT — commit
/// a fan claiming achieved G2 that measures within the commercial
/// gates on its live seams, including the TRANSVERSE curvature
/// residual (the honest G2 metric; the along-seam residual is
/// near-vacuous once G0+G1 hold — Meusnier). History this pins:
/// PK-122 recorded a "7.0e-2 G2 jump" that was actually the
/// along-seam echo of the then-large G1 error (2.5e-2); once the fan
/// rails were fixed the request committed — but the honest transverse
/// measurement then exposed a REAL 5.3e-1 curvature crease authored
/// by the rim solve's independently-fitted G2-row deltas, fixed by
/// pinning the k=2/k=7 columns to the shared corner translation
/// (setback_rim_reconcile). A regression to the decline path OR to a
/// dishonest claim fails this test.
///
/// Watertightness is pinned STRICT for both requests on both
/// fixtures — the PK-208 crack family (measured 236 welded open
/// edges on the pre-fix pentagon G1 control) is closed; see the
/// assertion comment for the three fixed authors.
fn blunt_g2_request_accepts(n: usize, fixture: &str) {
    let radii = vec![0.05; n];
    let result = fillet_n_gonal_pyramid_result_with_continuity(
        n,
        1.0,
        3.0,
        &radii,
        popkern_blend::ContinuityLevel::G2,
    );
    result
        .brep
        .validate()
        .unwrap_or_else(|e| panic!("[{fixture}] G2 request must leave a valid BRep: {e:?}"));
    let g2_claimed: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.n_sides == n && p.achieved_continuity == popkern_blend::ContinuityLevel::G2)
        .map(|p| p.face)
        .collect();
    assert!(
        !g2_claimed.is_empty(),
        "[{fixture}] the blunt fan must ACCEPT a G2 request (PK-187); it declined \
         ({} apex patches shipped, none claiming G2)",
        result
            .corner_patches
            .iter()
            .filter(|p| p.n_sides == n)
            .count()
    );
    let m = measure_live_cap_seams(&result.brep, &g2_claimed, N_SAMPLES);
    println!(
        "live-seam {fixture}(G2 request): seams={} g0={:.3e} g1={:.3e} g2={:.3e} g2t={:.3e}",
        m.seam_count, m.worst_g0_gap, m.worst_g1_normal, m.worst_g2_residual, m.worst_g2_transverse
    );
    assert!(m.seam_count > 0, "[{fixture}] fan seams must measure");
    helpers::assert_commercial_gates_live(fixture, &m);
    // Watertightness, STRICT on every blunt fixture (PK-208 closed):
    // the crack family's three authors are all fixed — displaced
    // pcurve seeds (cold-projection retry in the boundary-first ring
    // builder), staggered pre-boundary-first sampling (topology-scoped T-junction
    // zippering), and the rimsplit stub deformation (refit_cap_rim
    // now moves windowed-neighbour rim vertices onto the rows and
    // re-coheres their incident curves instead of ramping the exact
    // row segment onto stale vertex positions — the ramp left stub
    // curves 1.48e-3 off the cap and the tessellator rightly refused
    // them). The G1 control is pinned watertight too.
    let opens_g2 = welded_opens(&result.brep, fixture);
    assert_eq!(
        opens_g2, 0,
        "[{fixture}] accepted G2 fan must stay watertight; {opens_g2} open mesh edges"
    );
    let g1_control = fillet_n_gonal_pyramid_result_with_continuity(
        n,
        1.0,
        3.0,
        &radii,
        popkern_blend::ContinuityLevel::G1,
    );
    let opens_g1 = welded_opens(&g1_control.brep, fixture);
    assert_eq!(
        opens_g1, 0,
        "[{fixture}] the G1-control fan must stay watertight; {opens_g1} open mesh edges"
    );
}

fn welded_opens(brep: &popkern_core::topology::BRep, fixture: &str) -> usize {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .unwrap_or_else(|e| panic!("[{fixture}] fan must tessellate: {e:?}"));
    helpers::welded_open_edge_count(&mesh)
}

#[test]
fn pentagon5_blunt_g2_request_is_honest() {
    blunt_g2_request_accepts(5, "pentagon5_blunt");
}

/// Second blunt apex family (PK-187 acceptance 263: the improvement
/// must hold beyond the fixture it was developed against), pinned
/// fully watertight.
#[test]
fn hexagon6_blunt_g2_request_is_honest() {
    blunt_g2_request_accepts(6, "hexagon6_blunt");
}

#[test]
fn polar_passes_commercial_gates_on_pentagon5_uniform() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 0.8, &[0.05; 5], "pentagon5_uniform");
    assert_commercial_gates_live("pentagon5_uniform", &m);
}

#[test]
fn polar_passes_commercial_gates_on_pentagon5_mixed() {
    let m = run_n_gonal_pyramid_polar(
        5,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05],
        "pentagon5_mixed",
    );
    assert_commercial_gates_live("pentagon5_mixed", &m);
}

#[test]
fn pentagon5_mixed_no_longer_records_polar_dispatch_decline() {
    let result = fillet_n_gonal_pyramid_result(5, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    result
        .brep
        .validate()
        .expect("pentagon5_mixed polar route must leave a valid BRep");
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .all(|diagnostic| !matches!(
                diagnostic,
                popkern_blend::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "rational_gregory.dispatch.declined"
            )),
        "pentagon5_mixed should reach N>=5 Gregory/polar dispatch without fallback diagnostics"
    );
}

#[test]
fn polar_passes_commercial_gates_on_hexagon6_uniform() {
    let m = run_n_gonal_pyramid_polar(6, 1.0, 0.8, &[0.05; 6], "hexagon6_uniform");
    assert_commercial_gates_live("hexagon6_uniform", &m);
}

/// Re-pinned 2026-07-15 (PK-087): the support-cell topology this
/// test originally asserted was replaced by the setback vertex fan
/// (4 bi-9 pieces per sector); the commercial claim is now the live
/// seam gate, same as the sibling fixtures. (An interim "skew-trim
/// containment" theory pinned g1 at 1.2e-1 — FALSIFIED: the 8.0e-2
/// twist was an artifact of a projection-walk outlier-drop repair
/// that was itself reverted; the honest reading is g1 = 2.7e-4.)
#[test]
fn hexagon6_mixed_ships_commercial_setback_fan() {
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        "hexagon6_mixed",
    );
    assert_commercial_gates_live("hexagon6_mixed", &m);
}

#[test]
fn polar_passes_commercial_gates_on_pentagon5_variable() {
    // Variable-radius fixture. Wider 1:3 ratio spread
    // (0.02..0.06) than `pentagon5_mixed` (0.04..0.07, ratio ~1.75)
    // stresses the per-sector radius threading. Core-level tests
    // (`variable_radius_*` in `kp15_bi6_polar_bridge_tests`)
    // already verify the construction is homogeneous + finite on
    // asymmetric radii; this fixture verifies end-to-end through
    // the fillet pipeline + commercial G1/G2 gates.
    let m = run_n_gonal_pyramid_polar(
        5,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.04, 0.02],
        "pentagon5_variable",
    );
    assert_commercial_gates_live("pentagon5_variable", &m);
}

#[test]
fn polar_passes_commercial_gates_on_hexagon6_variable() {
    // Variable-radius — hexagonal apex with 1:3 spread.
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02],
        "hexagon6_variable",
    );
    assert_commercial_gates_live("hexagon6_variable", &m);
}

#[test]
fn polar_passes_commercial_gates_on_octagon8_mixed() {
    let m = run_n_gonal_pyramid_polar(
        8,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
        "octagon8_mixed",
    );
    assert_commercial_gates_live("octagon8_mixed", &m);
}

// ── Dispatch smoke test (flag flip must not crash the pipeline) ──
//
// Flipping `use_kp_polar_gregory = true` must not panic `fillet_edges`
// end-to-end. This proves the flag's plumbing is wired, NOT that the
// polar engine actually runs on the fixture — with the upstream
// N ≥ 5 gap (see module header) the polar promote function is never
// invoked, so this test is compatible with that gap. When the gap
// closes, this test will still pass and the commercial-gate fixtures
// above will start being flip-able.

#[test]
fn polar_flag_flip_does_not_crash_pipeline_at_n5() {
    let pyr = BRep::n_gonal_pyramid(5, 1.0, 0.8, Tolerance::default()).unwrap();
    let Some(vertex) = pick_vertex_with_valence(&pyr, 5) else {
        return;
    };
    let edges = incident_edges(&pyr, vertex);
    let (_brep, _faces) =
        fillet_corner_full(&pyr, &edges, &[0.04, 0.06, 0.05, 0.07, 0.05], true, true);
}

/// PK-187 acceptance 265 — Class-A inspection evidence across the
/// accepted blunt G2 fans' seams, in the measurable stripe-phase form
/// of the zebra render: a G2-continuous seam transitions stripes
/// smoothly (mismatch ratio near 0), a crease reads a large jump.
/// Same instrument and 0.25 bar as the PK-197 low-valence zebra gate.
#[test]
fn blunt_g2_fans_are_zebra_continuous_across_seams() {
    use popkern_core::geometry::Vec3;
    for (n, fixture) in [(5usize, "pentagon5_blunt"), (6, "hexagon6_blunt")] {
        let radii = vec![0.05; n];
        let result = fillet_n_gonal_pyramid_result_with_continuity(
            n,
            1.0,
            3.0,
            &radii,
            popkern_blend::ContinuityLevel::G2,
        );
        let caps: std::collections::BTreeSet<_> = result
            .corner_patches
            .iter()
            .filter(|p| p.n_sides == n)
            .map(|p| p.face)
            .collect();
        let brep = &result.brep;
        let mut measured = 0usize;
        let mut worst = 0.0_f64;
        for (_eh, edge) in brep.edges.iter() {
            let Some(he) = brep.half_edges.get(edge.half_edge) else {
                continue;
            };
            let Some(twin) = brep.half_edges.get(he.twin) else {
                continue;
            };
            let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_))
            else {
                continue;
            };
            let (fa, fb) = (la.face, lb.face);
            if fa == fb || !caps.contains(&fa) || !caps.contains(&fb) {
                continue;
            }
            if let Some(kink) = popkern_blend::zebra::zebra_boundary_kink(
                brep,
                fa,
                fb,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                0.02,
                24,
                1e-6,
            ) {
                measured += 1;
                worst = worst.max(kink);
            }
        }
        println!("[{fixture}] zebra: seams_measured={measured} worst_stripe_mismatch={worst:.4}");
        assert!(measured > 0, "[{fixture}] zebra measured no fan seams");
        assert!(
            worst < 0.25,
            "[{fixture}] accepted G2 fan shows a zebra stripe mismatch of {worst:.4} \
             across a seam — not Class-A continuous"
        );
    }
}
