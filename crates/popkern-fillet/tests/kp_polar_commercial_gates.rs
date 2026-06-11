//! Phase I — commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (Phase H, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "kp_polar_support_cell_diagnostics.rs"]
mod support_cell_diagnostics;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, collect_patches, fillet_corner_full,
    incident_edges, measure_ring, pick_vertex_with_valence,
};
use support_cell_diagnostics::assert_support_cell_quality_gate;

fn run_n_gonal_pyramid_polar(
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
    let (brep, faces) = fillet_corner_full(&pyr, &edges, radii, true, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

fn fillet_n_gonal_pyramid_result(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
) -> popkern_fillet::BlendResult {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
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
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None).unwrap()
}

/// Blunt-apex positive control (apex half-angle 18°, typical CAD).
/// Expected to take sphere-fit path and bypass polar entirely.
/// The `assert_all_binonic` assertion in `run_n_gonal_pyramid_polar`
/// will fail at the moment because the sphere-fit path produces a
/// single non-Gregory face, not N Gregory sub-patches. This is the
/// intended commercial behaviour for gentle corners — the test is
/// `#[ignore]` until the probe is refactored to assert "valid
/// commercial output" (sphere-fit OR polar-bi9) rather than
/// specifically "N Gregory bi-nonic patches."
///
/// Documented here as a reminder: pentagon5 at h=3.0 (18° apex) is
/// the REAL commercial-CAD test case; the `_uniform` / `_mixed` /
/// `_variable` sharp-apex tests below stress the polar path
/// specifically.
#[ignore]
#[test]
fn pentagon5_blunt_takes_sphere_fit_path() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 3.0, &[0.05; 5], "pentagon5_blunt");
    assert_commercial_gates("pentagon5_blunt", &m);
}
#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_uniform() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 0.8, &[0.05; 5], "pentagon5_uniform");
    assert_commercial_gates("pentagon5_uniform", &m);
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
    assert_commercial_gates("pentagon5_mixed", &m);
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
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "rational_gregory.dispatch.declined"
            )),
        "pentagon5_mixed should reach N>=5 Gregory/polar dispatch without fallback diagnostics"
    );
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_uniform() {
    let m = run_n_gonal_pyramid_polar(6, 1.0, 0.8, &[0.05; 6], "hexagon6_uniform");
    assert_commercial_gates("hexagon6_uniform", &m);
}

#[ignore]
#[test]
fn polar_support_cell_quality_gate_on_hexagon6_mixed() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    assert_support_cell_quality_gate("hexagon6_mixed", &result, 6);
}

#[test]
fn hexagon6_mixed_accepts_support_cell_topology() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    assert_support_cell_quality_gate("hexagon6_mixed", &result, 6);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_variable() {
    // Phase I.22 variable-radius fixture. Wider 1:3 ratio spread
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
    assert_commercial_gates("pentagon5_variable", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_variable() {
    // Phase I.22 variable-radius — hexagonal apex with 1:3 spread.
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02],
        "hexagon6_variable",
    );
    assert_commercial_gates("hexagon6_variable", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_octagon8_mixed() {
    let m = run_n_gonal_pyramid_polar(
        8,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
        "octagon8_mixed",
    );
    assert_commercial_gates("octagon8_mixed", &m);
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

/// Diagnostic probe — PK-085 G1 campaign measurement, run
/// manually with `--ignored`. Prints the constructive-cells
/// facts and the final shared-trim quality. Removal rule: delete
/// when the G1 campaign concludes (numbers recorded in the
/// PK-085 ledger at each step).
#[ignore]
#[test]
fn hexagon6_mixed_g1_campaign_probe() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    for diagnostic in &result.repair_report.apex_diagnostics {
        let popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
            continue;
        };
        let Some(detail) = cap.detail_code.as_deref() else {
            continue;
        };
        if detail.contains("chain_constructive")
            || (cap.stage == "kp20_polar_macro.dispatch.support_cell_shared_trim_quality"
                && detail.contains("measurement_scope=final_post_pipeline"))
        {
            println!("fact: {detail}");
        }
    }
}

/// Diagnostic probe — PK-085 G1 campaign worst-sample anatomy,
/// run manually with `--ignored`. Removal rule: delete when the
/// campaign concludes.
#[ignore]
#[test]
fn hexagon6_mixed_g1_worst_sample_anatomy_probe() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    let brep = &result.brep;
    let faces: std::collections::HashSet<_> = result
        .corner_patches
        .iter()
        .map(|patch| patch.face)
        .collect();
    let mut worst: Option<(f64, String)> = None;
    let mut visited = std::collections::HashSet::new();
    for &face in &faces {
        let Some(start) = brep
            .faces
            .get(face)
            .and_then(|f| brep.loops.get(f.outer_loop))
            .map(|l| l.half_edge)
        else {
            continue;
        };
        let mut cursor = start;
        for _ in 0..brep.half_edges.len() {
            let he = &brep.half_edges[cursor];
            if visited.insert(he.edge) {
                let twin = &brep.half_edges[he.twin];
                let twin_face = brep.loops[twin.loop_].face;
                if faces.contains(&twin_face)
                    && let (Some(pa), Some(pb)) = (he.pcurve.as_ref(), twin.pcurve.as_ref())
                    && let (Some(sa), Some(sb)) = (
                        brep.faces[face].surface.as_ref(),
                        brep.faces[twin_face].surface.as_ref(),
                    )
                {
                    for index in 0..=32 {
                        let f = index as f64 / 32.0;
                        let ta = pa.t_start + (pa.t_end - pa.t_start) * f;
                        let tb = pb.t_start + (pb.t_end - pb.t_start) * (1.0 - f);
                        let uva = pa.point_at(ta);
                        let uvb = pb.point_at(tb);
                        let (SurfaceGeometry::Gregory(ga), SurfaceGeometry::Gregory(gb)) =
                            (&sa.geometry, &sb.geometry)
                        else {
                            continue;
                        };
                        let sua = ga.partial_u(uva.x, uva.y);
                        let sva = ga.partial_v(uva.x, uva.y);
                        let sub_ = gb.partial_u(uvb.x, uvb.y);
                        let svb = gb.partial_v(uvb.x, uvb.y);
                        let (Some(na), Some(nb)) =
                            (sua.cross(&sva).normalized(), sub_.cross(&svb).normalized())
                        else {
                            continue;
                        };
                        let disagreement = 1.0 - na.dot(&nb).abs();
                        if worst.as_ref().is_none_or(|(w, _)| disagreement > *w) {
                            let detail = format!(
                                "g1={disagreement:.6e} f={f:.4} face_a={face:?} uv_a=({:.6},{:.6}) face_b={twin_face:?} uv_b=({:.6},{:.6})\\n  su_a={sua:?}\\n  sv_a={sva:?}\\n  su_b={sub_:?}\\n  sv_b={svb:?}",
                                uva.x, uva.y, uvb.x, uvb.y,
                            );
                            worst = Some((disagreement, detail));
                        }
                    }
                }
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
    if let Some((_, detail)) = worst {
        println!("worst_g1_sample: {detail}");
    }
}
