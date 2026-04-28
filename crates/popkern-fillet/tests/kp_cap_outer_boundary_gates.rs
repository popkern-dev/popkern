//! Cap-outer-boundary commercial-grade continuity gate.
//!
//! Sibling of `kp_binonic_commercial_gates.rs` — but where that
//! file measures continuity *between adjacent cap sub-patches*
//! at their shared inner spokes, this file measures continuity
//! between each cap sub-patch and the cylinder fillet roll it
//! sits on top of.
//!
//! # Why this exists
//!
//! Per Task #144 / 2026-04-27 failed-approaches entry: the
//! existing inter-sub-patch ring gates passed while the cap-outer
//! boundary silently drifted by 5–15% bbox in 3-D. FreeCAD shows
//! this as the "messy top" defect on every fillet zoo pyramid.
//! The user-visible quality bar is set here, not at the inner
//! spokes.
//!
//! # Status (2026-04-27)
//!
//! All fixtures `#[ignore]` because they're expected to fail
//! today — the gate is being introduced *before* the fix to
//! quantify the problem precisely. Once Task #144 closes, these
//! flip to passing and become regression guards.
//!
//! Threshold targets (mirroring `kp_binonic_commercial_gates`):
//! * G0 ≤ 1e-10 (position match at the seam)
//! * G1 ≤ 5e-2  (coplanarity at samples; ≈ 2.9° normal mismatch)

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    G0_GATE, G1_GATE, N_SAMPLES, fillet_corner, incident_edges, measure_cap_outer_boundary,
    pick_vertex_with_valence,
};

fn assert_cap_outer_gates(fixture: &str, m: &helpers::CapOuterMetrics) {
    println!(
        "cap_outer {fixture}: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
    assert!(
        m.edges_measured >= 1,
        "[{fixture}] no cap-vs-non-cap edges found — fixture broke"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "[{fixture}] cap-outer G0 gap {:.3e} exceeds {:.1e} gate",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_coplanarity < G1_GATE,
        "[{fixture}] cap-outer G1 coplanarity {:.3e} exceeds {:.1e} gate",
        m.worst_g1_coplanarity,
        G1_GATE
    );
}

fn run_cube(radii: &[f64]) -> helpers::CapOuterMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

fn run_n_pyramid(n: usize, h: f64, radii: &[f64]) -> helpers::CapOuterMetrics {
    let pyr = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

/// Phase 6c-C.1 win: cube_uniform G1 commercial-grade, locked
/// in as a regression guard. The 5e-2 gate is the
/// `kp_gate_helpers::G1_GATE` constant. Passed since commit
/// `016a780` (centroid-based interior-direction tie-break).
///
/// G0 still exceeds `G0_GATE` (1e-10) because the cap face's
/// 3 rolling-ball-footprint vertices aren't on the sphere; full
/// Phase 6c-C topology rebuild is needed to remove them. This
/// test gates only on G1 — graduates to the full
/// `assert_cap_outer_gates` once 6c-C.2/.3 land.
#[test]
fn cap_outer_boundary_g1_passes_on_cube_uniform() {
    let m = run_cube(&[0.10, 0.10, 0.10]);
    println!(
        "cap_outer cube_uniform: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
    assert!(m.edges_measured >= 1);
    assert!(
        m.worst_g1_coplanarity < helpers::G1_GATE,
        "cube_uniform G1 {:.3e} exceeds {:.1e} commercial gate",
        m.worst_g1_coplanarity,
        helpers::G1_GATE,
    );
}

#[ignore]
#[test]
fn cap_outer_boundary_passes_on_cube_mixed() {
    let m = run_cube(&[0.08, 0.10, 0.12]);
    assert_cap_outer_gates("cube_mixed", &m);
}

#[ignore]
#[test]
fn cap_outer_boundary_passes_on_cube_extreme() {
    let m = run_cube(&[0.02, 0.10, 0.20]);
    assert_cap_outer_gates("cube_extreme", &m);
}

#[ignore]
#[test]
fn cap_outer_boundary_passes_on_pyramid4_mixed() {
    let m = run_n_pyramid(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    assert_cap_outer_gates("pyramid4_mixed", &m);
}

/// N=5 — PnS path, not Gregory caps. Helper only counts cap edges
/// where one side is a Gregory face; PnS produces Sphere-surface
/// caps, so this fixture should report 0 edges measured. Kept here
/// so the matrix is documented; will be skipped or replaced when
/// the PnS sibling helper lands.
#[ignore]
#[test]
fn cap_outer_boundary_pyramid5_mixed_path_check() {
    let m = run_n_pyramid(5, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    println!(
        "cap_outer pyramid5_mixed: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

/// Probe-style: print numbers across the canonical fixtures with
/// no assertions so we can see the baseline without picking
/// thresholds first. Useful while developing the fix. Includes
/// `pyramid_n04_zoo` — the same geometry rendered by the
/// `fillet_zoo` showcase binary, where FreeCAD shows the
/// "messy top" defect.
///
/// G1 is reported but is *not* meaningful while G0 is large —
/// we're sampling each side's surface at 3-D points that are
/// 5-20cm apart, so the normals are essentially uncorrelated.
/// Once G0 drops below ~1e-3, G1 becomes informative.
/// Single-fixture probe for the cube_uniform path — used to trace
/// the per-vertex insert_vertex_blend invocation count without the
/// noise of the full 8-fixture matrix.
#[ignore]
#[test]
fn print_cube_uniform_only() {
    eprintln!("\n=== cube_uniform single-fixture probe ===");
    let m = run_cube(&[0.10, 0.10, 0.10]);
    eprintln!(
        "cube_uniform: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

#[ignore]
#[test]
fn print_pyramid3_uniform_only() {
    eprintln!("\n=== pyramid3_uniform single-fixture probe ===");
    let m = run_n_pyramid(3, 1.2, &[0.10, 0.10, 0.10]);
    eprintln!(
        "pyramid3_uniform: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

#[ignore]
#[test]
fn print_cap_outer_boundary_baseline() {
    eprintln!("\n=== cap-outer-boundary baseline (no assertions) ===");
    for (name, m) in [
        ("cube_uniform", run_cube(&[0.10, 0.10, 0.10])),
        ("cube_mixed", run_cube(&[0.08, 0.10, 0.12])),
        ("cube_extreme", run_cube(&[0.02, 0.10, 0.20])),
        (
            "pyramid3_uniform",
            run_n_pyramid(3, 1.2, &[0.10, 0.10, 0.10]),
        ),
        (
            "pyramid4_uniform",
            run_n_pyramid(4, 1.2, &[0.10, 0.10, 0.10, 0.10]),
        ),
        (
            "pyramid4_mixed",
            run_n_pyramid(4, 0.8, &[0.04, 0.06, 0.05, 0.07]),
        ),
        (
            "pyramid_n04_zoo",
            run_n_pyramid(4, 1.2, &[0.07, 0.09, 0.11, 0.13]),
        ),
        (
            "pyramid_n03_zoo",
            run_n_pyramid(3, 1.2, &[0.07, 0.10, 0.13]),
        ),
    ] {
        eprintln!(
            "  {name:18}: g0={:.3e}  g1={:.3e}  edges={}",
            m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
        );
    }
}
