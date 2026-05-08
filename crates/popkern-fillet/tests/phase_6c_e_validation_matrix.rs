//! Phase 6c-E validation matrix.
//!
//! Runs the apex-fillet fixture matrix under Stage 4.7e
//! (`POPKERN_SPHERE_FIT_BOUNDARY=1` + `POPKERN_STAGE4_WIDE_WIRE_TOL=1`)
//! and captures `cap_outer_g0` / `cap_outer_g1` for each fixture.
//!
//! Goal: confirm the Stage 4.7e win on `pyramid4_mixed`
//! (4.85e-10 m worst case) generalises to the rest of the matrix
//! before flipping the flag default in Phase 6c-F.
//!
//! Commercial gate: `cap_outer_g0` ≤ 5e-3 m.
//!
//! Each fixture runs as a separate test so failures localise.
//! `#[ignore]` keeps legacy/stage-specific diagnostics out of the
//! default gate when they are superseded by native A2 acceptance or
//! still need construction work.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    CapOuterMetrics, N_SAMPLES, fillet_corner, incident_edges,
    measure_cap_outer_boundary_including_pcurve_only, pick_vertex_with_valence,
};

/// Phase 6c-E commercial-grade gate. Stage 4.7e on
/// `pyramid4_mixed` reports 4.85e-10 m; 5e-3 m is the user-
/// visible quality bar.
const PHASE_6C_E_G0_GATE: f64 = 5e-3;

/// RAII env-var guard.
struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        unsafe {
            std::env::set_var(key, value);
        }
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.prior {
                Some(v) => std::env::set_var(self.key, v),
                None => std::env::remove_var(self.key),
            }
        }
    }
}

fn run_cube_with_stage4(radii: &[f64]) -> CapOuterMetrics {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    measure_cap_outer_boundary_including_pcurve_only(&brep, &faces, N_SAMPLES)
}

fn run_n_pyramid_with_stage4(n: usize, h: f64, radii: &[f64]) -> CapOuterMetrics {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    measure_cap_outer_boundary_including_pcurve_only(&brep, &faces, N_SAMPLES)
}

fn report(fixture: &str, m: &CapOuterMetrics) {
    eprintln!(
        "  {fixture:24}: cap_outer_g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

fn assert_g0(fixture: &str, m: &CapOuterMetrics) {
    assert!(
        m.edges_measured >= 1,
        "[{fixture}] no cap-vs-non-cap edges measured"
    );
    assert!(
        m.worst_g0_gap < PHASE_6C_E_G0_GATE,
        "[{fixture}] cap_outer_g0 {:.3e} exceeds {:.1e} commercial gate",
        m.worst_g0_gap,
        PHASE_6C_E_G0_GATE
    );
}

#[ignore]
#[test]
fn phase_6c_e_matrix_probe() {
    eprintln!(
        "\n=== Phase 6c-E validation matrix \
         (POPKERN_SPHERE_FIT_BOUNDARY=1, POPKERN_STAGE4_WIDE_WIRE_TOL=1) ==="
    );
    let cu = run_cube_with_stage4(&[0.10, 0.10, 0.10]);
    report("cube_uniform", &cu);
    let cm = run_cube_with_stage4(&[0.08, 0.10, 0.12]);
    report("cube_mixed", &cm);
    let ce = run_cube_with_stage4(&[0.02, 0.10, 0.20]);
    report("cube_extreme", &ce);
    let p3 = run_n_pyramid_with_stage4(3, 1.2, &[0.10, 0.10, 0.10]);
    report("pyramid3_uniform", &p3);
    let p4u = run_n_pyramid_with_stage4(4, 0.8, &[0.05, 0.05, 0.05, 0.05]);
    report("pyramid4_uniform", &p4u);
    let p4m = run_n_pyramid_with_stage4(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    report("pyramid4_mixed", &p4m);
}

#[test]
fn phase_6c_e_pyramid4_mixed_under_5mm() {
    let m = run_n_pyramid_with_stage4(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    report("pyramid4_mixed", &m);
    assert_g0("pyramid4_mixed", &m);
}

#[ignore]
#[test]
fn phase_6c_e_cube_uniform_under_5mm() {
    let m = run_cube_with_stage4(&[0.10, 0.10, 0.10]);
    report("cube_uniform", &m);
    assert_g0("cube_uniform", &m);
}

#[test]
fn phase_6c_e_cube_mixed_under_5mm() {
    let m = run_cube_with_stage4(&[0.08, 0.10, 0.12]);
    report("cube_mixed", &m);
    assert_g0("cube_mixed", &m);
}

/// Legacy Stage 4 diagnostic. The metric now measures pcurve-only
/// cap/support boundaries instead of silently skipping them; this
/// fixture exposes a real legacy cap_outer gap (~1.8e-1 m). The
/// commercial default gate for this fixture is native A2 in
/// `a2_5_cap_outer_g0_gregory_under_a2`.
#[ignore]
#[test]
fn phase_6c_e_cube_extreme_under_5mm() {
    let m = run_cube_with_stage4(&[0.02, 0.10, 0.20]);
    report("cube_extreme", &m);
    assert_g0("cube_extreme", &m);
}

#[ignore]
#[test]
fn phase_6c_e_pyramid3_uniform_under_5mm() {
    let m = run_n_pyramid_with_stage4(3, 1.2, &[0.10, 0.10, 0.10]);
    report("pyramid3_uniform", &m);
    assert_g0("pyramid3_uniform", &m);
}

#[ignore]
#[test]
fn phase_6c_e_pyramid4_uniform_under_5mm() {
    let m = run_n_pyramid_with_stage4(4, 0.8, &[0.05, 0.05, 0.05, 0.05]);
    report("pyramid4_uniform", &m);
    assert_g0("pyramid4_uniform", &m);
}

#[path = "support/phase_6c_e_breakdowns.rs"]
mod phase_6c_e_breakdowns;

#[path = "support/phase_6c_e_apex_diagnostics.rs"]
mod phase_6c_e_apex_diagnostics;

#[path = "support/phase_6c_e_vertex_coherence.rs"]
mod phase_6c_e_vertex_coherence;

#[path = "support/phase_6c_e_structural_mode_b.rs"]
mod phase_6c_e_structural_mode_b;
