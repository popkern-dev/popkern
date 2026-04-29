//! Phase 6c-D Stage 1 — boundary-injection probe.
//!
//! Validates that, with `POPKERN_SPHERE_FIT_BOUNDARY=1` set AND
//! `use_kp_binonic_gregory = true`, the K-P bi-nonic ring builder
//! consumes the sphere-fit-substituted hex-side boundary CPs
//! (`BlendBoundary::hex_boundary_cps_deg9`) instead of falling
//! back to `k_p_boundary` Hermite synthesis.
//!
//! # Stage 1 status: plumbing only
//!
//! Stage 1 *injects* the new boundary but does **not** re-fit the
//! patch's interior CPs to match. Empirically, the resulting
//! bi-nonic patch self-intersects (the row0/col9 boundary moved
//! but the strict-interior 64 CPs — still Hermite-blend output —
//! produced a fold). The promotion ladder rejects the fold via
//! `self_intersects_with_tolerance` and falls back to a lower
//! rung. Net effect: with injection on, fewer fixtures land
//! binonic than without.
//!
//! This is *expected*. The empirical fall-off proves the injection
//! is reaching the K-P guide builder — without injection the
//! patches don't fold. Stage 2 (re-fit strict-interior CPs against
//! the injected boundary) closes this; once it lands, binonic
//! lands again *and* the cap-outer boundary tracks the analytic
//! arc to the polynomial-fit floor.
//!
//! # Matrix probe
//!
//! `boundary_injection_matrix_probe` (ignored by default) prints
//! a side-by-side comparison of inject=false vs inject=true on
//! cube_mixed and pyramid4_mixed. Run with
//!
//! ```text
//! cargo test -p popkern-fillet --test kp_binonic_boundary_injection \
//!   -- --ignored --nocapture
//! ```
//!
//! to read the binonic-landing-rate and cap-outer continuity numbers.
//!
//! # Env-var safety
//!
//! Each test that sets `POPKERN_SPHERE_FIT_BOUNDARY` uses an
//! `EnvGuard` so the prior value is restored even if the test
//! panics — same-process other tests aren't perturbed.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, collect_patches, fillet_corner, incident_edges, measure_cap_outer_boundary,
    pick_vertex_with_valence,
};

/// RAII guard restoring the prior value of an env-var on drop.
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

struct ProbeResult {
    n_corner_faces: usize,
    n_gregory_patches: usize,
    n_binonic: usize,
    cap_outer_g0: f64,
    cap_outer_g1: f64,
    edges_measured: usize,
}

impl ProbeResult {
    fn print(&self, fixture: &str) {
        eprintln!(
            "  {fixture:18}: faces={} gregory={} binonic={} \
             cap_outer_g0={:.3e} cap_outer_g1={:.3e} edges={}",
            self.n_corner_faces,
            self.n_gregory_patches,
            self.n_binonic,
            self.cap_outer_g0,
            self.cap_outer_g1,
            self.edges_measured,
        );
    }
}

fn probe_cube(radii: &[f64], inject: bool) -> ProbeResult {
    let _guard = inject.then(|| EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1"));
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    finalize(&brep, &faces)
}

fn probe_square_pyramid(radii: &[f64], inject: bool) -> ProbeResult {
    let _guard = inject.then(|| EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1"));
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    finalize(&brep, &faces)
}

fn finalize(brep: &BRep, faces: &[popkern_core::topology::FaceHandle]) -> ProbeResult {
    let patches = collect_patches(brep, faces);
    let metrics = measure_cap_outer_boundary(brep, faces, N_SAMPLES);
    let n_binonic = patches
        .iter()
        .filter(|p| matches!(p, GregoryPatch::RectangularBiNonic(_)))
        .count();
    ProbeResult {
        n_corner_faces: faces.len(),
        n_gregory_patches: patches.len(),
        n_binonic,
        cap_outer_g0: metrics.worst_g0_gap,
        cap_outer_g1: metrics.worst_g1_coplanarity,
        edges_measured: metrics.edges_measured,
    }
}

/// Side-by-side comparison: cap-outer continuity and binonic
/// landing rate with vs without boundary injection. Diagnostic
/// only — Stage 1 expectation is that injection regresses the
/// binonic landing rate (folded interior, see module header).
/// Stage 2's interior re-fit will turn this into a real gate.
#[ignore]
#[test]
fn boundary_injection_matrix_probe() {
    eprintln!("\n=== Phase 6c-D Stage 1 boundary-injection matrix ===");
    eprintln!("(use_kp_binonic_gregory=true; ±POPKERN_SPHERE_FIT_BOUNDARY)");
    eprintln!("--- inject=false (baseline: Hermite-synth boundary) ---");
    let cm_off = probe_cube(&[0.08, 0.10, 0.12], false);
    cm_off.print("cube_mixed");
    let p4m_off = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], false);
    p4m_off.print("pyramid4_mixed");
    eprintln!("--- inject=true (Stage 1: substituted boundary) ---");
    let cm_on = probe_cube(&[0.08, 0.10, 0.12], true);
    cm_on.print("cube_mixed");
    let p4m_on = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], true);
    p4m_on.print("pyramid4_mixed");
}

/// Stage 1 wiring sanity: verify the *baseline* (no injection)
/// path still ships binonic on cube_mixed. Guards against the
/// Stage 1 plumbing accidentally regressing the inject-OFF path —
/// the `Option<&[...]>` parameter must thread `None` through
/// cleanly when the boundary's `hex_boundary_cps_deg9` is empty.
#[test]
fn cube_mixed_baseline_binonic_unaffected_by_stage1_plumbing() {
    let res = probe_cube(&[0.08, 0.10, 0.12], false);
    res.print("cube_mixed[no-inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 1 plumbing broke the inject=OFF baseline — cube_mixed should still ship binonic"
    );
}

/// Same baseline guard for pyramid4_mixed (the fixture that
/// already passes commercial gates per the H.7 matrix). With no
/// injection the K-P binonic dispatch must continue to land.
#[test]
fn pyramid4_mixed_baseline_binonic_unaffected_by_stage1_plumbing() {
    let res = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], false);
    res.print("pyramid4_mixed[no-inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 1 plumbing broke the inject=OFF baseline — pyramid4_mixed should still ship binonic"
    );
}

/// Phase 6c-D Stage 3 lock-in: under injection, both cube_mixed
/// and pyramid4_mixed must ship full binonic rings. The
/// row-9-based mid-edge formula (`row1[j] = (8/9)·row0[j] +
/// (1/9)·row9[j]`) + Coons interior fill produces fold-free
/// patches by construction, and the quality regression guard is
/// skipped under injection (the spoke-G1 metric is determined by
/// boundary geometry, not the L-M optimum, so the 1.2× threshold
/// doesn't apply).
///
/// Note: cap-outer continuity is **not** asserted yet. Stage 3
/// proves the dispatch ladder accepts injection-mode binonic; the
/// 3-D cap-vs-cylinder fidelity needs a follow-up (likely pcurve-
/// refit under the new boundary or wider edge-match tolerance in
/// `wire_substituted_boundary_to_edges`).
#[test]
fn cube_mixed_binonic_ships_under_injection() {
    let res = probe_cube(&[0.08, 0.10, 0.12], true);
    res.print("cube_mixed[inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 3 regression: cube_mixed under injection should ship at least one binonic patch"
    );
}

#[test]
fn pyramid4_mixed_binonic_ships_under_injection() {
    let res = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], true);
    res.print("pyramid4_mixed[inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 3 regression: pyramid4_mixed under injection should ship at least one binonic patch"
    );
}
