//! Phase 6c-X.A.7 perf sanity check.
//!
//! Compares Path A enabled vs disabled across the validation
//! matrix.  Asserts Path A's overhead is bounded — sanity gate
//! before flipping the default.  Not a microbenchmark; uses
//! `std::time::Instant` for ms-level coarseness.
//!
//! Activation guidance: ignored by default; run with
//! `cargo test phase_6c_a7_perf -- --ignored --nocapture`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::Instant;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{fillet_corner, incident_edges, pick_vertex_with_valence};

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

fn run_fixture(is_pyramid: bool, radii: &[f64], path_a: bool) -> std::time::Duration {
    // STAGE4_WIDE_WIRE_TOL and APEX_FACE_LIFT are default-on since
    // #207 cleanup.  APEX_PATH_A defaults on too but is the A/B knob
    // in this perf test.  SPHERE_FIT_BOUNDARY remains opt-in.
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let _g4 = EnvGuard::set("POPKERN_APEX_PATH_A", if path_a { "1" } else { "0" });

    let start = Instant::now();
    if is_pyramid {
        let pyr = BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, Tolerance::default()).unwrap();
        let vertex = pick_vertex_with_valence(&pyr, radii.len()).unwrap();
        let edges = incident_edges(&pyr, vertex);
        let _ = fillet_corner(&pyr, &edges, radii, true);
    } else {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = cube.vertices.keys().next().unwrap();
        let edges = incident_edges(&cube, vertex);
        let _ = fillet_corner(&cube, &edges, radii, true);
    }
    start.elapsed()
}

#[ignore]
#[test]
fn phase_6c_a7_path_a_overhead_bounded() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
        },
        Case {
            name: "cube_mixed",
            radii: &[0.08, 0.10, 0.12],
            is_pyramid: false,
        },
        Case {
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
        },
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    eprintln!("\n=== A.7 Path A perf overhead ===");
    eprintln!(
        "{:<18}  {:>10}  {:>10}  {:>10}",
        "fixture", "off (ms)", "on (ms)", "delta"
    );
    eprintln!("{:-<54}", "");
    let mut all_within_budget = true;
    for case in cases {
        // Warmup (allocator + JIT-ish caching).
        let _ = run_fixture(case.is_pyramid, case.radii, false);
        let _ = run_fixture(case.is_pyramid, case.radii, true);

        // Measured runs.
        let off = run_fixture(case.is_pyramid, case.radii, false);
        let on = run_fixture(case.is_pyramid, case.radii, true);
        let off_ms = off.as_secs_f64() * 1000.0;
        let on_ms = on.as_secs_f64() * 1000.0;
        let delta_pct = if off_ms > 0.0 {
            (on_ms - off_ms) / off_ms * 100.0
        } else {
            0.0
        };
        eprintln!(
            "{:<18}  {:>10.2}  {:>10.2}  {:>+9.1}%",
            case.name, off_ms, on_ms, delta_pct
        );
        // Budget: Path A overhead should not exceed 50% on any
        // single fixture.  Tighter than the slowest measurement
        // observed (cube_extreme at ~30% overhead from extra
        // diagnostic + relocation).  If overhead exceeds 50% on
        // a fixture, perf has regressed.
        if delta_pct > 50.0 {
            all_within_budget = false;
        }
    }
    assert!(
        all_within_budget,
        "Path A overhead exceeded 50% on at least one fixture",
    );
}
