//! Phase 6c-X.A.6 — validation matrix at full `validate_all`.
//!
//! Path A (apex topology rebuild) lands as A.3 + A.4 + A.5 behind
//! the `POPKERN_APEX_PATH_A=1` env flag.  A.6's gate is:
//! **with Path A enabled, every fixture in the validation matrix
//! must pass `BRep::validate_all()` cleanly** (zero validation
//! errors).
//!
//! `validate_all` is the expensive O(n²) self-intersection check
//! plus the structural and geometric checks that `validate()`
//! covers.  Path A claims to close the apex-region
//! vertex-to-curve-endpoint coherence gap (`validate_watertight` /
//! `validate_loop_closure`) for variable-radius cube and pyramid4
//! fixtures.  The diagnostic test below is the lock-in: it dumps
//! per-fixture error counts so progress can be tracked
//! incrementally without flipping the gate.
//!
//! When every fixture passes cleanly, A.7 flips
//! `POPKERN_APEX_PATH_A=1` to the default and removes the legacy
//! workarounds the new path supersedes.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

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

/// Wrap a fixture build under the Path A flag set.  Returns the
/// post-fillet BRep and the cap face handles.
fn build_with_path_a(
    is_pyramid: bool,
    radii: &[f64],
) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
    let _g4 = EnvGuard::set("POPKERN_APEX_PATH_A", "1");

    if is_pyramid {
        let pyr = BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, Tolerance::default()).unwrap();
        let vertex = pick_vertex_with_valence(&pyr, radii.len()).unwrap();
        let edges = incident_edges(&pyr, vertex);
        fillet_corner(&pyr, &edges, radii, true)
    } else {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = cube.vertices.keys().next().unwrap();
        let edges = incident_edges(&cube, vertex);
        fillet_corner(&cube, &edges, radii, true)
    }
}

/// A.6 baseline: dump `validate_all` error counts per fixture
/// under Path A.  Ignored by default — run with
/// `cargo test phase_6c_a6_validate_all_baseline -- --ignored
/// --nocapture` to see today's state.
#[ignore]
#[test]
fn phase_6c_a6_validate_all_baseline() {
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

    eprintln!("\n=== A.6 validate_all baseline (POPKERN_APEX_PATH_A=1) ===");
    eprintln!("tolerance.linear() = {:.3e}", Tolerance::default().linear());
    let mut totals: Vec<(String, usize)> = Vec::with_capacity(cases.len());
    for case in cases {
        let (brep, _faces) = build_with_path_a(case.is_pyramid, case.radii);
        let errors = brep.validate_all();
        eprintln!("  {:18}: {} validation errors", case.name, errors.len());
        for (i, err) in errors.iter().enumerate().take(5) {
            eprintln!("    [{i}] {err:?}");
        }
        if errors.len() > 5 {
            eprintln!("    ... ({} more)", errors.len() - 5);
        }
        totals.push((case.name.to_string(), errors.len()));
    }
    eprintln!("\nSummary:");
    for (name, n) in &totals {
        let mark = if *n == 0 { "✓" } else { "✗" };
        eprintln!("  {mark} {name:18} {n} errors");
    }
}

/// A.6 ratchet: per-fixture empirical-state regression guard.
/// As Path A closes more fixtures, the bounds tighten.  When
/// every count reaches zero, A.7 flips the default.
///
/// Today's snapshot (2026-05-01):
/// * pyramid4_mixed: 1 error (InconsistentOrientation on a face).
/// * pyramid4_uniform: 3 errors.
/// * cube_mixed / cube_extreme: 2 errors each.
/// * cube_uniform / pyramid3_uniform: 3 errors each.
///
/// This test asserts the current state of pyramid4_mixed (the
/// closest-to-clean fixture).  If a regression introduces
/// MORE errors, this trips; if work brings it to zero, the
/// assertion can tighten via the `≤ N` form below.
#[test]
fn phase_6c_a6_pyramid4_mixed_ratchet() {
    let (brep, _faces) = build_with_path_a(true, &[0.04, 0.06, 0.05, 0.07]);
    let errors = brep.validate_all();
    assert!(
        errors.len() <= 1,
        "pyramid4_mixed under Path A regression: {} error(s) (was 1):\n{}",
        errors.len(),
        errors
            .iter()
            .take(8)
            .map(|e| format!("  - {e:?}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );
}

// Empirical finding (2026-05-01): with or without
// `POPKERN_APEX_PATH_A=1`, `pyramid4_mixed`'s `validate_all`
// reports the same error count.  Path A's vertex-relocation +
// base-edge-repair work IS happening (POPKERN_APEX_PATH_A_TRACE
// shows 8 vertices moved, worst |V→T| = 5e-2 m), but the
// remaining `validate_all` error (`InconsistentOrientation` on
// a single face) is unrelated to the V→T coherence problem
// Path A targets.  The error reflects a degenerate surface
// normal at a face's parameter midpoint — distinct from the
// `validate_watertight` / `validate_loop_closure` checks Path
// A is designed to close.  A.7's flag-flip gate will need a
// richer regression metric than `validate_all().len()` —
// likely `cap_outer_g0` + a sweep of structural checks.
