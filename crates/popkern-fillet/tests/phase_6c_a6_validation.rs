//! Phase 6c-X.A.6 — composite validation gate for apex relocation.
//!
//! The apex topology rebuild lands as A.3 + A.4 + A.5 behind
//! the `POPKERN_APEX_PATH_A=1` env flag.  A.7's flag-flip gate
//! requires **a metric apex relocation actually moves** —
//! `validate_all`'s raw error count doesn't (the V↔T coherence it closes
//! doesn't surface there at the kernel's default linear tolerance,
//! and the remaining `InconsistentOrientation` errors are
//! degenerate-normal issues unrelated to V↔T).
//!
//! This module defines `MetricBreakdown`: a per-fixture decomposition
//! of `validate_all` errors into three buckets:
//!
//! * **Structural** — Euler-Poincaré, half-edge pairing/cycles,
//!   vertex consistency, manifold edges, vertex valence.  Must be
//!   zero on any healthy BRep.  Apex relocation doesn't change these.
//! * **Apex-relocation relevant** — `OpenLoop` (loop-closure gap from
//!   curve endpoints disagreeing with vertex positions) and
//!   `ShellGap` (adjacent half-edges' vertex positions diverging
//!   across an edge).  These are the gaps apex V↔T relocation
//!   is designed to close.
//! * **Other geometric** — `InconsistentOrientation` (degenerate
//!   surface normal at parameter midpoint), `CurveSurfaceDeviation`,
//!   `TwistedSurface`, `Precision`, `SelfIntersection`.  Apex relocation
//!   doesn't address these; they need separate work.
//!
//! Combined with the existing Phase 6c-E `cap_outer_g0` metric
//! (≤ 5e-3 m), the A.7 flip gate becomes:
//!   * structural == 0,
//!   * path_a_relevant == 0,
//!   * cap_outer_g0 ≤ 5e-3.
//!
//! `Other geometric` count is tracked but not gated — those go
//! through their own work items.

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/phase_6c_a6_validation_probes.rs"]
mod phase_6c_a6_validation_probes;
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

/// Wrap a fixture build under the apex-relocation flag set.  Returns the
/// post-fillet BRep and the cap face handles.
fn build_with_path_a(
    is_pyramid: bool,
    radii: &[f64],
) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
    build_with_flags(is_pyramid, radii, /* topology_rebuild */ false)
}

/// Variant that also enables `POPKERN_APEX_TOPOLOGY_REBUILD=1`
/// (Phase 6c-C 6→3 collapse for Sphere caps).  Used by the
/// auto-trigger probe to measure whether the rebuild closes the
/// path_a_relevant gaps apex relocation can't reach on collapsed-bridge
/// fixtures.
fn build_with_path_a_and_rebuild(
    is_pyramid: bool,
    radii: &[f64],
) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
    build_with_flags(is_pyramid, radii, /* topology_rebuild */ true)
}

fn build_with_flags(
    is_pyramid: bool,
    radii: &[f64],
    topology_rebuild: bool,
) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
    // POPKERN_STAGE4_WIDE_WIRE_TOL, POPKERN_APEX_FACE_LIFT, and
    // POPKERN_APEX_PATH_A are all default-on since #207 cleanup.
    // POPKERN_SPHERE_FIT_BOUNDARY remains opt-in (production
    // divergence on variable-radius corners).
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let _g5 = if topology_rebuild {
        Some(EnvGuard::set("POPKERN_APEX_TOPOLOGY_REBUILD", "1"))
    } else {
        None
    };

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

/// Decomposition of `validate_all` errors into apex-relocation concern
/// concern buckets.  See module docs for category definitions.
#[derive(Debug, Default)]
pub(crate) struct MetricBreakdown {
    pub structural: usize,
    pub path_a_relevant: usize,
    pub other_geometric: usize,
    /// Worst gap reported by an `OpenLoop` or `ShellGap` error,
    /// in metres.  0.0 when both are absent.  Useful for ratchets
    /// that watch for regressions in magnitude even when count
    /// stays flat.
    pub worst_path_a_gap: f64,
}

/// Bucket a `validate_all` error vector by category.
pub(crate) fn breakdown(errors: &[ValidationError]) -> MetricBreakdown {
    let mut out = MetricBreakdown::default();
    for err in errors {
        match err {
            ValidationError::EulerPoincare { .. }
            | ValidationError::BrokenTwin { .. }
            | ValidationError::BrokenCycle { .. }
            | ValidationError::InconsistentVertex { .. }
            | ValidationError::NonManifoldEdge { .. }
            | ValidationError::LowValence { .. } => out.structural += 1,
            ValidationError::OpenLoop { gap, .. } => {
                out.path_a_relevant += 1;
                if *gap > out.worst_path_a_gap {
                    out.worst_path_a_gap = *gap;
                }
            }
            ValidationError::ShellGap { gap, .. } => {
                out.path_a_relevant += 1;
                if *gap > out.worst_path_a_gap {
                    out.worst_path_a_gap = *gap;
                }
            }
            ValidationError::InconsistentOrientation { .. }
            | ValidationError::CurveSurfaceDeviation { .. }
            | ValidationError::TwistedSurface { .. }
            | ValidationError::Precision { .. }
            | ValidationError::SelfIntersection { .. } => out.other_geometric += 1,
        }
    }
    out
}

/// A.6 baseline: dump per-bucket validate_all decomposition.
/// Ignored — run with `--ignored --nocapture` for inspection.
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

    eprintln!("\n=== A.6 validate_all categorised baseline (POPKERN_APEX_PATH_A=1) ===");
    eprintln!("tolerance.linear() = {:.3e}", Tolerance::default().linear());
    eprintln!(
        "{:<18}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<48}", "");
    for case in cases {
        let (brep, _faces) = build_with_path_a(case.is_pyramid, case.radii);
        let errors = brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "{:<18}  {:>5}  {:>5}  {:>5}  {:>10}",
            case.name,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            if m.worst_path_a_gap > 0.0 {
                format!("{:.3e}", m.worst_path_a_gap)
            } else {
                "—".to_string()
            },
        );
        for (i, err) in errors.iter().enumerate().take(5) {
            eprintln!("    [{i}] {err:?}");
        }
        if errors.len() > 5 {
            eprintln!("    ... ({} more)", errors.len() - 5);
        }
    }
}

/// A.6 cap-kind dispatch lock-in: every uniform-radius fixture must
/// show `path_a_relevant == 0` under apex relocation.  After §6c-X.A.7's
/// auto-dispatch (suppress sphere-fit substitution on uniform
/// corners; let the original rolling-ball topology stand),
/// the uniform fixtures stay clean. Variable-radius and mixed
/// pyramid fixtures are tracked by diagnostic probes while their
/// repeated-apex support-continuation gaps remain open.
#[test]
fn phase_6c_a6_uniform_fixtures_path_a_clean() {
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
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
        },
    ];
    for case in cases {
        let (brep, _faces) = build_with_path_a(case.is_pyramid, case.radii);
        let errors = brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "[A.6] {}: structural={} pathA={} other={} worst_pathA_gap={:.3e}",
            case.name, m.structural, m.path_a_relevant, m.other_geometric, m.worst_path_a_gap,
        );
        assert_eq!(
            m.structural,
            0,
            "{}: structural must be 0; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        assert_eq!(
            m.path_a_relevant,
            0,
            "{}: path_a_relevant must be 0 under auto-dispatch; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        assert!(
            m.other_geometric <= 1,
            "{}: other_geometric regression {} (was ≤1):\n{}",
            case.name,
            m.other_geometric,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
}
