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
    // POPKERN_STAGE4_WIDE_WIRE_TOL, POPKERN_APEX_FACE_LIFT, and
    // POPKERN_APEX_PATH_A are all default-on since #207 cleanup.
    // POPKERN_SPHERE_FIT_BOUNDARY remains opt-in (production
    // divergence on variable-radius corners).
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");

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
            | ValidationError::DanglingReference { .. }
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
            | ValidationError::WindingDishonest { .. }
            | ValidationError::MissingPcurves { .. }
            | ValidationError::CurveSurfaceDeviation { .. }
            | ValidationError::TwistedSurface { .. }
            | ValidationError::Precision { .. }
            | ValidationError::SelfIntersection { .. } => out.other_geometric += 1,
        }
    }
    out
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
        // PK-239: zero. The PK-236 inventory these ceilings carried
        // (pyramid3_uniform = 9 base-corner end-termination chords
        // worst 1.063e-1 + 2 mitred-apex chords; pyramid4_uniform =
        // 12 substituted-boundary cap edges worst 1.5e-2) is closed
        // by the boundary-ring carve: the pyramid apexes and base
        // corners ride the false-sphere-claim gate arm (uniform
        // cylinder strips on mutually oblique axes) and carve exact
        // rings (measured on-surface deviation ≤ 3.8e-14), so every
        // fixture serves validate_all with zero geometric errors.
        max_other_geometric: usize,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
            max_other_geometric: 0,
        },
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
            max_other_geometric: 0,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
            max_other_geometric: 0,
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
            m.other_geometric <= case.max_other_geometric,
            "{}: other_geometric regression {} (measured ceiling {}):\n{}",
            case.name,
            m.other_geometric,
            case.max_other_geometric,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
}
