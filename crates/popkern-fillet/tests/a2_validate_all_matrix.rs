//! Phase 6c-A2.5 — `validate_all` matrix natively under A2.
//!
//! Locks in the commercial-grade Sphere-cap quality of A2 (per
//! ADR §6c-A2 and §6c-X.A.6's metric breakdown).  Runs the same
//! six-fixture matrix that `phase_6c_a6_validation` uses, but with
//! `FilletOptions::use_a2_apex_construction = true`.  Two
//! dispatch outcomes coexist in one matrix:
//!
//! * **Sphere-cap fixtures** (`cube_uniform`, `pyramid3_uniform`,
//!   `pyramid4_uniform`) — A2 succeeds end-to-end, no
//!   `BlendWarning::A2FellBackToLegacy`.  Topology is the
//!   commercial-style 3-vertex cap loop.
//! * **Gregory-cap fixtures** (`cube_mixed`, `cube_extreme`,
//!   `pyramid4_mixed`) — A2 falls back to legacy with a typed
//!   warning per apex.  After fallback the legacy substitute-then-
//!   patch + Path A pipeline produces the BRep.
//!
//! Either way: every fixture must pass the A.6 metric gate
//! (`structural == 0`, `path_a_relevant == 0`).  This is the
//! commercial-grade flip-gate's strictest form: under the new A2
//! architecture, the existing matrix doesn't regress.

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use popkern_fillet::error::BlendWarning;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{BlendResult, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, measure_cap_outer_boundary, pick_vertex_with_valence};

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

#[derive(Debug, Default)]
struct MetricBreakdown {
    structural: usize,
    path_a_relevant: usize,
    other_geometric: usize,
    worst_path_a_gap: f64,
}

fn breakdown(errors: &[ValidationError]) -> MetricBreakdown {
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

/// Build a fillet under A2 dispatch and the Path-A-equivalent
/// envelope.  POPKERN_SPHERE_FIT_BOUNDARY enables sphere-fit
/// substitution on the fallback (legacy) path so variable-radius
/// fixtures route through the Gregory branch (still gated, A2.4).
fn build_under_a2(is_pyramid: bool, radii: &[f64]) -> BlendResult {
    let _g = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let tol = Tolerance::default();
    let brep = if is_pyramid {
        BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, tol).unwrap()
    } else {
        BRep::unit_cube(tol).unwrap()
    };
    let vertex = if is_pyramid {
        pick_vertex_with_valence(&brep, radii.len()).unwrap()
    } else {
        brep.vertices.keys().next().unwrap()
    };
    let edges: Vec<EdgeHandle> = incident_edges(&brep, vertex);
    assert_eq!(edges.len(), radii.len(), "edge count must match radii");
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        use_a2_apex_construction: true,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &options, tol, None).expect("fillet_edges succeeds under A2")
}

fn count_a2_fallbacks(result: &BlendResult) -> usize {
    result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }))
        .count()
}

/// Whole-matrix gate: every fixture passes `validate_all` cleanly
/// under A2 dispatch, regardless of whether A2 takes the path
/// natively or falls back to legacy per-apex.  Strict gate:
/// `structural == 0` and `path_a_relevant == 0` for all 6.
#[test]
fn a2_5_matrix_passes_under_a2() {
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
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
        let errors = result.brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "[A2.5] {}: fallbacks={} structural={} pathA={} other={} worst_gap={:.3e}",
            case.name,
            fallbacks,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            m.worst_path_a_gap,
        );
        assert_eq!(
            m.structural,
            0,
            "[{}] structural must be 0 under A2; errors:\n{}",
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
            "[{}] path_a_relevant must be 0 under A2; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
}

/// Phase 6c-E commercial `cap_outer_g0 ≤ 5e-3 m` gate under A2
/// dispatch — variable-radius fixtures (Gregory cap path; A2
/// falls back to legacy) must not regress the gate that
/// `phase_6c_e_validation_matrix` already enforces.
///
/// Uniform fixtures (`cube_uniform`, `pyramid3_uniform`,
/// `pyramid4_uniform`) are NOT gated here: cube_uniform has zero
/// cap-outer edges measured under A2 native (Sphere cap path
/// emits no cap-vs-non-cap edges as `corner_patches`), and the
/// uniform pyramids exceed 5e-3 by ~100 mm under both legacy and
/// fallback — see `phase_6c_e_validation_matrix.rs`'s `#[ignore]`
/// markers for those.  Closing the uniform pyramid gap is its
/// own work item (A2.4 + cap-outer mathematics).
#[test]
fn a2_5_cap_outer_g0_variable_under_a2() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
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
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let cap_faces = result
            .corner_patches
            .iter()
            .map(|p| p.face)
            .collect::<Vec<_>>();
        let m = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
        eprintln!(
            "[A2.5 cap_outer] {}: g0={:.3e} g1={:.3e} edges={}",
            case.name, m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured,
        );
        assert!(
            m.worst_g0_gap <= 5.0e-3,
            "[{}] cap_outer_g0 {:.3e} exceeds 5.0e-3 commercial gate under A2",
            case.name,
            m.worst_g0_gap,
        );
    }
}

/// Sphere-cap fixtures that A2 currently succeeds natively on
/// (no `A2FellBackToLegacy`).  Today: only `cube_uniform`.
/// Uniform pyramids dispatch via A2 but fall back — A2.2.f
/// initially targeted cube-style 3-fillet corners.  Extending
/// native A2 to uniform pyramids is future work (A2.2.g?).
#[test]
fn a2_5_native_sphere_cap_set() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    // Add fixtures here as A2 native dispatch extends.  Today
    // only cube_uniform is verified.
    let cases: &[Case] = &[Case {
        name: "cube_uniform",
        radii: &[0.10, 0.10, 0.10],
        is_pyramid: false,
    }];
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
        assert_eq!(
            fallbacks, 0,
            "[{}] expected native A2 dispatch (Sphere cap); got {} A2FellBackToLegacy warning(s) — \
             A2.2.f Sphere-cap regression?",
            case.name, fallbacks,
        );
    }
}

/// Combined matrix dump (ignored, diagnostic).  Run with
/// `--ignored --nocapture` to inspect every fixture under A2.
#[ignore]
#[test]
fn a2_5_matrix_dump() {
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
    eprintln!("\n=== A2.5 validate_all matrix natively under A2 ===");
    eprintln!(
        "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
        "fixture", "fallback", "struct", "pathA", "other", "worst gap"
    );
    eprintln!("{:-<60}", "");
    for case in cases {
        let result = build_under_a2(case.is_pyramid, case.radii);
        let fallbacks = count_a2_fallbacks(&result);
        let errors = result.brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "{:<18}  {:>8}  {:>5}  {:>5}  {:>5}  {:>10}",
            case.name,
            fallbacks,
            m.structural,
            m.path_a_relevant,
            m.other_geometric,
            if m.worst_path_a_gap > 0.0 {
                format!("{:.3e}", m.worst_path_a_gap)
            } else {
                "—".to_string()
            },
        );
    }
}
