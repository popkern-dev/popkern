//! PK-222: conic / variable-profile cube-corner vertex blends and the
//! legacy N-form dispatch gate.
//!
//! Uniform conic and variable-radius cube corners collapse their
//! arc+bridge pairs during Phase 2a, so the blend face physically
//! carries the N = 3 half-edge loop. The legacy dispatch demanded the
//! 2N form there ("outer_half_edges=Some(3), expected 3 or 6") and
//! starved the whole rational route: these corners shipped the
//! polar-NURBS heuristic cap with NO record that the stronger path
//! was never attempted.
//!
//! The gate now accepts the N = 3 N-form (poke normalizes to 2N by
//! midpoint splits), and — because this arm is new dispatch reach
//! with no fixture history — its commits answer to the same
//! cap/support continuity acceptance the N >= 5 routes pass
//! (5e-3 / 5e-2). On today's mitred-ring topology the committed
//! sector complex measures cap_outer_g0 ~1e-1 and the dispatch
//! DECLINES with that fact recorded; the serve falls back to the
//! heuristic cap exactly as before, now with the full decline chain
//! in the diagnostics. The recorded next mechanism is the
//! mitred-ring -> hex-ring substrate (anchor vertices carry
//! junction-incompatible incidence: measured 6.4e-2 tear on the
//! second incident strip chord).

#[path = "kp_gate_helpers.rs"]
mod helpers;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_fillet::fillet_edges;
use popkern_fillet::types::{
    ApexConstructionDiagnostic, BlendProfile, FilletOptions, FilletRadius, FilletSpec,
};

use helpers::incident_edges;

fn cube_corner_result(
    profile: BlendProfile,
    radius_fn: impl Fn(usize) -> FilletRadius,
) -> popkern_fillet::BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let vertex = brep.vertices.keys().next().unwrap();
    let edges: Vec<EdgeHandle> = incident_edges(&brep, vertex);
    assert_eq!(edges.len(), 3);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .enumerate()
        .map(|(k, &e)| FilletSpec {
            edge: e,
            radius: radius_fn(k),
        })
        .collect();
    let options = FilletOptions {
        profile,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &options, tol, None).expect("cube corner fillet succeeds")
}

fn exact_boundary_details<'a>(
    result: &'a popkern_fillet::BlendResult,
) -> impl Iterator<Item = (&'a str, &'a str)> {
    result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::ExactBoundaryCap(d) => {
                Some((d.stage.as_str(), d.detail_code.as_deref().unwrap_or("")))
            }
            _ => None,
        })
}

/// The filed defect: the rational dispatch must never again decline a
/// conic cube corner as "not dispatch-ready" for carrying the
/// physical N = 3 loop. The fixture is the flat-rho small-radius
/// corner, which declines the (now-production, PK-225) setback rung
/// upstream and genuinely reaches the rational dispatch; there the
/// dispatch runs, measures, and either ships a complex that beats
/// the commercial bars or declines with the measured cap/support
/// fact. (The uniform rho=0.5 corner serves the setback fan before
/// the dispatch — pinned in setback_n3_corner_serve.rs.)
#[test]
fn conic_cube_corner_dispatch_is_measured_not_starved() {
    let result = cube_corner_result(BlendProfile::Conic { rho: 0.35 }, |_| {
        FilletRadius::Constant(0.1)
    });

    let starved = exact_boundary_details(&result).any(|(stage, detail)| {
        stage == "rational_gregory.dispatch.precondition_failed"
            && detail.contains("not dispatch-ready")
    });
    assert!(
        !starved,
        "N=3 N-form loop must be accepted by the dispatch gate, not reported dispatch-unready"
    );

    let measured_or_shipped = exact_boundary_details(&result)
        .any(|(stage, _)| stage == "rational_gregory.dispatch.cap_support_acceptance")
        || result.corner_patches.len() == 3
        || result
            .corner_patches
            .iter()
            .any(|cp| cp.sub_patch_provenance.is_some());
    assert!(
        measured_or_shipped,
        "the dispatch must either ship the sector complex or record the measured \
         cap/support decline; diagnostics: {:?}",
        result.repair_report.apex_diagnostics
    );
}

/// Serve contract for the conic / variable classes: structurally
/// sound, and the geometric ceiling is pinned so silent regressions
/// fail. Since PK-225/PK-226 these classes serve the KP20 setback fan
/// (validate-clean, deviation ~0); the 6e-2 ceiling remains as the
/// regression backstop for the legacy-fallback class.
#[test]
fn conic_and_variable_cube_corners_serve_structurally_sound() {
    for (name, profile, radius_fn) in [
        (
            "cube_uniform_conic",
            BlendProfile::Conic { rho: 0.5 },
            Box::new(|_: usize| FilletRadius::Constant(0.15)) as Box<dyn Fn(usize) -> FilletRadius>,
        ),
        (
            "cube_variable",
            BlendProfile::Circular,
            Box::new(|k: usize| FilletRadius::Variable {
                start: 0.10 + 0.02 * k as f64,
                end: 0.16 - 0.02 * k as f64,
            }),
        ),
        (
            "cube_variable_conic",
            BlendProfile::Conic { rho: 0.5 },
            Box::new(|k: usize| FilletRadius::Variable {
                start: 0.10 + 0.02 * k as f64,
                end: 0.16 - 0.02 * k as f64,
            }),
        ),
    ] {
        let result = cube_corner_result(profile, radius_fn);
        let mut structural = 0usize;
        let mut worst_deviation = 0.0f64;
        for error in result.brep.validate_all() {
            use popkern_core::error::ValidationError;
            match error {
                ValidationError::CurveSurfaceDeviation { max_deviation, .. } => {
                    worst_deviation = worst_deviation.max(max_deviation);
                }
                ValidationError::EulerPoincare { .. }
                | ValidationError::DanglingReference { .. }
                | ValidationError::BrokenTwin { .. }
                | ValidationError::BrokenCycle { .. }
                | ValidationError::InconsistentVertex { .. }
                | ValidationError::NonManifoldEdge { .. }
                | ValidationError::OpenLoop { .. }
                | ValidationError::ShellGap { .. } => structural += 1,
                _ => {}
            }
        }
        assert_eq!(structural, 0, "{name}: structural validation must be clean");
        assert!(
            worst_deviation <= 6.0e-2,
            "{name}: worst curve-surface deviation {worst_deviation:.3e} exceeds the recorded \
             mitred-ring class ceiling 6e-2"
        );
    }
}

/// The circular uniform corner keeps its exact rolling-ball sphere
/// serve — regression guard for the corner-fragment absorption and
/// contact-collapse generalizations landing with PK-222.
#[test]
fn circular_uniform_cube_corner_keeps_sphere_serve() {
    let result = cube_corner_result(BlendProfile::Circular, |_| FilletRadius::Constant(0.15));
    assert!(
        result
            .corner_patches
            .iter()
            .any(|cp| cp.method == popkern_fillet::types::CornerMethod::RollingBallSphere),
        "uniform circular corner must keep the rolling-ball sphere method"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "uniform circular corner must validate clean, got {errors:?}"
    );
}
