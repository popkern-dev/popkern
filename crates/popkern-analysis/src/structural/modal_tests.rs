use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};

const E: f64 = 70_000.0; // MPa
const NU: f64 = 0.33;
const RHO: f64 = 2700.0; // kg/m³

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        cg_max_iters: 2000,
        ..AnalysisOptions::default()
    }
}

/// A cantilever L = 40, b = 2, t = 1 mm clamped at x = 0, 4 cells
/// through t: the bending modes about the thin axis (z motion) against
/// Euler–Bernoulli, ω_n = (β_nL)² √(EI/(ρAL⁴)). Measured 2026-10-12:
/// 1.0026, 0.9998, 0.9955 of the first three (the wide-axis bending,
/// torsion and axial modes in between are told apart by their shapes).
#[test]
fn a_cantilever_bends_at_euler_bernoulli_frequencies() {
    let (l, b, t) = (40.0, 2.0, 1.0);
    let model = model_of(box_mesh([0.0; 3], [l, b, t]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let r =
        solve_modal(&model, &problem, RHO, 10, &options(0.25)).unwrap_or_else(|e| panic!("{e}"));
    // SI: E Pa, I m⁴, A m², L m.
    let (e_si, i_si, a_si, l_si) = (
        E * 1e6,
        b * t * t * t / 12.0 * 1e-12,
        b * t * 1e-6,
        l * 1e-3,
    );
    let base = (e_si * i_si / (RHO * a_si * l_si.powi(4))).sqrt() / (2.0 * std::f64::consts::PI);
    let exact: Vec<f64> = [1.875_104_1, 4.694_091_1, 7.854_757_4]
        .iter()
        .map(|bl| bl * bl * base)
        .collect();
    // Modes whose shape is mostly along z (bending about the thin axis).
    let z_modes: Vec<f64> = r
        .modes
        .iter()
        .zip(&r.shapes)
        .filter(|(_, shape)| {
            let sum = |c: usize| shape[0].iter().map(|d| d[c] * d[c]).sum::<f64>();
            sum(2) > 0.9 * (sum(0) + sum(1) + sum(2))
        })
        .map(|(m, _)| m.frequency_hz)
        .collect();
    let ratios: Vec<f64> = z_modes.iter().zip(&exact).map(|(f, e)| f / e).collect();
    assert_eq!(ratios.len(), 3, "{z_modes:?} vs {exact:?}");
    assert!(ratios.iter().all(|x| (x - 1.0).abs() <= 0.02), "{ratios:?}");
    assert!(
        r.modes.iter().all(|m| m.residual <= MODE_TOLERANCE),
        "{:?}",
        r.modes
    );
}

/// A free 10 mm cube: six rigid modes at (numerically) zero, the first
/// elastic mode converging, and the mass ρV exactly. Measured
/// 2026-10-12: 143.28 kHz at 10 cells across, 142.26 at 20 (0.7%).
#[test]
fn a_free_cube_has_six_rigid_modes_and_its_mass() {
    let model = model_of(box_mesh([0.0; 3], [10.0; 3]), 6);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces: vec![FaceMechanics::Free; 6],
    };
    let first = |h: f64| {
        let r =
            solve_modal(&model, &problem, RHO, 8, &options(h)).unwrap_or_else(|e| panic!("{e}"));
        assert!(
            (r.mass_kg / (RHO * 1e-6) - 1.0).abs() <= 1e-12,
            "mass {}",
            r.mass_kg
        );
        let elastic = r.modes[6].omega_squared;
        for m in &r.modes[..6] {
            assert!(
                m.omega_squared.abs() <= 1e-6 * elastic,
                "rigid {m:?} vs {elastic}"
            );
        }
        r.modes[6].frequency_hz
    };
    let (coarse, fine) = (first(1.0), first(0.5));
    assert!((coarse / fine - 1.0).abs() <= 0.02, "{coarse} vs {fine}");
}

/// The negative control for convergence: a step cap below what the
/// modes need declines rather than reporting unconverged modes.
#[test]
fn too_few_lanczos_steps_decline() {
    let model = model_of(box_mesh([0.0; 3], [10.0; 3]), 6);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces: vec![FaceMechanics::Free; 6],
    };
    let r = solve_modal_capped(
        &model,
        &problem,
        RHO,
        (
            8,
            ModalCaps {
                max_steps: 9,
                ..ModalCaps::default()
            },
        ),
        None,
        &options(2.5),
    );
    assert!(
        matches!(r, Err(AnalysisError::SolverNotConverged { .. })),
        "{r:?}"
    );
    let r = solve_modal(&model, &problem, RHO, 8, &options(2.5)).unwrap_or_else(|e| panic!("{e}"));
    assert!(r.report.lanczos_steps > 9, "{}", r.report.lanczos_steps);
}

/// A free 10 × 6 × 3 brick has six rigid modes. A single Lanczos run
/// sees one direction of a repeated eigenvalue (round-off aside) and
/// found five (measured 2026-10-12, 1 mm cells: five at ≤ 0.0053 and
/// the sixth slot taken by the first elastic mode, 5.4e11): the
/// negative control. Locking them (production) or deflated restarts
/// alone both recover all six.
#[test]
fn a_free_brick_has_six_rigid_modes_by_locking_or_by_restarts() {
    let model = model_of(box_mesh([0.0; 3], [10.0, 6.0, 3.0]), 6);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces: vec![FaceMechanics::Free; 6],
    };
    let rigid_count = |caps: ModalCaps| {
        let r = solve_modal_capped(&model, &problem, RHO, (7, caps), None, &options(1.0))
            .unwrap_or_else(|e| panic!("{e}"));
        let elastic = r.modes[6].omega_squared;
        r.modes
            .iter()
            .filter(|m| m.omega_squared.abs() <= 1e-6 * elastic)
            .count()
    };
    assert_eq!(rigid_count(ModalCaps::default()), 6);
    let restarts_only = ModalCaps {
        lock_rigid: false,
        ..ModalCaps::default()
    };
    assert_eq!(rigid_count(restarts_only), 6);
    let single_run = ModalCaps {
        lock_rigid: false,
        restarts: 0,
        ..ModalCaps::default()
    };
    assert_eq!(rigid_count(single_run), 5);
}

/// A square-section cantilever bends alike about both axes: its first
/// two modes are a repeated pair. Both are reported, equal, and the
/// verifying restart runs (and finds nothing new below the third).
#[test]
fn a_repeated_pair_is_reported_twice() {
    let model = model_of(box_mesh([0.0; 3], [20.0, 2.0, 2.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let r = solve_modal(&model, &problem, RHO, 3, &options(0.5)).unwrap_or_else(|e| panic!("{e}"));
    let w: Vec<f64> = r.modes.iter().map(|m| m.omega_squared).collect();
    assert!((w[1] / w[0] - 1.0).abs() <= 1e-6, "{w:?}");
    assert!(w[2] > 1.5 * w[1], "{w:?}");
    assert_eq!(r.report.lanczos_runs, 2);
}
