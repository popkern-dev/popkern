use super::*;
use crate::mesh_test_helpers::{box_mesh, model_of};

/// Aluminium 6061-T6 (the presets' cited values).
const K: f64 = 167.0;
const RHO: f64 = 2700.0;
const CP: f64 = 896.0;

fn erf(x: f64) -> f64 {
    // Maclaurin series; |x| ≤ 3 here.
    let mut term = x;
    let mut sum = x;
    for n in 1..80 {
        term *= -x * x / n as f64;
        sum += term / (2 * n + 1) as f64;
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

fn problem(
    faces: Vec<FaceCondition>,
    initial_c: f64,
    step_s: f64,
    steps: usize,
) -> TransientProblem {
    TransientProblem {
        conduction: ConductionProblem {
            conductivity: K,
            faces,
            source_w: 0.0,
            body_sources: Vec::new(),
            hdt_c: None,
        },
        density: RHO,
        specific_heat: CP,
        initial_c,
        step_s,
        steps,
    }
}

/// A bar x in [0, 40] mm at 20 °C whose x = 0 end is stepped to 100 °C:
/// T = Ts + (T0 − Ts) erf(x / 2√(αt)) while the heat has not reached the
/// far end. Measured 2026-10-12 at t = 0.5 s (h 0.5 mm, Δt 5 ms): within
/// 0.13% of the step at 1–8 mm depth.
#[test]
fn a_semi_infinite_bar_follows_erfc() {
    let model = model_of(box_mesh([0.0; 3], [40.0, 4.0, 4.0]), 6);
    let mut faces = vec![FaceCondition::Insulated; 6];
    faces[0] = FaceCondition::Fixed {
        t0_c: 100.0,
        gradient: [0.0; 3],
    };
    let (dt, steps) = (0.005, 100);
    let options = AnalysisOptions {
        target_cell_size: Some(0.5),
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let r = solve_transient(&model, &problem(faces, 20.0, dt, steps), &options)
        .unwrap_or_else(|e| panic!("{e}"));
    let alpha = K / (RHO * CP) * 1e6; // mm²/s
    let t = dt * steps as f64;
    let errors: Vec<(f64, f64)> = [1.0, 2.0, 4.0, 6.0, 8.0]
        .iter()
        .map(|&x| {
            let exact = 100.0 + (20.0 - 100.0) * erf(x / (2.0 * (alpha * t).sqrt()));
            let got = r
                .last
                .field
                .temperature_at([x, 2.0, 2.0])
                .ok()
                .flatten()
                .unwrap_or(f64::NAN);
            (x, (got - exact) / 80.0)
        })
        .collect();
    assert!(errors.iter().all(|e| e.1.abs() <= 0.02), "{errors:?}");
    let worst = r.steps.iter().map(|s| s.imbalance).fold(0.0, f64::max);
    assert!(worst <= 1e-6, "energy imbalance {worst}");
}

/// A 10 mm aluminium cube at 100 °C cooling by h = 50 W/(m²·K) to 20 °C
/// (Biot 0.003): T − T∞ = (T0 − T∞) exp(−t/τ), τ = ρcV/(hA). Measured
/// 2026-10-12 at t = τ in 200 steps: 0.37% (backward Euler's own lag is
/// 0.25%).
#[test]
fn a_lumped_cube_cools_exponentially() {
    let model = model_of(box_mesh([0.0; 3], [10.0; 3]), 6);
    let faces = vec![
        FaceCondition::Convection {
            h: 50.0,
            ambient_c: 20.0,
        };
        6
    ];
    let tau = RHO * CP * 1e-6 / (50.0 * 6e-4);
    let steps = 200;
    let options = AnalysisOptions {
        target_cell_size: Some(2.0),
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let r = solve_transient(
        &model,
        &problem(faces, 100.0, tau / steps as f64, steps),
        &options,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let got = r
        .last
        .field
        .temperature_at([5.0; 3])
        .ok()
        .flatten()
        .unwrap_or(f64::NAN);
    let exact = 20.0 + 80.0 * (-1.0f64).exp();
    assert!(
        ((got - 20.0) / (exact - 20.0) - 1.0).abs() <= 0.01,
        "{got} vs {exact}"
    );
    // Energy: what left equals what the cube lost.
    let lost = RHO * CP * 1e-6 * 80.0 * (1.0 - (-1.0f64).exp());
    let stored = r.steps.last().map_or(f64::NAN, |s| s.stored_j);
    assert!((-stored / lost - 1.0).abs() <= 0.01, "{stored} vs {lost}");
    let worst = r.steps.iter().map(|s| s.imbalance).fold(0.0, f64::max);
    assert!(worst <= 1e-6, "energy imbalance {worst}");
}

/// PH-077: a thin aluminium plate (20 × 20 × 1 mm, Biot ≈ 1e-4)
/// radiating from 300 °C (ε 0.9, no convection) to 20 °C against the
/// lumped balance ρcV dT/dt = −εσA(T⁴ − T∞⁴), integrated by RK4 in fine
/// steps. Measured 2026-10-12 after 30 s (150 steps of 0.2 s, radiation
/// relinearised each step): 206.82 °C against 206.60, 1.0012 of the
/// drop's remainder.
#[test]
fn a_radiating_plate_cools_as_the_lumped_balance() {
    use crate::thermal_setup::{ABSOLUTE_ZERO_C, STEFAN_BOLTZMANN};
    let (lx, ly, lz) = (20.0, 20.0, 1.0);
    let model = model_of(box_mesh([0.0; 3], [lx, ly, lz]), 6);
    let (emittance, ambient, start) = (0.9, 20.0, 300.0);
    let faces = vec![
        FaceCondition::Radiating {
            h: 0.0,
            emittance,
            ambient_c: ambient,
        };
        6
    ];
    let (dt, steps) = (0.2, 150);
    let options = AnalysisOptions {
        target_cell_size: Some(0.5),
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let r = solve_transient(&model, &problem(faces, start, dt, steps), &options)
        .unwrap_or_else(|e| panic!("{e}"));
    // The lumped reference (SI).
    let volume = lx * ly * lz * 1e-9;
    let area = 2.0 * (lx * ly + lx * lz + ly * lz) * 1e-6;
    let k = emittance * STEFAN_BOLTZMANN * area / (RHO * CP * volume);
    let a4 = (ambient - ABSOLUTE_ZERO_C).powi(4);
    let rate = |t: f64| -k * ((t - ABSOLUTE_ZERO_C).powi(4) - a4);
    let (mut t, fine) = (start, 1e-3);
    for _ in 0..((dt * steps as f64) / fine).round() as usize {
        let k1 = rate(t);
        let k2 = rate(t + 0.5 * fine * k1);
        let k3 = rate(t + 0.5 * fine * k2);
        let k4 = rate(t + fine * k3);
        t += fine / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    let got = r
        .last
        .field
        .temperature_at([10.0, 10.0, 0.5])
        .ok()
        .flatten()
        .unwrap_or(f64::NAN);
    let ratio = (got - ambient) / (t - ambient);
    assert!((ratio - 1.0).abs() <= 0.01, "{got} vs {t}");
    let worst = r.steps.iter().map(|s| s.imbalance).fold(0.0, f64::max);
    assert!(worst <= 1e-6, "energy imbalance {worst}");
}
