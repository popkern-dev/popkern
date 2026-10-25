use super::*;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::thermal::{ConductionProblem, solve_conduction};

/// A slab x in [0, 10] mm held at `t1` on x = 0 and losing heat from
/// x = 10 by `last`, every other face insulated.
fn slab(t1: f64, last: FaceCondition) -> ConductionProblem {
    let mut faces = vec![FaceCondition::Insulated; 6];
    faces[0] = FaceCondition::Fixed {
        t0_c: t1,
        gradient: [0.0; 3],
    };
    faces[1] = last;
    ConductionProblem {
        conductivity: 0.2,
        faces,
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    }
}

/// The far face's temperature from the exact balance
/// k (t1 − T)/L = h (T − T∞) + εσ(T⁴ − T∞⁴), by bisection.
fn exact(k: f64, l_m: f64, t1: f64, h: f64, emittance: f64, ambient: f64) -> f64 {
    let f = |t: f64| {
        let (tk, ak) = (t - ABSOLUTE_ZERO_C, ambient - ABSOLUTE_ZERO_C);
        k * (t1 - t) / l_m
            - h * (t - ambient)
            - emittance * STEFAN_BOLTZMANN * (tk.powi(4) - ak.powi(4))
    };
    let (mut lo, mut hi) = (ambient, t1);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

/// Measured 2026-10-12: the far face 219.16094 °C, 1.5e-12 °C from the
/// exact balance, in 5 passes (13 with the secant coefficient); one
/// linearisation at the ambient reads 272.13 °C.
#[test]
fn a_radiating_slab_reaches_the_exact_nonlinear_balance() {
    let model = model_of(box_mesh([0.0; 3], [10.0, 4.0, 4.0]), 6);
    let (t1, h, emittance, ambient) = (400.0, 5.0, 0.9, 20.0);
    let options = AnalysisOptions {
        target_cell_size: Some(1.0),
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let radiating = FaceCondition::Radiating {
        h,
        emittance,
        ambient_c: ambient,
    };
    let r =
        solve_conduction(&model, &slab(t1, radiating), &options).unwrap_or_else(|e| panic!("{e}"));
    let want = exact(0.2, 0.01, t1, h, emittance, ambient);
    let far = r
        .field
        .temperature_at([10.0, 2.0, 2.0])
        .ok()
        .flatten()
        .unwrap_or(f64::NAN);
    assert!((far - want).abs() <= 1e-9 * want, "{far} vs {want}");
    assert!(r.report.radiation_passes > 1 && r.report.radiation_passes < 20);
    let heat = &r.heat;
    assert!(heat.imbalance().abs() <= 1e-9 * heat.throughput());
    // Negative control: one linearisation at the ambient.
    let once = FaceCondition::Convection {
        h: h + tangent(ambient, ambient, emittance).0,
        ambient_c: ambient,
    };
    let r = solve_conduction(&model, &slab(t1, once), &options).unwrap_or_else(|e| panic!("{e}"));
    let far1 = r
        .field
        .temperature_at([10.0, 2.0, 2.0])
        .ok()
        .flatten()
        .unwrap_or(f64::NAN);
    assert!(
        (far1 - want).abs() > 0.05 * (want - ambient),
        "{far1} vs {want}"
    );
}
