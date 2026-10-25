use super::from_setup::solve_thermal;
use super::*;
use crate::face_ref::{Expect, FaceRef, Selector, Tolerance};
use crate::mesh_test_helpers::{box_mesh, model_of, rotated, tube};
use crate::report::ResolvedFace;
use crate::report::Warning;
use crate::setup::FaceDef;
use crate::thermal_setup::{FixedTemperature, ThermalSetup};

/// Gate options: a tight CG tolerance so exactness is visible.
fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 500,
        ..AnalysisOptions::default()
    }
}

fn solve(model: &AnalysisModel, problem: &ConductionProblem, h: f64) -> ThermalResult {
    let r = solve_conduction(model, problem, &options(h)).unwrap_or_else(|e| panic!("{e}"));
    let scale = r.heat.throughput();
    assert!(
        r.heat.imbalance().abs() <= 1e-9 * scale.max(1e-300),
        "heat imbalance {} of {scale}",
        r.heat.imbalance()
    );
    r
}

fn fixed(t0_c: f64) -> FaceCondition {
    FaceCondition::Fixed {
        t0_c,
        gradient: [0.0; 3],
    }
}

fn max_vertex_error(
    r: &ThermalResult,
    model: &AnalysisModel,
    exact: impl Fn([f64; 3]) -> f64,
) -> f64 {
    let mut worst = 0.0f64;
    for (s, temps) in r.vertex_temperature.iter().enumerate() {
        for (v, t) in model.solids[s].vertices.iter().zip(temps) {
            worst = worst.max((t - exact(*v)).abs());
        }
    }
    worst
}

#[test]
fn t1_a_slab_between_fixed_faces_is_exactly_linear() {
    let model = model_of(box_mesh([0.0; 3], [10.0, 4.0, 4.0]), 6);
    let mut faces = vec![FaceCondition::Insulated; 6];
    faces[0] = fixed(100.0);
    faces[1] = fixed(20.0);
    let problem = ConductionProblem {
        conductivity: 167.0,
        faces,
        source_w: 0.0,
        hdt_c: None,
        body_sources: Vec::new(),
    };
    let r = solve(&model, &problem, 0.7);
    let err = max_vertex_error(&r, &model, |x| 100.0 - 8.0 * x[0]);
    assert!(err <= 1e-8 * 80.0, "max error {err}");
    // Heat through the slab: k A ΔT / L, in W (mm → m).
    let q = 167.0 * (16e-6) * 80.0 / 10e-3;
    let into = -r.heat.face_out[0];
    let out = r.heat.face_out[1];
    assert!((into - q).abs() <= 1e-8 * q, "into face 0: {into} vs {q}");
    assert!((out - q).abs() <= 1e-8 * q, "out of face 1: {out} vs {q}");
    // Heat flux q = −k∇T = 167 W/(m·K) × 8 K/mm along +x, in W/m².
    let flux = 167.0 * 8.0 * 1e3;
    for v in r.vertex_heat_flux.iter().flatten() {
        let e = (v[0] - flux).abs() + v[1].abs() + v[2].abs();
        assert!(e <= 1e-8 * flux, "vertex flux {v:?} vs {flux}");
    }
    assert!((r.max_heat_flux_w_per_m2 - flux).abs() <= 1e-8 * flux);
}

#[test]
fn t2_a_slab_with_convection_is_exactly_linear() {
    let (k, h_c, t0, t_inf, l) = (20.0, 500.0, 90.0, 20.0, 8.0);
    let model = model_of(box_mesh([0.0; 3], [l, 3.0, 3.0]), 6);
    let mut faces = vec![FaceCondition::Insulated; 6];
    faces[0] = fixed(t0);
    faces[1] = FaceCondition::Convection {
        h: h_c,
        ambient_c: t_inf,
    };
    let problem = ConductionProblem {
        conductivity: k,
        faces,
        source_w: 0.0,
        hdt_c: None,
        body_sources: Vec::new(),
    };
    let r = solve(&model, &problem, 0.6);
    // k (t0 − ts)/L = h (ts − t∞), SI lengths.
    let lm = l * 1e-3;
    let ts = (k * t0 / lm + h_c * t_inf) / (k / lm + h_c);
    let err = max_vertex_error(&r, &model, |x| t0 + (ts - t0) * x[0] / l);
    assert!(err <= 1e-8 * (t0 - t_inf), "max error {err}");
}

#[test]
fn lf2_a_linear_field_on_all_faces_is_reproduced_at_five_offsets() {
    let gradient = [0.3, -0.2, 0.7];
    let field = |x: [f64; 3]| 5.0 + gradient[0] * x[0] + gradient[1] * x[1] + gradient[2] * x[2];
    let model = model_of(rotated(box_mesh([0.0; 3], [3.0; 3]), 17.0, [1.5; 3]), 6);
    let problem = ConductionProblem {
        conductivity: 50.0,
        faces: vec![
            FaceCondition::Fixed {
                t0_c: 5.0,
                gradient
            };
            6
        ],
        source_w: 0.0,
        hdt_c: None,
        body_sources: Vec::new(),
    };
    let mut dofs = Vec::new();
    for shift in [
        [0.0; 3],
        [0.13, 0.0, 0.0],
        [0.0, 0.29, 0.07],
        [0.31, 0.17, 0.23],
        [0.37; 3],
    ] {
        let o = AnalysisOptions {
            grid_shift: shift,
            ..options(0.4)
        };
        let r = solve_conduction(&model, &problem, &o).unwrap_or_else(|e| panic!("{e}"));
        assert!(r.heat.imbalance().abs() <= 1e-9 * r.heat.throughput());
        let err = max_vertex_error(&r, &model, field);
        assert!(err <= 1e-8 * 3.0, "shift {shift:?}: max error {err}");
        let inside = r
            .field
            .temperature_at([1.5; 3])
            .unwrap_or_else(|e| panic!("{e}"))
            .unwrap_or(f64::NAN);
        assert!((inside - field([1.5; 3])).abs() < 1e-7);
        // Identical discretisations solve bit-identically, so the final
        // residual's bits fingerprint the cut pattern.
        dofs.push((
            r.report.dof,
            r.report.cells_cut,
            r.report.residual_history.last().map(|v| v.to_bits()),
        ));
    }
    dofs.sort_unstable();
    dofs.dedup();
    assert!(
        dofs.len() >= 4,
        "the shifts must change the discretisation: {dofs:?}"
    );
}

#[test]
fn t3_a_tube_has_the_logarithmic_profile() {
    let (ri, ro, ti, to) = (2.0, 6.0, 100.0, 20.0);
    let model = model_of(tube(ri, ro, 4.0, 96), 3);
    let mut faces = vec![FaceCondition::Insulated; 3];
    faces[0] = fixed(ti);
    faces[1] = fixed(to);
    let k = 30.0;
    let problem = ConductionProblem {
        conductivity: k,
        faces,
        source_w: 0.0,
        hdt_c: None,
        body_sources: Vec::new(),
    };
    let r = solve(&model, &problem, 1.0);
    let exact = |rad: f64| ti + (to - ti) * (rad / ri).ln() / (ro / ri).ln();
    for rad in [3.0, 4.0, 5.0] {
        for angle in [0.3, 1.9, 4.0] {
            let p = [rad * f64::cos(angle), rad * f64::sin(angle), 2.0];
            let t = r
                .field
                .temperature_at(p)
                .unwrap_or_else(|e| panic!("{e}"))
                .unwrap_or(f64::NAN);
            assert!(
                (t - exact(rad)).abs() <= 0.02 * (ti - to),
                "r {rad}: {t} vs {}",
                exact(rad)
            );
        }
    }
    // Heat flow 2πkLΔT/ln(ro/ri), W.
    let q = 2.0 * std::f64::consts::PI * k * 4e-3 * (ti - to) / (ro / ri).ln();
    // Faceting (96 segments) and the h = 1 grid: within 2%.
    let into_inner = -r.heat.face_out[0];
    assert!((into_inner - q).abs() <= 0.02 * q, "{into_inner} vs {q}");
    // Radial heat flux at the bore vertices against k ΔT / (r ln(ro/ri)),
    // recovered like stresses. Measured: 0.789..0.907 at 4 cells through
    // the wall and 0.872..0.955 at 8 (the boundary under-read of a 1/r
    // field, as C5's hoop stress); max |q| 0.936 and 0.955 of the wall
    // value. Gated: within 25% at 4 cells, 15% at 8, and converging.
    let wall = 1e3 * k * (ti - to) / (ri * (ro / ri).ln());
    let bore = |r: &ThermalResult| -> (f64, f64) {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for (v, q) in model.solids[0].vertices.iter().zip(&r.vertex_heat_flux[0]) {
            let rad = v[0].hypot(v[1]);
            if (rad - ri).abs() < 1e-6 {
                let radial = (q[0] * v[0] + q[1] * v[1]) / rad / wall;
                lo = lo.min(radial);
                hi = hi.max(radial);
            }
        }
        (lo, hi)
    };
    let fine = solve(&model, &problem, 0.5);
    let ((lo4, hi4), (lo8, hi8)) = (bore(&r), bore(&fine));
    assert!(lo4 >= 0.75 && hi4 <= 1.05, "4 cells: {lo4}..{hi4}");
    assert!(lo8 >= 0.85 && hi8 <= 1.05, "8 cells: {lo8}..{hi8}");
    assert!(1.0 - lo8 < 0.7 * (1.0 - lo4), "{lo4} -> {lo8}");
}

#[test]
fn t4_a_pin_fin_matches_the_fin_equation() {
    // Al k = 200, 2 mm square, L = 100 mm, h = 10: m = 10 /m, mL = 1.
    let model = model_of(box_mesh([0.0; 3], [100.0, 2.0, 2.0]), 6);
    let mut faces = vec![
        FaceCondition::Convection {
            h: 10.0,
            ambient_c: 0.0
        };
        6
    ];
    faces[0] = fixed(100.0);
    faces[1] = FaceCondition::Insulated;
    let problem = ConductionProblem {
        conductivity: 200.0,
        faces,
        source_w: 0.0,
        hdt_c: None,
        body_sources: Vec::new(),
    };
    let r = solve(&model, &problem, 0.5);
    let tip = r
        .field
        .temperature_at([100.0, 1.0, 1.0])
        .unwrap_or_else(|e| panic!("{e}"))
        .unwrap_or(f64::NAN);
    let tip_exact = 100.0 / 1f64.cosh();
    assert!(
        (tip - tip_exact).abs() <= 0.03 * tip_exact,
        "tip {tip} vs {tip_exact}"
    );
    // Base heat flow √(kAhP) θ₀ tanh(mL), all into the fin through the base.
    let q_exact = (200.0f64 * 4e-6 * 10.0 * 8e-3).sqrt() * 100.0 * 1f64.tanh();
    let q_base = -r.heat.face_out[0];
    assert!(
        (q_base - q_exact).abs() <= 0.03 * q_exact,
        "base {q_base} vs {q_exact}"
    );
}

#[test]
fn hw1_a_hot_abs_block_warns_and_a_cool_one_does_not() {
    let model = model_of(box_mesh([0.0; 3], [6.0, 4.0, 4.0]), 6);
    let resolved = vec![ResolvedFace {
        name: "hot".to_string(),
        faces: vec![crate::model::FaceId(0)],
    }];
    for (t, warns) in [(90.0, true), (50.0, false)] {
        let setup = ThermalSetup {
            material: "ABS (FDM)".to_string(),
            tolerance: Tolerance::default(),
            faces: vec![FaceDef {
                name: "hot".to_string(),
                reference: FaceRef::Select(Selector::NearestTo([0.0, 2.0, 2.0])),
                expect: Expect::default(),
            }],
            fixed: vec![FixedTemperature {
                faces: vec!["hot".to_string()],
                temperature_c: t,
            }],
            flux: Vec::new(),
            convection: Vec::new(),
            sources: Vec::new(),
            transient: None,
            options: options(0.8),
        };
        let r = solve_thermal(&model, &setup, &resolved).unwrap_or_else(|e| panic!("{e}"));
        let warned = r
            .report
            .warnings
            .iter()
            .any(|w| matches!(w, Warning::ServiceTemperature { hdt_c, .. } if *hdt_c == 86.6));
        assert_eq!(warned, warns, "T = {t}");
        assert_eq!(r.report.resolved_faces, resolved);
    }
}
