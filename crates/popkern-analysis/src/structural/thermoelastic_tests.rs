use super::super::assemble::FaceMechanics;
use super::super::{ElasticProblem, StructuralResult, solve_thermoelastic};
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::options::AnalysisOptions;

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);
const E: f64 = 1000.0;
const NU: f64 = 0.3;
const ALPHA: f64 = 2.0e-5;

/// The bar [0, 10] × [0, 4] × [0, 3] on rollers at x = 0, y = 0, z = 0
/// (and at x = 10 when `held`), at temperature `t(x)` over 20 °C.
fn bar(
    held: bool,
    t: &dyn Fn([f64; 3]) -> Option<f64>,
) -> (crate::model::AnalysisModel, StructuralResult) {
    let model = model_of(box_mesh([0.0; 3], [10.0, 4.0, 3.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    for f in [0, 2, 4] {
        faces[f] = ROLLER;
    }
    if held {
        faces[1] = ROLLER;
    }
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.7),
        grid_shift: [0.13, 0.29, 0.41],
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let load = ThermalLoad {
        alpha: ALPHA,
        reference_c: 20.0,
        temperature: t,
    };
    let r =
        solve_thermoelastic(&model, &problem, &load, &options).unwrap_or_else(|e| panic!("{e}"));
    (model, r)
}

/// TE1: a uniform rise of 50 K expands the bar freely: u = αθ x exactly
/// and no stress.
#[test]
fn a_free_bar_expands_exactly_and_stress_free() {
    let theta = 50.0;
    let (model, r) = bar(false, &|_| Some(20.0 + theta));
    let worst = model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
        .map(|(p, u)| {
            (0..3)
                .map(|i| (u[i] - ALPHA * theta * p[i]).abs())
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    assert!(worst <= 1e-8 * ALPHA * theta * 10.0, "displacement {worst}");
    let s = E * ALPHA * theta;
    let peak = r.stresses[0]
        .surface
        .iter()
        .map(|x| x.stress.von_mises)
        .fold(0.0, f64::max);
    assert!(peak <= 1e-8 * s, "stress {peak} of {s}");
}

/// TE2: held between rollers at both ends, σ_xx = −Eαθ everywhere and
/// the sides expand by (1 + ν)αθ.
#[test]
fn a_bar_held_at_both_ends_reads_minus_e_alpha_theta() {
    let theta = 50.0;
    let (model, r) = bar(true, &|_| Some(20.0 + theta));
    let want = -E * ALPHA * theta;
    let worst = r.stresses[0]
        .surface
        .iter()
        .map(|x| (x.stress.stress[0] - want).abs())
        .fold(0.0, f64::max);
    assert!(worst <= 1e-8 * want.abs(), "σxx off by {worst} of {want}");
    let lateral = model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
        .map(|(p, u)| (u[1] - (1.0 + NU) * ALPHA * theta * p[1]).abs())
        .fold(0.0, f64::max);
    assert!(lateral <= 1e-8 * ALPHA * theta * 4.0, "{lateral}");
    // The held end carries the thrust Eαθ A (the face's force on its
    // support, outward).
    let thrust = r.face_reaction.get(1).map_or(f64::NAN, |f| f[0]);
    assert!(
        (thrust + want * 12.0).abs() <= 1e-9 * want.abs() * 12.0,
        "{thrust}"
    );
}

/// TE3: a temperature rising 5 K/mm along the bar, on its three
/// symmetry rollers. The x = 0 roller restrains the compatible free
/// field (its u_x varies over the face), so this is not stress-free and
/// its peak sits at that roller's edges (0.356, 0.401, 0.422 of Eα·50
/// K at h = 0.7, 0.35, 0.175: growing). The tip displacement converges
/// at second order. Measured 2026-10-12: tip u_x 4.17209e-3, 4.16688e-3,
/// 4.16551e-3 mm, successive changes in the ratio 3.8.
#[test]
fn a_temperature_gradient_converges() {
    let tip = |h: f64| {
        let model = model_of(box_mesh([0.0; 3], [10.0, 4.0, 3.0]), 6);
        let mut faces = vec![FaceMechanics::Free; 6];
        for f in [0, 2, 4] {
            faces[f] = ROLLER;
        }
        let problem = ElasticProblem {
            youngs: E,
            poisson: NU,
            faces,
        };
        let options = AnalysisOptions {
            target_cell_size: Some(h),
            grid_shift: [0.13, 0.29, 0.41].map(|x| x * h / 0.7),
            cg_rel_tol: 1e-12,
            ..AnalysisOptions::default()
        };
        let t = |p: [f64; 3]| Some(20.0 + 5.0 * p[0]);
        let load = ThermalLoad {
            alpha: ALPHA,
            reference_c: 20.0,
            temperature: &t,
        };
        let r = solve_thermoelastic(&model, &problem, &load, &options)
            .unwrap_or_else(|e| panic!("{e}"));
        r.field
            .displacement_at([10.0, 4.0, 3.0])
            .ok()
            .flatten()
            .map_or(f64::NAN, |u| u[0])
    };
    let (a, b, c) = (tip(0.7), tip(0.35), tip(0.175));
    let ratio = (a - b) / (b - c);
    assert!(ratio >= 3.0, "{a} {b} {c}: ratio {ratio}");
    assert!(((b - c) / c).abs() <= 1e-3, "{b} {c}");
}
