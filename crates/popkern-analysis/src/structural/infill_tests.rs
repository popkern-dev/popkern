use super::super::assemble::FaceMechanics;
use super::super::{ElasticProblem, StructuralResult, solve_elastic, solve_infilled};
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::options::AnalysisOptions;

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);
const E: f64 = 2300.0;
const NU: f64 = 0.35;

/// Bean's worked values (notes 11): E/E_m at D = 0.15 … 0.75, and the
/// solid at D = 1.
#[test]
fn bean_fits_match_the_worked_values() {
    for (d, want) in [
        (0.15, 0.059),
        (0.20, 0.078),
        (0.25, 0.098),
        (0.50, 0.245),
        (0.75, 0.507),
    ] {
        let (e, _) = gyroid(d, 1.0, 0.35);
        assert!((e - want).abs() <= 0.0015, "D {d}: {e} vs {want}");
    }
    let (e, nu) = gyroid(1.0, E, NU);
    assert!(
        (e - E).abs() <= 1e-12 * E && (nu - NU).abs() <= 1e-12,
        "{e} {nu}"
    );
    // ν = E/2G − 1 stays a valid isotropic ratio where Bean's cubic ν
    // passes 0.5.
    for d in [0.1, 0.25, 0.5, 0.75] {
        let (_, nu) = gyroid(d, E, NU);
        assert!(nu > 0.0 && nu < 0.5, "D {d}: ν {nu}");
    }
}

/// The bar [0, L] × [0, 10]² on rollers at x = 0, y = 0, z = 0, pulled
/// by `traction` (MPa) at x = L.
fn bar(infill: Option<Infill>) -> (crate::model::AnalysisModel, StructuralResult) {
    let model = model_of(box_mesh([0.0; 3], [40.0, 10.0, 10.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    for f in [0, 2, 4] {
        faces[f] = ROLLER;
    }
    faces[1] = FaceMechanics::Traction([1.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.5),
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let r = match infill {
        Some(f) => solve_infilled(&model, &problem, &f, &options),
        None => solve_elastic(&model, &problem, &options),
    }
    .unwrap_or_else(|e| panic!("{e}"));
    (model, r)
}

/// The loaded face's centre's displacement.
fn end_displacement(r: &StructuralResult) -> f64 {
    r.field
        .displacement_at([40.0, 5.0, 5.0])
        .ok()
        .flatten()
        .map_or(f64::NAN, |u| u[0])
}

/// A full infill is the solid.
#[test]
fn a_full_infill_is_the_solid() {
    let (_, solid) = bar(None);
    let (_, full) = bar(Some(Infill {
        density: 1.0,
        walls_mm: 1.0,
        line_width_mm: 0.4,
    }));
    let (a, b) = (end_displacement(&solid), end_displacement(&full));
    assert!((a - b).abs() <= 1e-10 * a.abs(), "{a} vs {b}");
}

/// Axial stiffness against the walls-and-core composite: away from the
/// ends, walls and core share the strain, σ_nom / E_par with E_par the
/// area-weighted modulus. Measured 2026-10-12 (D 0.25, walls 1 mm, h
/// 0.5): the strain between x = 10 and 30 reads 1.012–1.018 of it
/// across the section. (The loaded face's centre moves 24% more than a
/// series-and-parallel estimate: the 1 mm solid end cap bends over the
/// soft core, which plane sections do not see.)
#[test]
fn an_infilled_bar_is_as_stiff_as_its_walls_and_core() {
    let (d, w) = (0.25, 1.0);
    let (_, r) = bar(Some(Infill {
        density: d,
        walls_mm: w,
        line_width_mm: 0.4,
    }));
    let (e_core, _) = gyroid(d, E, NU);
    let a = 100.0;
    let a_core = (10.0 - 2.0 * w) * (10.0 - 2.0 * w);
    let e_par = ((a - a_core) * E + a_core * e_core) / a;
    let at = |r: &StructuralResult, p: [f64; 3]| {
        r.field
            .displacement_at(p)
            .ok()
            .flatten()
            .map_or(f64::NAN, |u| u[0])
    };
    for (y, z) in [(5.0, 5.0), (0.5, 0.5), (0.5, 5.0), (2.0, 5.0)] {
        let strain = (at(&r, [30.0, y, z]) - at(&r, [10.0, y, z])) / 20.0;
        let ratio = strain * e_par;
        assert!((ratio - 1.0).abs() <= 0.05, "({y}, {z}): {ratio}");
    }
    // Negative control: the solid bar's strain is σ/E, far below.
    let (_, solid) = bar(None);
    let strain = (at(&solid, [30.0, 5.0, 5.0]) - at(&solid, [10.0, 5.0, 5.0])) / 20.0;
    assert!(strain * e_par < 0.5, "{}", strain * e_par);
}

/// PH-078: walls and core are one polymer, so a uniform rise expands an
/// infilled bar freely: u = αθx and no stress, as in the solid.
#[test]
fn a_uniform_rise_expands_an_infilled_bar_freely() {
    use super::super::solve_with;
    use super::super::thermoelastic::ThermalLoad;
    let model = model_of(box_mesh([0.0; 3], [40.0, 10.0, 10.0]), 6);
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
        target_cell_size: Some(0.5),
        cg_rel_tol: 1e-12,
        ..AnalysisOptions::default()
    };
    let (alpha, theta) = (6.8e-5, 40.0);
    let t = |_: [f64; 3]| Some(20.0 + theta);
    let load = ThermalLoad {
        alpha,
        reference_c: 20.0,
        temperature: &t,
    };
    let infill = Infill {
        density: 0.25,
        walls_mm: 1.0,
        line_width_mm: 0.4,
    };
    let (model, r) = (
        &model,
        solve_with(&model, &problem, Some(&load), Some(&infill), &options)
            .unwrap_or_else(|e| panic!("{e}")),
    );
    let worst = model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
        .map(|(p, u)| {
            (0..3)
                .map(|i| (u[i] - alpha * theta * p[i]).abs())
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    assert!(worst <= 1e-8 * alpha * theta * 40.0, "{worst}");
    let peak = r.stresses[0]
        .surface
        .iter()
        .map(|x| x.stress.von_mises)
        .fold(0.0, f64::max);
    assert!(peak <= 1e-8 * E * alpha * theta, "{peak}");
}

/// PH-078: an infilled cantilever (60 × 4 × 4 mm, walls 1 mm, D 0.25)
/// against the composite Euler–Bernoulli beam: EI = E_m I_walls + E* I_core,
/// ρA = ρ (A_walls + D A_core). Measured 2026-10-12 (h 0.5): 261.47 Hz
/// against 263.44, 0.9925.
#[test]
fn an_infilled_cantilever_rings_as_its_composite_beam() {
    use super::super::modal::solve_modal_with;
    let (l, b, w, d) = (60.0, 4.0, 1.0, 0.25);
    let model = model_of(box_mesh([0.0; 3], [l, b, b]), 6);
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
    let options = AnalysisOptions {
        target_cell_size: Some(0.5),
        cg_rel_tol: 1e-12,
        cg_max_iters: 2000,
        ..AnalysisOptions::default()
    };
    let rho = 1240.0;
    let infill = Infill {
        density: d,
        walls_mm: w,
        line_width_mm: 0.4,
    };
    let r = solve_modal_with(&model, &problem, rho, 2, Some(&infill), &options)
        .unwrap_or_else(|e| panic!("{e}"));
    let (e_core, _) = gyroid(d, E, NU);
    let c = b - 2.0 * w;
    let (i_all, i_core) = (b * b * b * b / 12.0, c * c * c * c / 12.0);
    let ei = (E * (i_all - i_core) + e_core * i_core) * 1e6 * 1e-12;
    let rho_a = rho * ((b * b - c * c) + d * c * c) * 1e-6;
    let f1 = 1.875_104_1f64.powi(2) / (2.0 * std::f64::consts::PI * (l * 1e-3).powi(2))
        * (ei / rho_a).sqrt();
    let got = r.modes[0].frequency_hz;
    assert!((got / f1 - 1.0).abs() <= 0.03, "{got} vs {f1}");
}
