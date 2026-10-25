//! Vector multigrid gates (PH-025): C2 iterations and accuracy, C5
//! (Lamé, displacement), and E2/X8 (AgFEM γ over sliver configurations,
//! with AgFEM off as the negative control).

use super::assemble::FaceMechanics;
use super::thin_wall_tests::{rotate_vector, tip};
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of, quarter_tube};

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);

fn clamp() -> FaceMechanics {
    FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    })
}

/// PCG iterations to 1e-8 for the 10 × 1 × 1 cantilever.
fn iterations(cells: f64, degrees: f64) -> usize {
    let centre = [0.0, 0.5, 0.5];
    let model = model_of(
        crate::mesh_test_helpers::rotated(box_mesh([0.0; 3], [10.0, 1.0, 1.0]), degrees, centre),
        6,
    );
    let down = rotate_vector([0.0, 0.0, -1.0], degrees);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = clamp();
    faces[1] = FaceMechanics::Traction(down.map(|d| 1e-3 * d));
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(1.0 / cells),
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    solve_elastic(&model, &problem, &options)
        .unwrap_or_else(|e| panic!("{e}"))
        .report
        .cg_iterations
}

#[test]
fn c2_cantilever_accuracy_and_mesh_independent_iterations() {
    // Timoshenko: PL³/3EI + PL/κGA, E = 1, ν = 0.3, 1 × 1 section.
    let exact = 1e-3 * 1000.0 / (3.0 / 12.0) + 1e-3 * 10.0 / (5.0 / 6.0 / 2.6);
    let tip2 = tip(2.0, 0.3, 0.0);
    assert!((tip2 - exact).abs() <= 0.05 * exact, "{tip2} vs {exact}");
    for degrees in [0.0, 17.0] {
        let its: Vec<usize> = [2.0, 4.0, 8.0]
            .iter()
            .map(|&c| iterations(c, degrees))
            .collect();
        assert!(its.iter().all(|&i| i < 60), "{degrees}: {its:?}");
    }
}

/// Worst relative error of the inner-surface radial displacement of
/// the plane-strain quarter tube against Lamé, and the PCG iterations.
fn lame_inner_error(h: f64) -> (f64, usize) {
    let (ri, ro, p, e, nu) = (2.0, 6.0, 1e-3, 1.0, 0.3);
    let model = model_of(quarter_tube(ri, ro, 2.0, 48), 5);
    // Plane strain (rollers on both ends), symmetry rollers on the cut
    // planes, pressure inside.
    let faces = vec![
        FaceMechanics::Pressure(p),
        FaceMechanics::Free,
        ROLLER,
        ROLLER,
        ROLLER,
    ];
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-10,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
    let a = p * ri * ri / (ro * ro - ri * ri);
    let u_r = (1.0 + nu) / e * ((1.0 - 2.0 * nu) * a * ri + a * ro * ro / ri);
    let mut worst = 0.0f64;
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        let rad = (v[0] * v[0] + v[1] * v[1]).sqrt();
        if (rad - ri).abs() < 1e-9 {
            worst = worst.max(((d[0] * v[0] + d[1] * v[1]) / rad - u_r).abs() / u_r);
        }
    }
    (worst, r.report.cg_iterations)
}

/// C5, displacement part: the inner radial displacement converges to
/// Lamé at second order (measured 7.4% / 1.5% / 0.4% at 4 / 8 / 16
/// cells through the wall). The stress criterion (σ_θθ(R_i) within 5%
/// at 4 cells) needs recovery and is gated in PH-026.
#[test]
fn c5_a_thick_cylinder_under_pressure_converges_to_lame() {
    let (coarse, _) = lame_inner_error(1.0);
    let (fine, its) = lame_inner_error(0.5);
    assert!(fine <= 0.02, "8 cells: {fine}");
    assert!(coarse / fine >= 3.0, "order: {coarse} -> {fine}");
    assert!(its < 60, "{its}");
}

/// Tip deflection and iterations of a 4 × 1 × 1 cantilever on cells of
/// size `h` with AgFEM threshold `gamma`.
fn sliver_case(h: f64, gamma: f64) -> Result<(f64, usize), AnalysisError> {
    let model = model_of(box_mesh([0.0; 3], [4.0, 1.0, 1.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = clamp();
    faces[1] = FaceMechanics::Traction([0.0, 0.0, -1e-3]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        aggregation_threshold: gamma,
        cg_max_iters: 300,
        coarsest_dofs: 600,
        // The far end and the sides land at a different fraction of a
        // cell for every h: slivers of every depth.
        grid_shift: [0.0, 0.37 * h, 0.61 * h],
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&model, &problem, &options)?;
    let tip = r
        .field
        .displacement_at([4.0, 0.5, 0.5])?
        .map_or(f64::NAN, |d| -d[2]);
    Ok((tip, r.report.cg_iterations))
}

#[test]
fn e2_x8_agfem_keeps_iterations_and_answers_steady_over_slivers() {
    let hs: Vec<f64> = (0..50).map(|i| 0.25 + 0.012 * i as f64 / 49.0).collect();
    for gamma in [0.1, 0.25, 0.5] {
        let mut tips = Vec::new();
        for &h in &hs {
            let (tip, its) =
                sliver_case(h, gamma).unwrap_or_else(|e| panic!("gamma {gamma} h {h}: {e}"));
            assert!(its < 60, "gamma {gamma} h {h}: {its} iterations");
            tips.push(tip);
        }
        let (lo, hi) = tips
            .iter()
            .fold((f64::INFINITY, 0.0f64), |(l, u), &t| (l.min(t), u.max(t)));
        assert!(hi - lo <= 0.05 * lo, "gamma {gamma}: tips {lo}..{hi}");
    }
    // Negative control: without aggregation some sliver case fails to
    // converge within the budget, or needs far more iterations.
    let worst = hs
        .iter()
        .map(|&h| match sliver_case(h, 0.0) {
            Ok((_, its)) => its,
            Err(_) => usize::MAX,
        })
        .max()
        .unwrap_or(0);
    assert!(
        worst >= 60,
        "AgFEM off still converged in {worst} iterations"
    );
}
