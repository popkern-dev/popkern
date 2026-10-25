use super::super::assemble::FaceMechanics;
use super::super::recovery_tests::{plate, plate_faces};
use super::super::{ElasticProblem, solve_elastic};
use super::*;
use crate::element::vector_boundary::PieceSupport;

/// C6's converged hole peak (PH-058: 13 to 20 cells per radius).
const CONVERGED: f64 = 3.42334579850262e-3;

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-10,
        cg_max_iters: 400,
        max_active_cells: 2_000_000,
        ..AnalysisOptions::default()
    }
}

/// The C6 plate in tension: model and problem.
fn c6() -> (crate::model::AnalysisModel, ElasticProblem) {
    let (model, segments) = plate();
    let mut faces = plate_faces(model.faces.len());
    for (f, _, n) in &segments {
        if n[0] > 0.5 {
            faces[*f as usize] = FaceMechanics::Traction([1e-3, 0.0, 0.0]);
        }
    }
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    (model, problem)
}

#[test]
fn a_linear_field_is_reproduced_through_the_caps() {
    let (model, _) = plate();
    let g = [[1e-3, 2e-3, 0.0], [-1e-3, 5e-4, 1e-3], [0.0, 1e-3, -2e-3]];
    let fixed = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: g,
    });
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces: vec![fixed; model.faces.len()],
    };
    let parent = solve_elastic(&model, &problem, &options(1.0)).unwrap_or_else(|e| panic!("{e}"));
    let sub = solve_submodel(
        &model,
        &problem,
        &parent,
        SolidId(0),
        [-1.0, -1.0, -1.0],
        [6.0, 4.0, 3.0],
        &options(0.25),
        Caps::Displacement,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    // The x and y cuts are capped; z is beyond the plate.
    assert_eq!(sub.caps.len(), 2, "{:?}", sub.caps);
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for (v, d) in sub.model.solids[0]
        .vertices
        .iter()
        .zip(&sub.result.vertex_displacement[0])
    {
        for i in 0..3 {
            let exact: f64 = (0..3).map(|j| g[i][j] * v[j]).sum();
            worst = worst.max((d[i] - exact).abs());
            scale = scale.max(exact.abs());
        }
    }
    assert!(worst / scale <= 1e-8, "{}", worst / scale);
}

/// The coarse plate (2 cells per hole radius), its submodel around the
/// hole at cell size `h` reaching `reach` along x, and both peaks over
/// the converged one.
fn refined(reach: f64, h: f64) -> (f64, f64) {
    let (model, problem) = c6();
    let parent = solve_elastic(&model, &problem, &options(1.0)).unwrap_or_else(|e| panic!("{e}"));
    let coarse = parent.stresses[0]
        .surface
        .iter()
        .filter(|s| s.face.0 == 4)
        .map(|s| s.stress.von_mises)
        .fold(0.0, f64::max);
    let sub = solve_submodel(
        &model,
        &problem,
        &parent,
        SolidId(0),
        [-1.0, -1.0, -1.0],
        [reach, 7.0, 3.0],
        &options(h),
        Caps::Displacement,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let (peak, _) = sub
        .peak_on(FaceId(4), 0.5)
        .unwrap_or_else(|| panic!("no hole samples"));
    (coarse / CONVERGED, peak / CONVERGED)
}

/// Measured 2026-10-11 against C6's converged peak: the coarse plate
/// (2 cells per hole radius) 0.541; its submodel at 12 cells per radius
/// 0.981 with the box to x = 6 (3 radii) and 0.985 to x = 10; the whole
/// plate at the same cell size 0.992. The coarse field on the caps costs
/// 1.2% and 0.7%.
#[test]
fn a_submodel_recovers_the_hole_peak_and_does_not_depend_on_its_box() {
    let h = 1.0 / 6.0;
    let (coarse, near) = refined(6.0, h);
    let (_, far) = refined(10.0, h);
    let (model, problem) = c6();
    let full = solve_elastic(&model, &problem, &options(h)).unwrap_or_else(|e| panic!("{e}"));
    let full_peak = full.stresses[0]
        .surface
        .iter()
        .filter(|s| s.face.0 == 4)
        .map(|s| s.stress.von_mises)
        .fold(0.0, f64::max)
        / CONVERGED;
    assert!(coarse < 0.6, "coarse {coarse}");
    assert!((near - 1.0).abs() <= 0.03, "{near}");
    assert!(
        (near - full_peak).abs() <= 0.02,
        "{near} vs whole plate {full_peak}"
    );
    assert!((near - far).abs() <= 0.01, "{near} vs {far}");
}

/// The box reaches 4 radii, widened to two coarse cells (PH-104).
#[test]
fn a_submodel_box_reaches_four_radii_or_two_coarse_cells() {
    assert_eq!(box_half(4.0, 0.4, 0.1), 1.6);
    assert_eq!(box_half(4.0, 0.4, 2.6), 5.2);
}
