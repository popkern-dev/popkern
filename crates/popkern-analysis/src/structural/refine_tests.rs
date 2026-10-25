//! PH-036: stress raisers re-solved in submodels (`--refine`) against
//! the closed forms, from the coarse resolutions linear SPR misses at:
//! C5 (Lamé bore) at 4 cells through the wall, C6 (Roark plate with a
//! hole) and C7 (Kirsch) at 3 cells per radius, over four grid shifts.

use super::assemble::FaceMechanics;
use super::recovery_tests::{kirsch_problem, plate, plate_faces};
use super::submodel::refine_raiser;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{model_of, quarter_tube};
use crate::model::FaceId;

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);

/// Grid shifts as fractions of the cell size.
const SHIFTS: [[f64; 3]; 4] = [
    [0.0; 3],
    [0.31, 0.17, 0.0],
    [0.53, 0.71, 0.0],
    [0.11, 0.43, 0.0],
];

/// (coarse, refined) peak principal stress of `face` (radius `radius`)
/// over `exact`, at cell size `h` and each shift.
fn refined(
    model: &crate::model::AnalysisModel,
    problem: &ElasticProblem,
    (face, radius): (u32, f64),
    h: f64,
    exact: f64,
) -> Vec<(f64, f64)> {
    SHIFTS
        .iter()
        .map(|s| {
            let options = AnalysisOptions {
                target_cell_size: Some(h),
                grid_shift: s.map(|x| x * h),
                cg_rel_tol: 1e-10,
                cg_max_iters: 400,
                ..AnalysisOptions::default()
            };
            let r = solve_elastic(model, problem, &options).unwrap_or_else(|e| panic!("{e}"));
            let coarse = r.stresses[0]
                .surface
                .iter()
                .filter(|x| x.face == FaceId(face))
                .map(|x| x.stress.max_principal)
                .fold(0.0, f64::max);
            let sub = refine_raiser(model, problem, &r, FaceId(face), radius, &options)
                .unwrap_or_else(|e| panic!("{e}"));
            (coarse / exact, sub.peak_principal / exact)
        })
        .collect()
}

/// Measured 2026-10-11, over the Lamé hoop stress at the bore: refined
/// 1.042 at every shift (the box holds the whole quarter tube); linear
/// SPR at the vertices averages 34% off here.
#[test]
fn c5_the_bore_refines_to_lame_from_4_cells() {
    let (ri, ro, p) = (2.0, 6.0, 1e-3);
    let model = model_of(quarter_tube(ri, ro, 2.0, 48), 5);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces: vec![
            FaceMechanics::Pressure(p),
            FaceMechanics::Free,
            ROLLER,
            ROLLER,
            ROLLER,
        ],
    };
    let hoop = p * ri * ri / (ro * ro - ri * ri) * (1.0 + ro * ro / (ri * ri));
    let r = refined(&model, &problem, (0, ri), 1.0, hoop);
    assert!(r.iter().all(|x| (x.1 - 1.0).abs() <= 0.05), "{r:?}");
}

/// Measured 2026-10-11, over Roark's K_t σ_nom (plane stress): coarse
/// face peak 0.812–0.950, refined 1.004–1.005.
#[test]
fn c6_the_hole_refines_to_roark_from_3_cells() {
    let (model, segments) = plate();
    let sigma = 1e-3;
    let mut faces = plate_faces(model.faces.len());
    for (f, _, n) in &segments {
        if n[0] > 0.5 {
            faces[*f as usize] = FaceMechanics::Traction([sigma, 0.0, 0.0]);
        }
    }
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let q: f64 = 1.0 / 3.0;
    let kt = 3.0 - 3.13 * q + 3.66 * q * q - 1.53 * q * q * q;
    let r = refined(
        &model,
        &problem,
        (4, 2.0),
        2.0 / 3.0,
        kt * sigma * 12.0 / 8.0,
    );
    assert!(r.iter().all(|x| (x.1 - 1.0).abs() <= 0.10), "{r:?}");
    // Negative control: the coarse peak misses by more than 10% at some
    // shift.
    assert!(r.iter().any(|x| x.0 < 0.9), "{r:?}");
}

/// Measured 2026-10-11, over Kirsch's 3t at the hole: coarse face peak
/// 0.822–0.960, refined 1.003–1.004.
#[test]
fn c7_the_hole_refines_to_kirsch_from_3_cells() {
    let (model, segments) = plate();
    let t = 1e-3;
    let problem = kirsch_problem(&model, &segments, t);
    let r = refined(&model, &problem, (4, 2.0), 2.0 / 3.0, 3.0 * t);
    assert!(r.iter().all(|x| (x.1 - 1.0).abs() <= 0.10), "{r:?}");
    assert!(r.iter().any(|x| x.0 < 0.9), "{r:?}");
}
