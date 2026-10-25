//! PH-058: the bore under-read table's second source, C6 (a plate with
//! a central hole in tension, Roark 17.1 case 7a).

use super::assemble::FaceMechanics;
use super::recovery_tests::{plate, plate_faces};
use super::*;

/// Grid shifts as fractions of h.
const SHIFTS: [[f64; 3]; 4] = [
    [0.0; 3],
    [0.21, 0.37, 0.0],
    [0.43, 0.11, 0.29],
    [0.67, 0.53, 0.13],
];

/// The C6 plate in tension at cell size `h` and grid shift `shift`
/// (fractions of h): the model (face 4 the hole) and the result.
fn solve_c6(h: f64, shift: [f64; 3]) -> (crate::model::AnalysisModel, StructuralResult) {
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
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        grid_shift: shift.map(|s| s * h),
        cg_rel_tol: 1e-10,
        cg_max_iters: 400,
        max_active_cells: 2_000_000,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
    (model, r)
}

/// C6's converged hole peak: 13 to 20 cells per radius give 3.4096e-3
/// to 3.4233e-3 over two grid shifts each (0.985-0.989 of Roark's
/// plane-stress K_t σ_nom); this is the h = 0.1 value.
const CONVERGED: f64 = 3.42334579850262e-3;

/// Measured (2026-10-11): the hole peak over CONVERGED at 4 shifts:
/// 1 cell 0.448-0.763, 2 cells 0.540-0.639, 3 cells 0.742-0.864, 4
/// cells 0.825-0.923, 5 cells 0.879-0.956, 6 cells 0.917-0.977, 8 cells
/// 0.956-0.989. The bore check's expected under-read covers every shift
/// at 2 and 3 cells, rows C6 sets.
#[test]
fn the_bore_table_covers_c6_at_every_shift() {
    use crate::assess::fillet::FilletStatus;
    use crate::assess::{AssessOptions, bore_checks};
    use crate::model::Surface;
    for cells in [2.0, 3.0] {
        let h = 2.0 / cells;
        for shift in SHIFTS {
            let (mut model, r) = solve_c6(h, shift);
            model.faces[4].surface = Surface::Cylinder {
                origin: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                radius: 2.0,
            };
            let checks = bore_checks(&model, &r, &AssessOptions::default())
                .unwrap_or_else(|e| panic!("{e}"));
            let [check] = checks.as_slice() else {
                panic!("{checks:?}");
            };
            let FilletStatus::UnderResolved { expected_low_by } = check.status else {
                panic!("{check:?}");
            };
            let shortfall = 1.0 - check.peak_von_mises / CONVERGED;
            assert!(
                expected_low_by >= shortfall,
                "{cells} cells {shift:?}: {expected_low_by} < {shortfall}"
            );
        }
    }
}
