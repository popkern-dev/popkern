//! PH-061: a sliver whose bottom lies on a cell wall, tessellated so
//! that some of its bottom triangles touch the sliver's side pieces only
//! at a vertex. Before, those pieces stayed with the material across
//! the slot while the sliver's moments counted the volume above them,
//! and a linear field read off by 7e-8 (planetary_gearbox's tooth tips).

use crate::mesh_test_helpers::{extrude_faces, model_of};
use crate::model::AnalysisModel;
use crate::options::AnalysisOptions;
use crate::thermal::{ConductionProblem, FaceCondition, solve_conduction};

/// The slotted plate (a 6 × 2 × 1 plate, slot of width 0.075 centred at
/// y = 1 from x = 1) with the upper prong's slot face split at x = 3.5
/// and its bottom fanned so the triangles touching that face only at
/// (3.5, b) come first.
fn plate() -> AnalysisModel {
    let w = 0.075;
    let (a, b) = (1.0 - 0.5 * w, 1.0 + 0.5 * w);
    let outline = [
        [0.0, 0.0],
        [6.0, 0.0],
        [6.0, a],
        [1.0, a],
        [1.0, b],
        [3.5, b],
        [6.0, b],
        [6.0, 2.0],
        [3.5, 2.0],
        [0.0, 2.0],
    ];
    let tris = [
        [5, 8, 9],
        [5, 7, 8],
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 9],
        [4, 5, 9],
        [5, 6, 7],
    ];
    let mut mesh = extrude_faces(&outline, &tris, 0.0, 1.0, &[2, 3, 2, 2, 2, 2, 4, 2, 2, 5]);
    mesh.exact_volume = Some(12.0 - 5.0 * w);
    model_of(mesh, 6)
}

/// The slot face 0.005 under a cell wall, the bottom on one.
fn options() -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(0.25),
        grid_shift: [0.0, 0.2075, 0.0],
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    }
}

/// Measured 2026-10-11: worst vertex error 5.1e-7 before, exact (1e-9
/// gate) now.
#[test]
fn a_linear_temperature_is_exact_beside_the_sliver() {
    let model = plate();
    let g = [0.7, -0.3, 0.45];
    let problem = ConductionProblem {
        conductivity: 20.0,
        faces: vec![
            FaceCondition::Fixed {
                t0_c: 1.0,
                gradient: g,
            };
            6
        ],
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let r = solve_conduction(&model, &problem, &options()).unwrap_or_else(|e| panic!("{e}"));
    assert!(r.report.split_cells > 0);
    let exact = |p: &[f64; 3]| 1.0 + (0..3).map(|a| g[a] * p[a]).sum::<f64>();
    let worst = model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_temperature[0])
        .map(|(p, t)| (t - exact(p)).abs())
        .fold(0.0, f64::max);
    assert!(worst <= 1e-9, "{worst}");
}
