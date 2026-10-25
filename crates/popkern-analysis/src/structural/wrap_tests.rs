//! PH-062: a sliver wrapping around a corner of the material across a
//! slot (an L-shaped strip thinner than the lattice spacing).

use super::assemble::FaceMechanics;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{extrude_faces, model_of};
use crate::moments::Splitting;

const G: f64 = 0.075;
const T: f64 = 1.0;

/// A 6 × 4 × 1 plate with an L-shaped slot of width `G`: up from the
/// y = 0 edge at x in [x0, x0 + G] to y0 + G, then left to x = 1. Block
/// A (x < x0, y < y0) hangs from the root x < 1; the rest (B) wraps
/// round A's corner. Faces: 0/1 z sides, 2 A's bottom (y = 0), 3 the
/// slot, 4 B's bottom, 5 other sides, 6 the top (y = 4).
fn hooked(x0: f64, y0: f64, g: f64) -> crate::model::AnalysisModel {
    let outline = [
        [0.0, 0.0],
        [x0, 0.0],
        [x0, y0],
        [1.0, y0],
        [1.0, y0 + g],
        [x0 + g, y0 + g],
        [x0 + g, 0.0],
        [6.0, 0.0],
        [6.0, 4.0],
        [0.0, 4.0],
    ];
    let tris = [
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 9],
        [4, 5, 9],
        [5, 8, 9],
        [5, 6, 7],
        [5, 7, 8],
    ];
    let mut mesh = extrude_faces(&outline, &tris, 0.0, T, &[2, 3, 3, 3, 3, 3, 4, 5, 6, 5]);
    let slot = g * y0 + g * (x0 - 1.0) + g * g;
    mesh.exact_volume = Some(T * (24.0 - slot));
    model_of(mesh, 7)
}

/// A's corner for B's slot faces 0.005 inside the cell [3.75, 4] ×
/// [2.25, 2.5] at h = 0.25, unshifted, with a slot of width `g`.
fn corner(g: f64) -> (f64, f64) {
    (4.0 - 0.005 - g, 2.5 - 0.005 - g)
}

/// A's mean downward displacement at its bottom when pulled down with
/// the top clamped, at cell size `h`.
fn sag(g: f64, h: f64, split: Splitting) -> (f64, StructuralResult) {
    let (x0, y0) = corner(g);
    let model = hooked(x0, y0, g);
    let mut faces = vec![FaceMechanics::Free; 7];
    faces[6] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[2] = FaceMechanics::Traction([0.0, -1e-3, 0.0]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        max_active_cells: 400_000,
        coarsest_dofs: 600,
        cg_max_iters: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic_with(
        &model,
        &problem,
        &options,
        crate::recover::Fallback::InCell,
        None,
        split,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let points = [0.5, 1.5, 2.5, 3.5];
    let sum: f64 = points
        .iter()
        .map(|&x| {
            r.field
                .displacement_at([x, 0.05, 0.5 * T])
                .unwrap_or_else(|e| panic!("{e}"))
                .unwrap_or_else(|| panic!("outside"))[1]
        })
        .sum();
    (-sum / points.len() as f64, r)
}

/// Measured 2026-10-11, A's sag over each plate's own reference at
/// h = 0.1: the 0.075 slot at h = 0.25 reads 0.915 (0.063 before, when
/// the L-shaped strip of B in the corner cell, rejected as a sliver
/// because its bounding box held A's corner, welded A to B); the
/// resolved 0.5 slot 0.937; welded 0.03.
#[test]
fn a_strip_wrapping_a_corner_across_a_slot_does_not_weld() {
    let reference = sag(G, 0.1, Splitting::Exact).0;
    let (slit, r) = sag(G, 0.25, Splitting::Exact);
    let resolved = sag(0.5, 0.25, Splitting::Exact).0 / sag(0.5, 0.1, Splitting::Exact).0;
    let slit = slit / reference;
    assert!(r.report.split_cells > 0);
    assert!((slit - resolved).abs() <= 0.05, "{slit} vs {resolved}");
    // Negative control: welded, A hangs from B's corner too.
    let welded = sag(G, 0.25, Splitting::Welded).0 / reference;
    assert!(welded < 0.1, "welded {welded}");
}

/// The corner cell holds A's corner and B's L-shaped strip as two
/// fragments, and a linear temperature is exact through them.
#[test]
fn the_wrapping_strip_is_its_own_fragment_and_exact() {
    use crate::clip::clip_solid;
    use crate::grid::Grid;
    use crate::moments::rasterise_with;
    use crate::thermal::{ConductionProblem, FaceCondition, solve_conduction};
    let (x0, y0) = corner(G);
    let model = hooked(x0, y0, G);
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let grid = Grid::for_model([0.0; 3], [6.0, 4.0, T], &options).unwrap_or_else(|e| panic!("{e}"));
    let mesh = &model.solids[0];
    let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
    let m = rasterise_with(&grid, mesh.id, &clipped, 2, &|_| false, Splitting::Exact)
        .unwrap_or_else(|e| panic!("{e}"));
    let cell = grid
        .cell_of_point([3.8, 2.3, 0.1])
        .and_then(|ijk| grid.cell_index(ijk))
        .unwrap_or_else(|| panic!("corner cell"));
    let sc = m
        .fragments
        .split
        .iter()
        .find(|s| s.cell == cell)
        .unwrap_or_else(|| panic!("corner cell not split"));
    assert_eq!(m.entries(cell).len(), 2);
    assert_eq!(sc.slivers.len(), 1);
    // A's corner block (0.17 square) and the strip (0.005 thick).
    let volumes: Vec<f64> = m.entries(cell).map(|e| m.cells[e].m[0]).collect();
    let strip = 0.005 * (0.25 + 0.245) * 0.25;
    assert!(
        (volumes[0] - 0.17 * 0.17 * 0.25).abs() <= 1e-12,
        "{volumes:?}"
    );
    assert!((volumes[1] - strip).abs() <= 1e-12, "{volumes:?}");
    let g = [0.7, -0.3, 0.45];
    let problem = ConductionProblem {
        conductivity: 20.0,
        faces: vec![
            FaceCondition::Fixed {
                t0_c: 1.0,
                gradient: g,
            };
            7
        ],
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let r = solve_conduction(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
    let worst = mesh
        .vertices
        .iter()
        .zip(&r.vertex_temperature[0])
        .map(|(p, t)| (t - 1.0 - (0..3).map(|a| g[a] * p[a]).sum::<f64>()).abs())
        .fold(0.0, f64::max);
    assert!(worst <= 1e-9, "{worst}");
}
