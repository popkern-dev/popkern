use super::*;
use crate::clip::clip_solid;
use crate::mesh_test_helpers::{bounds, tube, u_prism};
use crate::model::SolidMesh;
use crate::options::AnalysisOptions;

/// Whether the surface point `p` of `mesh` counts as flat at h = 0.25.
fn flat(mesh: &SolidMesh, p: [f64; 3]) -> bool {
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        ..AnalysisOptions::default()
    };
    let (lo, hi) = bounds(mesh);
    let grid = Grid::for_model(lo, hi, &options).unwrap_or_else(|e| panic!("{e}"));
    let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
    let cell = grid
        .cell_of_point(p)
        .and_then(|ijk| grid.cell_index(ijk))
        .unwrap_or_else(|| panic!("{p:?} outside the grid"));
    FlatCells::new(&grid, &clipped).is_flat(cell, p)
}

#[test]
fn flat_and_convex_surfaces_extrapolate_and_concave_ones_do_not() {
    // The U prism: x 0..3, y 0..2, z 0..1, a slot x 1..2 down to y = 0.5.
    let u = u_prism(1.0);
    // A face, with the faces two cells away behind it.
    assert!(flat(&u, [0.2, 0.0, 0.4]));
    // Two cells under the slot, its walls' planes are within reach and
    // in front: planes count unbounded, which only errs toward
    // interpolation.
    assert!(!flat(&u, [1.4, 0.0, 0.4]));
    // A convex edge.
    assert!(flat(&u, [2.9, 0.0, 0.4]));
    // Beside a re-entrant corner (the slot's floor next to its wall).
    assert!(!flat(&u, [1.1, 0.5, 0.4]));
    // A bore: each facet lies in front of its neighbours' planes.
    let t = tube(2.0, 6.0, 4.0, 16);
    // The middle of the first facet (between angles 0 and π/8) of a
    // circle of radius `r`.
    let on_facet = |r: f64| {
        let a = std::f64::consts::PI / 16.0;
        let rho = r * a.cos();
        [rho * a.cos(), rho * a.sin(), 1.9]
    };
    assert!(!flat(&t, on_facet(2.0)));
    // The same tube's outside is convex.
    assert!(flat(&t, on_facet(6.0)));
}
