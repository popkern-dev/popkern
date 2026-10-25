use super::*;

/// A grid with exactly 5×3×7 cells (3×1×5 inside the padding).
fn grid_5_3_7() -> Grid {
    match Grid::covering([0.0, 0.0, 0.0], [3.0, 1.0, 5.0], 1.0) {
        Ok(g) => g,
        Err(e) => panic!("{e}"),
    }
}

#[test]
fn cell_and_node_indices_are_bijections() {
    let g = grid_5_3_7();
    assert_eq!(g.dims(), [5, 3, 7]);
    for index in 0..g.cell_count() {
        let ijk = g.cell_ijk(index).unwrap_or_else(|| panic!("cell {index}"));
        assert_eq!(g.cell_index(ijk), Some(index));
    }
    for (k, j, i) in (0..7).flat_map(|k| (0..3).flat_map(move |j| (0..5).map(move |i| (k, j, i)))) {
        let index = g
            .cell_index([i, j, k])
            .unwrap_or_else(|| panic!("{i},{j},{k}"));
        assert_eq!(g.cell_ijk(index), Some([i, j, k]));
    }
    assert_eq!(g.node_dims(), [6, 4, 8]);
    for index in 0..g.node_count() {
        let ijk = g.node_ijk(index).unwrap_or_else(|| panic!("node {index}"));
        assert_eq!(g.node_index(ijk), Some(index));
    }
    assert_eq!(g.cell_index([5, 0, 0]), None);
    assert_eq!(g.cell_ijk(g.cell_count()), None);
    assert_eq!(g.node_index([0, 4, 0]), None);
    assert_eq!(g.node_ijk(g.node_count()), None);
}

#[test]
fn iteration_is_lexicographic_with_x_fastest() {
    let g = grid_5_3_7();
    let cells: Vec<[usize; 3]> = g.cells().collect();
    assert_eq!(cells.len(), 105);
    assert_eq!(
        &cells[..6],
        &[
            [0, 0, 0],
            [1, 0, 0],
            [2, 0, 0],
            [3, 0, 0],
            [4, 0, 0],
            [0, 1, 0]
        ]
    );
    assert_eq!(cells[15], [0, 0, 1]);
    for (index, ijk) in cells.iter().enumerate() {
        assert_eq!(g.cell_index(*ijk), Some(index));
    }
}

#[test]
fn cell_nodes_follow_the_corner_order() {
    let g = grid_5_3_7();
    let ijk = [2, 1, 3];
    let nodes = g.cell_nodes(ijk).unwrap_or_else(|| panic!("in grid"));
    for (corner, node) in CELL_CORNERS.iter().zip(nodes) {
        let expected = [0, 1, 2].map(|a| ijk[a] + corner[a]);
        assert_eq!(g.node_ijk(node), Some(expected));
    }
    assert_eq!(g.cell_nodes([4, 2, 7]), None);
}

/// Every bbox lies inside the grid with a full cell of margin.
fn assert_padded(min: [f64; 3], max: [f64; 3], h: f64) {
    let g = Grid::covering(min, max, h).unwrap_or_else(|e| panic!("{e}"));
    let lo = g.origin();
    for a in 0..3 {
        let hi = lo[a] + h * g.dims()[a] as f64;
        // Coordinates round at their own magnitude, so a margin may fall
        // short of h by a few ulps of the coordinate.
        let slack = 8.0 * f64::EPSILON * min[a].abs().max(max[a].abs()).max(h);
        assert!(
            min[a] - lo[a] >= h - slack,
            "low margin on {a}: {min:?} {lo:?} h {h}"
        );
        assert!(
            hi - max[a] >= h - slack,
            "high margin on {a}: {max:?} {hi} h {h}"
        );
        // and no more than two cells of slack, so the padding is minimal
        assert!(hi - max[a] < 2.0 * h + 1e-9, "excess on {a}");
    }
}

#[test]
fn the_box_has_a_full_cell_of_margin_on_every_side() {
    let offsets = [
        [0.0, 0.0, 0.0],
        [-13.37, 2.5e3, 0.1],
        [1e6, -1e6, 7.0 / 3.0],
    ];
    let shapes = [[1.0, 1.0, 1.0], [40.0, 2.0, 0.6], [0.3, 17.0, 120.0]];
    for offset in offsets {
        for shape in shapes {
            let max = [0, 1, 2].map(|a| offset[a] + shape[a]);
            for h in [0.1, 0.3, 1.0 / 3.0, 0.7] {
                assert_padded(offset, max, h);
            }
        }
    }
}

#[test]
fn a_derived_cell_size_keeps_within_the_budget() {
    let shapes = [
        [1.0, 1.0, 1.0],
        [40.0, 2.0, 0.6],
        [0.3, 17.0, 120.0],
        [100.0, 1.0, 1.0],
    ];
    for shape in shapes {
        for budget in [1, 7, 1000, 30_000] {
            let options = AnalysisOptions {
                max_active_cells: budget,
                ..AnalysisOptions::default()
            };
            let g = Grid::for_model([0.0; 3], shape, &options).unwrap_or_else(|e| panic!("{e}"));
            let [a, b, c] = g.dims().map(|d| d - 2);
            assert!(
                a * b * c <= budget,
                "{shape:?} budget {budget}: {a}x{b}x{c}"
            );
            // Minimal: either the cube-root start already fit, or one
            // search step smaller would exceed the budget.
            let start = (shape[0] * shape[1] * shape[2] / budget as f64).pcbrt();
            let [a, b, c] = cover(shape, g.cell_size() / DERIVE_GROWTH);
            assert!(
                g.cell_size() == start || a * b * c > budget,
                "{shape:?} budget {budget}: h {} is not minimal",
                g.cell_size()
            );
        }
    }
    let fixed = AnalysisOptions {
        target_cell_size: Some(0.25),
        ..AnalysisOptions::default()
    };
    let g = Grid::for_model([0.0; 3], [1.0; 3], &fixed).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(g.cell_size(), 0.25);
    assert_eq!(g.dims(), [6, 6, 6]);
}

#[test]
fn points_map_to_their_cells() {
    let g = grid_5_3_7();
    assert_eq!(g.cell_of_point([-0.5, -0.5, -0.5]), Some([0, 0, 0]));
    assert_eq!(g.cell_of_point([0.0, 0.0, 0.0]), Some([1, 1, 1]));
    assert_eq!(g.cell_of_point([4.0, 2.0, 6.0]), Some([4, 2, 6]));
    assert_eq!(g.cell_of_point([4.01, 0.0, 0.0]), None);
    assert_eq!(g.cell_of_point([f64::NAN, 0.0, 0.0]), None);
    assert_eq!(g.node_position([1, 1, 1]), [0.0, 0.0, 0.0]);
}

#[test]
fn degenerate_inputs_are_typed_errors() {
    assert_eq!(
        Grid::covering([0.0; 3], [1.0, 0.0, 1.0], 0.1),
        Err(GridError::ZeroExtent { axis: 1 })
    );
    assert_eq!(
        Grid::covering([0.0; 3], [1.0, 1.0, f64::INFINITY], 0.1),
        Err(GridError::NonFiniteBounds)
    );
    assert_eq!(
        Grid::covering([0.0; 3], [1.0; 3], 0.0),
        Err(GridError::BadCellSize(0.0))
    );
    assert!(matches!(
        Grid::covering([0.0; 3], [1.0; 3], f64::NAN),
        Err(GridError::BadCellSize(h)) if h.is_nan()
    ));
    assert!(matches!(
        Grid::covering([0.0; 3], [1.0; 3], 1e-9),
        Err(GridError::TooLarge { .. })
    ));
    assert!(matches!(
        Grid::covering([0.0; 3], [1e300; 3], 1e-300),
        Err(GridError::TooLarge { .. })
    ));
    let none = AnalysisOptions {
        max_active_cells: 0,
        ..AnalysisOptions::default()
    };
    assert_eq!(
        Grid::for_model([0.0; 3], [1.0; 3], &none),
        Err(GridError::NoBudget)
    );
}
