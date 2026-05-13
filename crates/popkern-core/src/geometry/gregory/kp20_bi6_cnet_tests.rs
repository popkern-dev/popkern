use super::*;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.approx_eq(&b, 1e-12),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

fn affine_grid(offset: f64) -> [[Point3; 4]; 4] {
    let mut grid = [[p(0.0, 0.0, 0.0); 4]; 4];
    for (row, cps) in grid.iter_mut().enumerate() {
        for (col, cp) in cps.iter_mut().enumerate() {
            *cp = p(offset + row as f64, col as f64, row as f64 - col as f64);
        }
    }
    grid
}

#[test]
fn cubic_bspline_to_bezier_preserves_constant_grid() {
    let grid = [[p(1.0, -2.0, 3.0); 4]; 4];
    let bezier = cubic_bspline_grid_to_bezier(&grid);

    for row in bezier {
        for cp in row {
            assert_close(cp, p(1.0, -2.0, 3.0));
        }
    }
}

#[test]
fn cubic_bspline_to_bezier_matches_tensor_product_matrix_entry() {
    let grid = affine_grid(0.25);
    let bezier = cubic_bspline_grid_to_bezier(&grid);
    let i = 2;
    let j = 3;
    let mut expected = p(0.0, 0.0, 0.0);
    for row in 0..4 {
        for col in 0..4 {
            let w = CUBIC_BSPLINE_TO_BEZIER[i][row] * CUBIC_BSPLINE_TO_BEZIER[j][col];
            expected.x += w * grid[row][col].x;
            expected.y += w * grid[row][col].y;
            expected.z += w * grid[row][col].z;
        }
    }

    assert_close(bezier[i][j], expected);
}

#[test]
fn regular_cnet_extracts_depth_two_tensor_border() {
    let cnet = Kp20RegularCNet {
        sectors: vec![affine_grid(0.0), affine_grid(10.0), affine_grid(20.0)],
    };

    let tensor = kp20_tensor_border_from_regular_cnet(&cnet).unwrap();

    assert_eq!(tensor.sectors.len(), 3);
    let bezier0 = cubic_bspline_grid_to_bezier(&cnet.sectors[0]);
    for row in 0..3 {
        for col in 0..4 {
            assert_close(tensor.sectors[0][row][col], bezier0[row][col]);
        }
    }
}

#[test]
fn regular_cnet_rejects_short_or_non_finite_input() {
    let short = Kp20RegularCNet {
        sectors: vec![affine_grid(0.0), affine_grid(1.0)],
    };
    assert!(kp20_tensor_border_from_regular_cnet(&short).is_err());

    let mut bad_grid = affine_grid(0.0);
    bad_grid[1][2].z = f64::NAN;
    let bad = Kp20RegularCNet {
        sectors: vec![bad_grid, affine_grid(1.0), affine_grid(2.0)],
    };
    assert!(kp20_tensor_border_from_regular_cnet(&bad).is_err());
}
