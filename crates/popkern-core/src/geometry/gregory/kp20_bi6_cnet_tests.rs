use super::*;
use crate::geometry::gregory::catmull_clark_star_point;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn cap_net(n: usize) -> Kp20CapNet {
    let mut sectors = Vec::new();
    for k in 0..n {
        let mut sector = [p(0.0, 0.0, 0.0); 6];
        for (i, cp) in sector.iter_mut().enumerate() {
            *cp = p(k as f64, i as f64 + 1.0, 0.1 * (k + i) as f64);
        }
        sectors.push(sector);
    }
    Kp20CapNet {
        center: p(0.25, -0.5, 0.75),
        sectors,
    }
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
fn cap_net_star_point_matches_shared_catmull_clark_limit() {
    let cap_net = cap_net(5);
    let c5s: Vec<Point3> = cap_net.sectors.iter().map(|sector| sector[4]).collect();
    let c6s: Vec<Point3> = cap_net.sectors.iter().map(|sector| sector[5]).collect();
    let expected = catmull_clark_star_point(5, cap_net.center, &c5s, &c6s);
    let actual = kp20_star_point_from_cap_net(&cap_net).unwrap();

    assert_close(actual, expected);
}

#[test]
fn cap_net_star_point_uses_n3_corrected_limit_weights() {
    let cap_net = cap_net(3);
    let actual = kp20_star_point_from_cap_net(&cap_net).unwrap();
    let (gamma5, gamma6) = kp20_limit_weights(3).unwrap();
    let center_weight = 1.0 - 3.0 * gamma5 - 3.0 * gamma6;
    let mut expected = p(
        center_weight * cap_net.center.x,
        center_weight * cap_net.center.y,
        center_weight * cap_net.center.z,
    );
    for sector in &cap_net.sectors {
        expected.x += gamma5 * sector[4].x + gamma6 * sector[5].x;
        expected.y += gamma5 * sector[4].y + gamma6 * sector[5].y;
        expected.z += gamma5 * sector[4].z + gamma6 * sector[5].z;
    }

    assert_close(actual, expected);
}

#[test]
fn cap_net_star_point_rejects_invalid_or_unpublished_valence() {
    assert!(kp20_star_point_from_cap_net(&cap_net(2)).is_err());
    assert!(kp20_star_point_from_cap_net(&cap_net(4)).is_err());
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
