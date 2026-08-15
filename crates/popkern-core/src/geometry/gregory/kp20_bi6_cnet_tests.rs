// Paper-notation test module: loop indices, helper arities, and tuple shapes
// mirror KP15/KP20/PR08 subscripts and equation forms; keep the index form.
#![allow(clippy::needless_range_loop)]

use super::super::kp20_bi6_cap::Kp20MacroPatch;
use super::super::kp20_bi6_tau_solve_assemble::{
    Kp20PlanarTauPaperSolveInput, solve_kp20_planar_tau_paper_sector_f3_diagnosed,
};
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

#[test]
fn cap_net_figure_labels_follow_kp15_fig1a_radial_columns() {
    let cap_net = cap_net(5);
    let sector = 2;
    let columns = kp20_cap_net_figure_sector_columns(&cap_net, sector).unwrap();

    assert_eq!(
        KP20_CAP_NET_FIGURE_SECTOR_COLUMNS,
        [
            [Kp20CapNetFigureLabel::C4, Kp20CapNetFigureLabel::C1],
            [Kp20CapNetFigureLabel::C5, Kp20CapNetFigureLabel::C2],
            [Kp20CapNetFigureLabel::C6, Kp20CapNetFigureLabel::C3],
        ]
    );
    assert_close(columns[0][0], cap_net.sectors[sector][3]);
    assert_close(columns[0][1], cap_net.sectors[sector][0]);
    assert_close(columns[1][0], cap_net.sectors[sector][4]);
    assert_close(columns[1][1], cap_net.sectors[sector][1]);
    assert_close(columns[2][0], cap_net.sectors[sector][5]);
    assert_close(columns[2][1], cap_net.sectors[sector][2]);
    assert_close(
        kp20_cap_net_figure_point(&cap_net, sector, Kp20CapNetFigureLabel::C7).unwrap(),
        cap_net.center,
    );
}

#[test]
fn cap_net_figure_labels_reject_out_of_range_sector() {
    let cap_net = cap_net(3);
    assert!(kp20_cap_net_figure_sector_columns(&cap_net, 3).is_err());
    assert!(kp20_cap_net_figure_point(&cap_net, 3, Kp20CapNetFigureLabel::C1).is_err());
    assert!(kp20_cap_net_figure_point(&cap_net, 3, Kp20CapNetFigureLabel::C7).is_ok());
}

#[test]
fn characteristic_cap_net_lifts_stam_eigenvector_by_fourier_mode_one() {
    let n = 5;
    let scalars = Kp20CharacteristicScalars::for_valence(n).unwrap();
    let cap_net = kp20_characteristic_cap_net(n).unwrap();

    assert_eq!(cap_net.sectors.len(), n);
    assert_close(cap_net.center, p(0.0, 0.0, 0.0));
    for k in 0..n {
        let theta = core::f64::consts::TAU * k as f64 / n as f64;
        let (s, c) = theta.sin_cos();
        for m in 0..6 {
            let value = scalars.eigenvector[m];
            assert_close(cap_net.sectors[k][m], p(value * c, value * s, 0.0));
        }
    }
}

#[test]
fn characteristic_cap_net_is_rotationally_covariant() {
    let n = 7;
    let cap_net = kp20_characteristic_cap_net(n).unwrap();
    let theta = core::f64::consts::TAU / n as f64;
    let (s, c) = theta.sin_cos();

    for k in 0..n {
        let next = (k + 1) % n;
        for m in 0..6 {
            let current = cap_net.sectors[k][m];
            let expected = p(
                c * current.x - s * current.y,
                s * current.x + c * current.y,
                0.0,
            );
            assert_close(cap_net.sectors[next][m], expected);
        }
    }
}

#[test]
fn characteristic_cap_net_rejects_invalid_valence() {
    assert!(kp20_characteristic_cap_net(2).is_err());
}

#[test]
fn figure_tensor_border_rows_follow_rotational_characteristic_cap_net() {
    let n = 5;
    let cap_net = kp20_characteristic_cap_net(n).unwrap();
    let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).unwrap();
    let theta = core::f64::consts::TAU / n as f64;
    let (s, c) = theta.sin_cos();

    assert_eq!(rows.sectors.len(), n);
    for k in 0..n {
        let next = (k + 1) % n;
        for row in 0..3 {
            for col in 0..4 {
                let current = rows.sectors[k][row][col];
                let expected = p(
                    c * current.x - s * current.y,
                    s * current.x + c * current.y,
                    0.0,
                );
                assert_close(rows.sectors[next][row][col], expected);
            }
        }
    }
}

#[test]
fn figure_tensor_border_rows_are_rejected_by_paper_tau_constraints() {
    let cap_net = kp20_characteristic_cap_net(5).unwrap();
    let rows = kp20_cap_net_figure_tensor_border_rows(&cap_net).unwrap();
    let tensor_border = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();
    let solved = solve_kp20_planar_tau_paper_sector_f3_diagnosed(&Kp20PlanarTauPaperSolveInput {
        valence: 5,
        sector: 0,
        characteristic_tensor_border: &tensor_border,
        initial: Kp20MacroPatch::zero(),
        regularization_weight: 1e-10,
    })
    .unwrap();

    assert!(solved.kkt_condition_number.is_some_and(f64::is_finite));
    assert!(
        solved.max_constraint_residual > 1.0,
        "diagnostic Fig. 9 tensor-border rows should remain rejected until the regular Catmull-Clark ring is derived; residual {}",
        solved.max_constraint_residual
    );
}

#[test]
fn figure_tensor_border_rows_reject_short_cap_net() {
    let cap_net = cap_net(2);
    assert!(kp20_cap_net_figure_tensor_border_rows(&cap_net).is_err());
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

fn tensor_border_rows_from_grid(grid: &[[Point3; 4]; 4]) -> [[Point3; 4]; 3] {
    [grid[0], grid[1], grid[2]]
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
fn cubic_bspline_rows_to_tensor_border_matches_full_grid_conversion() {
    let grid = affine_grid(0.25);
    let rows = tensor_border_rows_from_grid(&grid);
    let partial = cubic_bspline_rows_to_tensor_border_rows(&rows);
    let full = cubic_bspline_grid_to_bezier(&grid);

    for row in 0..3 {
        for col in 0..4 {
            assert_close(partial[row][col], full[row][col]);
        }
    }
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
fn regular_cnet_rows_extract_same_depth_two_tensor_border_as_full_grid() {
    let full = Kp20RegularCNet {
        sectors: vec![affine_grid(0.0), affine_grid(10.0), affine_grid(20.0)],
    };
    let rows = Kp20RegularTensorBorderCNet {
        sectors: full
            .sectors
            .iter()
            .map(tensor_border_rows_from_grid)
            .collect(),
    };

    let full_tensor = kp20_tensor_border_from_regular_cnet(&full).unwrap();
    let row_tensor = kp20_tensor_border_from_regular_cnet_rows(&rows).unwrap();

    for sector in 0..full_tensor.sectors.len() {
        for row in 0..3 {
            for col in 0..4 {
                assert_close(
                    row_tensor.sectors[sector][row][col],
                    full_tensor.sectors[sector][row][col],
                );
            }
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

#[test]
fn regular_cnet_rows_reject_short_or_non_finite_input() {
    let short = Kp20RegularTensorBorderCNet {
        sectors: vec![
            tensor_border_rows_from_grid(&affine_grid(0.0)),
            tensor_border_rows_from_grid(&affine_grid(1.0)),
        ],
    };
    assert!(kp20_tensor_border_from_regular_cnet_rows(&short).is_err());

    let mut bad_rows = tensor_border_rows_from_grid(&affine_grid(0.0));
    bad_rows[1][2].z = f64::NAN;
    let bad = Kp20RegularTensorBorderCNet {
        sectors: vec![
            bad_rows,
            tensor_border_rows_from_grid(&affine_grid(1.0)),
            tensor_border_rows_from_grid(&affine_grid(2.0)),
        ],
    };
    assert!(kp20_tensor_border_from_regular_cnet_rows(&bad).is_err());
}
