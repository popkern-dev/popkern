//! KP20 c-net to bi-3 tensor-border conversion.
//!
//! KP20 can consume an already extracted bi-3 tensor-border, but the
//! paper-strict path starts from the surrounding regular cubic c-net.
//! This module implements the uniform cubic B-spline to Bezier
//! conversion used to obtain that tensor-border, side-by-side with the
//! existing KP15 bridge.

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp_polar_shared::weighted_sum;
use super::kp20_bi6_math::kp20_limit_weights;
use crate::geometry::Point3;

/// Uniform cubic B-spline to Bezier conversion matrix for one knot span.
pub const CUBIC_BSPLINE_TO_BEZIER: [[f64; 4]; 4] = [
    [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0, 0.0],
    [0.0, 4.0 / 6.0, 2.0 / 6.0, 0.0],
    [0.0, 2.0 / 6.0, 4.0 / 6.0, 0.0],
    [0.0, 1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0],
];

/// Regular 4x4 c-net stencil per KP20 sector.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20RegularCNet {
    /// One regular 4x4 cubic B-spline stencil per sector.
    pub sectors: Vec<[[Point3; 4]; 4]>,
}

/// KP18/KP20 extraordinary cap-net with `6n + 1` unique points.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20CapNet {
    /// Extraordinary input node `c_7`.
    pub center: Point3,
    /// Per-sector nodes `c^k_1..c^k_6`.
    pub sectors: Vec<[Point3; 6]>,
}

/// Compute KP18/KP20's Catmull-Clark limit point from the `6n + 1`
/// cap-net.
pub fn kp20_star_point_from_cap_net(cap_net: &Kp20CapNet) -> Result<Point3, String> {
    let n = cap_net.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP20 cap-net star point needs at least 3 sectors, got {n}"
        ));
    }
    let (gamma5, gamma6) = kp20_limit_weights(n)?;
    let center_weight = 1.0 - n as f64 * gamma5 - n as f64 * gamma6;
    let mut pairs = Vec::with_capacity(1 + 2 * n);
    pairs.push((center_weight, cap_net.center));
    for sector in &cap_net.sectors {
        pairs.push((gamma5, sector[4]));
        pairs.push((gamma6, sector[5]));
    }
    Ok(weighted_sum(&pairs))
}

/// Convert a regular cubic c-net ring into KP20's depth-2 bi-3
/// tensor-border.
pub fn kp20_tensor_border_from_regular_cnet(
    cnet: &Kp20RegularCNet,
) -> Result<TensorBorderBi3, String> {
    if cnet.sectors.len() < 3 {
        return Err(format!(
            "KP20 c-net tensor-border conversion needs at least 3 sectors, got {}",
            cnet.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(cnet.sectors.len());
    for (sector_index, grid) in cnet.sectors.iter().enumerate() {
        validate_grid(sector_index, grid)?;
        let bezier = cubic_bspline_grid_to_bezier(grid);
        sectors.push([bezier[0], bezier[1], bezier[2]]);
    }
    Ok(TensorBorderBi3 { sectors })
}

/// Apply the cubic B-spline to Bezier matrix in both tensor directions.
pub fn cubic_bspline_grid_to_bezier(grid: &[[Point3; 4]; 4]) -> [[Point3; 4]; 4] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = tensor_weighted_point(i, j, grid);
        }
    }
    out
}

fn tensor_weighted_point(i: usize, j: usize, grid: &[[Point3; 4]; 4]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for row in 0..4 {
        for col in 0..4 {
            let w = CUBIC_BSPLINE_TO_BEZIER[i][row] * CUBIC_BSPLINE_TO_BEZIER[j][col];
            out.x += w * grid[row][col].x;
            out.y += w * grid[row][col].y;
            out.z += w * grid[row][col].z;
        }
    }
    out
}

fn validate_grid(sector_index: usize, grid: &[[Point3; 4]; 4]) -> Result<(), String> {
    for (row, cps) in grid.iter().enumerate() {
        for (col, p) in cps.iter().enumerate() {
            if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
                return Err(format!(
                    "KP20 c-net sector {sector_index} grid[{row}][{col}] is non-finite"
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "kp20_bi6_cnet_tests.rs"]
mod tests;
