//! KP20 c-net to bi-3 tensor-border conversion.
//!
//! KP20 can consume an already extracted bi-3 tensor-border, but the
//! paper-strict path starts from the surrounding regular cubic c-net.
//! This module implements the uniform cubic B-spline to Bezier
//! conversion used to obtain that tensor-border, side-by-side with the
//! existing KP15 bridge.

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp_polar_shared::weighted_sum;
use super::kp20_bi6_math::{Kp20CharacteristicScalars, kp20_limit_weights};
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

/// Regular 4x3 B-spline control rows that determine a depth-2 tensor-border.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20RegularTensorBorderCNet {
    /// One partial regular stencil per sector. Rows are depth rows and columns
    /// run along the boundary curve.
    pub sectors: Vec<[[Point3; 4]; 3]>,
}

/// KP18/KP20 extraordinary cap-net with `6n + 1` unique points.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20CapNet {
    /// Extraordinary input node `c_7`.
    pub center: Point3,
    /// Per-sector nodes `c^k_1..c^k_6`.
    pub sectors: Vec<[Point3; 6]>,
}

/// KP15 Fig. 1a / KP18 c-net label in one sector of the `6n + 1`
/// cap-net.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20CapNetFigureLabel {
    /// `c^k_1`, the outer node of the left radial column.
    C1,
    /// `c^k_2`, the outer node of the middle radial column.
    C2,
    /// `c^k_3`, the outer node of the right radial column.
    C3,
    /// `c^k_4`, the inner node of the left radial column.
    C4,
    /// `c^k_5`, the inner node of the middle radial column.
    C5,
    /// `c^k_6`, the inner node of the right radial column.
    C6,
    /// The shared extraordinary node `c_7`.
    C7,
}

impl Kp20CapNetFigureLabel {
    fn sector_index(self) -> Option<usize> {
        match self {
            Self::C1 => Some(0),
            Self::C2 => Some(1),
            Self::C3 => Some(2),
            Self::C4 => Some(3),
            Self::C5 => Some(4),
            Self::C6 => Some(5),
            Self::C7 => None,
        }
    }
}

/// The KP15 Fig. 1a sector label layout as radial columns.
///
/// Each column is `[inner, outer]` relative to the extraordinary node.
/// This does not by itself define the surrounding regular 4x4 B-spline
/// stencil; it only records the verified figure labels used by that
/// later derivation.
pub const KP20_CAP_NET_FIGURE_SECTOR_COLUMNS: [[Kp20CapNetFigureLabel; 2]; 3] = [
    [Kp20CapNetFigureLabel::C4, Kp20CapNetFigureLabel::C1],
    [Kp20CapNetFigureLabel::C5, Kp20CapNetFigureLabel::C2],
    [Kp20CapNetFigureLabel::C6, Kp20CapNetFigureLabel::C3],
];

/// Read one KP15 Fig. 1a label from a [`Kp20CapNet`] sector.
pub fn kp20_cap_net_figure_point(
    cap_net: &Kp20CapNet,
    sector: usize,
    label: Kp20CapNetFigureLabel,
) -> Result<Point3, String> {
    if let Some(index) = label.sector_index() {
        Ok(*cap_net
            .sectors
            .get(sector)
            .ok_or_else(|| {
                format!(
                    "KP20 cap-net figure sector {sector} outside valence {}",
                    cap_net.sectors.len()
                )
            })?
            .get(index)
            .expect("figure label index is in c1..c6 range"))
    } else {
        Ok(cap_net.center)
    }
}

/// Return the three radial columns of one sector in the verified
/// KP15 Fig. 1a label order.
pub fn kp20_cap_net_figure_sector_columns(
    cap_net: &Kp20CapNet,
    sector: usize,
) -> Result<[[Point3; 2]; 3], String> {
    Ok([
        [
            kp20_cap_net_figure_point(cap_net, sector, Kp20CapNetFigureLabel::C4)?,
            kp20_cap_net_figure_point(cap_net, sector, Kp20CapNetFigureLabel::C1)?,
        ],
        [
            kp20_cap_net_figure_point(cap_net, sector, Kp20CapNetFigureLabel::C5)?,
            kp20_cap_net_figure_point(cap_net, sector, Kp20CapNetFigureLabel::C2)?,
        ],
        [
            kp20_cap_net_figure_point(cap_net, sector, Kp20CapNetFigureLabel::C6)?,
            kp20_cap_net_figure_point(cap_net, sector, Kp20CapNetFigureLabel::C3)?,
        ],
    ])
}

/// Build the planar Catmull-Clark characteristic cap-net used by
/// KP18/KP20.
///
/// The six sector values are Stam's frequency-one right eigenvector
/// entries, lifted around the extraordinary vertex by the inverse DFT:
/// `c^k_m := e_m (cos(2πk/n), sin(2πk/n), 0)`. The center value is
/// zero because the subdominant characteristic eigenfunctions vanish
/// at the extraordinary vertex. This creates the c-net-level source
/// for the characteristic ring; conversion to a sector tensor-border
/// still depends on the surrounding 4x4 stencil layout.
pub fn kp20_characteristic_cap_net(n: usize) -> Result<Kp20CapNet, String> {
    let scalars = Kp20CharacteristicScalars::for_valence(n)?;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = core::f64::consts::TAU * k as f64 / n as f64;
        let (s, c) = theta.sin_cos();
        let mut sector = [Point3::new(0.0, 0.0, 0.0); 6];
        for (m, cp) in sector.iter_mut().enumerate() {
            let value = scalars.eigenvector[m];
            *cp = Point3::new(value * c, value * s, 0.0);
        }
        sectors.push(sector);
    }
    Ok(Kp20CapNet {
        center: Point3::new(0.0, 0.0, 0.0),
        sectors,
    })
}

/// Diagnostic Fig. 9 cap-net row extraction for the characteristic tensor-border.
///
/// This helper deliberately encodes the direct visual inference from KP18
/// Fig. 9: the outer sector row is `c1,c2,c3`, the next row is `c4,c5,c6`,
/// and the innermost regular row collapses to the extraordinary node. It is
/// useful as a measured candidate for the KP20 characteristic-border campaign,
/// but it is not the accepted paper source until the regular-ring extraction
/// is derived from the Catmull-Clark subdivision matrix.
#[cfg(test)]
pub fn kp20_cap_net_figure_tensor_border_rows(
    cap_net: &Kp20CapNet,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let n = cap_net.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP20 cap-net tensor-border row extraction needs at least 3 sectors, got {n}"
        ));
    }
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let next = (k + 1) % n;
        let current = cap_net
            .sectors
            .get(k)
            .expect("sector index is bounded by cap-net valence");
        let next_sector = cap_net
            .sectors
            .get(next)
            .expect("next sector index is bounded by cap-net valence");
        sectors.push([
            [current[0], current[1], current[2], next_sector[0]],
            [current[3], current[4], current[5], next_sector[3]],
            [
                cap_net.center,
                cap_net.center,
                cap_net.center,
                cap_net.center,
            ],
        ]);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
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

/// Convert regular 4x3 B-spline control rows into KP20's depth-2 bi-3
/// tensor-border.
///
/// KP20 Fig. 3b uses exactly the first three BB rows of the surrounding
/// bi-3 spline ring. Under the uniform cubic B-spline to Bezier matrix,
/// those rows depend only on the first three source rows, so a full 4x4
/// regular stencil is not required for tensor-border extraction.
pub fn kp20_tensor_border_from_regular_cnet_rows(
    cnet: &Kp20RegularTensorBorderCNet,
) -> Result<TensorBorderBi3, String> {
    if cnet.sectors.len() < 3 {
        return Err(format!(
            "KP20 c-net tensor-border row conversion needs at least 3 sectors, got {}",
            cnet.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(cnet.sectors.len());
    for (sector_index, rows) in cnet.sectors.iter().enumerate() {
        validate_border_rows(sector_index, rows)?;
        sectors.push(cubic_bspline_rows_to_tensor_border_rows(rows));
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

/// Apply the cubic B-spline to Bezier matrix to the first three tensor-border
/// rows.
pub fn cubic_bspline_rows_to_tensor_border_rows(rows: &[[Point3; 4]; 3]) -> [[Point3; 4]; 3] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = tensor_border_weighted_point(i, j, rows);
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

fn tensor_border_weighted_point(i: usize, j: usize, rows: &[[Point3; 4]; 3]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for row in 0..3 {
        for col in 0..4 {
            let w = CUBIC_BSPLINE_TO_BEZIER[i][row] * CUBIC_BSPLINE_TO_BEZIER[j][col];
            out.x += w * rows[row][col].x;
            out.y += w * rows[row][col].y;
            out.z += w * rows[row][col].z;
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

fn validate_border_rows(sector_index: usize, rows: &[[Point3; 4]; 3]) -> Result<(), String> {
    for (row, cps) in rows.iter().enumerate() {
        for (col, p) in cps.iter().enumerate() {
            if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
                return Err(format!(
                    "KP20 c-net tensor-border sector {sector_index} rows[{row}][{col}] is non-finite"
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "kp20_bi6_cnet_tests.rs"]
mod tests;
