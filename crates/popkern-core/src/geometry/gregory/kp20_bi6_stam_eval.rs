//! Stam/PR08 Catmull-Clark characteristic evaluator substrate.
//!
//! This is a diagnostic scaffold for the KP20 characteristic tensor-border
//! campaign. It encodes the source-backed parts of Stam 1998 §3: the three
//! regular bicubic B-spline patches, the dyadic tile dispatch, and the
//! characteristic-ring self-similarity. The remaining open step is mapping
//! PR08 Fig. 6.3's q-labels into these three regular patch grids.

#[cfg(test)]
use super::kp20_bi6_cc_ring::Kp20CatmullClarkCharacteristicRing;
#[cfg(test)]
use super::kp20_bi6_cnet::Kp20RegularTensorBorderCNet;
use super::kp20_bi6_cnet::cubic_bspline_grid_to_bezier;
use crate::geometry::Point3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20StamPatch {
    LowerRight,
    CenterRight,
    UpperLeft,
}

impl Kp20StamPatch {
    const fn index(self) -> usize {
        match self {
            Self::LowerRight => 0,
            Self::CenterRight => 1,
            Self::UpperLeft => 2,
        }
    }
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20StamGridSide {
    Top,
    Bottom,
    Left,
    Right,
}

#[cfg(test)]
impl Kp20StamGridSide {
    pub const fn all() -> [Self; 4] {
        [Self::Top, Self::Bottom, Self::Left, Self::Right]
    }

    pub const fn name(self) -> &'static str {
        match self {
            Self::Top => "top",
            Self::Bottom => "bottom",
            Self::Left => "left",
            Self::Right => "right",
        }
    }
}

#[cfg(test)]
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20StamGridTransform {
    Identity,
    Rotate90,
    Rotate180,
    Rotate270,
    Transpose,
    AntiTranspose,
    FlipRows,
    FlipColumns,
}

#[cfg(test)]
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20Pr08StamProjectedGrids {
    pub sectors: Vec<[[[Point3; 4]; 4]; 3]>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20StamTile {
    pub patch: Kp20StamPatch,
    pub level: u32,
    pub u: f64,
    pub v: f64,
}

/// Dispatch `(u,v)` to Stam's dyadic tiles and local patch coordinates.
///
/// The returned `level` is the `n` in Stam equations (5)-(9), so the
/// characteristic scaling factor is `lambda.powi(level - 1)`.
pub fn kp20_stam_tile_for(u: f64, v: f64) -> Result<Kp20StamTile, String> {
    if !u.is_finite() || !v.is_finite() || u < 0.0 || v < 0.0 || u > 1.0 || v > 1.0 {
        return Err(format!(
            "Stam tile coordinates must be finite in [0,1]^2, got ({u}, {v})"
        ));
    }
    if u == 0.0 && v == 0.0 {
        return Ok(Kp20StamTile {
            patch: Kp20StamPatch::LowerRight,
            level: 1,
            u: 0.0,
            v: 0.0,
        });
    }

    let mut level = 1_u32;
    let mut scale = 2.0_f64;
    loop {
        let lower = 1.0 / scale;
        let upper = 2.0 / scale;
        if u >= lower && u <= upper && v <= lower {
            return Ok(Kp20StamTile {
                patch: Kp20StamPatch::LowerRight,
                level,
                u: scale * u - 1.0,
                v: scale * v,
            });
        }
        if u >= lower && u <= upper && v >= lower && v <= upper {
            return Ok(Kp20StamTile {
                patch: Kp20StamPatch::CenterRight,
                level,
                u: scale * u - 1.0,
                v: scale * v - 1.0,
            });
        }
        if u <= lower && v >= lower && v <= upper {
            return Ok(Kp20StamTile {
                patch: Kp20StamPatch::UpperLeft,
                level,
                u: scale * u,
                v: scale * v - 1.0,
            });
        }
        level += 1;
        scale *= 2.0;
        if level > 60 {
            return Err(format!(
                "Stam tile dispatch did not converge for ({u}, {v})"
            ));
        }
    }
}

pub fn kp20_stam_eval_characteristic(
    grids: &[[[Point3; 4]; 4]; 3],
    lambda: f64,
    u: f64,
    v: f64,
) -> Result<Point3, String> {
    if !lambda.is_finite() || lambda <= 0.0 || lambda >= 1.0 {
        return Err(format!(
            "Stam characteristic evaluator needs lambda in (0,1), got {lambda}"
        ));
    }
    let tile = kp20_stam_tile_for(u, v)?;
    let value = eval_cubic_bspline_grid(&grids[tile.patch.index()], tile.u, tile.v)?;
    Ok(scale_point(value, lambda.powi(tile.level as i32 - 1)))
}

pub fn eval_cubic_bspline_grid(grid: &[[Point3; 4]; 4], u: f64, v: f64) -> Result<Point3, String> {
    if !u.is_finite()
        || !v.is_finite()
        || u < -1e-12
        || v < -1e-12
        || u > 1.0 + 1e-12
        || v > 1.0 + 1e-12
    {
        return Err(format!(
            "cubic B-spline patch coordinates must be finite in [0,1]^2, got ({u}, {v})"
        ));
    }
    let bezier = cubic_bspline_grid_to_bezier(grid);
    let bu = cubic_bernstein(u.clamp(0.0, 1.0));
    let bv = cubic_bernstein(v.clamp(0.0, 1.0));
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (i, bu_i) in bu.iter().enumerate() {
        for (j, bv_j) in bv.iter().enumerate() {
            let w = bu_i * bv_j;
            out.x += w * bezier[i][j].x;
            out.y += w * bezier[i][j].y;
            out.z += w * bezier[i][j].z;
        }
    }
    Ok(out)
}

/// Stam Fig. 6's three 4x4 picking grids, using the one-based control
/// vertex numbers printed in the paper for a valence-`n` patch.
pub fn kp20_stam_figure6_picking_grids(n: usize) -> Result<[[[usize; 4]; 4]; 3], String> {
    if n < 3 {
        return Err(format!("Stam Fig. 6 picking grids need n >= 3, got {n}"));
    }
    let b = 2 * n;
    Ok([
        [
            [b + 7, b + 6, b + 2, b + 10],
            [4, 5, b + 3, b + 11],
            [1, 6, b + 4, b + 12],
            [8, 7, b + 5, b + 13],
        ],
        [
            [b + 16, b + 15, b + 14, b + 9],
            [b + 7, b + 6, b + 2, b + 10],
            [4, 5, b + 3, b + 11],
            [1, 6, b + 4, b + 12],
        ],
        [
            [b + 17, b + 16, b + 15, b + 14],
            [b + 8, b + 7, b + 6, b + 2],
            [3, 4, 5, b + 3],
            [2, 1, 6, b + 4],
        ],
    ])
}

#[cfg(test)]
pub fn kp20_pr08_stam_visual_projected_grids(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Result<Kp20Pr08StamProjectedGrids, String> {
    if ring.sectors.len() < 3 {
        return Err(format!(
            "PR08-to-Stam visual projection needs at least 3 sectors, got {}",
            ring.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(ring.sectors.len());
    for sector_index in 0..ring.sectors.len() {
        let previous = &ring.sectors[(sector_index + ring.sectors.len() - 1) % ring.sectors.len()];
        let sector = &ring.sectors[sector_index];
        let next = &ring.sectors[(sector_index + 1) % ring.sectors.len()];

        let lower_right = [
            [sector[3], sector[6], sector[11], sector[12]],
            [sector[2], sector[5], sector[9], sector[10]],
            [sector[0], sector[1], sector[4], sector[7]],
            [previous[0], previous[2], previous[3], previous[12]],
        ];
        let center_right = [
            [next[1], next[4], next[7], next[8]],
            lower_right[0],
            lower_right[1],
            lower_right[2],
        ];
        let upper_left = [
            [next[3], next[6], next[11], next[12]],
            [next[2], sector[3], sector[6], sector[11]],
            [sector[2], sector[3], sector[5], sector[9]],
            [next[0], sector[0], sector[1], sector[4]],
        ];
        sectors.push([lower_right, center_right, upper_left]);
    }
    Ok(Kp20Pr08StamProjectedGrids { sectors })
}

#[cfg(test)]
pub fn kp20_pr08_stam_appendix_tensor_projected_grids(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Result<Kp20Pr08StamProjectedGrids, String> {
    let visual = kp20_pr08_stam_visual_projected_grids(ring)?;
    let mut sectors = Vec::with_capacity(visual.sectors.len());
    for grids in &visual.sectors {
        sectors.push([
            visual_grid_to_appendix_tensor_order(&grids[Kp20StamPatch::LowerRight.index()]),
            visual_grid_to_appendix_tensor_order(&grids[Kp20StamPatch::CenterRight.index()]),
            visual_grid_to_appendix_tensor_order(&grids[Kp20StamPatch::UpperLeft.index()]),
        ]);
    }
    Ok(Kp20Pr08StamProjectedGrids { sectors })
}

#[cfg(test)]
pub fn kp20_pr08_middle_normalized_projected_grids(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Result<Kp20Pr08StamProjectedGrids, String> {
    if ring.sectors.len() < 3 {
        return Err(format!(
            "PR08 middle-normalized projection needs at least 3 sectors, got {}",
            ring.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(ring.sectors.len());
    for sector_index in 0..ring.sectors.len() {
        let sector = &ring.sectors[sector_index];
        let next = &ring.sectors[(sector_index + 1) % ring.sectors.len()];

        // PR08 p.114 gives the middle patch normalization stencil as the
        // lower-right 3x3 tensor-product block
        //   [q3 q4 q13; q6 q7 q12; q9 q10 q13].
        // The first row/column are the adjacent visible Fig. 6.3 labels that
        // complete that block into one regular 4x4 B-spline grid.
        let middle_normalized = [
            [next[1], next[2], next[5], next[8]],
            [sector[1], sector[2], sector[3], sector[12]],
            [sector[4], sector[5], sector[6], sector[11]],
            [sector[7], sector[8], sector[9], sector[12]],
        ];
        sectors.push([middle_normalized, middle_normalized, middle_normalized]);
    }
    Ok(Kp20Pr08StamProjectedGrids { sectors })
}

#[cfg(test)]
pub fn kp20_pr08_middle_patch_normalization_value(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_index: usize,
) -> Result<Point3, String> {
    let sector = ring.sectors.get(sector_index).ok_or_else(|| {
        format!(
            "PR08 normalization sector {sector_index} outside valence {}",
            ring.sectors.len()
        )
    })?;
    let weights = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
    let stencil = [
        [sector[2], sector[3], sector[12]],
        [sector[5], sector[6], sector[11]],
        [sector[8], sector[9], sector[12]],
    ];
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (i, wi) in weights.iter().enumerate() {
        for (j, wj) in weights.iter().enumerate() {
            let weight = wi * wj;
            out.x += weight * stencil[i][j].x;
            out.y += weight * stencil[i][j].y;
            out.z += weight * stencil[i][j].z;
        }
    }
    Ok(out)
}

#[cfg(test)]
pub fn kp20_pr08_stam_projected_tensor_border_side_rows(
    projection: &Kp20Pr08StamProjectedGrids,
    patch: Kp20StamPatch,
    side: Kp20StamGridSide,
    reversed_depth: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    if projection.sectors.len() < 3 {
        return Err(format!(
            "PR08-to-Stam projected tensor-border rows need at least 3 sectors, got {}",
            projection.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(projection.sectors.len());
    for grids in &projection.sectors {
        let mut rows = side_rows(&grids[patch.index()], side);
        if reversed_depth {
            rows = [rows[2], rows[1], rows[0]];
        }
        sectors.push(rows);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

#[cfg(test)]
pub fn kp20_pr08_stam_projected_transformed_tensor_border_side_rows(
    projection: &Kp20Pr08StamProjectedGrids,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
    side: Kp20StamGridSide,
    reversed_depth: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    if projection.sectors.len() < 3 {
        return Err(format!(
            "PR08-to-Stam transformed tensor-border rows need at least 3 sectors, got {}",
            projection.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(projection.sectors.len());
    for grids in &projection.sectors {
        let transformed = transform_grid(&grids[patch.index()], transform);
        let mut rows = side_rows(&transformed, side);
        if reversed_depth {
            rows = [rows[2], rows[1], rows[0]];
        }
        sectors.push(rows);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

#[cfg(test)]
pub fn kp20_pr08_stam_projected_transformed_grid(
    projection: &Kp20Pr08StamProjectedGrids,
    sector: usize,
    patch: Kp20StamPatch,
    transform: Kp20StamGridTransform,
) -> Result<[[Point3; 4]; 4], String> {
    let grids = projection.sectors.get(sector).ok_or_else(|| {
        format!(
            "PR08-to-Stam transformed grid sector {sector} outside valence {}",
            projection.sectors.len()
        )
    })?;
    Ok(transform_grid(&grids[patch.index()], transform))
}

fn scale_point(p: Point3, scale: f64) -> Point3 {
    Point3::new(scale * p.x, scale * p.y, scale * p.z)
}

#[cfg(test)]
fn side_rows(grid: &[[Point3; 4]; 4], side: Kp20StamGridSide) -> [[Point3; 4]; 3] {
    match side {
        Kp20StamGridSide::Top => [grid[0], grid[1], grid[2]],
        Kp20StamGridSide::Bottom => [grid[3], grid[2], grid[1]],
        Kp20StamGridSide::Left => [
            [grid[0][0], grid[1][0], grid[2][0], grid[3][0]],
            [grid[0][1], grid[1][1], grid[2][1], grid[3][1]],
            [grid[0][2], grid[1][2], grid[2][2], grid[3][2]],
        ],
        Kp20StamGridSide::Right => [
            [grid[0][3], grid[1][3], grid[2][3], grid[3][3]],
            [grid[0][2], grid[1][2], grid[2][2], grid[3][2]],
            [grid[0][1], grid[1][1], grid[2][1], grid[3][1]],
        ],
    }
}

#[cfg(test)]
// Index form: each `(row, col)` destination maps through a symmetry-dependent
// `(source_row, source_col)`; the index arithmetic IS the grid transform.
#[allow(clippy::needless_range_loop)]
fn transform_grid(grid: &[[Point3; 4]; 4], transform: Kp20StamGridTransform) -> [[Point3; 4]; 4] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for row in 0..4 {
        for col in 0..4 {
            let (source_row, source_col) = match transform {
                Kp20StamGridTransform::Identity => (row, col),
                Kp20StamGridTransform::Rotate90 => (3 - col, row),
                Kp20StamGridTransform::Rotate180 => (3 - row, 3 - col),
                Kp20StamGridTransform::Rotate270 => (col, 3 - row),
                Kp20StamGridTransform::Transpose => (col, row),
                Kp20StamGridTransform::AntiTranspose => (3 - col, 3 - row),
                Kp20StamGridTransform::FlipRows => (3 - row, col),
                Kp20StamGridTransform::FlipColumns => (row, 3 - col),
            };
            out[row][col] = grid[source_row][source_col];
        }
    }
    out
}

#[cfg(test)]
fn visual_grid_to_appendix_tensor_order(grid: &[[Point3; 4]; 4]) -> [[Point3; 4]; 4] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for u in 0..4 {
        for v in 0..4 {
            out[u][v] = grid[3 - v][u];
        }
    }
    out
}

fn cubic_bernstein(t: f64) -> [f64; 4] {
    let s = 1.0 - t;
    [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t]
}

#[cfg(test)]
#[path = "kp20_bi6_stam_eval_tests.rs"]
mod tests;
