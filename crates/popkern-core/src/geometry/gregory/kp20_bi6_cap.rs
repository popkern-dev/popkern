//! KP20 side-by-side bi-6 macro-patch infrastructure.
//!
//! This is the first production-shaped slice of the KP20 backend. It
//! does not replace the KP15 single-patch cap. It owns the 2x2 bi-6
//! macro-patch data model and the paper-level boundary operations that
//! can be verified before guide sampling and S0 optimisation are wired.

#![allow(clippy::needless_range_loop)]

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp_polar_shared::weighted_sum;
use super::kp20_bi6_math::kp20_g2_strip_weights;
use super::kp20_border_degree::{elevate_bi3_rows_to_sextic, scalar_cubic_times_sextic_reduced};
use crate::geometry::Point3;

/// KP20 sector quadrant in canonical paper order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20PieceId {
    /// Piece 1: closer to the extraordinary point, left of `o -> m`.
    NearOriginLeft,
    /// Piece 2: closer to the extraordinary point, right of `o -> m`.
    NearOriginRight,
    /// Piece 3: closer to the outer midpoint, right of `o -> m`.
    NearMidpointRight,
    /// Piece 4: closer to the outer midpoint, left of `o -> m`.
    NearMidpointLeft,
}

impl Kp20PieceId {
    /// All pieces in the stable sector-local storage order.
    pub const fn all() -> [Self; 4] {
        [
            Self::NearOriginLeft,
            Self::NearOriginRight,
            Self::NearMidpointRight,
            Self::NearMidpointLeft,
        ]
    }

    /// Stable zero-based index into [`Kp20MacroPatch::pieces`].
    pub const fn index(self) -> usize {
        match self {
            Self::NearOriginLeft => 0,
            Self::NearOriginRight => 1,
            Self::NearMidpointRight => 2,
            Self::NearMidpointLeft => 3,
        }
    }
}

/// One bi-6 tensor-product patch in a KP20 macro sector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20Bi6Patch {
    /// Tensor-product BB control grid, indexed `[u][v]`.
    pub cps: [[Point3; 7]; 7],
}

impl Kp20Bi6Patch {
    /// Construct a patch with all BB coefficients at the origin.
    pub fn zero() -> Self {
        Self {
            cps: [[Point3::new(0.0, 0.0, 0.0); 7]; 7],
        }
    }
}

/// The four bi-6 pieces for one sector of the KP20 2x2 macro-patch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20MacroPatch {
    /// Four sector quadrants in [`Kp20PieceId::all`] order.
    pub pieces: [Kp20Bi6Patch; 4],
}

impl Kp20MacroPatch {
    /// Construct a macro-patch whose four pieces are zero patches.
    pub fn zero() -> Self {
        Self {
            pieces: [Kp20Bi6Patch::zero(); 4],
        }
    }

    /// Immutable access to a sector quadrant by semantic piece id.
    pub fn piece(&self, id: Kp20PieceId) -> &Kp20Bi6Patch {
        &self.pieces[id.index()]
    }

    /// Mutable access to a sector quadrant by semantic piece id.
    pub fn piece_mut(&mut self, id: Kp20PieceId) -> &mut Kp20Bi6Patch {
        &mut self.pieces[id.index()]
    }
}

/// Degree-6 rows produced from one half of a split bi-3 tensor-border.
///
/// Rows `0..=2` are the position, first-cross, and second-cross
/// boundary rows after KP20 equation (14)'s reparameterisation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20BoundaryHalf {
    /// Degree-6 boundary rows `0..=2`, indexed `[depth][u]`.
    pub rows: [[Point3; 7]; 3],
}

/// Per-sector pair of reparameterised tensor-border halves.
///
/// `near_c` is the first half of the sector-local tensor-border span
/// and `near_m` is the second half. A sector-separating strip uses
/// `near_c` from the two adjacent sectors; it does not use both halves
/// from the same sector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20SplitTensorBorder {
    /// Half closer to the cap-net node `c`.
    pub near_c: Kp20BoundaryHalf,
    /// Half closer to the midpoint `m`.
    pub near_m: Kp20BoundaryHalf,
}

/// Split and reparameterise all bi-3 tensor-border sectors for KP20
/// equation (14).
pub fn split_kp20_tensor_border(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP20 tensor-border needs at least 3 sectors, got {n}"
        ));
    }
    let c = (core::f64::consts::TAU / n as f64).cos();
    split_kp20_tensor_border_with_d1_3(input, 2.0 * c / 3.0)
}

/// Variant of [`split_kp20_tensor_border`] with an explicit KP20
/// equation (14) `d¹_3` value.
pub fn split_kp20_tensor_border_with_d1_3(
    input: &TensorBorderBi3,
    d1_3: f64,
) -> Result<Vec<Kp20SplitTensorBorder>, String> {
    let n = input.sectors.len();
    if n < 3 {
        return Err(format!(
            "KP20 tensor-border needs at least 3 sectors, got {n}"
        ));
    }
    if !d1_3.is_finite() {
        return Err(format!(
            "KP20 tensor-border d1_3 must be finite, got {d1_3}"
        ));
    }
    let c = (core::f64::consts::TAU / n as f64).cos();
    let (d0, d1) = kp20_eq14_d_profiles(c, d1_3);

    let mut out = Vec::with_capacity(n);
    for (sector_index, sector) in input.sectors.iter().enumerate() {
        validate_bi3_sector(sector_index, sector)?;
        let mut left = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        let mut right = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for row in 0..3 {
            let (l, r) = split_cubic_half(sector[row]);
            left[row] = l;
            right[row] = r;
        }
        out.push(Kp20SplitTensorBorder {
            near_c: reparameterise_half(elevate_bi3_rows_to_sextic(left), d0),
            near_m: reparameterise_half(elevate_bi3_rows_to_sextic(right), d1),
        });
    }
    Ok(out)
}

fn kp20_eq14_d_profiles(c: f64, d1_3: f64) -> ([f64; 4], [f64; 4]) {
    let delta = d1_3 / 4.0 - c / 12.0;
    (
        [0.0, 0.0, 0.0, delta],
        [delta, 2.0 * delta, 4.0 * delta, d1_3],
    )
}

/// Initialise one KP20 macro-patch per tensor-border sector and assign
/// the m-side boundary rows from the split/reparameterised border.
pub fn initialise_kp20_macro_patches_from_tensor_border(
    input: &TensorBorderBi3,
) -> Result<Vec<Kp20MacroPatch>, String> {
    let split = split_kp20_tensor_border(input)?;
    let mut sectors = vec![Kp20MacroPatch::zero(); split.len()];
    for k in 0..split.len() {
        let previous = (k + split.len() - 1) % split.len();
        sectors[k].piece_mut(Kp20PieceId::NearMidpointLeft).cps[0..3]
            .copy_from_slice(&split[previous].near_c.rows);
        sectors[k].piece_mut(Kp20PieceId::NearMidpointRight).cps[0..3]
            .copy_from_slice(&split[k].near_c.rows);
    }
    Ok(sectors)
}

/// Apply KP20 equation (12): derive rows 0, 1, and 5 of an upper
/// piece from rows 2, 4, and 6 of the lower piece.
pub fn kp20_propagate_c2_upper_from_lower(lower: &Kp20Bi6Patch, upper: &mut Kp20Bi6Patch) {
    upper.cps[0] = lower.cps[6];
    for j in 0..7 {
        let half_delta = weighted_sum(&[(0.25, lower.cps[4][j]), (-0.25, lower.cps[2][j])]);
        upper.cps[5][j] = lower.cps[6][j] + (half_delta - Point3::new(0.0, 0.0, 0.0));
        upper.cps[1][j] = lower.cps[6][j] - (half_delta - Point3::new(0.0, 0.0, 0.0));
    }
}

/// Apply KP20 equation (13) to the sector-separating G2 strip rows.
///
/// The two patches are the paired pieces on opposite sides of the
/// sector-separating curve. Their `j = 1` rows are overwritten for
/// `i = 4, 5, 6`.
pub fn kp20_apply_g2_strip(p: &mut Kp20Bi6Patch, q: &mut Kp20Bi6Patch, b0: f64) {
    let (w1, w2, w3) = kp20_g2_strip_weights(b0);
    for i in 5..=6 {
        let quarter = weighted_sum(&[(0.25, p.cps[i][2]), (-0.25, q.cps[i][2])]);
        p.cps[i][1] = p.cps[i][0] + (quarter - Point3::new(0.0, 0.0, 0.0));
        q.cps[i][1] = p.cps[i][0] - (quarter - Point3::new(0.0, 0.0, 0.0));
    }

    let quarter_42 = weighted_sum(&[(0.25, p.cps[4][2]), (-0.25, q.cps[4][2])]);
    let correction = weighted_sum(&[
        (-w3, p.cps[5][0]),
        (w3, p.cps[6][0]),
        (w2, p.cps[5][2]),
        (-w2, q.cps[5][2]),
        (w1, p.cps[6][2]),
        (-w1, q.cps[6][2]),
    ]);
    p.cps[4][1] = p.cps[4][0]
        + (quarter_42 - Point3::new(0.0, 0.0, 0.0))
        + (correction - Point3::new(0.0, 0.0, 0.0));
    q.cps[4][1] = p.cps[4][0] - (quarter_42 - Point3::new(0.0, 0.0, 0.0))
        + (correction - Point3::new(0.0, 0.0, 0.0));
}

fn split_cubic_half(row: [Point3; 4]) -> ([Point3; 4], [Point3; 4]) {
    let p10 = midpoint(row[0], row[1]);
    let p11 = midpoint(row[1], row[2]);
    let p12 = midpoint(row[2], row[3]);
    let p20 = midpoint(p10, p11);
    let p21 = midpoint(p11, p12);
    let p30 = midpoint(p20, p21);
    ([row[0], p10, p20, p30], [p30, p21, p12, row[3]])
}

/// Reparameterise one raw tensor-border piece (degree 6 along, three
/// depth rows) by KP20 equation (14)'s ρ with an explicit cubic `d`
/// profile (BB coefficients along the piece's own parameter).
///
/// Public entry for callers whose border pieces are natively split
/// (e.g. the characteristic ring, whose sector boundary has a genuine
/// parametric corner at `c`), so the de Casteljau halving of
/// [`split_kp20_tensor_border`] does not apply. Cubic pieces elevate
/// exactly first (`elevate_bi3_rows_to_sextic`).
pub fn kp20_reparameterise_border_piece(rows: [[Point3; 7]; 3], d: [f64; 4]) -> Kp20BoundaryHalf {
    reparameterise_half(rows, d)
}

/// KP20 equation (14) `d` profiles `(d⁰, d¹)` for the two border
/// pieces of one sector, given `c = cos(2π/n)` and the free parameter
/// `d¹₃` (paper default `2c/3`).
pub fn kp20_eq14_border_d_profiles(c: f64, d1_3: f64) -> ([f64; 4], [f64; 4]) {
    kp20_eq14_d_profiles(c, d1_3)
}

/// Equation (14) on degree-6 rows: `row1 = ½(row0 + row1)`,
/// `row2 = 0.2·row0 + 0.6·row1 + 0.2·row2 + 0.1·d·(row1 − row0)`. The
/// product is degree 9 and is reduced to degree 6 with both end
/// 2-jets preserved exactly (see `kp20_border_degree`); for
/// elevated-cubic input the reduction is exact and reproduces the
/// historical cubic pipeline.
fn reparameterise_half(rows: [[Point3; 7]; 3], d: [f64; 4]) -> Kp20BoundaryHalf {
    let source_row0 = rows[0];
    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 7];
    let mut row2 = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0_usize..7 {
        row1[i] = weighted_sum(&[(0.5, rows[0][i]), (0.5, rows[1][i])]);
        row2[i] = weighted_sum(&[(0.2, rows[0][i]), (0.6, rows[1][i]), (0.2, rows[2][i])]);
    }
    let diff: [crate::geometry::Vec3; 7] = core::array::from_fn(|i| rows[1][i] - rows[0][i]);
    let d_times_diff = scalar_cubic_times_sextic_reduced(d, diff);
    for i in 0_usize..7 {
        row2[i] = row2[i] + d_times_diff[i] * 0.1;
    }
    Kp20BoundaryHalf {
        rows: [source_row0, row1, row2],
    }
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    weighted_sum(&[(0.5, a), (0.5, b)])
}

fn validate_bi3_sector(sector_index: usize, sector: &[[Point3; 4]; 3]) -> Result<(), String> {
    for row in 0..3 {
        for col in 0..4 {
            let p = sector[row][col];
            if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
                return Err(format!(
                    "KP20 tensor-border sector {sector_index} row {row} CP {col} is non-finite"
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "kp20_bi6_cap_tests.rs"]
mod tests;
