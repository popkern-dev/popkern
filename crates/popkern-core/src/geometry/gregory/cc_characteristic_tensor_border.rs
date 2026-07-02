//! Characteristic tensor-border in the KP20 sector chart.
//!
//! KP20 §4.3 constructs the reparameterization τ so that each sector is
//! "G²-connected using ρ with the tensor-border of the characteristic
//! ring of Catmull-Clark subdivision". This module extracts that
//! tensor-border from the reference-grade `cc_characteristic_ring` and
//! expresses it in the chart the KP20 τ/cap substrate already uses:
//!
//! - Sector `k` spans the wedge between spoke `k` (angle `2πk/n`) and
//!   spoke `k+1`, counterclockwise, in the unit n-gon frame: the spoke
//!   endpoints `m^k` (where the sector-separating curves meet the
//!   ring) are normalized onto the unit vertices `e^{2πik/n}`.
//! - The sector's outer boundary consists of TWO bi-3 pieces meeting
//!   at the sector corner `c` (a genuine parametric corner of the
//!   characteristic ring's inner boundary): `near_c` runs from `m^k`
//!   (along index 0) to `c` (along index 3), `near_m` from `c` to
//!   `m^{k+1}`. This matches the operational convention of
//!   `split_kp20_tensor_border`, whose de Casteljau halving is the
//!   smooth-rim production adaptation of the same two-piece layout.
//! - Rows are indexed `[depth][along]` with depth 0 on the boundary
//!   and depths 1..2 the cap-side C² prolongation of the ring surface
//!   (Bernstein normal reflection `r₁ = 2q₀ − q₁`,
//!   `r₂ = 4q₀ − 4q₁ + q₂` of the ring-side BB rows), exactly the
//!   convention of the planar regular fixture that closes the τ rows.
//!
//! Parametric scales match by construction: the level-1 ring knot
//! interval equals the cap's 2×2 macro-piece size (the dyadic
//! consistency that motivates KP20's 2×2 layout), so no additional
//! derivative scaling is applied.

use super::cc_characteristic_ring::{CcCharacteristicRing, cc_characteristic_ring};
use crate::geometry::Point3;

/// One bi-3 tensor-border piece: rows `[depth][along]`, depth 0 on the
/// cap boundary, along running with the sector's counterclockwise
/// orientation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CcBorderPiece {
    /// Bi-3 BB rows (planar: `z = 0`).
    pub rows: [[Point3; 4]; 3],
}

/// The two-piece characteristic tensor-border of one sector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CcSectorBorder {
    /// Piece from `m^k` (along 0) to the sector corner `c` (along 3).
    pub near_c: CcBorderPiece,
    /// Piece from the sector corner `c` (along 0) to `m^{k+1}`
    /// (along 3).
    pub near_m: CcBorderPiece,
}

/// Characteristic tensor-border for all sectors, in the KP20 unit
/// n-gon frame.
#[derive(Debug, Clone, PartialEq)]
pub struct CcCharacteristicTensorBorder {
    /// Valence.
    pub n: usize,
    /// Catmull-Clark subdominant eigenvalue.
    pub lambda: f64,
    /// Per-sector border pieces, sector `k` between spokes `k`, `k+1`.
    pub sectors: Vec<CcSectorBorder>,
    /// The normalized sector-corner value `c` of sector 0 (the other
    /// sectors' corners are its rotations).
    pub corner: [f64; 2],
}

/// Extract the characteristic tensor-border for valence `n`, in the
/// KP20 unit n-gon frame (spoke endpoints on unit vertices,
/// counterclockwise).
pub fn cc_characteristic_tensor_border(n: usize) -> Result<CcCharacteristicTensorBorder, String> {
    let ring = cc_characteristic_ring(n)?;

    // Ring-side 4x3 B-spline windows of the two inner-boundary edges
    // of one sector (one Stam patch quadrant [0, 1/2]^2):
    //
    // - E_u = {u = 1/2} x [0, 1/2]: the s = 0 edge of subpatch
    //   Omega_1, whose picked grid is grids[0][t_row][s_col]. Depth
    //   direction = s (columns), along = t (rows), t = 0 at the v = 0
    //   spoke and t = 1 at the corner c.
    // - E_v = [0, 1/2] x {v = 1/2}: the t = 0 edge of Omega_3
    //   (grids[2]). Depth = t (rows), along = s (columns), s = 0 at
    //   the u = 0 spoke and s = 1 at c.
    let e_u = ring_side_bb_block(&ring, EdgeSelect::OmegaOne);
    let e_v = ring_side_bb_block(&ring, EdgeSelect::OmegaThree);

    // Normalization: send the v = 0 spoke endpoint m = chi(1/2, 0) to
    // the unit vertex 1 + 0i (KP20's tau boundary ends within ~1% of
    // the unit n-gon vertex: tau_03 = mu0 * v-dot).
    let m_u = bspline_corner_value(&ring.grids[0]);
    let m_v = bspline_corner_value(&ring.grids[2]);
    let m_u_norm_sq = m_u[0] * m_u[0] + m_u[1] * m_u[1];
    if m_u_norm_sq <= 1.0e-30 {
        return Err(format!(
            "characteristic border spoke endpoint degenerate for n = {n}"
        ));
    }
    let zeta = [m_u[0] / m_u_norm_sq, -m_u[1] / m_u_norm_sq];

    // Orientation: after normalization the other spoke endpoint must
    // sit at +2*pi/n (counterclockwise). If it lands at -2*pi/n,
    // conjugate the whole ring (the seed's phase convention is
    // orientation-agnostic; KP20's chart is counterclockwise).
    let m_v_normalized = complex_mul(zeta, m_v);
    let conjugate = m_v_normalized[1] < 0.0;

    let normalize = |value: [f64; 2]| -> [f64; 2] {
        let scaled = complex_mul(zeta, value);
        if conjugate {
            [scaled[0], -scaled[1]]
        } else {
            scaled
        }
    };

    // Sector 0's pieces in the counterclockwise chart:
    // near_c runs from m^0 = 1 (the normalized E_u start) to c along
    // E_u; near_m runs from c to m^1 = e^{2*pi*i/n} along E_v
    // REVERSED (E_v's natural along direction starts at the u = 0
    // spoke).
    let near_c_rows = map_block(&e_u, normalize, false);
    let near_m_rows = map_block(&e_v, normalize, true);

    let corner = [near_c_rows[0][3].x, near_c_rows[0][3].y];
    let alpha = core::f64::consts::TAU / n as f64;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let phase = k as f64 * alpha;
        let rot = [phase.cos(), phase.sin()];
        sectors.push(CcSectorBorder {
            near_c: CcBorderPiece {
                rows: rotate_block(&near_c_rows, rot),
            },
            near_m: CcBorderPiece {
                rows: rotate_block(&near_m_rows, rot),
            },
        });
    }

    Ok(CcCharacteristicTensorBorder {
        n,
        lambda: ring.lambda,
        sectors,
        corner,
    })
}

/// Tensor-border of a genuine Catmull-Clark ring on a node-wise
/// lifted characteristic net: nodes at the normalized planar
/// characteristic positions with `z = lift(x, y)` applied PER NODE,
/// subdivided once through the verified Stam machinery and extracted
/// with the same window/prolongation pipeline. Because the lift is
/// node-wise, the resulting ring is a real bicubic CC spline whose
/// height is NOT a polynomial function of the plane — a
/// production-grade, exactly C²-compatible curved fixture. The lift
/// must be rotation-symmetric (`lift(R·p) = lift(p)`) so sector 0's
/// pieces rotate to the other sectors.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn cc_lifted_ring_border_pieces(
    n: usize,
    lift: impl Fn(f64, f64) -> f64,
) -> Result<Vec<[[[Point3; 4]; 3]; 2]>, String> {
    use super::cc_characteristic_ring::{cc_extended_subdivision_rows, cc_picking_vectors};
    let ring = cc_characteristic_ring(n)?;

    // Normalization, exactly as the planar extractor derives it.
    let m_u = bspline_corner_value(&ring.grids[0]);
    let m_v = bspline_corner_value(&ring.grids[2]);
    let norm_sq = m_u[0] * m_u[0] + m_u[1] * m_u[1];
    let zeta = [m_u[0] / norm_sq, -m_u[1] / norm_sq];
    let conjugate = complex_mul(zeta, m_v)[1] < 0.0;
    let normalize = |value: [f64; 2]| -> [f64; 2] {
        let scaled = complex_mul(zeta, value);
        if conjugate {
            [scaled[0], -scaled[1]]
        } else {
            scaled
        }
    };

    // Node-wise lifted 3D net, one CC subdivision, picked grids.
    let net3: Vec<[f64; 3]> = ring
        .extended_net
        .iter()
        .map(|value| {
            let p = normalize(*value);
            [p[0], p[1], lift(p[0], p[1])]
        })
        .collect();
    let a_bar = cc_extended_subdivision_rows(n);
    let refined3: Vec<[f64; 3]> = a_bar
        .iter()
        .map(|row| {
            let mut acc = [0.0_f64; 3];
            for (w, v) in row.iter().zip(&net3) {
                acc[0] += w * v[0];
                acc[1] += w * v[1];
                acc[2] += w * v[2];
            }
            acc
        })
        .collect();
    let picks = cc_picking_vectors(n);
    let grid = |k: usize, r: usize, c: usize| -> [f64; 3] { refined3[picks[k][r * 4 + c] - 1] };

    // Per-edge 4x3 windows, along/depth B-to-BB, cap-side reflection —
    // the planar pipeline, run on three components.
    let block = |omega_one: bool, reverse_along: bool| -> [[Point3; 4]; 3] {
        let window = |a: usize, d: usize| -> [f64; 3] {
            if omega_one {
                grid(0, a, d)
            } else {
                grid(2, d, a)
            }
        };
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for axis in 0..3_usize {
            let mut along_bb = [[0.0_f64; 4]; 3];
            for (d, column) in along_bb.iter_mut().enumerate() {
                let p = [
                    window(0, d)[axis],
                    window(1, d)[axis],
                    window(2, d)[axis],
                    window(3, d)[axis],
                ];
                *column = [
                    (p[0] + 4.0 * p[1] + p[2]) / 6.0,
                    (2.0 * p[1] + p[2]) / 3.0,
                    (p[1] + 2.0 * p[2]) / 3.0,
                    (p[1] + 4.0 * p[2] + p[3]) / 6.0,
                ];
            }
            for (a, _) in (0..4_usize).enumerate() {
                let source_a = if reverse_along { 3 - a } else { a };
                let q0 =
                    (along_bb[0][source_a] + 4.0 * along_bb[1][source_a] + along_bb[2][source_a])
                        / 6.0;
                let q1 = (2.0 * along_bb[1][source_a] + along_bb[2][source_a]) / 3.0;
                let q2 = (along_bb[1][source_a] + 2.0 * along_bb[2][source_a]) / 3.0;
                let values = [q0, 2.0 * q0 - q1, 4.0 * q0 - 4.0 * q1 + q2];
                for (depth, value) in values.iter().enumerate() {
                    match axis {
                        0 => rows[depth][a].x = *value,
                        1 => rows[depth][a].y = *value,
                        _ => rows[depth][a].z = *value,
                    }
                }
            }
        }
        rows
    };
    let near_c_rows = block(true, false);
    let near_m_rows = block(false, true);

    let alpha = core::f64::consts::TAU / n as f64;
    let rotate3 = |rows: &[[Point3; 4]; 3], k: usize| -> [[Point3; 4]; 3] {
        let rot = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
        let mut out = *rows;
        for row in out.iter_mut() {
            for cp in row.iter_mut() {
                *cp = Point3::new(
                    cp.x * rot[0] - cp.y * rot[1],
                    cp.x * rot[1] + cp.y * rot[0],
                    cp.z,
                );
            }
        }
        out
    };
    Ok((0..n)
        .map(|k| [rotate3(&near_c_rows, k), rotate3(&near_m_rows, k)])
        .collect())
}

enum EdgeSelect {
    OmegaOne,
    OmegaThree,
}

/// Ring-side BB block of one inner-boundary edge: `[depth][along]`
/// complex values, depth 0 on the boundary, depth increasing INTO the
/// ring, along in the edge's natural grid direction.
fn ring_side_bb_block(ring: &CcCharacteristicRing, edge: EdgeSelect) -> [[[f64; 2]; 4]; 3] {
    // 4x3 B-spline window: along index a = 0..3, depth index d = 0..2.
    let window = |a: usize, d: usize| -> [f64; 2] {
        match edge {
            // Omega_1 grid[r][c]: r = t (along), c = s (depth).
            EdgeSelect::OmegaOne => ring.grids[0][a][d],
            // Omega_3 grid[r][c]: r = t (depth), c = s (along).
            EdgeSelect::OmegaThree => ring.grids[2][d][a],
        }
    };

    // Along-direction B-to-BB conversion per depth column (one full
    // uniform-cubic span: standard span formulas).
    let mut along_bb = [[[0.0_f64; 2]; 4]; 3];
    for (d, column) in along_bb.iter_mut().enumerate() {
        let p: [[f64; 2]; 4] = [window(0, d), window(1, d), window(2, d), window(3, d)];
        *column = [
            combine(&[(1.0 / 6.0, p[0]), (4.0 / 6.0, p[1]), (1.0 / 6.0, p[2])]),
            combine(&[(2.0 / 3.0, p[1]), (1.0 / 3.0, p[2])]),
            combine(&[(1.0 / 3.0, p[1]), (2.0 / 3.0, p[2])]),
            combine(&[(1.0 / 6.0, p[1]), (4.0 / 6.0, p[2]), (1.0 / 6.0, p[3])]),
        ];
    }

    // Depth-direction B-to-BB at the boundary span end: BB depth rows
    // 0..2 from B-spline depth columns 0..2 (the fourth column does
    // not influence depth <= 2).
    let mut out = [[[0.0_f64; 2]; 4]; 3];
    for a in 0..4 {
        let q0 = combine(&[
            (1.0 / 6.0, along_bb[0][a]),
            (4.0 / 6.0, along_bb[1][a]),
            (1.0 / 6.0, along_bb[2][a]),
        ]);
        let q1 = combine(&[(2.0 / 3.0, along_bb[1][a]), (1.0 / 3.0, along_bb[2][a])]);
        let q2 = combine(&[(1.0 / 3.0, along_bb[1][a]), (2.0 / 3.0, along_bb[2][a])]);
        out[0][a] = q0;
        out[1][a] = q1;
        out[2][a] = q2;
    }
    out
}

/// Convert a ring-side block to the cap-side chart: apply the complex
/// normalization, optionally reverse the along direction, and replace
/// the ring-side depth rows with the cap-side C² prolongation
/// (Bernstein normal reflection).
fn map_block(
    block: &[[[f64; 2]; 4]; 3],
    normalize: impl Fn([f64; 2]) -> [f64; 2],
    reverse_along: bool,
) -> [[Point3; 4]; 3] {
    let pick = |d: usize, a: usize| -> [f64; 2] {
        let a = if reverse_along { 3 - a } else { a };
        normalize(block[d][a])
    };
    let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
    for a in 0..4_usize {
        let q0 = pick(0, a);
        let q1 = pick(1, a);
        let q2 = pick(2, a);
        // Cap-side prolongation: r0 = q0, r1 = 2 q0 - q1,
        // r2 = 4 q0 - 4 q1 + q2.
        let r1 = [2.0 * q0[0] - q1[0], 2.0 * q0[1] - q1[1]];
        let r2 = [
            4.0 * q0[0] - 4.0 * q1[0] + q2[0],
            4.0 * q0[1] - 4.0 * q1[1] + q2[1],
        ];
        let [row0, row1, row2] = &mut rows;
        row0[a] = Point3::new(q0[0], q0[1], 0.0);
        row1[a] = Point3::new(r1[0], r1[1], 0.0);
        row2[a] = Point3::new(r2[0], r2[1], 0.0);
    }
    rows
}

fn rotate_block(rows: &[[Point3; 4]; 3], rot: [f64; 2]) -> [[Point3; 4]; 3] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
    for (out_row, row) in out.iter_mut().zip(rows) {
        for (out_cp, cp) in out_row.iter_mut().zip(row) {
            *out_cp = Point3::new(
                cp.x * rot[0] - cp.y * rot[1],
                cp.x * rot[1] + cp.y * rot[0],
                0.0,
            );
        }
    }
    out
}

/// Bicubic B-spline patch value at local corner `(0, 0)` of a picked
/// 4x4 grid: tensor product of the uniform-cubic endpoint mask
/// `[1/6, 2/3, 1/6]` over the 3x3 control points nearest the corner.
fn bspline_corner_value(grid: &[[[f64; 2]; 4]; 4]) -> [f64; 2] {
    let mask = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
    let mut acc = [0.0_f64; 2];
    for (r, wr) in mask.iter().enumerate() {
        for (c, wc) in mask.iter().enumerate() {
            acc[0] += wr * wc * grid[r][c][0];
            acc[1] += wr * wc * grid[r][c][1];
        }
    }
    acc
}

fn combine(terms: &[(f64, [f64; 2])]) -> [f64; 2] {
    let mut acc = [0.0_f64; 2];
    for (w, v) in terms {
        acc[0] += w * v[0];
        acc[1] += w * v[1];
    }
    acc
}

fn complex_mul(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

#[cfg(test)]
#[path = "cc_characteristic_tensor_border_tests.rs"]
mod tests;
