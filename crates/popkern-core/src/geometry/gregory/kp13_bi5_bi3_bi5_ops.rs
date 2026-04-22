//! Bi-3 ↔ bi-5 BB-grid operations for the KP13 cap assembly
//! pipeline: 1D / 2D degree raise and bi-5 subdivision via
//! de Casteljau.
//!
//! These are classical primitives, collected here so the I.14.b
//! outer-layer transfer and related cap-assembly modules stay
//! focused on the KP13 geometry.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use crate::geometry::Point3;

/// 1-D degree raise from degree 3 to degree 4.
///
/// For a univariate BB polynomial `p(t) = Σ_{k=0..3} C_k · B^3_k(t)`,
/// the same polynomial in degree-4 BB form has coefficients
/// `D_k = (k / 4) · C_{k-1} + ((4 - k) / 4) · C_k` (with out-of-
/// range `C` terms treated as zero).
pub fn degree_raise_1d_3_to_4(c: &[Point3; 4]) -> [Point3; 5] {
    let mut d = [Point3::new(0.0, 0.0, 0.0); 5];
    for k in 0..5 {
        let a = k as f64 / 4.0;
        let b = (4 - k) as f64 / 4.0;
        let prev = if k > 0 {
            c[k - 1]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        let cur = if k < 4 {
            c[k]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        d[k] = Point3::new(
            a * prev.x + b * cur.x,
            a * prev.y + b * cur.y,
            a * prev.z + b * cur.z,
        );
    }
    d
}

/// 1-D degree raise from degree 4 to degree 5.
pub fn degree_raise_1d_4_to_5(c: &[Point3; 5]) -> [Point3; 6] {
    let mut d = [Point3::new(0.0, 0.0, 0.0); 6];
    for k in 0..6 {
        let a = k as f64 / 5.0;
        let b = (5 - k) as f64 / 5.0;
        let prev = if k > 0 {
            c[k - 1]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        let cur = if k < 5 {
            c[k]
        } else {
            Point3::new(0.0, 0.0, 0.0)
        };
        d[k] = Point3::new(
            a * prev.x + b * cur.x,
            a * prev.y + b * cur.y,
            a * prev.z + b * cur.z,
        );
    }
    d
}

/// 1-D degree raise 3 → 5 (two steps).
#[inline]
pub fn degree_raise_1d_3_to_5(c: &[Point3; 4]) -> [Point3; 6] {
    let d4 = degree_raise_1d_3_to_4(c);
    degree_raise_1d_4_to_5(&d4)
}

/// 2-D tensor-product degree raise from bi-3 (4 × 4) to bi-5
/// (6 × 6). Applies the 1-D raise along each row, then along each
/// column of the intermediate grid.
pub fn bi3_to_bi5_grid(bi3: &[[Point3; 4]; 4]) -> [[Point3; 6]; 6] {
    // Step 1: raise each fixed-v row (indexed by j) in the u
    // direction → intermediate 6 × 4 grid.
    let mut inter = [[Point3::new(0.0, 0.0, 0.0); 4]; 6];
    for j in 0..4 {
        let row = [bi3[0][j], bi3[1][j], bi3[2][j], bi3[3][j]];
        let raised = degree_raise_1d_3_to_5(&row);
        for i in 0..6 {
            inter[i][j] = raised[i];
        }
    }
    // Step 2: raise each fixed-u column in the v direction →
    // final 6 × 6 grid.
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for i in 0..6 {
        let col = [inter[i][0], inter[i][1], inter[i][2], inter[i][3]];
        let raised = degree_raise_1d_3_to_5(&col);
        out[i].copy_from_slice(&raised);
    }
    out
}

/// 1-D de Casteljau subdivision of a bi-5 BB curve at `t = 1/2`.
///
/// Returns `(left, right)` where the `left` CPs describe the
/// curve on `t ∈ [0, 1/2]` re-parameterised to `[0, 1]`, and
/// `right` describes `t ∈ [1/2, 1]`.
pub fn decasteljau_bi5_at_half(c: &[Point3; 6]) -> ([Point3; 6], [Point3; 6]) {
    // level[k] holds the degree-(5-k) sub-curves' CPs; we reduce
    // iteratively.
    let mut level: [[Point3; 6]; 6] = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    level[0].copy_from_slice(c);
    for k in 1..6 {
        for i in 0..(6 - k) {
            level[k][i] = midpoint(level[k - 1][i], level[k - 1][i + 1]);
        }
    }
    // Left sub-curve CPs: level[k][0] for k = 0..=5.
    // Right sub-curve CPs: level[5 - k][k] for k = 0..=5 (reading
    // along the "diagonal" of the de Casteljau triangle from
    // the final value to the original endpoint).
    let mut left = [Point3::new(0.0, 0.0, 0.0); 6];
    let mut right = [Point3::new(0.0, 0.0, 0.0); 6];
    for k in 0..6 {
        left[k] = level[k][0];
        right[k] = level[5 - k][k];
    }
    (left, right)
}

#[inline]
fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}

/// Subdivide a bi-5 6 × 6 CP grid at the midpoint `(u, v) =
/// (1/2, 1/2)` into four 6 × 6 sub-patches, arranged in
/// `out[iu][jv]` order matching the KP13 cap layout:
///
/// - `out[0][0]` — `(u, v) ∈ [0, 1/2]² → O` region.
/// - `out[0][1]` — `(u, v) ∈ [0, 1/2] × [1/2, 1] → L` region.
/// - `out[1][0]` — `(u, v) ∈ [1/2, 1] × [0, 1/2] → R` region.
/// - `out[1][1]` — `(u, v) ∈ [1/2, 1]² → M` region.
///
/// Each sub-patch is re-parameterised to its own `[0, 1]²`.
pub fn subdivide_bi5_2x2(cps: &[[Point3; 6]; 6]) -> [[[[Point3; 6]; 6]; 2]; 2] {
    // Step 1: for each fixed-v row j, split in u → two sub-rows
    // (left and right in u).
    let mut u_split: [[[Point3; 6]; 6]; 2] = [[[Point3::new(0.0, 0.0, 0.0); 6]; 6]; 2];
    for j in 0..6 {
        let row = [
            cps[0][j], cps[1][j], cps[2][j], cps[3][j], cps[4][j], cps[5][j],
        ];
        let (left, right) = decasteljau_bi5_at_half(&row);
        for i in 0..6 {
            u_split[0][i][j] = left[i];
            u_split[1][i][j] = right[i];
        }
    }

    // Step 2: for each u-sub-grid, split in v → the four final
    // sub-patches.
    let mut out: [[[[Point3; 6]; 6]; 2]; 2] = [[[[Point3::new(0.0, 0.0, 0.0); 6]; 6]; 2]; 2];
    for iu in 0..2 {
        for i in 0..6 {
            let col = [
                u_split[iu][i][0],
                u_split[iu][i][1],
                u_split[iu][i][2],
                u_split[iu][i][3],
                u_split[iu][i][4],
                u_split[iu][i][5],
            ];
            let (bottom, top) = decasteljau_bi5_at_half(&col);
            out[iu][0][i].copy_from_slice(&bottom);
            out[iu][1][i].copy_from_slice(&top);
        }
    }
    out
}

#[cfg(test)]
#[path = "kp13_bi5_bi3_bi5_ops_tests.rs"]
mod tests;
