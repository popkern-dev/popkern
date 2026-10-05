//! PnS cap — G¹ cap-cap continuity constraints (eq 2, eq 2³).
//!
//! Between adjacent sectors `cur = p^k` and `nxt = p^{k+1}` sharing
//! the boundary `p(u, 0)`, the G¹ join equation (eq 1) is
//! ```text
//!   nxt_v(u, 0) + cur_v(u, 0) = 2c(1 − u) · cur_u(u, 0),
//!   c := cos(2π/n)
//! ```
//!
//! For **bi-4** patches (degree 4 in u along the boundary), this
//! unfolds into eq (2)'s five explicit assignments that pin down
//! the inner-most row of `cur` (the j = 0 row except the end
//! corners, which are shared) from the input tensor-border /
//! β-reparameterized data + the central point. See `14bi234.pdf`
//! Fig. 8c for labelling.
//!
//! For **bi-3** patches (degree 3, used when n ∈ {3, 5}), eq (2³)
//! replaces (2) with `cur_{20}` only; the other j=0 coefs come
//! directly from b ∘ β since there are fewer of them. The n=3/n=5
//! bi-3 paths ALSO require a separate (vertex-enclosure) system
//! for the inner twist CPs `cur_{11}, nxt_{11}`, given by eq (7) —
//! that's a separate solver (coming in a later commit).
//!
//! All formulas here are **closed-form and local** (per-pair of
//! adjacent sectors); no global solve needed.

#![allow(dead_code)]

use crate::geometry::Point3;
use popkern_math::PortableMath;

/// Compute the G¹ -enforced j = 0 row of a bi-4 sector from its
/// j = 1 row and the central point (eq 2 of `14bi234.pdf`).
///
/// Inputs:
/// - `cur00`: the shared central point (all N patches meet here).
/// - `cur_row1`: current sector's j=1 layer `[cur_{01}, cur_{11},
///   cur_{21}, cur_{31}, cur_{41}]` (5 Point3).
/// - `nxt_row1`: next sector's j=1 layer `[nxt_{01}, nxt_{11},
///   nxt_{21}, nxt_{31}, nxt_{41}]` (5 Point3).
/// - `n`: valence.
///
/// Outputs a 5-tuple `(cur10, cur20, nxt21, cur30, cur40)` — the
/// CPs that eq (2) determines. Note that `nxt_{21}` belongs to the
/// *next* sector's j=1 row (in-place update), and the caller must
/// write it there.
///
/// Eq (2):
/// ```text
///   cur_{10} := (2(c−1)·cur_{00} + cur_{01} + nxt_{01}) / (2c)
///   cur_{20} := [4c(cur_{11} + nxt_{11})
///                + (3c − 4)(cur_{01} + nxt_{01} + 2(c − 1)·cur_{00})] / (6c²)
///   nxt_{21} := −cur_{21} + 2·cur_{20} + c·(cur_{30} − cur_{20})
///   cur_{30} := [4(cur_{31} + nxt_{31}) − c·(cur_{41} + nxt_{41})] / (2(4 − c))
///   cur_{40} := (cur_{41} + nxt_{41}) / 2
/// ```
///
/// Evaluation order: cur10, cur20, cur30, cur40, nxt21 (each uses
/// only previously-computed values).
pub fn bi4_g1_row_from_row1(
    cur00: Point3,
    cur_row1: [Point3; 5],
    nxt_row1: [Point3; 5],
    n: usize,
) -> Result<Bi4G1Row, Bi4G1SolveError> {
    if n < 3 {
        return Err(Bi4G1SolveError::ValenceTooSmall(n));
    }
    // At n = 4, c = cos(π/2) = 0 and the *solved* presentation of
    // eq (2) below divides by 2c and 6c². That is a property of the
    // presentation, not of the geometry: the underlying G¹ relations
    // are division-free (see [`bi4_g1_residuals`]) and hold at n = 4
    // to machine zero. What genuinely changes at n = 4 is the rank —
    // the system stops *determining* `cur_{10}` and `cur_{20}` and
    // instead *constrains* them, leaving them as free DOF fixed by
    // the central a∘σ expansion. So there is nothing to solve here,
    // and returning a value would be a fabrication.
    if n == 4 {
        return Err(Bi4G1SolveError::UnderdeterminedAtN4);
    }
    let c = (2.0 * std::f64::consts::PI / (n as f64)).pcos();
    let [cur01, cur11, cur21, cur31, cur41] = cur_row1;
    let [nxt01, nxt11, _nxt21, nxt31, nxt41] = nxt_row1;

    // cur_{10}
    let two_c_minus_1 = 2.0 * (c - 1.0);
    let cur10 = scale_add(
        &[(two_c_minus_1, cur00), (1.0, cur01), (1.0, nxt01)],
        1.0 / (2.0 * c),
    );

    // cur_{20}
    let three_c_minus_4 = 3.0 * c - 4.0;
    let inner = scale_add(&[(1.0, cur01), (1.0, nxt01), (2.0 * (c - 1.0), cur00)], 1.0);
    let cur20 = scale_add(
        &[(4.0 * c, cur11), (4.0 * c, nxt11), (three_c_minus_4, inner)],
        1.0 / (6.0 * c * c),
    );

    // cur_{30}
    let cur30 = scale_add(
        &[(4.0, cur31), (4.0, nxt31), (-c, cur41), (-c, nxt41)],
        1.0 / (2.0 * (4.0 - c)),
    );

    // cur_{40}
    let cur40 = Point3::new(
        0.5 * (cur41.x + nxt41.x),
        0.5 * (cur41.y + nxt41.y),
        0.5 * (cur41.z + nxt41.z),
    );

    // nxt_{21}: depends on cur21 (input), cur20 + cur30 (just
    // computed). Eq (2):
    //   nxt_{21} := -cur_{21} + 2·cur_{20} + c·(cur_{30} − cur_{20})
    //            = -cur_{21} + (2 - c)·cur_{20} + c·cur_{30}
    let nxt21 = scale_add(&[(-1.0, cur21), (2.0 - c, cur20), (c, cur30)], 1.0);

    Ok(Bi4G1Row {
        cur10,
        cur20,
        cur30,
        cur40,
        nxt21,
    })
}

/// Why the solved form of eq (2) could not produce a j = 0 row.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum Bi4G1SolveError {
    /// The sector construction needs at least a triangle.
    #[error("PnS bi-4 G¹: n must be ≥ 3, got {0}")]
    ValenceTooSmall(usize),
    /// `n = 4` makes `c = 0`, and eq (2)'s solved form divides by
    /// `2c` / `6c²`. There the G¹ system does not determine
    /// `cur_{10}` / `cur_{20}`; it constrains them. Verify such caps
    /// with [`bi4_g1_residuals`] instead, which is division-free and
    /// valid at every valence.
    #[error(
        "PnS bi-4 G¹: at n=4 (c=0) eq (2) constrains rather than determines the j=0 row; use bi4_g1_residuals to verify instead of solving"
    )]
    UnderdeterminedAtN4,
}

/// Residuals of the five G¹ relations of eq (2), in **division-free**
/// form — valid at every valence including `n = 4`.
///
/// # Why this exists
///
/// Eq (2) is printed as five assignments that *solve* for the j = 0
/// row, and that solve divides by `2c` and `6c²` with
/// `c = cos(2π/n)`. At `n = 4`, `c = 0` and the printed form blows
/// up — which this kernel long read as "the bi-4 cap does not support
/// n = 4". It does. Clearing the denominators (and substituting the
/// first relation into the second to drop one power of `c`) gives an
/// equivalent system with no division at all:
///
/// ```text
///   E1:  cur_{01} + nxt_{01}   = 2c·cur_{10} + 2(1−c)·cur_{00}
///   E2:  2(cur_{11} + nxt_{11}) = 3c·cur_{20} + (4−3c)·cur_{10}
///   E3:  cur_{21} + nxt_{21}   = (2−c)·cur_{20} + c·cur_{30}
///   E4:  2(4−c)·cur_{30}       = 4(cur_{31} + nxt_{31}) − c(cur_{41} + nxt_{41})
///   E5:  2·cur_{40}            = cur_{41} + nxt_{41}
/// ```
///
/// Measured against caps built by [`super::pns_cap_assemble`], every
/// residual sits at machine zero (≤ 1.4e-15) for
/// `n ∈ {4, 5, 6, 8}` — n = 4 included, and no worse there than
/// anywhere else. This is the honest G¹ gate for the bi-4 cap;
/// [`bi4_g1_row_from_row1`] remains the *constructor* for valences
/// where the system is determined.
///
/// Returns the five residual vectors in E1..E5 order. Each is a
/// difference of `Point3` combinations, so its magnitude carries the
/// same length units as the input control points — callers compare it
/// against their own explicit tolerance.
pub fn bi4_g1_residuals(
    cur00: Point3,
    cur_row0: [Point3; 5],
    cur_row1: [Point3; 5],
    nxt_row1: [Point3; 5],
    n: usize,
) -> [[f64; 3]; 5] {
    let c = (2.0 * std::f64::consts::PI / (n as f64)).pcos();
    let [_cur00_dup, cur10, cur20, cur30, cur40] = cur_row0;
    let [cur01, cur11, cur21, cur31, cur41] = cur_row1;
    let [nxt01, nxt11, nxt21, nxt31, nxt41] = nxt_row1;

    let comb = |terms: &[(f64, Point3)]| -> [f64; 3] {
        let mut acc = [0.0_f64; 3];
        for (s, p) in terms {
            acc[0] += s * p.x;
            acc[1] += s * p.y;
            acc[2] += s * p.z;
        }
        acc
    };

    [
        comb(&[
            (1.0, cur01),
            (1.0, nxt01),
            (-2.0 * c, cur10),
            (-2.0 * (1.0 - c), cur00),
        ]),
        comb(&[
            (2.0, cur11),
            (2.0, nxt11),
            (-3.0 * c, cur20),
            (-(4.0 - 3.0 * c), cur10),
        ]),
        comb(&[(1.0, cur21), (1.0, nxt21), (-(2.0 - c), cur20), (-c, cur30)]),
        comb(&[
            (2.0 * (4.0 - c), cur30),
            (-4.0, cur31),
            (-4.0, nxt31),
            (c, cur41),
            (c, nxt41),
        ]),
        comb(&[(2.0, cur40), (-1.0, cur41), (-1.0, nxt41)]),
    ]
}

/// Largest Euclidean magnitude among the five [`bi4_g1_residuals`]
/// vectors — the single number a G¹ gate asserts on.
pub fn bi4_g1_worst_residual(
    cur00: Point3,
    cur_row0: [Point3; 5],
    cur_row1: [Point3; 5],
    nxt_row1: [Point3; 5],
    n: usize,
) -> f64 {
    bi4_g1_residuals(cur00, cur_row0, cur_row1, nxt_row1, n)
        .iter()
        .map(|r| (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt())
        .fold(0.0_f64, f64::max)
}

/// Result of [`bi4_g1_row_from_row1`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bi4G1Row {
    /// Current sector's `p_{10}`.
    pub cur10: Point3,
    /// Current sector's `p_{20}`.
    pub cur20: Point3,
    /// Current sector's `p_{30}`.
    pub cur30: Point3,
    /// Current sector's `p_{40}`.
    pub cur40: Point3,
    /// Next sector's `p_{21}` (consistency constraint; the caller
    /// writes this back into the next sector's row 1).
    pub nxt21: Point3,
}

/// Compute the G¹ -enforced `p_{20}` of a bi-3 sector from its
/// j = 1 row and the central point (eq 2³ of `14bi234.pdf` §4.2).
///
/// Only `cur_{20}` is derived here — the bi-3 case requires a
/// separate `n × n` vertex-enclosure system (eq 7) to determine
/// the twist CPs `cur_{11}, nxt_{11}`; that's a solver, shipped
/// separately.
///
/// Eq (2³):
/// ```text
///   cur_{20} := [3(cur_{21} + nxt_{21}) − c(cur_{31} + nxt_{31})] / (2(3 − c))
/// ```
pub fn bi3_g1_row20(
    cur21: Point3,
    cur31: Point3,
    nxt21: Point3,
    nxt31: Point3,
    n: usize,
) -> Point3 {
    assert!(
        n == 3 || n == 5,
        "PnS bi-3 G¹: only n ∈ {{3, 5}} supported, got {n}"
    );
    let c = (2.0 * std::f64::consts::PI / (n as f64)).pcos();
    scale_add(
        &[(3.0, cur21), (3.0, nxt21), (-c, cur31), (-c, nxt31)],
        1.0 / (2.0 * (3.0 - c)),
    )
}

/// Linear combination `scale * Σ w_i · p_i` over a list of
/// `(weight, point)` pairs.
fn scale_add(pairs: &[(f64, Point3)], scale: f64) -> Point3 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (w, p) in pairs {
        x += w * p.x;
        y += w * p.y;
        z += w * p.z;
    }
    Point3::new(scale * x, scale * y, scale * z)
}

#[cfg(test)]
#[path = "pns_g1_constraints_tests.rs"]
mod tests;
