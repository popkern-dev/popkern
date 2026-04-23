//! Phase IV PnS — G¹ cap-cap continuity constraints (eq 2, eq 2³).
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
) -> Bi4G1Row {
    assert!(n >= 3, "PnS bi-4 G¹: n must be ≥ 3, got {n}");
    let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();
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

    Bi4G1Row {
        cur10,
        cur20,
        cur30,
        cur40,
        nxt21,
    }
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
    let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();
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
