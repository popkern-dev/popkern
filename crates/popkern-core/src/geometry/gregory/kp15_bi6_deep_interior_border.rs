//! KP15 bi-6 tensor-border reparameterisation, script convention.
//!
//! Faithful port of the two border stages of the KP15 companion
//! Maple script (`docs/papers/n6.txt`, reproduced by the sympy
//! oracle `tools/reference/kp15-bi6-cap/n6_oracle.py`):
//!
//! - [`bspline_net_to_bi3_bb`] — the script's `BtoBB`: converts a
//!   `4 × 4` bicubic B-spline net into its bi-3 BB (Bézier) form.
//! - [`reparam_tensor_border_row`] — the script's `TBdata`: applies
//!   the `β` tensor-border reparameterisation, emitting the three
//!   depth rows (`j = 0, 1, 2`) of the reparameterised bi-6 border
//!   at the seven along-edge positions `i = 0..=6`.
//!
//! # Indexing convention (critical)
//!
//! All indices here are in the **script** convention: script
//! `p[i, j]` equals paper `p[6 − i, 6 − j]`, and script `[6, 6]` is
//! the extraordinary point. See
//! `tools/reference/kp15-bi6-cap/README.md`. The oracle's fixture
//! (`n6_cap_fixture.json`) is emitted in this same convention, so
//! this port is validated bit-close against it.

use super::super::kp_polar_shared::weighted_sum;
use crate::geometry::Point3;

/// Convert a `4 × 4` bicubic B-spline control net into its bi-3
/// BB-form coefficients `b[i][j]`, `i, j ∈ 0..=3`.
///
/// Port of the script's `BtoBB` (oracle `n6_oracle.py::BtoBB`).
/// The net is indexed `net[row][col]`; the script's accessor is
/// `p(i, j) = net[j][i]`.
pub fn bspline_net_to_bi3_bb(net: &[[Point3; 4]; 4]) -> [[Point3; 4]; 4] {
    let p = |i: usize, j: usize| net[j][i];
    let mut b = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];

    b[0][0] = weighted_sum(&[
        (16.0 / 36.0, p(1, 1)),
        (4.0 / 36.0, p(1, 0)),
        (4.0 / 36.0, p(2, 1)),
        (4.0 / 36.0, p(1, 2)),
        (4.0 / 36.0, p(0, 1)),
        (1.0 / 36.0, p(0, 0)),
        (1.0 / 36.0, p(2, 0)),
        (1.0 / 36.0, p(2, 2)),
        (1.0 / 36.0, p(0, 2)),
    ]);
    b[3][0] = weighted_sum(&[
        (16.0 / 36.0, p(2, 1)),
        (4.0 / 36.0, p(2, 0)),
        (4.0 / 36.0, p(3, 1)),
        (4.0 / 36.0, p(2, 2)),
        (4.0 / 36.0, p(1, 1)),
        (1.0 / 36.0, p(1, 0)),
        (1.0 / 36.0, p(3, 0)),
        (1.0 / 36.0, p(3, 2)),
        (1.0 / 36.0, p(1, 2)),
    ]);
    b[0][3] = weighted_sum(&[
        (16.0 / 36.0, p(1, 2)),
        (4.0 / 36.0, p(1, 1)),
        (4.0 / 36.0, p(2, 2)),
        (4.0 / 36.0, p(1, 3)),
        (4.0 / 36.0, p(0, 2)),
        (1.0 / 36.0, p(0, 1)),
        (1.0 / 36.0, p(2, 1)),
        (1.0 / 36.0, p(2, 3)),
        (1.0 / 36.0, p(0, 3)),
    ]);
    b[3][3] = weighted_sum(&[
        (16.0 / 36.0, p(2, 2)),
        (4.0 / 36.0, p(2, 1)),
        (4.0 / 36.0, p(3, 2)),
        (4.0 / 36.0, p(2, 3)),
        (4.0 / 36.0, p(1, 2)),
        (1.0 / 36.0, p(1, 1)),
        (1.0 / 36.0, p(3, 1)),
        (1.0 / 36.0, p(3, 3)),
        (1.0 / 36.0, p(1, 3)),
    ]);

    let (r49, r29, r19, r118) = (4.0 / 9.0, 2.0 / 9.0, 1.0 / 9.0, 1.0 / 18.0);
    b[1][0] = weighted_sum(&[
        (r49, p(1, 1)),
        (r29, p(2, 1)),
        (r19, p(1, 0)),
        (r19, p(1, 2)),
        (r118, p(2, 0)),
        (r118, p(2, 2)),
    ]);
    b[2][0] = weighted_sum(&[
        (r29, p(1, 1)),
        (r49, p(2, 1)),
        (r118, p(1, 0)),
        (r118, p(1, 2)),
        (r19, p(2, 0)),
        (r19, p(2, 2)),
    ]);
    b[1][3] = weighted_sum(&[
        (r49, p(1, 2)),
        (r29, p(2, 2)),
        (r19, p(1, 1)),
        (r19, p(1, 3)),
        (r118, p(2, 1)),
        (r118, p(2, 3)),
    ]);
    b[2][3] = weighted_sum(&[
        (r29, p(1, 2)),
        (r49, p(2, 2)),
        (r118, p(1, 1)),
        (r118, p(1, 3)),
        (r19, p(2, 1)),
        (r19, p(2, 3)),
    ]);
    b[0][1] = weighted_sum(&[
        (r49, p(1, 1)),
        (r29, p(1, 2)),
        (r19, p(0, 1)),
        (r19, p(2, 1)),
        (r118, p(0, 2)),
        (r118, p(2, 2)),
    ]);
    b[0][2] = weighted_sum(&[
        (r29, p(1, 1)),
        (r49, p(1, 2)),
        (r118, p(0, 1)),
        (r118, p(2, 1)),
        (r19, p(0, 2)),
        (r19, p(2, 2)),
    ]);
    b[3][1] = weighted_sum(&[
        (r49, p(2, 1)),
        (r29, p(2, 2)),
        (r19, p(1, 1)),
        (r19, p(3, 1)),
        (r118, p(1, 2)),
        (r118, p(3, 2)),
    ]);
    b[3][2] = weighted_sum(&[
        (r29, p(2, 1)),
        (r49, p(2, 2)),
        (r118, p(1, 1)),
        (r118, p(3, 1)),
        (r19, p(1, 2)),
        (r19, p(3, 2)),
    ]);
    b[1][1] = weighted_sum(&[
        (r49, p(1, 1)),
        (r29, p(2, 1)),
        (r29, p(1, 2)),
        (r19, p(2, 2)),
    ]);
    b[2][1] = weighted_sum(&[
        (r49, p(2, 1)),
        (r29, p(1, 1)),
        (r29, p(2, 2)),
        (r19, p(1, 2)),
    ]);
    b[1][2] = weighted_sum(&[
        (r49, p(1, 2)),
        (r29, p(1, 1)),
        (r29, p(2, 2)),
        (r19, p(2, 1)),
    ]);
    b[2][2] = weighted_sum(&[
        (r49, p(2, 2)),
        (r29, p(2, 1)),
        (r29, p(1, 2)),
        (r19, p(1, 1)),
    ]);

    b
}

/// Reparameterise a bi-3 tensor border into the three depth rows
/// (`j = 0, 1, 2`) of the bi-6 border, at along-edge positions
/// `i = 0..=6`.
///
/// Port of the script's `TBdata` (oracle `n6_oracle.py::TBdata`).
/// `p(i, j)` reads the bi-3 border coefficient at along-edge index
/// `i ∈ 0..=3` and depth `j ∈ 0..=2`. `c = cos(2π/n)`, `d2` and
/// `tau` are the [`super::kp15_bi6_reparam::Bi6Reparams`] shape
/// scalars. Returns `rows[j][i]`.
pub fn reparam_tensor_border_row(
    p: impl Fn(usize, usize) -> Point3,
    c: f64,
    d2: f64,
    tau: f64,
) -> [[Point3; 7]; 3] {
    let r = |num: f64, den: f64| num / den;
    let mut rows = [[Point3::new(0.0, 0.0, 0.0); 7]; 3];

    // Row 0 (position) — pure bi-3 → bi-6 degree raise of the edge.
    rows[0] = [
        p(0, 0),
        weighted_sum(&[(0.5, p(0, 0)), (0.5, p(1, 0))]),
        weighted_sum(&[(0.2, p(0, 0)), (0.2, p(2, 0)), (0.6, p(1, 0))]),
        weighted_sum(&[
            (r(1.0, 20.0), p(0, 0)),
            (r(9.0, 20.0), p(1, 0)),
            (r(9.0, 20.0), p(2, 0)),
            (r(1.0, 20.0), p(3, 0)),
        ]),
        weighted_sum(&[(0.2, p(3, 0)), (0.2, p(1, 0)), (0.6, p(2, 0))]),
        weighted_sum(&[(0.5, p(2, 0)), (0.5, p(3, 0))]),
        p(3, 0),
    ];

    // Row 1 (first cross-derivative).
    rows[1] = [
        weighted_sum(&[(0.5, p(0, 0)), (0.5, p(0, 1))]),
        weighted_sum(&[
            (0.25, p(0, 0)),
            (0.25, p(1, 0)),
            (0.25, p(0, 1)),
            (0.25, p(1, 1)),
        ]),
        weighted_sum(&[
            (r(1.0, 10.0), p(0, 0)),
            (r(1.0, 10.0), p(0, 1)),
            (r(3.0, 10.0), p(1, 0)),
            (r(1.0, 10.0), p(2, 0)),
            (r(3.0, 10.0), p(1, 1)),
            (r(1.0, 10.0), p(2, 1)),
        ]),
        weighted_sum(&[
            (r(1.0, 40.0), p(0, 0)),
            (r(9.0, 40.0), p(1, 0)),
            (r(9.0, 40.0), p(2, 0)),
            (r(1.0, 40.0), p(3, 0)),
            (r(1.0, 40.0), p(0, 1)),
            (r(9.0, 40.0), p(1, 1)),
            (r(9.0, 40.0), p(2, 1)),
            (r(1.0, 40.0), p(3, 1)),
        ]),
        weighted_sum(&[
            (r(1.0, 10.0), p(1, 0)),
            (r(3.0, 10.0), p(2, 0)),
            (r(1.0, 10.0), p(3, 0)),
            (r(1.0, 10.0), p(1, 1)),
            (r(3.0, 10.0), p(2, 1)),
            (r(1.0, 10.0), p(3, 1)),
        ]),
        weighted_sum(&[
            (0.25, p(2, 0)),
            (0.25, p(3, 0)),
            (0.25, p(2, 1)),
            (0.25, p(3, 1)),
        ]),
        weighted_sum(&[(0.5, p(3, 0)), (0.5, p(3, 1))]),
    ];

    // Row 2 (second cross-derivative) — reshaped by `β` via d2/c/τ.
    rows[2][0] = weighted_sum(&[(0.6, p(0, 1)), (0.2, p(0, 0)), (0.2, p(0, 2))]);
    rows[2][1] = weighted_sum(&[
        (r(1.0, 10.0), p(1, 2)),
        (r(1.0, 10.0), p(0, 0)),
        (r(1.0, 10.0), p(1, 0)),
        (r(3.0, 10.0), p(0, 1)),
        (r(3.0, 10.0), p(1, 1)),
        (r(1.0, 10.0), p(0, 2)),
    ]);
    rows[2][2] = weighted_sum(&[
        (r(3.0, 25.0), p(1, 2)),
        (r(1.0, 25.0), p(2, 2)),
        (r(1.0, 25.0), p(0, 0)),
        (r(3.0, 25.0), p(1, 0)),
        (r(1.0, 25.0), p(2, 0)),
        (r(3.0, 25.0), p(0, 1)),
        (r(9.0, 25.0), p(1, 1)),
        (r(3.0, 25.0), p(2, 1)),
        (r(1.0, 25.0), p(0, 2)),
        (d2 / 50.0, p(0, 1)),
        (-d2 / 25.0, p(0, 0)),
        (d2 / 50.0, p(1, 0)),
    ]);
    rows[2][3] = weighted_sum(&[
        ((18.0 * d2 + 54.0) / 600.0, p(2, 0)),
        ((3.0 * d2 + 2.0 * c + 18.0) / 600.0, p(0, 1)),
        (
            (-12.0 * c + 54.0 - 42.0 * d2) / 600.0 + 12.0 * c / (600.0 * tau),
            p(1, 0),
        ),
        (
            (10.0 * c - 6.0 * d2 + 6.0) / 600.0 - 12.0 * c / (600.0 * tau),
            p(0, 0),
        ),
        ((27.0 * d2 + 162.0) / 600.0, p(1, 1)),
        (54.0 / 600.0, p(1, 2)),
        (162.0 / 600.0, p(2, 1)),
        (6.0 / 600.0, p(3, 0)),
        (54.0 / 600.0, p(2, 2)),
        (6.0 / 600.0, p(3, 2)),
        (18.0 / 600.0, p(3, 1)),
        (6.0 / 600.0, p(0, 2)),
    ]);
    rows[2][4] = weighted_sum(&[
        (6.0 / 150.0, p(1, 2)),
        (18.0 / 150.0, p(2, 2)),
        (6.0 / 150.0, p(3, 2)),
        ((6.0 + 6.0 * c - 5.0 * d2 - 8.0 * c / tau) / 150.0, p(1, 0)),
        (
            (18.0 - 10.0 * d2 - 8.0 * c + 8.0 * c / tau) / 150.0,
            p(2, 0),
        ),
        ((6.0 + 3.0 * d2) / 150.0, p(3, 0)),
        ((18.0 + 2.0 * c + 3.0 * d2) / 150.0, p(1, 1)),
        ((54.0 + 9.0 * d2) / 150.0, p(2, 1)),
        (18.0 / 150.0, p(3, 1)),
    ]);
    rows[2][5] = weighted_sum(&[
        (6.0 / 60.0, p(2, 2)),
        (6.0 / 60.0, p(3, 2)),
        ((6.0 - 4.0 * d2 + 2.0 * c - 4.0 * c / tau) / 60.0, p(2, 0)),
        ((6.0 - 2.0 * d2 - 4.0 * c + 4.0 * c / tau) / 60.0, p(3, 0)),
        ((18.0 + 2.0 * c + 3.0 * d2) / 60.0, p(2, 1)),
        ((18.0 + 3.0 * d2) / 60.0, p(3, 1)),
    ]);
    rows[2][6] = weighted_sum(&[
        (0.2, p(3, 2)),
        (0.2, p(3, 0)),
        (0.6, p(3, 1)),
        (c / 15.0 + d2 / 10.0, p(3, 1)),
        (-c / 15.0 - d2 / 10.0, p(3, 0)),
    ]);

    rows
}

#[cfg(test)]
#[path = "kp15_bi6_deep_interior_border_tests.rs"]
mod tests;
