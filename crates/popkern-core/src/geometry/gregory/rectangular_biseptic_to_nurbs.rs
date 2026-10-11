//! Exact bidegree-(11, 11) rational-NURBS conversion for
//! [`RectangularGregoryBiseptic`] — the degree-7 analogue of
//! [`super::rectangular_biquintic_full_to_nurbs`].
//!
//! Structure is identical to the biquintic-full converter (the
//! 4-corner-pair rational-blend denominator clearing) with
//! bidegree raised from `(5+4, 5+4) = (9, 9)` to
//! `(7+4, 7+4) = (11, 11)`, producing a 12×12 = 144-CP clamped
//! NURBS on `[0, 1]²` that matches `point_at` to floating-point
//! precision. ε-regularised at the four corners.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::rectangular_biseptic::RectangularGregoryBiseptic;
use super::rectangular_to_nurbs::{
    Bivar, accumulate_scaled, bernstein_basis_77, denom_bernstein, linear_bernstein,
};
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

impl RectangularGregoryBiseptic {
    /// Convert to an exact bidegree-(11, 11) rational Bezier
    /// surface expressed as a clamped [`NurbsSurface`] on
    /// `[0, 1]²`.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let denoms = denom_bernstein();
        let d_full = denoms[0]
            .multiply(&denoms[1])
            .multiply(&denoms[2])
            .multiply(&denoms[3]);

        let mut n_x = Bivar::zeros(11, 11);
        let mut n_y = Bivar::zeros(11, 11);
        let mut n_z = Bivar::zeros(11, 11);

        // Non-corner-adjacent positions (28 boundary + 32 free
        // interior = 60) contribute as plain `B_i^7 B_j^7 · P · D_full`.
        for i in 0..8 {
            for j in 0..8 {
                let is_corner_adjacent = matches!((i, j), (1, 1) | (1, 6) | (6, 1) | (6, 6));
                if is_corner_adjacent {
                    continue;
                }
                let p = self.grid[i][j];
                let term = bernstein_basis_77(i, j).multiply(&d_full);
                accumulate_scaled(&mut n_x, &term, p.x);
                accumulate_scaled(&mut n_y, &term, p.y);
                accumulate_scaled(&mut n_z, &term, p.z);
            }
        }

        let d_over_dk: [Bivar; 4] = [
            denoms[1].multiply(&denoms[2]).multiply(&denoms[3]),
            denoms[0].multiply(&denoms[2]).multiply(&denoms[3]),
            denoms[0].multiply(&denoms[1]).multiply(&denoms[3]),
            denoms[0].multiply(&denoms[1]).multiply(&denoms[2]),
        ];

        // Per-corner (wv, wu) Bernstein-(1, 1) tables. Same as
        // biquintic-full; the corner layout (0 at (0,0), 1 at
        // (0, last), 2 at (last, last), 3 at (last, 0)) is the
        // same modulo the `last` index.
        let wv_bern: [Bivar; 4] = [
            linear_bernstein(0.0, 1.0, 0.0, 1.0),
            linear_bernstein(1.0, 0.0, 1.0, 0.0),
            linear_bernstein(1.0, 0.0, 1.0, 0.0),
            linear_bernstein(0.0, 1.0, 0.0, 1.0),
        ];
        let wu_bern: [Bivar; 4] = [
            linear_bernstein(0.0, 0.0, 1.0, 1.0),
            linear_bernstein(0.0, 0.0, 1.0, 1.0),
            linear_bernstein(1.0, 1.0, 0.0, 0.0),
            linear_bernstein(1.0, 1.0, 0.0, 0.0),
        ];

        let corner_adjacent = [
            ((1usize, 1usize), 0usize),
            ((1, 6), 1),
            ((6, 6), 2),
            ((6, 1), 3),
        ];
        for ((i, j), k) in corner_adjacent {
            let (p_u, p_v) = self.corner_inner_cps[k];
            let base = bernstein_basis_77(i, j).multiply(&d_over_dk[k]);
            let term_u = wv_bern[k].multiply(&base);
            let term_v = wu_bern[k].multiply(&base);
            accumulate_scaled(&mut n_x, &term_u, p_u.x);
            accumulate_scaled(&mut n_y, &term_u, p_u.y);
            accumulate_scaled(&mut n_z, &term_u, p_u.z);
            accumulate_scaled(&mut n_x, &term_v, p_v.x);
            accumulate_scaled(&mut n_y, &term_v, p_v.y);
            accumulate_scaled(&mut n_z, &term_v, p_v.z);
        }

        let d_elev = d_full.elevate_to(11, 11);

        let mut cps: Vec<Point3> = Vec::with_capacity(144);
        let mut weights: Vec<f64> = Vec::with_capacity(144);
        for i_u in 0..12 {
            for i_v in 0..12 {
                weights.push(d_elev.get(i_u, i_v));
                cps.push(Point3::new(
                    n_x.get(i_u, i_v),
                    n_y.get(i_u, i_v),
                    n_z.get(i_u, i_v),
                ));
            }
        }

        let epsilon = 1e-14;
        let corners = [
            (0usize, 0usize, self.grid[0][0]),
            (0, 11, self.grid[0][7]),
            (11, 0, self.grid[7][0]),
            (11, 11, self.grid[7][7]),
        ];
        for (i_u, i_v, p) in corners {
            let idx = i_u * 12 + i_v;
            weights[idx] = epsilon;
            cps[idx] = Point3::new(epsilon * p.x, epsilon * p.y, epsilon * p.z);
        }

        let control_points: Vec<Point3> = cps
            .iter()
            .zip(weights.iter())
            .map(|(wcp, &w)| Point3::new(wcp.x / w, wcp.y / w, wcp.z / w))
            .collect();

        let mut knots = Vec::with_capacity(24);
        knots.extend(std::iter::repeat_n(0.0, 12));
        knots.extend(std::iter::repeat_n(1.0, 12));

        NurbsSurface {
            degree_u: 11,
            degree_v: 11,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: 12,
            num_v: 12,
            weights,
        }
    }
}
