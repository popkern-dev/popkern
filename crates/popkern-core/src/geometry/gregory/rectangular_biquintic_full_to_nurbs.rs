//! Exact bidegree-(9, 9) rational-NURBS conversion for
//! [`RectangularGregoryBiquinticFull`].
//!
//! Mirrors the deg5 conversion in `rectangular_deg5_to_nurbs.rs`
//! with one structural difference: in deg5, ALL 16 interior
//! Bernstein terms `(i, j) ∈ {1..=4}²` share a quadrant's
//! rational inner-CP pair; in biquintic-full only 4 of the 16
//! (the corner-adjacent positions `(1,1), (1,4), (4,1), (4,4)`)
//! are rational-blend driven. The other 12 free interior CPs are
//! plain Bezier, handled identically to boundary CPs in the
//! polynomial-times-denominator accumulation.
//!
//! Reuses the `Bivar` arithmetic helpers from
//! [`super::rectangular_to_nurbs`].
//!
//! # Algorithm
//!
//! Denominator `D(u, v) = D_0 · D_1 · D_2 · D_3` (bidegree `(4, 4)`),
//! same as cubic/deg5 forms. Clearing the corner denominators
//! gives a bidegree-`(5 + 4, 5 + 4) = (9, 9)` rational Bezier
//! with 100 control points on a clamped `[0, 1]²` NURBS surface.
//! ε-regularization at the four corners reproduces the analytic
//! Gregory corner position exactly.
//!
//! The "free interior" positions (12 of them) contribute as
//! `B_i^5 B_j^5 · P_{ij} · D`, identical to boundary positions.
//! The 4 corner-adjacent positions contribute as
//! `B_i^5 B_j^5 · (w_v P_{k,u} + w_u P_{k,v}) · (D / D_k)`,
//! identical to the deg5 form.

use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
use super::rectangular_to_nurbs::{
    Bivar, accumulate_scaled, bernstein_basis_55, denom_bernstein, linear_bernstein,
};
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

impl RectangularGregoryBiquinticFull {
    /// Convert to an exact bidegree-(9, 9) rational Bezier
    /// surface expressed as a clamped [`NurbsSurface`] on
    /// `[0, 1]²`.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let denoms = denom_bernstein();
        let d_full = denoms[0]
            .multiply(&denoms[1])
            .multiply(&denoms[2])
            .multiply(&denoms[3]);

        let mut n_x = Bivar::zeros(9, 9);
        let mut n_y = Bivar::zeros(9, 9);
        let mut n_z = Bivar::zeros(9, 9);

        for i in 0..6 {
            for j in 0..6 {
                let is_corner_adjacent = matches!((i, j), (1, 1) | (1, 4) | (4, 1) | (4, 4));
                if is_corner_adjacent {
                    continue;
                }
                let p = self.grid[i][j];
                let term = bernstein_basis_55(i, j).multiply(&d_full);
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

        // Per-corner (wv, wu) Bernstein-(1, 1) coefficient tables.
        // Corner 0 at (1,1): wv = v,     wu = u
        // Corner 1 at (1,4): wv = 1 - v, wu = u
        // Corner 2 at (4,4): wv = 1 - v, wu = 1 - u
        // Corner 3 at (4,1): wv = v,     wu = 1 - u
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

        // The 4 corner-adjacent Bernstein indices → corner k.
        let corner_adjacent = [
            ((1usize, 1usize), 0usize),
            ((1, 4), 1),
            ((4, 4), 2),
            ((4, 1), 3),
        ];
        for ((i, j), k) in corner_adjacent {
            let (p_u, p_v) = self.corner_inner_cps[k];
            let base = bernstein_basis_55(i, j).multiply(&d_over_dk[k]);
            let term_u = wv_bern[k].multiply(&base);
            let term_v = wu_bern[k].multiply(&base);
            accumulate_scaled(&mut n_x, &term_u, p_u.x);
            accumulate_scaled(&mut n_y, &term_u, p_u.y);
            accumulate_scaled(&mut n_z, &term_u, p_u.z);
            accumulate_scaled(&mut n_x, &term_v, p_v.x);
            accumulate_scaled(&mut n_y, &term_v, p_v.y);
            accumulate_scaled(&mut n_z, &term_v, p_v.z);
        }

        let d_elev = d_full.elevate_to(9, 9);

        let mut cps: Vec<Point3> = Vec::with_capacity(100);
        let mut weights: Vec<f64> = Vec::with_capacity(100);
        for i_u in 0..10 {
            for i_v in 0..10 {
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
            (0, 9, self.grid[0][5]),
            (9, 0, self.grid[5][0]),
            (9, 9, self.grid[5][5]),
        ];
        for (i_u, i_v, p) in corners {
            let idx = i_u * 10 + i_v;
            weights[idx] = epsilon;
            cps[idx] = Point3::new(epsilon * p.x, epsilon * p.y, epsilon * p.z);
        }

        let control_points: Vec<Point3> = cps
            .iter()
            .zip(weights.iter())
            .map(|(wcp, &w)| Point3::new(wcp.x / w, wcp.y / w, wcp.z / w))
            .collect();

        let mut knots = Vec::with_capacity(20);
        knots.extend(std::iter::repeat_n(0.0, 10));
        knots.extend(std::iter::repeat_n(1.0, 10));

        NurbsSurface {
            degree_u: 9,
            degree_v: 9,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: 10,
            num_v: 10,
            weights,
        }
    }
}
