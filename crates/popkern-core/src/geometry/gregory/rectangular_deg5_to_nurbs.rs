//! Exact bidegree-(9, 9) rational-NURBS conversion for
//! [`RectangularGregoryDeg5`]. Mirrors the cubic conversion in
//! `rectangular_to_nurbs.rs` with bidegree bumped from (5, 5) ×
//! (4, 4) = (9, 9) and 16 interior Bernstein terms per sub-patch
//! instead of 4.
//!
//! Shares the `Bivar` arithmetic helpers with the cubic conversion
//! via the `rectangular_to_nurbs` module.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::rectangular_deg5::RectangularGregoryDeg5;
use super::rectangular_to_nurbs::{
    Bivar, accumulate_scaled, bernstein_basis_55, denom_bernstein, linear_bernstein,
};
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

impl RectangularGregoryDeg5 {
    /// Convert the quintic Gregory patch to an exact bidegree-(9, 9)
    /// rational Bezier surface, represented as a clamped
    /// [`NurbsSurface`] on `[0, 1]²`.
    ///
    /// Uses the same denominator-clearing algorithm as the cubic case
    /// — 4 corner denominators `D_0..D_3` cleared to give a bidegree-
    /// `(5 + 4, 5 + 4) = (9, 9)` rational Bezier. ε-regularized at
    /// the four corners. The 4-inner-CP-pair structure (one pair per
    /// corner-quadrant) means all 16 interior bidegree-5 Bernstein
    /// terms share their quadrant's pair.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let denoms = denom_bernstein();
        let d_full = denoms[0]
            .multiply(&denoms[1])
            .multiply(&denoms[2])
            .multiply(&denoms[3]);

        let mut n_x = Bivar::zeros(9, 9);
        let mut n_y = Bivar::zeros(9, 9);
        let mut n_z = Bivar::zeros(9, 9);

        // Non-inner contribution: for each non-inner `(i, j)` in the
        // 6×6 grid, accumulate `B_i^5 B_j^5 · P_{ij} · D`. Inner
        // positions are `(i, j) ∈ {1, 2, 3, 4} × {1, 2, 3, 4}`.
        for i in 0..6 {
            for j in 0..6 {
                let is_inner = (1..=4).contains(&i) && (1..=4).contains(&j);
                if is_inner {
                    continue;
                }
                let p = self.boundary_cps[i][j];
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

        // Per-quadrant `(wv, wu)` tables matching
        // `RectangularGregoryDeg5::quadrant_weights`.
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

        for i in 1..=4 {
            for j in 1..=4 {
                let k = match (i, j) {
                    (1..=2, 1..=2) => 0,
                    (1..=2, 3..=4) => 1,
                    (3..=4, 3..=4) => 2,
                    (3..=4, 1..=2) => 3,
                    _ => unreachable!(),
                };
                let (p_u, p_v) = self.inner_cps[k];
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
            (0usize, 0usize, self.boundary_cps[0][0]),
            (0, 9, self.boundary_cps[0][5]),
            (9, 0, self.boundary_cps[5][0]),
            (9, 9, self.boundary_cps[5][5]),
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
