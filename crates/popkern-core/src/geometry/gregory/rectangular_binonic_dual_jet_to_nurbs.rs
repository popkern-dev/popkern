//! Exact rational-NURBS conversion for
//! [`RectangularGregoryBiNonicDualJet`].
//!
//! Denominator clearing over the activated corners' squared-weight
//! denominators `D_k = wu² + wv²` (bidegree (2, 2) each), via the
//! shared [`Bivar`] helpers — the same algorithm as the cubic and
//! quintic Gregory conversions with the denominator degree bumped.
//! With `A` activated corners the result is a bidegree
//! `(9 + 2A, 9 + 2A)` rational Bezier; the setback fan activates one
//! corner per corner piece, giving (11, 11).
//!
//! The conversion is polynomial-identity exact away from the
//! activated corners; *at* an activated corner the denominator
//! vanishes, so that grid entry is ε-regularized (weight `1e-14`,
//! control point pinned to the analytic corner position) — the
//! established convention of `rectangular_to_nurbs.rs`.

use super::rectangular_binonic_dual_jet::{
    DUAL_JET_CORNER_POSITIONS, DualJetCell, RectangularGregoryBiNonicDualJet, dual_jet_cell_of,
    dual_jet_cell_position,
};
use super::rectangular_to_nurbs::{Bivar, accumulate_scaled};
use crate::geometry::Point3;
use crate::geometry::gregory::RectangularGregoryBiNonic;
use crate::geometry::nurbs_surface::NurbsSurface;

/// `wu²` for a corner as a bidegree-(2, 2) Bernstein bivariate:
/// `u² = B_2^2(u)` (or `(1−u)² = B_0^2(u)`), constant in `v`.
fn wu_squared(corner: usize) -> Bivar {
    let (ic, _) = DUAL_JET_CORNER_POSITIONS[corner];
    let iu = if ic == 0 { 2 } else { 0 };
    let mut b = Bivar::zeros(2, 2);
    for j in 0..=2 {
        b.coefs[iu * 3 + j] = 1.0;
    }
    b
}

/// `wv²` for a corner, constant in `u`.
fn wv_squared(corner: usize) -> Bivar {
    let (_, jc) = DUAL_JET_CORNER_POSITIONS[corner];
    let jv = if jc == 0 { 2 } else { 0 };
    let mut b = Bivar::zeros(2, 2);
    for i in 0..=2 {
        b.coefs[i * 3 + jv] = 1.0;
    }
    b
}

fn add_bivar(a: &Bivar, b: &Bivar) -> Bivar {
    let mut out = a.clone();
    for (o, s) in out.coefs.iter_mut().zip(b.coefs.iter()) {
        *o += s;
    }
    out
}

fn bernstein_basis_99(i: usize, j: usize) -> Bivar {
    let mut b = Bivar::zeros(9, 9);
    b.coefs[i * 10 + j] = 1.0;
    b
}

impl RectangularGregoryBiNonicDualJet {
    /// Convert to an exact rational Bezier surface (as a clamped
    /// [`NurbsSurface`] on `[0, 1]²`). With no activated corner the
    /// patch is a plain bi-nonic and converts exactly at (9, 9).
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let active: Vec<usize> = (0..4).filter(|&k| self.corners[k].is_some()).collect();
        if active.is_empty() {
            return RectangularGregoryBiNonic { grid: self.grid }.to_nurbs_exact();
        }

        // D_k = wu² + wv² per activated corner; D = Π D_k.
        let denoms: Vec<Bivar> = active
            .iter()
            .map(|&k| add_bivar(&wu_squared(k), &wv_squared(k)))
            .collect();
        let mut d_full = denoms[0].clone();
        for d in denoms.iter().skip(1) {
            d_full = d_full.multiply(d);
        }
        let out_deg = 9 + 2 * active.len();

        let mut n_x = Bivar::zeros(out_deg, out_deg);
        let mut n_y = Bivar::zeros(out_deg, out_deg);
        let mut n_z = Bivar::zeros(out_deg, out_deg);

        // Plain cells: B_i^9 B_j^9 · P · D.
        for i in 0..10 {
            for j in 0..10 {
                let blended = dual_jet_cell_of(i, j)
                    .map(|(corner, _)| self.corners[corner].is_some())
                    .unwrap_or(false);
                if blended {
                    continue;
                }
                let p = self.grid[i][j];
                let term = bernstein_basis_99(i, j).multiply(&d_full);
                accumulate_scaled(&mut n_x, &term, p.x);
                accumulate_scaled(&mut n_y, &term, p.y);
                accumulate_scaled(&mut n_z, &term, p.z);
            }
        }

        // Blended cells: B_i^9 B_j^9 · (wv²·P_u + wu²·P_v) · (D / D_k).
        for (slot, &k) in active.iter().enumerate() {
            let mut d_over_dk: Option<Bivar> = None;
            for (other_slot, d) in denoms.iter().enumerate() {
                if other_slot == slot {
                    continue;
                }
                d_over_dk = Some(match d_over_dk {
                    Some(acc) => acc.multiply(d),
                    None => d.clone(),
                });
            }
            let block = match self.corners[k].as_ref() {
                Some(block) => block,
                None => continue,
            };
            let (t_wv, t_wu) = (wv_squared(k), wu_squared(k));
            for (cell, pair) in [
                (DualJetCell::Cell11, block.cell_11),
                (DualJetCell::Cell12, block.cell_12),
                (DualJetCell::Cell21, block.cell_21),
                (DualJetCell::Cell22, block.cell_22),
            ] {
                let (i, j) = dual_jet_cell_position(k, cell);
                let base = match &d_over_dk {
                    Some(rest) => bernstein_basis_99(i, j).multiply(rest),
                    None => bernstein_basis_99(i, j),
                };
                let (p_u, p_v) = pair;
                let term_u = t_wv.multiply(&base);
                let term_v = t_wu.multiply(&base);
                accumulate_scaled(&mut n_x, &term_u, p_u.x);
                accumulate_scaled(&mut n_y, &term_u, p_u.y);
                accumulate_scaled(&mut n_z, &term_u, p_u.z);
                accumulate_scaled(&mut n_x, &term_v, p_v.x);
                accumulate_scaled(&mut n_y, &term_v, p_v.y);
                accumulate_scaled(&mut n_z, &term_v, p_v.z);
            }
        }

        let d_elev = d_full.elevate_to(out_deg, out_deg);
        let count = out_deg + 1;
        let mut cps: Vec<Point3> = Vec::with_capacity(count * count);
        let mut weights: Vec<f64> = Vec::with_capacity(count * count);
        for i_u in 0..count {
            for i_v in 0..count {
                weights.push(d_elev.get(i_u, i_v));
                cps.push(Point3::new(
                    n_x.get(i_u, i_v),
                    n_y.get(i_u, i_v),
                    n_z.get(i_u, i_v),
                ));
            }
        }

        // ε-regularize every vanishing-weight slot. The squared
        // denominator vanishes QUADRATICALLY at an activated corner,
        // so after elevation a whole corner-adjacent neighbourhood of
        // Bernstein weight coefficients is exactly zero (measured:
        // slots (n,n), (n,n−1), (n−1,n), (n−1,n−1) at one active
        // corner) with matching zero numerators — a 0/0 slot, not a
        // pole. Each such slot takes weight ε with its control point
        // pinned to the analytic surface at the slot's Greville
        // parameters, so the O(ε) perturbation stays on the surface.
        let epsilon = 1e-14;
        let max_weight = weights.iter().fold(0.0_f64, |a, &w| a.max(w.abs()));
        let zero_bar = 1e-12 * max_weight.max(f64::MIN_POSITIVE);
        for i_u in 0..count {
            for i_v in 0..count {
                let idx = i_u * count + i_v;
                if weights[idx].abs() > zero_bar {
                    continue;
                }
                let p = self.point_at(i_u as f64 / out_deg as f64, i_v as f64 / out_deg as f64);
                weights[idx] = epsilon;
                cps[idx] = Point3::new(epsilon * p.x, epsilon * p.y, epsilon * p.z);
            }
        }

        let control_points: Vec<Point3> = cps
            .iter()
            .zip(weights.iter())
            .map(|(wcp, &w)| Point3::new(wcp.x / w, wcp.y / w, wcp.z / w))
            .collect();

        let mut knots = Vec::with_capacity(2 * count);
        knots.extend(std::iter::repeat_n(0.0, count));
        knots.extend(std::iter::repeat_n(1.0, count));

        NurbsSurface {
            degree_u: out_deg,
            degree_v: out_deg,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: count,
            num_v: count,
            weights,
        }
    }
}
