//! One sector-separating piece pair of the KP20 bi-9 resampling: the
//! previous sector's piece `p` (seam on its `v = 0` side, along `u`)
//! and the current sector's `p̃` (seam on its `u = 0` side, along `v`),
//! both parameterised outward from `o`.
//!
//! Unknowns per coordinate (20): the shared value row `V`, and each
//! side's first and second differences `F`, `S` (middles 3–6). With
//! `∂_v p = 9 F_p`, `∂²_vv p = 72 S_p`, `∂_u p = V'`:
//!
//! ```text
//! eq (8): 9 F_c + 9 F_p − b V'                          = 0
//! eq (9): 72 S_c − 72 S_p + 18 b F_p' − b² V'' − b b' V' = 0
//! ```
//!
//! collocated at 16 points (eq (9) has degree ≤ 11). The fit targets are
//! the target nets' band jets on each side (minimal movement in jet
//! space).

// Coefficient nets are addressed by (i, j) index pairs throughout.
#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::{MIDDLE, Net, Side, b9, samples, slot};
use crate::geometry::variational_cap::linalg;
use crate::geometry::{Point3, Vec3};

const UNKNOWNS: usize = 20;
/// Offsets of the five unknown blocks.
const V: usize = 0;
const FP: usize = 4;
const SP: usize = 8;
const FC: usize = 12;
const SC: usize = 16;

fn v3(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}

/// Jet coefficient `(V, F, S)` of a net's band at `along`.
fn jet_at(net: &Net, side: Side, along: usize) -> [Vec3; 3] {
    let r = |d: usize| {
        let (i, j) = slot(side, d, along);
        v3(net[i][j])
    };
    let (a, b, c) = (r(0), r(1), r(2));
    [a, b - a, c - b * 2.0 + a]
}

fn quad(coeffs: [f64; 3], u: f64) -> (f64, f64) {
    let (a, b, c) = (coeffs[0], coeffs[1], coeffs[2]);
    let value = a * (1.0 - u) * (1.0 - u) + 2.0 * b * u * (1.0 - u) + c * u * u;
    let slope = 2.0 * ((b - a) * (1.0 - u) + (c - b) * u);
    (value, slope)
}

/// Rows of eq (8)/(9) at `u`: coefficients over the 20 unknowns, and
/// the held part (per coordinate) moved to the right-hand side.
fn relation_rows(prev: &Net, curr: &Net, b: [f64; 3], u: f64) -> [([f64; UNKNOWNS], Vec3); 2] {
    let (bu, bp) = quad(b, u);
    let mut eq8 = [0.0; UNKNOWNS];
    let mut eq9 = [0.0; UNKNOWNS];
    let (mut h8, mut h9) = (Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
    for along in 0..10 {
        let [w, w1, w2] = b9(along, u);
        let c8 = [(-bu * w1, V), (9.0 * w, FP), (9.0 * w, FC)];
        let c9 = [
            (-bu * bu * w2 - bu * bp * w1, V),
            (18.0 * bu * w1, FP),
            (-72.0 * w, SP),
            (72.0 * w, SC),
        ];
        match MIDDLE.iter().position(|&m| m == along) {
            Some(k) => {
                for (c, off) in c8 {
                    eq8[off + k] += c;
                }
                for (c, off) in c9 {
                    eq9[off + k] += c;
                }
            }
            None => {
                let [vp, fp, sp] = jet_at(prev, Side::J0, along);
                let [_, fc, sc] = jet_at(curr, Side::I0, along);
                h8 = h8 + vp * (-bu * w1) + fp * (9.0 * w) + fc * (9.0 * w);
                h9 = h9
                    + vp * (-bu * bu * w2 - bu * bp * w1)
                    + fp * (18.0 * bu * w1)
                    + sp * (-72.0 * w)
                    + sc * (72.0 * w);
            }
        }
    }
    [(eq8, h8 * -1.0), (eq9, h9 * -1.0)]
}

/// The worst eq (8)/(9) residual and C⁰ gap along the pair's seam, from
/// the full nets.
pub(crate) fn sector_seam_residual(prev: &Net, curr: &Net, b: [f64; 3]) -> [f64; 3] {
    let mut worst = [0.0f64; 3];
    for k in 0..=40 {
        let u = k as f64 / 40.0;
        let (bu, bp) = quad(b, u);
        let (mut vp, mut vp1, mut vp2, mut fp, mut fp1, mut sp, mut fc, mut sc, mut vc) = (
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
        );
        for along in 0..10 {
            let [w, w1, w2] = b9(along, u);
            let [a, f, s] = jet_at(prev, Side::J0, along);
            let [c0, g, t] = jet_at(curr, Side::I0, along);
            vp = vp + a * w;
            vp1 = vp1 + a * w1;
            vp2 = vp2 + a * w2;
            fp = fp + f * w;
            fp1 = fp1 + f * w1;
            sp = sp + s * w;
            fc = fc + g * w;
            sc = sc + t * w;
            vc = vc + c0 * w;
        }
        let eq8 = fc * 9.0 + fp * 9.0 - vp1 * bu;
        let eq9 = sc * 72.0 - sp * 72.0 + fp1 * (18.0 * bu) - vp2 * (bu * bu) - vp1 * (bu * bp);
        worst[0] = worst[0].max((vp - vc).length());
        worst[1] = worst[1].max(eq8.length());
        worst[2] = worst[2].max(eq9.length());
    }
    worst
}

/// Re-fit the pair's band middles: least squares to the target nets'
/// band jets over the null space of the collocated eq (8)/(9) rows. Returns `None`
/// when the rows are inconsistent at the start (the start's own middles
/// must satisfy them) or a solve fails.
pub(crate) fn fit_sector_seam(
    prev: &mut Net,
    curr: &mut Net,
    b: [f64; 3],
    target_prev: &Net,
    target_curr: &Net,
) -> Option<()> {
    // Constraint rows C x = d (same C for every coordinate).
    let colloc = samples(16);
    let mut c_rows: Vec<[f64; UNKNOWNS]> = Vec::new();
    let mut d_rows: Vec<Vec3> = Vec::new();
    for &u in &colloc {
        for (row, rhs) in relation_rows(prev, curr, b, u) {
            c_rows.push(row);
            d_rows.push(rhs);
        }
    }
    let m = c_rows.len();
    let mut c = Mat::<f64>::zeros(m, UNKNOWNS);
    for (r, row) in c_rows.iter().enumerate() {
        for (k, w) in row.iter().enumerate() {
            c.write(r, k, *w);
        }
    }
    let svd = linalg::svd(&c)?;
    let smax = svd.s.first().copied().unwrap_or(0.0);
    let rank = svd.s.iter().filter(|&&s| s > 1e-10 * smax).count();
    // Particular solution per coordinate: V_r Σ⁻¹ U_rᵀ d.
    let mut xp = [[0.0f64; UNKNOWNS]; 3];
    for axis in 0..3 {
        let d: Vec<f64> = d_rows.iter().map(|v| [v.x, v.y, v.z][axis]).collect();
        for k in 0..rank {
            let proj: f64 = (0..m).map(|r| svd.u.read(r, k) * d[r]).sum::<f64>() / svd.s[k];
            for x in 0..UNKNOWNS {
                xp[axis][x] += svd.v.read(x, k) * proj;
            }
        }
        // Consistency: the rows must be satisfiable (the start satisfies them).
        let resid: f64 = (0..m)
            .map(|r| {
                ((0..UNKNOWNS)
                    .map(|x| c.read(r, x) * xp[axis][x])
                    .sum::<f64>()
                    - d[r])
                    .abs()
            })
            .fold(0.0, f64::max);
        let scale = d.iter().fold(1.0f64, |a, v| a.max(v.abs()));
        if resid > 1e-9 * scale {
            return None;
        }
    }
    let dim = UNKNOWNS - rank;
    // Objective rows: V ≈ value, F/S ≈ exact jets per side.
    let fit = samples(24);
    let mut a_rows: Vec<[f64; UNKNOWNS]> = Vec::new();
    let mut b_rows: Vec<Vec3> = Vec::new();
    // Band jets of the targets along the seam (value, first and second
    // difference), evaluated as degree-9 curves.
    let band = |net: &Net, side: Side, u: f64| {
        let mut out = [Vec3::new(0.0, 0.0, 0.0); 3];
        for along in 0..10 {
            let w = b9(along, u)[0];
            let jet = jet_at(net, side, along);
            for q in 0..3 {
                out[q] = out[q] + jet[q] * w;
            }
        }
        out
    };
    for &u in &fit {
        let (jp, jc) = (
            band(target_prev, Side::J0, u),
            band(target_curr, Side::I0, u),
        );
        let goals = [
            (V, (jp[0] + jc[0]) * 0.5, Side::J0, 0usize, &*prev),
            (FP, jp[1], Side::J0, 1, &*prev),
            (SP, jp[2], Side::J0, 2, &*prev),
            (FC, jc[1], Side::I0, 1, &*curr),
            (SC, jc[2], Side::I0, 2, &*curr),
        ];
        for (off, goal, side, q, net) in goals {
            let mut row = [0.0; UNKNOWNS];
            let mut known = Vec3::new(0.0, 0.0, 0.0);
            for along in 0..10 {
                let w = b9(along, u)[0];
                match MIDDLE.iter().position(|&mm| mm == along) {
                    Some(k) => row[off + k] = w,
                    None => known = known + jet_at(net, side, along)[q] * w,
                }
            }
            a_rows.push(row);
            b_rows.push(goal - known);
        }
    }
    let mut solution = xp;
    if dim > 0 {
        let rows = a_rows.len();
        let mut az = Mat::<f64>::zeros(rows, dim);
        for (r, row) in a_rows.iter().enumerate() {
            for k in 0..dim {
                let col = rank + k;
                az.write(
                    r,
                    k,
                    (0..UNKNOWNS).map(|x| row[x] * svd.v.read(x, col)).sum(),
                );
            }
        }
        let mut rhs = Mat::<f64>::zeros(rows, 3);
        for (r, (row, goal)) in a_rows.iter().zip(&b_rows).enumerate() {
            for axis in 0..3 {
                let ax: f64 = (0..UNKNOWNS).map(|x| row[x] * xp[axis][x]).sum();
                rhs.write(r, axis, [goal.x, goal.y, goal.z][axis] - ax);
            }
        }
        let y = linalg::lstsq(az, rhs, 1e-12)?;
        for axis in 0..3 {
            for k in 0..dim {
                let col = rank + k;
                for x in 0..UNKNOWNS {
                    solution[axis][x] += svd.v.read(x, col) * y.read(k, axis);
                }
            }
        }
    }
    let at = |off: usize, k: usize| {
        Vec3::new(
            solution[0][off + k],
            solution[1][off + k],
            solution[2][off + k],
        )
    };
    for (k, &along) in MIDDLE.iter().enumerate() {
        let (v, fp, sp, fc, sc) = (at(V, k), at(FP, k), at(SP, k), at(FC, k), at(SC, k));
        let writes = [
            (Side::J0, 0, v, true),
            (Side::J0, 1, v + fp, true),
            (Side::J0, 2, v + fp * 2.0 + sp, true),
            (Side::I0, 0, v, false),
            (Side::I0, 1, v + fc, false),
            (Side::I0, 2, v + fc * 2.0 + sc, false),
        ];
        for (side, depth, value, on_prev) in writes {
            let (i, j) = slot(side, depth, along);
            let net = if on_prev { &mut *prev } else { &mut *curr };
            net[i][j] = Point3::new(value.x, value.y, value.z);
        }
    }
    Some(())
}
