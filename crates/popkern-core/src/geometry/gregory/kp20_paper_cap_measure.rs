//! Cap-level gate measurements for the KP20 paper assembly.

use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_local_constraint::{kp20_eq8_jacobian_dense, kp20_eq9_jacobian_dense};
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use super::kp20_paper_cap::Kp20PaperCapReport;
use crate::geometry::Point3;

const P1: Kp20PieceId = Kp20PieceId::NearOriginLeft;
const P2: Kp20PieceId = Kp20PieceId::NearOriginRight;
const P4: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

/// Evaluate the cap-level gates on the assembled sectors.
pub(crate) fn measure_paper_cap(
    n: usize,
    sectors: &[Kp20MacroPatch],
    placements_per_sector: &[Vec<(Kp20PieceId, usize, usize, Point3)>],
    guide_fit_residual: f64,
    max_strip_movement: f64,
) -> Kp20PaperCapReport {
    let c = (core::f64::consts::TAU / n as f64).cos();
    let mut max_c0 = 0.0_f64;
    let mut max_g1 = 0.0_f64;
    let mut max_g2 = 0.0_f64;

    for curve in 0..n {
        let prev = (curve + n - 1) % n;
        let curr = curve;
        // C0 along the two curve halves.
        for i in 0..=6_usize {
            let near_o = delta(
                sectors[prev].piece(P1).cps[i][0],
                sectors[curr].piece(P1).cps[0][i],
            );
            let near_m = delta(
                sectors[prev].piece(P2).cps[i][0],
                sectors[curr].piece(P4).cps[0][i],
            );
            max_c0 = max_c0.max(near_o).max(near_m);
        }
        // eq (8)/(9) residual BB coefficients per half pair.
        for (prev_piece, curr_piece, b) in [
            (P1, P1, kp20_b_coeffs_near_o(c)),
            (P2, P4, kp20_b_coeffs_near_m(c)),
        ] {
            for (row_index, jac) in kp20_eq8_jacobian_dense(b)
                .into_iter()
                .chain(kp20_eq9_jacobian_dense(b))
                .enumerate()
            {
                let mut acc = [0.0_f64; 3];
                for side in 0..2_usize {
                    for i in 0..7_usize {
                        for j in 0..7_usize {
                            let w = jac[side * 49 + i * 7 + j];
                            if w == 0.0 {
                                continue;
                            }
                            let cp = if side == 0 {
                                sectors[prev].piece(prev_piece).cps[i][j]
                            } else {
                                sectors[curr].piece(curr_piece).cps[j][i]
                            };
                            acc[0] += w * cp.x;
                            acc[1] += w * cp.y;
                            acc[2] += w * cp.z;
                        }
                    }
                }
                let magnitude = acc[0].abs().max(acc[1].abs()).max(acc[2].abs());
                if row_index < 7 {
                    max_g1 = max_g1.max(magnitude);
                } else {
                    max_g2 = max_g2.max(magnitude);
                }
            }
        }
    }

    // Gray rows must equal the eq (14) bands (step iv must not touch
    // them).
    let mut max_border = 0.0_f64;
    for (s, placements) in placements_per_sector.iter().enumerate() {
        for &(piece, i, j, value) in placements {
            max_border = max_border.max(delta(sectors[s].piece(piece).cps[i][j], value));
        }
    }

    // Lemma-4 apex consistency: sector s's spoke-(s+1) edge cubic run
    // equals sector s+1's spoke-(s+1) edge run (orders 0..=3 as edge
    // BB coefficients).
    let mut lemma4 = 0.0_f64;
    for s in 0..n {
        let next = (s + 1) % n;
        for i in 0..=3_usize {
            lemma4 = lemma4.max(delta(
                sectors[s].piece(P1).cps[i][0],
                sectors[next].piece(P1).cps[0][i],
            ));
        }
    }

    // Surface regularity: minimum |du x dv| over corner-region
    // samples of every piece.
    let mut min_regularity = f64::INFINITY;
    for sector in sectors {
        for piece_id in Kp20PieceId::all() {
            let piece = sector.piece(piece_id);
            for &(u, v) in &[
                (0.0, 0.0),
                (0.5, 0.0),
                (1.0, 0.0),
                (0.0, 0.5),
                (0.5, 0.5),
                (1.0, 0.5),
                (0.0, 1.0),
                (0.5, 1.0),
                (1.0, 1.0),
            ] {
                let (du, dv) = bi6_first_partials(piece, u, v);
                let cross = [
                    du[1] * dv[2] - du[2] * dv[1],
                    du[2] * dv[0] - du[0] * dv[2],
                    du[0] * dv[1] - du[1] * dv[0],
                ];
                let norm = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
                min_regularity = min_regularity.min(norm);
            }
        }
    }

    Kp20PaperCapReport {
        n,
        sectors: sectors.to_vec(),
        guide_fit_residual,
        max_c0_residual: max_c0,
        max_g1_residual: max_g1,
        max_g2_residual: max_g2,
        max_border_residual: max_border,
        lemma4_delta: lemma4,
        max_strip_movement,
        input_corner_mismatch: 0.0,
        min_regularity,
    }
}

fn delta(a: Point3, b: Point3) -> f64 {
    (a.x - b.x)
        .abs()
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

fn bernstein6(t: f64) -> ([f64; 7], [f64; 7]) {
    let mut b = [0.0_f64; 7];
    let mut db = [0.0_f64; 7];
    let s = 1.0 - t;
    let choose = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    for (k, coeff) in choose.iter().enumerate() {
        b[k] = coeff * t.powi(k as i32) * s.powi((6 - k) as i32);
        let dk = k as f64;
        db[k] = coeff
            * (dk * t.powi((k as i32 - 1).max(0)) * s.powi((6 - k) as i32)
                - (6.0 - dk) * t.powi(k as i32) * s.powi((5 - k as i32).max(0)));
    }
    (b, db)
}

fn bi6_first_partials(
    piece: &super::kp20_bi6_cap::Kp20Bi6Patch,
    u: f64,
    v: f64,
) -> ([f64; 3], [f64; 3]) {
    let (bu, dbu) = bernstein6(u);
    let (bv, dbv) = bernstein6(v);
    let mut du = [0.0_f64; 3];
    let mut dv = [0.0_f64; 3];
    for i in 0..7 {
        for j in 0..7 {
            let cp = piece.cps[i][j];
            du[0] += dbu[i] * bv[j] * cp.x;
            du[1] += dbu[i] * bv[j] * cp.y;
            du[2] += dbu[i] * bv[j] * cp.z;
            dv[0] += bu[i] * dbv[j] * cp.x;
            dv[1] += bu[i] * dbv[j] * cp.y;
            dv[2] += bu[i] * dbv[j] * cp.z;
        }
    }
    (du, dv)
}
