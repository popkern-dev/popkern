//! Measurement side of the KP20 paper τ solve: affine seed, sampled
//! Jacobian summary, and the printed eq (21)/(22) cubic-expansion
//! comparison. Companion to `kp20_paper_tau`.

use super::kp20_bi6_cap::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_tau::Kp20TauMap;
use super::kp20_bi6_tau_solve_assemble::Kp20PlanarTauPaperSolveError;
use crate::geometry::Point3;

const P1: Kp20PieceId = Kp20PieceId::NearOriginLeft;
const P2: Kp20PieceId = Kp20PieceId::NearOriginRight;
const P3: Kp20PieceId = Kp20PieceId::NearMidpointRight;
const P4: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

// ---- Seed, orientation, jet gate ------------------------------------------

/// Affine seed: the identity-like map placing each quadrant at its
/// domain position in the global frame.
pub(crate) fn kp20_paper_tau_affine_seed(n: usize) -> Kp20MacroPatch {
    let alpha = core::f64::consts::TAU / n as f64;
    let spoke0 = [1.0, 0.0];
    let spoke1 = [alpha.cos(), alpha.sin()];
    let mut patch = Kp20MacroPatch::zero();
    let offsets: [(Kp20PieceId, f64, f64); 4] = [
        (P1, 0.0, 0.0),
        (P2, 0.5, 0.0),
        (P3, 0.5, 0.5),
        (P4, 0.0, 0.5),
    ];
    for (piece_id, a0, b0) in offsets {
        let piece = patch.piece_mut(piece_id);
        for i in 0..=6_usize {
            for j in 0..=6_usize {
                let a = a0 + 0.5 * i as f64 / 6.0;
                let b = b0 + 0.5 * j as f64 / 6.0;
                piece.cps[i][j] = Point3::new(
                    a * spoke1[0] + b * spoke0[0],
                    a * spoke1[1] + b * spoke0[1],
                    0.0,
                );
            }
        }
    }
    patch
}

/// Sampled Jacobian summary (negatives, near-zeros, minimum) over all
/// four pieces on a 5×5 grid, with all pieces sharing the global
/// domain orientation.
pub(crate) fn kp20_paper_tau_orientation_summary(tau: &Kp20MacroPatch) -> (usize, usize, f64) {
    let samples = [0.0, 0.25, 0.5, 0.75, 1.0];
    let mut negatives = 0;
    let mut near_zero = 0;
    let mut min_det = f64::INFINITY;
    for piece_id in Kp20PieceId::all() {
        let piece = tau.piece(piece_id);
        for &u in &samples {
            for &v in &samples {
                let (du, dv) = bi6_first_partials(piece, u, v);
                // The chart's (a, b) = (spoke₁, spoke₀) frame is
                // clockwise, so a correctly-oriented map has NEGATIVE
                // raw determinant; report the chart-signed value.
                let det = -(du[0] * dv[1] - du[1] * dv[0]);
                min_det = min_det.min(det);
                if det < -1.0e-10 {
                    negatives += 1;
                } else if det.abs() <= 1.0e-10 {
                    near_zero += 1;
                }
            }
        }
    }
    (negatives, near_zero, min_det)
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

fn bi6_first_partials(piece: &Kp20Bi6Patch, u: f64, v: f64) -> ([f64; 2], [f64; 2]) {
    let (bu, dbu) = bernstein6(u);
    let (bv, dbv) = bernstein6(v);
    let mut du = [0.0_f64; 2];
    let mut dv = [0.0_f64; 2];
    for i in 0..7 {
        for j in 0..7 {
            let cp = piece.cps[i][j];
            du[0] += dbu[i] * bv[j] * cp.x;
            du[1] += dbu[i] * bv[j] * cp.y;
            dv[0] += bu[i] * dbv[j] * cp.x;
            dv[1] += bu[i] * dbv[j] * cp.y;
        }
    }
    (du, dv)
}

/// Corner partials `∂^a_u ∂^b_v P(0,0)` (a, b ≤ 3) of one bi-6 piece
/// in its own parameters, from BB differences.
fn bi6_corner_partials(piece: &Kp20Bi6Patch) -> [[[f64; 2]; 4]; 4] {
    let falling = |order: usize| -> f64 { (0..order).map(|k| (6 - k) as f64).product() };
    let binom = |order: usize, k: usize| -> f64 {
        let mut value = 1.0;
        for step in 0..k {
            value = value * (order - step) as f64 / (step + 1) as f64;
        }
        value
    };
    let mut out = [[[0.0_f64; 2]; 4]; 4];
    for (a, out_a) in out.iter_mut().enumerate() {
        for (b, slot) in out_a.iter_mut().enumerate() {
            let mut acc = [0.0_f64; 2];
            for p in 0..=a {
                for q in 0..=b {
                    let sign = if (a - p + b - q) % 2 == 0 { 1.0 } else { -1.0 };
                    let w = sign * binom(a, p) * binom(b, q);
                    let cp = piece.cps[p][q];
                    acc[0] += w * cp.x;
                    acc[1] += w * cp.y;
                }
            }
            let scale = falling(a) * falling(b);
            slot[0] = acc[0] * scale;
            slot[1] = acc[1] * scale;
        }
    }
    out
}

/// Worst per-order deviation between the solved τ's cubic expansion at
/// `o` and the printed eq (21)/(22) table, compared at the
/// triangle-BB coefficient level.
///
/// The solved corner jet (sector-scale: one macro piece spans half the
/// `o → m` parameter) is converted to the ten total-degree-3 BB
/// coefficients of `Kp20TauMap`'s chart, and compared against the
/// printed table. KP20 eq (21) OVERWRITES τ₀₁/τ₀₂ (the presentation
/// gauge of the boundary run), so those two coefficients (and their
/// transposes) are excluded; the µ-carrying coefficients τ₀₃, τ₃₀,
/// τ₁₁, τ₁₂, τ₂₁ are the published cross-check. Deltas are reported
/// per total order: `delta[3]` bounds |µ₀|, |µ₂|, |µ₃| deviations and
/// `delta[2]` the τ₁₁ (µ₁ = 1/3) deviation; τ₀₀ sits in `delta[0]`.
pub(crate) fn kp20_paper_tau_jet_delta(
    n: usize,
    tau: &Kp20MacroPatch,
) -> Result<[f64; 4], Kp20PlanarTauPaperSolveError> {
    use super::direct_linalg::solve_dense;

    let printed =
        Kp20TauMap::for_sector(n, 0).map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let partials = bi6_corner_partials(tau.piece(P1));

    // Sector-scale jet as a 10-vector in the (a, b) partial order used
    // below.
    let index_pairs: [(usize, usize); 10] = [
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
        (1, 2),
        (2, 0),
        (2, 1),
        (3, 0),
    ];

    // Basis matrix: jet entries of each unit triangle-BB coefficient.
    let mut matrix = vec![vec![0.0_f64; 10]; 10];
    for (col, &(ci, cj)) in index_pairs.iter().enumerate() {
        let mut unit = Kp20TauMap::from_cps([[0.0; 2]; 10]);
        unit.cps[super::kp20_bi6_tau::flat_index(3, ci, cj)] = [1.0, 0.0];
        let jet = unit.jet3(0.0, 0.0);
        for (row, &(a, b)) in index_pairs.iter().enumerate() {
            matrix[row][col] = jet.d[a][b][0];
        }
    }

    // Solve for the triangle-BB coefficients of the solved jet, both
    // axes, both (u, v) index assignments (the chart here may
    // transpose the printed map's convention; the accepted assignment
    // is pinned by the regression tests).
    let mut best: Option<[f64; 4]> = None;
    for swap in [false, true] {
        let mut solved_cps = [[0.0_f64; 2]; 10];
        let mut ok = true;
        for axis in 0..2_usize {
            let rhs: Vec<f64> = index_pairs
                .iter()
                .map(|&(a, b)| {
                    let (pa, pb) = if swap { (b, a) } else { (a, b) };
                    partials[pa][pb][axis] * 2.0_f64.powi((a + b) as i32)
                })
                .collect();
            match solve_dense(matrix.clone(), rhs) {
                Some(solution) => {
                    for (slot, value) in solution.into_iter().enumerate() {
                        let (ci, cj) = index_pairs[slot];
                        solved_cps[super::kp20_bi6_tau::flat_index(3, ci, cj)][axis] = value;
                    }
                }
                None => ok = false,
            }
        }
        if !ok {
            continue;
        }
        let mut delta = [0.0_f64; 4];
        for &(ci, cj) in &index_pairs {
            // Exclude the eq (21)-overwritten presentation
            // coefficients (first two boundary CPs on both edges).
            if (ci == 0 && (cj == 1 || cj == 2)) || (cj == 0 && (ci == 1 || ci == 2)) {
                continue;
            }
            let order = ci + cj;
            let idx = super::kp20_bi6_tau::flat_index(3, ci, cj);
            let solved = solved_cps[idx];
            let target = printed.cps[idx];
            delta[order] = delta[order]
                .max((solved[0] - target[0]).abs())
                .max((solved[1] - target[1]).abs());
        }
        if best
            .map(|current| delta.iter().sum::<f64>() < current.iter().sum::<f64>())
            .unwrap_or(true)
        {
            best = Some(delta);
        }
    }
    best.ok_or_else(|| {
        Kp20PlanarTauPaperSolveError::TensorBorder(
            "cubic jet-to-BB conversion is singular".to_string(),
        )
    })
}
