//! KP20 constrained degree-9 resampling (PK-599).
//!
//! KP20 step (iii) samples `h = g ∘ τ` into bi-6 pieces from four
//! corner (3,3)-jet blocks, averaged where they overlap (paper Fig. 9c).
//! `h` is not bi-6: on the equal-radius cube corner the averaged middle
//! rows carry a 4× curvature spike (43 against the composition's 10.7)
//! that survives every later stage — the setback cap's ripple. A bi-9
//! piece reproduces `h` (bound measured: 10.70, deviation 1.5e-5).
//!
//! Start from the elevated bi-6 cap `X₀` after step (iv): it satisfies
//! every relation this module keeps (internal C², sector-seam C⁰ and
//! eq (8)/(9), the placed rim bands, the o-jets), and elevation keeps
//! them exact. Each net partitions into owned sets — corner 3×3 blocks
//! and rim-band middles held at `X₀`; internal-seam band middles from
//! one shared jet fit per seam; sector-seam band middles from one
//! constrained fit per piece pair; the 4×4 interior by least squares to
//! `h`. Every relation is re-measured at the end; any failure declines
//! and the caller keeps the bi-6 cap. Design:
//! `docs/architecture/kp20-constrained-bi9-sampling.md`.

use super::kp_degree_elevate::elevate_bi6_grid_to_bi9;
use super::kp20_bi6_cap::Kp20MacroPatch;
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use super::kp20_bi9_exact::ExactPiece;
use super::kp20_bi9_seams::{
    Net, Side, fit_band_side, fit_internal_seam, fit_sector_seam, sector_seam_residual, slot,
};
use crate::geometry::Point3;

/// Piece indices in [`super::Kp20PieceId::all`] order.
const P1: usize = 0;
const P2: usize = 1;
const P3: usize = 2;
const P4: usize = 3;

/// Relation gates, relative to the cap's size: C⁰ and eq (8) (first
/// order) and eq (9) (second order, larger magnitudes).
const FIRST_ORDER_GATE: f64 = 1e-8;
const SECOND_ORDER_GATE: f64 = 1e-6;

/// A resampled cap.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20Bi9Resample {
    /// Bi-9 piece nets `[i][j]`, sector-major, pieces in
    /// [`super::Kp20PieceId::all`] order (the bi-6 cap's flattening).
    pub pieces: Vec<[[Point3; 10]; 10]>,
    /// Worst movement of the target `X₀ + (S₉ − S)` from the elevated
    /// start `X₀` (the sampling artefact removed).
    pub deviation_before: f64,
    /// Worst movement of the resampled cap from `X₀`.
    pub deviation_after: f64,
    /// Worst relation residuals at the end: internal C², sector C⁰,
    /// eq (8), eq (9).
    pub relation_residual: [f64; 4],
}

/// Why the resampling declined (the bi-6 cap serves).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Kp20Bi9Decline {
    /// Step (iv) ran `SoftG`: the start does not satisfy eq (8)/(9).
    #[error("bi-9 resampling needs the HardG2 start")]
    SoftEnforcement,
    /// The start's relations exceed the gates (a misread relation).
    #[error("start violates {relation} on seam {seam}: {residual:.3e} (gate {gate:.1e})")]
    StartInfeasible {
        /// Relation name.
        relation: &'static str,
        /// Seam index (sector of the previous side).
        seam: usize,
        /// Measured residual.
        residual: f64,
        /// Gate.
        gate: f64,
    },
    /// A seam fit failed (inconsistent rows or a singular solve).
    #[error("{kind} seam fit failed at sector {sector}, pair {pair}")]
    SeamFit {
        /// `internal` or `sector`.
        kind: &'static str,
        /// Sector.
        sector: usize,
        /// Pair index within the sector.
        pair: usize,
    },
    /// The interior fit of a piece failed.
    #[error("interior fit failed at sector {sector}, piece {piece}")]
    InteriorFit {
        /// Sector.
        sector: usize,
        /// Piece.
        piece: usize,
    },
    /// A relation exceeds its gate after resampling.
    #[error("resampled cap violates {relation}: {residual:.3e} (gate {gate:.1e})")]
    Residual {
        /// Relation name.
        relation: &'static str,
        /// Measured residual.
        residual: f64,
        /// Gate.
        gate: f64,
    },
}

fn relation_gates(scale: f64) -> [f64; 3] {
    [
        FIRST_ORDER_GATE * scale,
        FIRST_ORDER_GATE * scale,
        SECOND_ORDER_GATE * scale,
    ]
}

/// The worst sector-seam residuals `[C⁰, eq (8), eq (9)]` over all pairs.
fn seam_residuals(n: usize, nets: &[[Net; 4]]) -> Vec<(usize, [f64; 3])> {
    let c = (core::f64::consts::TAU / n as f64).cos();
    (0..n)
        .map(|s| {
            let (prev, curr) = (&nets[s], &nets[(s + 1) % n]);
            let near_o = sector_seam_residual(&prev[P1], &curr[P1], kp20_b_coeffs_near_o(c));
            let near_m = sector_seam_residual(&prev[P2], &curr[P4], kp20_b_coeffs_near_m(c));
            (s, core::array::from_fn(|k| near_o[k].max(near_m[k])))
        })
        .collect()
}

/// Internal seams of a sector: (minus piece, minus side, plus piece, plus side).
const INTERNAL: [(usize, Side, usize, Side); 4] = [
    (P1, Side::I9, P2, Side::I0),
    (P4, Side::I9, P3, Side::I0),
    (P1, Side::J9, P4, Side::J0),
    (P2, Side::J9, P3, Side::J0),
];

fn internal_c2_residual(nets: &[[Net; 4]]) -> f64 {
    let mut worst = 0.0f64;
    for sector in nets {
        for &(m, sm, p, sp) in &INTERNAL {
            for along in 0..10 {
                let get = |net: &Net, side: Side, d: usize| {
                    let (i, j) = slot(side, d, along);
                    net[i][j]
                };
                let v = |q: Point3| crate::geometry::Vec3::new(q.x, q.y, q.z);
                let (m0, m1, m2) = (
                    v(get(&sector[m], sm, 0)),
                    v(get(&sector[m], sm, 1)),
                    v(get(&sector[m], sm, 2)),
                );
                let (p0, p1, p2) = (
                    v(get(&sector[p], sp, 0)),
                    v(get(&sector[p], sp, 1)),
                    v(get(&sector[p], sp, 2)),
                );
                let c0 = (m0 - p0).length();
                let c1 = ((p1 - p0) - (m0 - m1)).length();
                let c2 = ((p2 - p1 * 2.0 + p0) - (m2 - m1 * 2.0 + m0)).length();
                worst = worst.max(c0).max(c1).max(c2);
            }
        }
    }
    worst
}

/// Worst distance between two sets of nets (sampled on a grid).
fn deviation(a: &[[Net; 4]], b: &[[Net; 4]]) -> f64 {
    let mut worst = 0.0f64;
    for (sa, sb) in a.iter().zip(b) {
        for (na, nb) in sa.iter().zip(sb) {
            for i in 0..=12 {
                for j in 0..=12 {
                    let (u, v) = (i as f64 / 12.0, j as f64 / 12.0);
                    worst = worst.max(eval(na, u, v).distance_to(&eval(nb, u, v)));
                }
            }
        }
    }
    worst
}

/// Evaluate a bi-9 net.
pub(super) fn eval(net: &Net, u: f64, v: f64) -> Point3 {
    let (bu, bv): (Vec<f64>, Vec<f64>) = (
        (0..10)
            .map(|i| super::kp20_bi9_seams::b9(i, u)[0])
            .collect(),
        (0..10)
            .map(|j| super::kp20_bi9_seams::b9(j, v)[0])
            .collect(),
    );
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for i in 0..10 {
        for j in 0..10 {
            let w = bu[i] * bv[j];
            x += w * net[i][j].x;
            y += w * net[i][j].y;
            z += w * net[i][j].z;
        }
    }
    Point3::new(x, y, z)
}

/// Resample the cap. `sampled` is the step (iii) bi-6 sampling and
/// `start` the bi-6 cap after step (iv) with `HardG2` enforcement; `tau`
/// is the unrotated sector-0 τ. The target is `start + (h − sampled)`:
/// the composition with the later stages' displacement kept.
pub fn resample_kp20_cap_bi9(
    n: usize,
    guide: &[Kp20Tri6GuideSector],
    tau: &Kp20MacroPatch,
    sampled: &[Kp20MacroPatch],
    start: &[Kp20MacroPatch],
) -> Result<Kp20Bi9Resample, Kp20Bi9Decline> {
    let elevate = |caps: &[Kp20MacroPatch]| -> Vec<[Net; 4]> {
        caps.iter()
            .map(|m| core::array::from_fn(|k| elevate_bi6_grid_to_bi9(m.pieces[k].cps)))
            .collect()
    };
    let (start9, sampled9) = (elevate(start), elevate(sampled));
    let exacts: Vec<Vec<ExactPiece<'_>>> = (0..n)
        .map(|s| {
            tau.pieces
                .iter()
                .map(|piece| ExactPiece::new(n, s, &guide[s], piece))
                .collect()
        })
        .collect();
    let scale = cap_scale(&start9);
    let gates = relation_gates(scale);
    let names = ["sector C0", "eq (8)", "eq (9)"];
    for (seam, res) in seam_residuals(n, &start9) {
        for k in 0..3 {
            if res[k] > gates[k] {
                return Err(Kp20Bi9Decline::StartInfeasible {
                    relation: names[k],
                    seam,
                    residual: res[k],
                    gate: gates[k],
                });
            }
        }
    }
    // 1. Accurate bi-9 sampling S₉ of h from the elevated step (iii)
    //    sampling: corner blocks (exact 2-jets) held; every band middle
    //    and interior re-fitted to h (internal seams shared, so C² exact).
    let mut s9 = sampled9.clone();
    for (s, sector) in s9.iter_mut().enumerate() {
        for (pair, &(m, sm, p, sp)) in INTERNAL.iter().enumerate() {
            let (mut minus, mut plus) = (sector[m], sector[p]);
            fit_internal_seam(&mut minus, &mut plus, sm, sp, &exacts[s][m], &exacts[s][p]).ok_or(
                Kp20Bi9Decline::SeamFit {
                    kind: "internal",
                    sector: s,
                    pair,
                },
            )?;
            sector[m] = minus;
            sector[p] = plus;
        }
        let sides = [
            (P1, Side::J0),
            (P2, Side::J0),
            (P1, Side::I0),
            (P4, Side::I0),
            (P2, Side::I9),
            (P3, Side::I9),
            (P3, Side::J9),
            (P4, Side::J9),
        ];
        for (pair, (piece, side)) in sides.into_iter().enumerate() {
            fit_band_side(&mut sector[piece], side, &exacts[s][piece]).ok_or(
                Kp20Bi9Decline::SeamFit {
                    kind: "sector side",
                    sector: s,
                    pair,
                },
            )?;
        }
        for k in 0..4 {
            super::kp20_bi9_interior::fit_interior(&mut sector[k], &exacts[s][k]).ok_or(
                Kp20Bi9Decline::InteriorFit {
                    sector: s,
                    piece: k,
                },
            )?;
        }
    }
    // 2–3. Target T₉ = X₀ + (S₉ − S): the sampling artefact removed,
    //      every later-stage change kept. Placement owns the rim bands, so
    //      the correction is zero there.
    let mut nets = start9.clone();
    for s in 0..n {
        for k in 0..4 {
            for i in 0..10 {
                for j in 0..10 {
                    if rim_band(k, i, j) {
                        continue;
                    }
                    let d = s9[s][k][i][j] - sampled9[s][k][i][j];
                    nets[s][k][i][j] = nets[s][k][i][j] + d;
                }
            }
        }
    }
    let deviation_before = deviation(&start9, &nets);
    // 4. Sector seams: minimal movement from T₉ onto C⁰ + eq (8)/(9).
    let target = nets.clone();
    let c = (core::f64::consts::TAU / n as f64).cos();
    for s in 0..n {
        let t = (s + 1) % n;
        for (pair, (a, b, coeffs)) in [
            (P1, P1, kp20_b_coeffs_near_o(c)),
            (P2, P4, kp20_b_coeffs_near_m(c)),
        ]
        .into_iter()
        .enumerate()
        {
            let (mut prev, mut curr) = (nets[s][a], nets[t][b]);
            fit_sector_seam(&mut prev, &mut curr, coeffs, &target[s][a], &target[t][b]).ok_or(
                Kp20Bi9Decline::SeamFit {
                    kind: "sector",
                    sector: s,
                    pair,
                },
            )?;
            nets[s][a] = prev;
            nets[t][b] = curr;
        }
    }
    // Every relation, re-measured.
    let c2 = internal_c2_residual(&nets);
    if c2 > gates[0] {
        return Err(Kp20Bi9Decline::Residual {
            relation: "internal C2",
            residual: c2,
            gate: gates[0],
        });
    }
    let mut worst = [0.0f64; 3];
    for (_, res) in seam_residuals(n, &nets) {
        for k in 0..3 {
            worst[k] = worst[k].max(res[k]);
        }
    }
    for k in 0..3 {
        if worst[k] > gates[k] {
            return Err(Kp20Bi9Decline::Residual {
                relation: names[k],
                residual: worst[k],
                gate: gates[k],
            });
        }
    }
    let deviation_after = deviation(&start9, &nets);
    Ok(Kp20Bi9Resample {
        pieces: nets
            .iter()
            .flat_map(|sector| sector.iter().copied())
            .collect(),
        deviation_before,
        deviation_after,
        relation_residual: [c2, worst[0], worst[1], worst[2]],
    })
}

/// Whether coefficient `(i, j)` of piece `k` lies in a placement-owned
/// rim band (depth 0–2 of a rim edge: P2 `u = 1`, P3 `u = 1` and
/// `v = 1`, P4 `v = 1`).
fn rim_band(k: usize, i: usize, j: usize) -> bool {
    match k {
        P2 => i >= 7,
        P3 => i >= 7 || j >= 7,
        P4 => j >= 7,
        _ => false,
    }
}

/// Bounding-box diagonal of the nets (floor 1e-9).
fn cap_scale(nets: &[[Net; 4]]) -> f64 {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for net in nets.iter().flatten() {
        for p in net.iter().flatten() {
            for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[k] = lo[k].min(c);
                hi[k] = hi[k].max(c);
            }
        }
    }
    (0..3)
        .map(|k| (hi[k] - lo[k]).powi(2))
        .sum::<f64>()
        .sqrt()
        .max(1e-9)
}

#[cfg(test)]
#[path = "kp20_bi9_resample_tests.rs"]
mod tests;
