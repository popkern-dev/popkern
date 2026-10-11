//! KP20 §4.3 paper τ solved from the derived characteristic
//! tensor-border.
//!
//! WP2 of the recovery plan
//! (`docs/architecture/kp20-paper-faithful-recovery-plan.md`): τ is a
//! planar bi-6 cap of 2×2 quadrant pieces per sector whose input
//! tensor-border is the Catmull-Clark characteristic ring's
//! tensor-border (WP1b), constrained per §4.2/§4.3 and F₃-minimized.
//! The row system here is built from first principles in one explicit
//! chart (below), reusing only geometry-independent substrate: the F₃
//! KKT solver, the internal-C² quadrant lattice, the eq (8)/(9)
//! residual Jacobians, and the eq (14) border reparameterizer.
//!
//! # Chart
//!
//! Canonical sector 0, solved in GLOBAL plane coordinates. Domain
//! coordinates `(a, b)` with global position `a·spoke₁ + b·spoke₀`
//! (`spoke₀` at angle 0, `spoke₁` at `2π/n`): `o` at `(0,0)`, `m⁰`
//! (spoke-0 rim point) at `(0,1)`, `m¹` at `(1,0)`, the sector corner
//! `c` near `(1,1)`. Quadrant pieces with globally-aligned BB indices
//! `(i ↔ a, j ↔ b)`:
//!
//! - `NearOriginLeft`  = P1 = `[0,½]²` (owns `o` at `cps[0][0]`),
//! - `NearOriginRight` = P2 = `[½,1]×[0,½]` (rim edge at `i = 6`),
//! - `NearMidpointRight` = P3 = `[½,1]²` (corner piece, rim edges at
//!   `i = 6` and `j = 6`),
//! - `NearMidpointLeft`  = P4 = `[0,½]×[½,1]` (rim edge at `j = 6`).
//!
//! Diagonal symmetry maps P1, P3 to themselves and swaps P2 ↔ P4,
//! matching KP20 Fig. 8b.
//!
//! # Border pins (eq 14)
//!
//! Each rim edge is parameterized from `c` toward `m` and split at the
//! macro-piece boundary: t⁰ = the c-half with `d⁰ ~ [0,0,0,δ]`
//! (vanishing to second order at `c`), t¹ = the m-half with
//! `d¹ ~ [δ,2δ,4δ,d¹₃]` reaching the free parameter `d¹₃ = 2c/3` at
//! `m` — the end where §4.2's "at m" consistency with eq (13) lives.
//! This is mirror-symmetric within the sector and rotationally
//! consistent across sectors (`d = d¹₃` on both sides of every spoke
//! point `m`, `d = 0` at every corner `c`).

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::cc_characteristic_tensor_border::cc_characteristic_tensor_border;
use super::kp20_bi6_cap::{
    Kp20BoundaryHalf, Kp20MacroPatch, Kp20PieceId, kp20_eq14_border_d_profiles,
    kp20_reparameterise_border_piece,
};
use super::kp20_bi6_local_constraint::{kp20_eq8_jacobian_dense, kp20_eq9_jacobian_dense};
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use super::kp20_bi6_tau_f3_solve::{
    Kp20PlanarTauCoeff, Kp20PlanarTauF3ScalarSolveInput, Kp20PlanarTauF3SolveOutput,
    Kp20PlanarTauLinearConstraint, Kp20PlanarTauScalarCoeff, Kp20PlanarTauScalarLinearConstraint,
    kp20_planar_tau_internal_c2_constraints, solve_kp20_planar_tau_f3_with_scalar_constraints,
};
use super::kp20_bi6_tau_solve_assemble::Kp20PlanarTauPaperSolveError;
use super::kp20_border_degree::{
    elevate_bi3_rows_to_sextic, reverse_sextic_rows, split_sextic_rows,
};
use super::kp20_paper_tau_measure::{
    kp20_paper_tau_affine_seed, kp20_paper_tau_jet_delta, kp20_paper_tau_orientation_summary,
};
use crate::geometry::Point3;
use popkern_math::PortableMath;

const P1: Kp20PieceId = Kp20PieceId::NearOriginLeft;
const P2: Kp20PieceId = Kp20PieceId::NearOriginRight;
const P3: Kp20PieceId = Kp20PieceId::NearMidpointRight;
const P4: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

/// Solved paper τ for one valence, with its acceptance measurements.
#[derive(Debug, Clone)]
pub struct Kp20PaperTauReport {
    /// Valence.
    pub n: usize,
    /// Solved planar τ (canonical sector 0, global plane coords).
    pub tau: Kp20MacroPatch,
    /// Maximum linear row residual of the F₃ KKT solve.
    pub max_row_residual: f64,
    /// Sampled Jacobian summary over all four pieces.
    pub negative_jacobian_samples: usize,
    /// Near-zero Jacobian samples.
    pub near_zero_jacobian_samples: usize,
    /// Minimum sampled Jacobian determinant.
    pub min_jacobian: f64,
    /// Worst per-order deviation of the solved cubic jet at `o` from
    /// the printed eq (21)/(22) table (orders 0..=3, sector-scale
    /// parameters). `None` at a valence the paper does not tabulate
    /// (n = 4, the regular valence): the τ solve itself is
    /// valence-generic, only this cross-check needs the table.
    pub jet_delta: Option<[f64; 4]>,
}

/// Solve the KP20 paper τ for one valence from the derived
/// characteristic tensor-border (canonical sector 0; the full τ is its
/// rotations).
///
/// τ depends only on `n` (the characteristic map of the valence), so the
/// result is memoised per valence for the process: every KP20 corner of
/// a valence reuses one solve (PK-604 measured 0.085 s per solve at
/// n = 8). The cached value is the exact result of the first solve.
pub fn solve_kp20_paper_tau(n: usize) -> Result<Kp20PaperTauReport, Kp20PlanarTauPaperSolveError> {
    use std::collections::HashMap;
    use std::sync::{Mutex, OnceLock};
    type Cache = Mutex<HashMap<usize, Result<Kp20PaperTauReport, Kp20PlanarTauPaperSolveError>>>;
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(hit) = cache
        .lock()
        .unwrap_or_else(|poison| poison.into_inner())
        .get(&n)
    {
        return hit.clone();
    }
    let solved = compute_kp20_paper_tau(n);
    cache
        .lock()
        .unwrap_or_else(|poison| poison.into_inner())
        .insert(n, solved.clone());
    solved
}

/// The uncached τ solve behind [`solve_kp20_paper_tau`].
fn compute_kp20_paper_tau(n: usize) -> Result<Kp20PaperTauReport, Kp20PlanarTauPaperSolveError> {
    let solved = solve_kp20_paper_tau_raw(n)?;
    let orientation = kp20_paper_tau_orientation_summary(&solved.tau_macro_patch);
    let jet_delta = match super::kp20_bi6_math::kappa_for_valence(n) {
        Some(_) => Some(kp20_paper_tau_jet_delta(n, &solved.tau_macro_patch)?),
        None => None,
    };
    Ok(Kp20PaperTauReport {
        n,
        tau: solved.tau_macro_patch,
        max_row_residual: solved.max_constraint_residual,
        negative_jacobian_samples: orientation.0,
        near_zero_jacobian_samples: orientation.1,
        min_jacobian: orientation.2,
        jet_delta,
    })
}

/// Diagnostic variant returning the raw solve output (no gates).
pub(crate) fn solve_kp20_paper_tau_raw(
    n: usize,
) -> Result<Kp20PlanarTauF3SolveOutput, Kp20PlanarTauPaperSolveError> {
    let border =
        cc_characteristic_tensor_border(n).map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let mut vector_constraints = kp20_planar_tau_internal_c2_constraints();
    vector_constraints.extend(kp20_paper_tau_border_pin_rows(n, &border.sectors[0]));
    vector_constraints.extend(kp20_paper_tau_origin_rows());
    let mut scalar_constraints = kp20_paper_tau_rotational_g1_g2_rows(n);
    scalar_constraints.extend(kp20_paper_tau_diagonal_symmetry_rows(n));
    let seed = kp20_paper_tau_affine_seed(n);
    solve_kp20_planar_tau_f3_with_scalar_constraints(&Kp20PlanarTauF3ScalarSolveInput {
        initial: seed,
        vector_constraints,
        scalar_constraints,
        regularization_weight: 1.0e-10,
    })
    .map_err(Kp20PlanarTauPaperSolveError::from)
}

// ---- Border pins --------------------------------------------------------

/// Equation (14) border pins for the canonical sector under the chosen
/// layout, pinned onto P4/P3 (the `m⁰` edge) and P3/P2 (the `m¹`
/// edge).
fn kp20_paper_tau_border_pin_rows(
    n: usize,
    sector: &super::cc_characteristic_tensor_border::CcSectorBorder,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    let c = (core::f64::consts::TAU / n as f64).pcos();
    kp20_paper_tau_border_pin_rows_with_d1_3(n, sector, 2.0 * c / 3.0)
}

/// Equation (14) border pins for the canonical sector: each rim edge
/// parameterized c→m and split at the macro seam, t⁰ = c-half with
/// d⁰, t¹ = m-half with d¹ (d¹₃ at m), pinned onto P4/P3 (the m⁰
/// edge) and P3/P2 (the m¹ edge).
///
/// Rejected alternatives (recorded falsifications): parameterizing
/// m→c (d¹₃ at the corner) breaks the otherwise exactly-closed
/// rotational G¹/G² rows at ~0.58 and is d¹₃-independent; treating
/// t⁰/t¹ as the two whole edges breaks them at ~0.4–1.9; flipping the
/// d sign on one edge breaks them the same way. The accepted reading
/// closes every family at solver tolerance.
pub(crate) fn kp20_paper_tau_border_pin_rows_with_d1_3(
    n: usize,
    sector: &super::cc_characteristic_tensor_border::CcSectorBorder,
    d1_3: f64,
) -> Vec<Kp20PlanarTauLinearConstraint> {
    kp20_paper_border_band_placements(
        n,
        &elevate_bi3_rows_to_sextic(sector.near_c.rows),
        &elevate_bi3_rows_to_sextic(sector.near_m.rows),
        d1_3,
    )
    .into_iter()
    .map(|(piece, i, j, value)| {
        Kp20PlanarTauLinearConstraint::pin(Kp20PlanarTauCoeff::new(piece, i, j), [value.x, value.y])
    })
    .collect()
}

/// Equation (14) border-band placements for one sector, shared by the
/// τ solve (as pins) and the §4.4 cap assembly (as overwrites).
///
/// Accepted layout (WP2): each rim edge parameterized c→m, split at
/// the macro seam, d⁰ on the c-half, d¹ reaching `d¹₃` at m; the
/// corner-overlap 3×3 of the corner piece carries the two strips'
/// AVERAGE (the unique diagonally-symmetric choice on cornered
/// characteristic data; exact for smooth production borders where the
/// two strips agree). Rejected alternatives are recorded
/// falsifications: m→c parameterization (d¹₃ at the corner) and the
/// two-whole-edges reading both break the otherwise-exact rotational
/// G¹/G² closure; a per-edge d sign flip does the same.
///
/// Rows are degree 6 along (PK-398); cubic sources elevate exactly
/// and reproduce the historical placements.
pub(crate) fn kp20_paper_border_band_placements(
    n: usize,
    near_c_rows: &[[Point3; 7]; 3],
    near_m_rows: &[[Point3; 7]; 3],
    d1_3: f64,
) -> Vec<(Kp20PieceId, usize, usize, Point3)> {
    let c = (core::f64::consts::TAU / n as f64).pcos();
    let (d0, d1) = kp20_eq14_border_d_profiles(c, d1_3);

    // Edge A: [c -> m^0] = reverse of the near_c piece ([m^0 -> c]);
    // domain b = 1 line, a descending from c toward m^0.
    let edge_a = reverse_sextic_rows(near_c_rows);
    let (a_c_half, a_m_half) = split_sextic_rows(&edge_a);
    let band_a_c = kp20_reparameterise_border_piece(a_c_half, d0);
    let band_a_m = kp20_reparameterise_border_piece(a_m_half, d1);
    // Edge B: [c -> m^1] = the near_m piece as stored.
    let (b_c_half, b_m_half) = split_sextic_rows(near_m_rows);
    let band_b_c = kp20_reparameterise_border_piece(b_c_half, d0);
    let band_b_m = kp20_reparameterise_border_piece(b_m_half, d1);

    // Bands run c->m while piece indices ascend toward c: reversed.
    let mut a_rows = Vec::new();
    band_placements(&mut a_rows, P3, &band_a_c, true, true);
    let mut b_rows = Vec::new();
    band_placements(&mut b_rows, P3, &band_b_c, false, true);
    average_overlapping_placements(&mut a_rows, &mut b_rows);
    let mut rows = Vec::with_capacity(4 * 21);
    rows.extend(a_rows);
    rows.extend(b_rows);
    band_placements(&mut rows, P4, &band_a_m, true, true);
    band_placements(&mut rows, P2, &band_b_m, false, true);
    rows
}

fn band_placements(
    rows: &mut Vec<(Kp20PieceId, usize, usize, Point3)>,
    piece: Kp20PieceId,
    band: &Kp20BoundaryHalf,
    j_band: bool,
    reverse_along: bool,
) {
    for depth in 0..3_usize {
        for along in 0..7_usize {
            let source = band.rows[depth][if reverse_along { 6 - along } else { along }];
            let (i, j) = if j_band {
                (along, 6 - depth)
            } else {
                (6 - depth, along)
            };
            rows.push((piece, i, j, source));
        }
    }
}

/// Corner-compatibility measure for one sector's two border pieces:
/// the maximum disagreement between the two eq (14) bands on the
/// corner piece's overlapping 3×3 block at `c`. Zero (to roundoff)
/// exactly when the two pieces are Hermite data of ONE surface patch
/// cornered at `c` (e.g., a subdivision-ring extraction); large
/// values flag invalid input such as a smooth cubic split at its
/// parameter midpoint, whose averaged corner bands kink the cap rim.
pub fn kp20_paper_border_corner_mismatch(
    n: usize,
    near_c_rows: &[[Point3; 7]; 3],
    near_m_rows: &[[Point3; 7]; 3],
    d1_3: f64,
) -> f64 {
    let c = (core::f64::consts::TAU / n as f64).pcos();
    let (d0, _) = kp20_eq14_border_d_profiles(c, d1_3);
    let edge_a = reverse_sextic_rows(near_c_rows);
    let (a_c_half, _) = split_sextic_rows(&edge_a);
    let band_a_c = kp20_reparameterise_border_piece(a_c_half, d0);
    let (b_c_half, _) = split_sextic_rows(near_m_rows);
    let band_b_c = kp20_reparameterise_border_piece(b_c_half, d0);
    let mut a_rows = Vec::new();
    band_placements(&mut a_rows, P3, &band_a_c, true, true);
    let mut b_rows = Vec::new();
    band_placements(&mut b_rows, P3, &band_b_c, false, true);
    let mut worst = 0.0_f64;
    let mut worst_slot = (0_usize, 0_usize);
    for &(piece, i, j, a_value) in &a_rows {
        for &(bp, bi, bj, b_value) in &b_rows {
            if piece == bp && i == bi && j == bj {
                let d = (a_value.x - b_value.x)
                    .abs()
                    .max((a_value.y - b_value.y).abs())
                    .max((a_value.z - b_value.z).abs());
                if d > worst {
                    worst = d;
                    worst_slot = (i, j);
                }
            }
        }
    }
    if crate::trace::trace_enabled("kp-corner-mismatch") {
        eprintln!("    corner mismatch worst at slot {worst_slot:?}: {worst:e}");
    }
    worst
}

fn average_overlapping_placements(
    primary: &mut [(Kp20PieceId, usize, usize, Point3)],
    secondary: &mut Vec<(Kp20PieceId, usize, usize, Point3)>,
) {
    secondary.retain(|&(piece, i, j, value)| {
        if let Some(existing) = primary
            .iter_mut()
            .find(|(p, pi, pj, _)| *p == piece && *pi == i && *pj == j)
        {
            existing.3 = Point3::new(
                (existing.3.x + value.x) / 2.0,
                (existing.3.y + value.y) / 2.0,
                (existing.3.z + value.z) / 2.0,
            );
            false
        } else {
            true
        }
    });
}

// ---- Origin rows (Fig 8b pin + rotational fixpoint) ---------------------

fn kp20_paper_tau_origin_rows() -> Vec<Kp20PlanarTauLinearConstraint> {
    // Rotational C0 at the extraordinary point: tau(o) is a fixed
    // point of the sector rotation, hence the origin.
    let mut rows = vec![Kp20PlanarTauLinearConstraint::pin(
        Kp20PlanarTauCoeff::new(P1, 0, 0),
        [0.0, 0.0],
    )];
    // Fig. 8b indicated averaging on the two spoke boundary runs at o:
    // the second boundary coefficient is the midpoint of its
    // neighbors (KP20 §4.3, black dot := ½ red + ½ blue).
    rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
        (Kp20PlanarTauCoeff::new(P1, 0, 1), 2.0),
        (Kp20PlanarTauCoeff::new(P1, 0, 0), -1.0),
        (Kp20PlanarTauCoeff::new(P1, 0, 2), -1.0),
    ]));
    rows.push(Kp20PlanarTauLinearConstraint::homogeneous(vec![
        (Kp20PlanarTauCoeff::new(P1, 1, 0), 2.0),
        (Kp20PlanarTauCoeff::new(P1, 0, 0), -1.0),
        (Kp20PlanarTauCoeff::new(P1, 2, 0), -1.0),
    ]));
    rows
}

// ---- Rotational G1/G2 rows across the spoke ------------------------------

fn scalar(piece: Kp20PieceId, i: usize, j: usize, axis: usize) -> Kp20PlanarTauScalarCoeff {
    Kp20PlanarTauScalarCoeff::new(Kp20PlanarTauCoeff::new(piece, i, j), axis)
}

/// KP20 eq (8)/(9) rows across spoke 0 between the previous sector
/// (the canonical patch rotated by −2π/n) and the canonical sector.
///
/// Previous side: pieces P1, P2 in their natural frame (spoke at
/// `v = 0`, `u` ascending o→m along the domain `b = 0` edge of the
/// previous sector, which is global spoke 0). Current side: pieces
/// P1, P4 with the spoke at their `a = 0` edge, so their Jacobian
/// indices enter transposed. Near-o pair (P1_prev, P1_curr) with
/// `b ~ c·[2,2,1]`; near-m pair (P2_prev, P4_curr) with
/// `b̲ ~ [c,0,0]` (value `c` at the macro seam, 0 at `m`).
fn kp20_paper_tau_rotational_g1_g2_rows(n: usize) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let alpha = core::f64::consts::TAU / n as f64;
    let c = alpha.pcos();
    // R^{-1}: rotation by −alpha applied to previous-side values.
    let (rc, rs) = (alpha.pcos(), -alpha.psin());
    let mut rows = Vec::new();

    let mut add_pair = |prev_piece: Kp20PieceId, curr_piece: Kp20PieceId, b: [f64; 3]| {
        let jacobians: Vec<[f64; 98]> = kp20_eq8_jacobian_dense(b)
            .into_iter()
            .chain(kp20_eq9_jacobian_dense(b))
            .collect();
        for jac in jacobians {
            for axis in 0..2_usize {
                let mut terms = Vec::new();
                for side in 0..2_usize {
                    for i in 0..7_usize {
                        for j in 0..7_usize {
                            let w = jac[side * 49 + i * 7 + j];
                            if w == 0.0 {
                                continue;
                            }
                            if side == 0 {
                                // Previous side, rotated by R^{-1}:
                                // (R^{-1} q)_x = rc·qx − rs·qy,
                                // (R^{-1} q)_y = rs·qx + rc·qy.
                                let (wx, wy) = if axis == 0 {
                                    (w * rc, -w * rs)
                                } else {
                                    (w * rs, w * rc)
                                };
                                terms.push((scalar(prev_piece, i, j, 0), wx));
                                terms.push((scalar(prev_piece, i, j, 1), wy));
                            } else {
                                // Current side, transposed indices.
                                terms.push((scalar(curr_piece, j, i, axis), w));
                            }
                        }
                    }
                }
                rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(terms));
            }
        }
    };

    add_pair(P1, P1, kp20_b_coeffs_near_o(c));
    add_pair(P2, P4, kp20_b_coeffs_near_m(c));
    // Cross-sector C0 along the separating curve: the previous
    // sector's (i, 0) row, rotated by R^{-1}, is the current sector's
    // (0, i) row (near o: P1 ↔ P1; near m: P2 ↔ P4). At n >= 5 the
    // G rows above and the border pins imply it (measured 8e-13); at
    // the regular valence n = 4 (c = 0, b(u) ≡ 0) they do not, and
    // the solve left the spoke edge 6.6e-3 off its spoke — every
    // sampled cap then broke C0 across sectors by 1.3e-2 (PK-518).
    for (prev_piece, curr_piece) in [(P1, P1), (P2, P4)] {
        for i in 0..7_usize {
            for axis in 0..2_usize {
                let (wx, wy) = if axis == 0 { (rc, -rs) } else { (rs, rc) };
                rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                    (scalar(prev_piece, i, 0, 0), wx),
                    (scalar(prev_piece, i, 0, 1), wy),
                    (scalar(curr_piece, 0, i, axis), -1.0),
                ]));
            }
        }
    }
    rows
}

// ---- Diagonal symmetry ----------------------------------------------------

/// KP20 §4.3: τ is diagonally symmetric. In the canonical global
/// frame the sector-bisector reflection is `S(z) = w·conj(z)` with
/// `w = e^{2πi/n}`; the domain swap exchanges `(i, j)` and maps
/// P1→P1, P3→P3, P2↔P4.
fn kp20_paper_tau_diagonal_symmetry_rows(n: usize) -> Vec<Kp20PlanarTauScalarLinearConstraint> {
    let alpha = core::f64::consts::TAU / n as f64;
    let (wx, wy) = (alpha.pcos(), alpha.psin());
    let mut rows = Vec::new();
    // S(x, y) = (wx·x + wy·y, wy·x − wx·y).
    let mut add = |target: Kp20PieceId, source: Kp20PieceId| {
        for i in 0..7_usize {
            for j in 0..7_usize {
                if target == source && (i, j) > (j, i) {
                    continue; // each unordered pair once
                }
                for axis in 0..2_usize {
                    let (sx, sy) = if axis == 0 { (wx, wy) } else { (wy, -wx) };
                    rows.push(Kp20PlanarTauScalarLinearConstraint::homogeneous(vec![
                        (scalar(target, i, j, axis), 1.0),
                        (scalar(source, j, i, 0), -sx),
                        (scalar(source, j, i, 1), -sy),
                    ]));
                }
            }
        }
    };
    add(P1, P1);
    add(P3, P3);
    add(P2, P4);
    rows
}

#[cfg(test)]
#[path = "kp20_paper_tau_tests.rs"]
mod tests;
