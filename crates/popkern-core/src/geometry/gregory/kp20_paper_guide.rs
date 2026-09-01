//! KP18-faithful guide fit via corner jets of `g ∘ χ̃` (WP3).
//!
//! KP18 §4.2 (the construction KP20 §4.3 delegates to) determines the
//! total-degree-6 guide `g` by least squares between the input
//! tensor-border and the tensor-border assembled from corner jets of
//! `g ∘ χ̃` — the guide composed with the Catmull-Clark characteristic
//! tensor-border. The earlier `kp20_bi6_guide_lsq` fit compared guide
//! rows directly in the sector chart (a structural surrogate,
//! recovery-plan RC2); this module keeps its verified KP18 constraint
//! rows and per-axis KKT solver, and replaces only the fit target.
//!
//! Both sides of the fit run through the same pipeline: order-2
//! corner jets → `hermite_to_bi5_corner_block` → per-piece 6×3 bi-5
//! border rows. For a bi-3 input piece this reconstruction equals its
//! exact degree-raise (a cubic is determined by its end 2-jets), so
//! the identity gates in the module tests are exact.
//!
//! Chart: the guide sector's local coordinates follow the kernel
//! convention (first coordinate along spoke `k+1`, second along spoke
//! `k`); border pieces are `[m^k → c]` and `[c → m^{k+1}]` per sector
//! with rows `[depth][along]`, depth 0 on the boundary (the WP1b
//! chart). At the sector corner `c` the two reflected prolongations
//! of a cornered ring agree through order 2 — which is exactly why
//! KP18's corner-jet constructor stops at order 2 — so the c-corner
//! jets are well-defined.

use super::cc_characteristic_tensor_border::{
    CcCharacteristicTensorBorder, cc_characteristic_tensor_border,
};
use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp20_bi6_guide_lsq::{
    Kp20GuideLsqError, LinearRow, ResidualRows, independent_rows, solve_axis, var_index,
};
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use super::kp20_bi6_math::hermite_to_bi5_corner_block;
use super::kp20_border_degree::elevate_bi3_rows_to_sextic;
pub use super::kp20_paper_border_piece::Kp20PaperGuideBorderPiece;
use crate::geometry::Point3;

const TRI6_CP_COUNT: usize = 28;

/// Two-piece-per-sector border input for the paper guide fit
/// (`[m^k → c]` then `[c → m^{k+1}]`).
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20PaperGuideBorder {
    /// Per-sector border pieces.
    pub sectors: Vec<[Kp20PaperGuideBorderPiece; 2]>,
}

impl Kp20PaperGuideBorder {
    /// Split a smooth production tensor-border (one cubic per sector,
    /// the sector corner at its midpoint) into the two-piece form.
    pub fn from_tensor_border(border: &TensorBorderBi3) -> Self {
        let sectors = border
            .sectors
            .iter()
            .map(|rows| {
                let mut first = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
                let mut second = first;
                for depth in 0..3 {
                    let (a, b) = split_cubic(rows[depth]);
                    first[depth] = a;
                    second[depth] = b;
                }
                [
                    Kp20PaperGuideBorderPiece::from_bi3(first),
                    Kp20PaperGuideBorderPiece::from_bi3(second),
                ]
            })
            .collect();
        Self { sectors }
    }

    /// Lift the planar characteristic tensor-border (z = 0). Used by
    /// the identity gates and available for planar diagnostics.
    pub fn from_characteristic(border: &CcCharacteristicTensorBorder) -> Self {
        let sectors = border
            .sectors
            .iter()
            .map(|sector| {
                [
                    Kp20PaperGuideBorderPiece::from_bi3(sector.near_c.rows),
                    Kp20PaperGuideBorderPiece::from_bi3(sector.near_m.rows),
                ]
            })
            .collect();
        Self { sectors }
    }
}

/// Solved paper guide with measured postconditions.
#[derive(Debug, Clone)]
pub struct Kp20PaperGuideOutput {
    /// One total-degree-6 guide sector per valence sector.
    pub sectors: Vec<Kp20Tri6GuideSector>,
    /// Maximum absolute CP residual of the fitted bi-5 border rows.
    pub max_fit_residual: f64,
    /// Maximum absolute residual of the KP18 equality constraints.
    pub max_constraint_residual: f64,
    /// KKT condition number of the constrained solve.
    pub kkt_condition_number: f64,
}

/// Fit the KP18/KP20 guide to `border` through corner jets of
/// `g ∘ χ̃`, with the extraordinary-point value pinned to
/// `star_point`.
pub fn solve_kp20_paper_guide(
    border: &Kp20PaperGuideBorder,
    star_point: Point3,
) -> Result<Kp20PaperGuideOutput, Kp20GuideLsqError> {
    solve_kp20_paper_guide_with_fairing(border, star_point, 0.0)
}

/// [`solve_kp20_paper_guide`] with an explicit fairing weight: a
/// least-squares penalty on the discrete Laplacian of each sector's
/// control net (zero for affine maps, so the identity gates are
/// unaffected). This is recorded shape conditioning, not part of the
/// paper's construction: on inputs the degree-6 guide cannot fit
/// exactly, the plain LSQ spreads the misfit as interior oscillation
/// that the Hermite sampling copies into the cap; a small weight
/// trades a slightly larger fit residual for a fair interior. `0.0`
/// is the paper-faithful default.
pub fn solve_kp20_paper_guide_with_fairing(
    border: &Kp20PaperGuideBorder,
    star_point: Point3,
    fairing_weight: f64,
) -> Result<Kp20PaperGuideOutput, Kp20GuideLsqError> {
    let n = border.sectors.len();
    if n < 3 || n == 4 || n > 10 {
        return Err(Kp20GuideLsqError::UnsupportedValence(n));
    }
    let chi = cc_characteristic_tensor_border(n).map_err(Kp20GuideLsqError::InvalidTensorBorder)?;

    // Fit rows: for each sector and piece, the 18 bi-5 border CPs of
    // g∘χ̃ (linear in g's coefficients) match the input piece's CPs.
    let mut rows = Vec::new();
    let mut targets: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for (sector, pieces) in border.sectors.iter().enumerate() {
        let chi_pieces = [
            &chi.sectors[sector].near_c.rows,
            &chi.sectors[sector].near_m.rows,
        ];
        for (piece_index, piece) in pieces.iter().enumerate() {
            let chi_rows = planar_rows(chi_pieces[piece_index]);
            for end in 0..2_usize {
                // g∘χ̃ corner block, linear in sector coefficients.
                let chi_jet = border_corner_jet_planar(&chi_rows, end);
                let weights = composed_corner_block_weights(n, sector, &chi_jet);
                // Input corner block (exact degree-raise for cubics).
                let input_jet = border_corner_jet(&piece.rows, end);
                let input_block = jet_to_bi5_block(&input_jet);
                for depth in 0..3_usize {
                    for along in 0..3_usize {
                        let mut terms = Vec::new();
                        for (cp, weight) in weights[along][depth].iter().enumerate() {
                            if *weight != 0.0 {
                                terms.push((var_index(sector, cp), *weight));
                            }
                        }
                        rows.push(LinearRow { terms });
                        let target = input_block[along][depth];
                        targets[0].push(target.x);
                        targets[1].push(target.y);
                        targets[2].push(target.z);
                    }
                }
            }
        }
    }

    if fairing_weight > 0.0 {
        for sector in 0..n {
            for i in 0..=6_usize {
                for j in 0..=(6 - i) {
                    let neighbors: [(i64, i64); 6] =
                        [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)];
                    let mut terms =
                        vec![(var_index(sector, super::tri6_flat_index(6, i, j)), 0.0_f64)];
                    let mut count = 0.0;
                    for (di, dj) in neighbors {
                        let ni = i as i64 + di;
                        let nj = j as i64 + dj;
                        if ni < 0 || nj < 0 || ni + nj > 6 {
                            continue;
                        }
                        terms.push((
                            var_index(sector, super::tri6_flat_index(6, ni as usize, nj as usize)),
                            -fairing_weight,
                        ));
                        count += 1.0;
                    }
                    if count < 2.0 {
                        continue;
                    }
                    terms[0].1 = fairing_weight * count;
                    rows.push(LinearRow { terms });
                    targets[0].push(0.0);
                    targets[1].push(0.0);
                    targets[2].push(0.0);
                }
            }
        }
    }
    let residuals = ResidualRows { rows, targets };
    let all_constraints = kp18_guide_constraint_rows(n);
    let vars = n * TRI6_CP_COUNT;
    let constraints = independent_rows(vars, &all_constraints);

    let mut solved = [
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
        vec![0.0_f64; vars],
    ];
    let mut kkt_condition_number = 0.0_f64;
    for (axis, solved_axis) in solved.iter_mut().enumerate() {
        let center = match axis {
            0 => star_point.x,
            1 => star_point.y,
            _ => star_point.z,
        };
        let axis_solution = solve_axis(
            vars,
            &residuals,
            &residuals.targets[axis],
            &constraints,
            center,
        )
        .ok_or(Kp20GuideLsqError::SingularSystem(axis))?;
        kkt_condition_number = kkt_condition_number.max(axis_solution.kkt_condition_number);
        *solved_axis = axis_solution.values;
    }

    let mut sectors = vec![Kp20Tri6GuideSector::zero(); n];
    for (s, sector) in sectors.iter_mut().enumerate() {
        for (cp, slot) in sector.cps.iter_mut().enumerate() {
            *slot = Point3::new(
                solved[0][var_index(s, cp)],
                solved[1][var_index(s, cp)],
                solved[2][var_index(s, cp)],
            );
        }
    }

    let max_fit_residual = max_row_residual(&residuals, &solved);
    let max_constraint_residual = {
        let dummy_targets: [Vec<f64>; 3] = [
            vec![0.0; all_constraints.len()],
            vec![0.0; all_constraints.len()],
            vec![0.0; all_constraints.len()],
        ];
        let constraint_residuals = ResidualRows {
            rows: all_constraints,
            targets: dummy_targets,
        };
        max_row_residual(&constraint_residuals, &solved)
    };

    Ok(Kp20PaperGuideOutput {
        sectors,
        max_fit_residual,
        max_constraint_residual,
        kkt_condition_number,
    })
}

/// KP18 §3 guide constraint rows, exactly as printed: per sector
/// boundary the C⁰ join (eq 2, i = 0..6) and C¹ join (eq 3,
/// i = 1..6) with weights `w₀ = −1, w₁ = 2c, w₂ = 2(1−c)`, plus the
/// center-expansion single-valuedness rows: the C² join (eq 4) at
/// i = 2, 3 and the C³ join (eq 5) at i = 3. Nothing else — the
/// guide is C¹ with a unique cubic expansion at the center, and the
/// identity map satisfies every row (asserted by the module tests).
///
/// Chart: guide storage `cps[tri6_flat_index(6, q, r)]` maps to the
/// paper's barycentric `b_{p,q,r}` with `p = 6 − q − r` the
/// center-corner power, `q` the own-spoke power, `r` the next-spoke
/// power; the shared boundary between sectors `s` and `s+1` is
/// sector s's `q = 0` edge and sector s+1's `r = 0` edge.
pub(crate) fn kp18_guide_constraint_rows(n: usize) -> Vec<LinearRow> {
    let c = (core::f64::consts::TAU / n as f64).cos();
    let (w0, w1, w2) = (-1.0, 2.0 * c, 2.0 * (1.0 - c));
    let flat = |q: usize, r: usize| super::tri6_flat_index(6, q, r);
    let mut rows = Vec::new();
    for s in 0..n {
        let next = (s + 1) % n;
        let this_var = |q: usize, r: usize| var_index(s, flat(q, r));
        let next_var = |q: usize, r: usize| var_index(next, flat(q, r));
        // Eq (2), C0: b^{s+1}_{d-i,i,0} = b^s_{d-i,0,i}, i = 0..6.
        for i in 0..=6_usize {
            rows.push(LinearRow {
                terms: vec![(next_var(i, 0), 1.0), (this_var(0, i), -1.0)],
            });
        }
        // Eq (3), C1: b^{s+1}_{d-i,i-1,1} = w0 b^s_{d-i,1,i-1}
        //   + w1 b^s_{d-i,0,i} + w2 b^s_{d-i+1,0,i-1}, i = 1..6.
        for i in 1..=6_usize {
            rows.push(LinearRow {
                terms: vec![
                    (next_var(i - 1, 1), 1.0),
                    (this_var(1, i - 1), -w0),
                    (this_var(0, i), -w1),
                    (this_var(0, i - 1), -w2),
                ],
            });
        }
        // Eq (4), C2 at i = 2, 3 (center single-valuedness only).
        for i in 2..=3_usize {
            rows.push(LinearRow {
                terms: vec![
                    (next_var(i - 2, 2), 1.0),
                    (this_var(2, i - 2), -w0 * w0),
                    (this_var(1, i - 1), -2.0 * w0 * w1),
                    (this_var(0, i), -w1 * w1),
                    (this_var(1, i - 2), -2.0 * w0 * w2),
                    (this_var(0, i - 1), -2.0 * w1 * w2),
                    (this_var(0, i - 2), -w2 * w2),
                ],
            });
        }
        // Eq (5), C3 at i = 3.
        rows.push(LinearRow {
            terms: vec![
                (next_var(0, 3), 1.0),
                (this_var(3, 0), -w0 * w0 * w0),
                (this_var(2, 1), -3.0 * w0 * w0 * w1),
                (this_var(1, 2), -3.0 * w0 * w1 * w1),
                (this_var(0, 3), -w1 * w1 * w1),
                (this_var(2, 0), -3.0 * w0 * w0 * w2),
                (this_var(1, 1), -6.0 * w0 * w1 * w2),
                (this_var(0, 2), -3.0 * w1 * w1 * w2),
                (this_var(1, 0), -3.0 * w0 * w2 * w2),
                (this_var(0, 1), -3.0 * w1 * w2 * w2),
                (this_var(0, 0), -w2 * w2 * w2),
            ],
        });
    }
    rows
}

/// Test hook: the assembled fit rows and targets for a border.
#[cfg(test)]
pub(crate) fn kp20_paper_guide_fit_rows_for_test(
    border: &Kp20PaperGuideBorder,
) -> (Vec<LinearRow>, [Vec<f64>; 3]) {
    let n = border.sectors.len();
    let chi = cc_characteristic_tensor_border(n).expect("chi");
    let mut rows = Vec::new();
    let mut targets: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for (sector, pieces) in border.sectors.iter().enumerate() {
        let chi_pieces = [
            &chi.sectors[sector].near_c.rows,
            &chi.sectors[sector].near_m.rows,
        ];
        for (piece_index, piece) in pieces.iter().enumerate() {
            let chi_rows = planar_rows(chi_pieces[piece_index]);
            for end in 0..2_usize {
                let chi_jet = border_corner_jet_planar(&chi_rows, end);
                let weights = composed_corner_block_weights(n, sector, &chi_jet);
                let input_jet = border_corner_jet(&piece.rows, end);
                let input_block = jet_to_bi5_block(&input_jet);
                for depth in 0..3_usize {
                    for along in 0..3_usize {
                        let mut terms = Vec::new();
                        for (cp, weight) in weights[along][depth].iter().enumerate() {
                            if *weight != 0.0 {
                                terms.push((var_index(sector, cp), *weight));
                            }
                        }
                        rows.push(LinearRow { terms });
                        let target = input_block[along][depth];
                        targets[0].push(target.x);
                        targets[1].push(target.y);
                        targets[2].push(target.z);
                    }
                }
            }
        }
    }
    (rows, targets)
}

fn max_row_residual(residuals: &ResidualRows, solved: &[Vec<f64>; 3]) -> f64 {
    let mut worst = 0.0_f64;
    for (row_index, row) in residuals.rows.iter().enumerate() {
        for (axis_solution, axis_targets) in solved.iter().zip(&residuals.targets) {
            let mut acc = 0.0;
            for &(var, weight) in &row.terms {
                acc += weight * axis_solution[var];
            }
            worst = worst.max((acc - axis_targets[row_index]).abs());
        }
    }
    worst
}

fn split_cubic(row: [Point3; 4]) -> ([Point3; 4], [Point3; 4]) {
    let mid =
        |a: Point3, b: Point3| Point3::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, (a.z + b.z) / 2.0);
    let p01 = mid(row[0], row[1]);
    let p12 = mid(row[1], row[2]);
    let p23 = mid(row[2], row[3]);
    let p012 = mid(p01, p12);
    let p123 = mid(p12, p23);
    let center = mid(p012, p123);
    ([row[0], p01, p012, center], [center, p123, p23, row[3]])
}

fn planar_rows(rows: &[[Point3; 4]; 3]) -> [[[f64; 2]; 4]; 3] {
    let mut out = [[[0.0_f64; 2]; 4]; 3];
    for (depth, row) in rows.iter().enumerate() {
        for (along, cp) in row.iter().enumerate() {
            out[depth][along] = [cp.x, cp.y];
        }
    }
    out
}

// ---- Corner jets ----------------------------------------------------------

/// Order-2 corner jet `J[p][q] = ∂^p_s ∂^q_t` of a degree-6-along,
/// depth-2 border piece at its start (`end = 0`) or end (`end = 1`)
/// corner, with the along direction `s` pointing INTO the piece at
/// both corners (frame flip at the end corner) and `t` the depth
/// direction.
fn border_corner_jet(rows: &[[Point3; 7]; 3], end: usize) -> [[Point3; 3]; 3] {
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut jet = [[zero; 3]; 3];
    // Depth-derivative row curves of the bi-3 (depth degree 3) data:
    // T_t = 3(r1 − r0), T_tt = 6(r2 − 2 r1 + r0).
    let mut depth_curves = [[zero; 7]; 3];
    for along in 0..7 {
        let r0 = rows[0][along];
        let r1 = rows[1][along];
        let r2 = rows[2][along];
        depth_curves[0][along] = r0;
        depth_curves[1][along] = Point3::new(
            3.0 * (r1.x - r0.x),
            3.0 * (r1.y - r0.y),
            3.0 * (r1.z - r0.z),
        );
        depth_curves[2][along] = Point3::new(
            6.0 * (r2.x - 2.0 * r1.x + r0.x),
            6.0 * (r2.y - 2.0 * r1.y + r0.y),
            6.0 * (r2.z - 2.0 * r1.z + r0.z),
        );
    }
    for (q, curve) in depth_curves.iter().enumerate() {
        let (value, d1, d2) = sextic_end_jet(curve, end);
        let sign = |p: usize| if end == 1 && p % 2 == 1 { -1.0 } else { 1.0 };
        jet[0][q] = value;
        jet[1][q] = scale(d1, sign(1));
        jet[2][q] = scale(d2, sign(2));
    }
    jet
}

/// Planar variant for the cubic characteristic rows, returning
/// `[f64; 2]` jets (the rows are elevated exactly; end 2-jets are
/// degree-invariant).
fn border_corner_jet_planar(rows: &[[[f64; 2]; 4]; 3], end: usize) -> [[[f64; 2]; 3]; 3] {
    let mut lifted = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
    for depth in 0..3 {
        for along in 0..4 {
            let v = rows[depth][along];
            lifted[depth][along] = Point3::new(v[0], v[1], 0.0);
        }
    }
    let jet = border_corner_jet(&elevate_bi3_rows_to_sextic(lifted), end);
    let mut out = [[[0.0_f64; 2]; 3]; 3];
    for p in 0..3 {
        for q in 0..3 {
            out[p][q] = [jet[p][q].x, jet[p][q].y];
        }
    }
    out
}

/// End 2-jet of a degree-6 row: `d1 = 6(b − a)`, `d2 = 30(c − 2b + a)`.
fn sextic_end_jet(cps: &[Point3; 7], end: usize) -> (Point3, Point3, Point3) {
    let (a, b, c) = if end == 0 {
        (cps[0], cps[1], cps[2])
    } else {
        (cps[6], cps[5], cps[4])
    };
    let value = a;
    let d1 = Point3::new(6.0 * (b.x - a.x), 6.0 * (b.y - a.y), 6.0 * (b.z - a.z));
    let d2 = Point3::new(
        30.0 * (c.x - 2.0 * b.x + a.x),
        30.0 * (c.y - 2.0 * b.y + a.y),
        30.0 * (c.z - 2.0 * b.z + a.z),
    );
    (value, d1, d2)
}

fn scale(p: Point3, s: f64) -> Point3 {
    Point3::new(p.x * s, p.y * s, p.z * s)
}

/// Bi-5 corner block from an order-2 corner jet (per axis via the
/// verified `hermite_to_bi5_corner_block` weights). Index order:
/// `block[along][depth]`.
fn jet_to_bi5_block(jet: &[[Point3; 3]; 3]) -> [[Point3; 3]; 3] {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 3]; 3];
    for axis in 0..3_usize {
        let mut h = [[0.0_f64; 3]; 3];
        for (p, row) in h.iter_mut().enumerate() {
            for (q, slot) in row.iter_mut().enumerate() {
                let v = jet[p][q];
                *slot = match axis {
                    0 => v.x,
                    1 => v.y,
                    _ => v.z,
                };
            }
        }
        let block = hermite_to_bi5_corner_block(h);
        for (p, row) in block.iter().enumerate() {
            for (q, value) in row.iter().enumerate() {
                match axis {
                    0 => out[p][q].x = *value,
                    1 => out[p][q].y = *value,
                    _ => out[p][q].z = *value,
                }
            }
        }
    }
    out
}

// ---- g ∘ χ̃ chain rule ----------------------------------------------------

/// Truncated bivariate Taylor jet through order (2, 2): entry
/// `[p][q]` is the Taylor COEFFICIENT of `s^p t^q` (the derivative
/// divided by `p! q!`). Products truncate exactly — the (2,2)-jet of
/// a product depends only on the factors' (2,2)-jets — so composing
/// the degree-6 guide with the border's corner Taylor through this
/// algebra yields the exact composed corner jet, including the mixed
/// third- and fourth-order entries the KP18 corner-jet constructor
/// requires.
#[derive(Debug, Clone, Copy)]
struct Jet22([[f64; 3]; 3]);

impl Jet22 {
    fn constant(value: f64) -> Self {
        let mut out = [[0.0_f64; 3]; 3];
        out[0][0] = value;
        Jet22(out)
    }

    fn mul(&self, other: &Jet22) -> Jet22 {
        let mut out = [[0.0_f64; 3]; 3];
        for p1 in 0..3 {
            for q1 in 0..3 {
                if self.0[p1][q1] == 0.0 {
                    continue;
                }
                for p2 in 0..(3 - p1) {
                    for q2 in 0..(3 - q1) {
                        out[p1 + p2][q1 + q2] += self.0[p1][q1] * other.0[p2][q2];
                    }
                }
            }
        }
        Jet22(out)
    }

    fn add_scaled(&mut self, other: &Jet22, scale: f64) {
        for p in 0..3 {
            for q in 0..3 {
                self.0[p][q] += scale * other.0[p][q];
            }
        }
    }

    /// Derivatives `h[p][q] = ∂^p_s ∂^q_t` from the coefficients.
    fn derivatives(&self) -> [[f64; 3]; 3] {
        let factorial = [1.0, 1.0, 2.0];
        let mut out = [[0.0_f64; 3]; 3];
        for p in 0..3 {
            for q in 0..3 {
                out[p][q] = self.0[p][q] * factorial[p] * factorial[q];
            }
        }
        out
    }
}

/// Per-CP weights of the composed corner block: for guide sector
/// `sector`, returns `weights[along][depth][cp]` such that the bi-5
/// block of `g ∘ χ̃` at this corner is `Σ_cp weight · g.cps[cp]`
/// (identical weights on every axis, since χ̃ is planar).
fn composed_corner_block_weights(
    n: usize,
    sector: usize,
    chi_jet: &[[[f64; 2]; 3]; 3],
) -> [[Vec<f64>; 3]; 3] {
    // Convert the χ̃ jet from the global plane to the guide sector's
    // local chart. The guide's own convention (pinned by the
    // cross-sector shared-spoke probes: sector s's `jet3(0, t)` runs
    // along the same spoke as sector s+1's `jet3(t, 0)`) is
    // guide-u along spoke k and guide-v along spoke k+1 — the
    // TRANSPOSE of the τ-local convention.
    let alpha = core::f64::consts::TAU / n as f64;
    let a0 = sector as f64 * alpha;
    let a1 = (sector + 1) as f64 * alpha;
    let spoke0 = [a0.cos(), a0.sin()];
    let spoke1 = [a1.cos(), a1.sin()];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    // (guide-u, guide-v) = (component along spoke0, along spoke1).
    let to_local = |g: [f64; 2]| -> [f64; 2] {
        [
            (spoke1[0] * g[1] - spoke1[1] * g[0]) / det,
            (g[0] * spoke0[1] - g[1] * spoke0[0]) / det,
        ]
    };

    // Border Taylor in the guide chart as (2,2)-jets.
    let factorial = [1.0, 1.0, 2.0];
    let mut u_jet = Jet22([[0.0; 3]; 3]);
    let mut v_jet = Jet22([[0.0; 3]; 3]);
    for p in 0..3 {
        for q in 0..3 {
            let local = to_local(chi_jet[p][q]);
            u_jet.0[p][q] = local[0] / (factorial[p] * factorial[q]);
            v_jet.0[p][q] = local[1] / (factorial[p] * factorial[q]);
        }
    }
    let mut w_jet = Jet22::constant(1.0);
    w_jet.add_scaled(&u_jet, -1.0);
    w_jet.add_scaled(&v_jet, -1.0);

    // Power tables through degree 6.
    let powers = |base: &Jet22| -> [Jet22; 7] {
        let mut out = [Jet22::constant(1.0); 7];
        for k in 1..7 {
            out[k] = out[k - 1].mul(base);
        }
        out
    };
    let u_pow = powers(&u_jet);
    let v_pow = powers(&v_jet);
    let w_pow = powers(&w_jet);

    let mut weights: [[Vec<f64>; 3]; 3] = Default::default();
    for row in weights.iter_mut() {
        for slot in row.iter_mut() {
            *slot = vec![0.0_f64; TRI6_CP_COUNT];
        }
    }
    // Basis B_{q_idx, r_idx} = C · u^{q_idx} v^{r_idx} w^{6-q_idx-r_idx}
    // with the multinomial coefficient C.
    let multinomial = |a: usize, b: usize, c: usize| -> f64 {
        let factorial = |k: usize| -> f64 { (1..=k).map(|x| x as f64).product::<f64>().max(1.0) };
        factorial(6) / (factorial(a) * factorial(b) * factorial(c))
    };
    for q_idx in 0..=6_usize {
        for r_idx in 0..=(6 - q_idx) {
            let cp = super::tri6_flat_index(6, q_idx, r_idx);
            let coeff = multinomial(q_idx, r_idx, 6 - q_idx - r_idx);
            let composed = u_pow[q_idx]
                .mul(&v_pow[r_idx])
                .mul(&w_pow[6 - q_idx - r_idx]);
            let mut h = composed.derivatives();
            for row in h.iter_mut() {
                for value in row.iter_mut() {
                    *value *= coeff;
                }
            }
            let block = hermite_to_bi5_corner_block(h);
            for p in 0..3 {
                for q in 0..3 {
                    weights[p][q][cp] = block[p][q];
                }
            }
        }
    }
    weights
}

#[cfg(test)]
#[path = "kp20_paper_guide_tests.rs"]
mod tests;
