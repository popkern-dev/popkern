//! §5.3 Phase-10 polynomial-identity G1 solver for a **single**
//! shared spoke between two bidegree-7 two-stripe Gregory patches.
//! # The identity (Karciauskas-Peters form)
//!
//! For patches `A`, `B` sharing a boundary spoke `C(t)`:
//!
//! ```text
//!   xtan_A(t) + α · xtan_B(t) + β(t) · T(t) ≡ 0     ∀ t ∈ [0, 1]
//! ```
//!
//! where `xtan_A`, `xtan_B` are the cross-boundary tangent
//! polynomials (degree 7 for two-stripe), `T(t) = C'(t)` is the
//! boundary tangent (degree 6), and `β(t)` is a scalar blending
//! polynomial (degree 1 by default — matches the `xtan` degree
//! after multiplication by `T`). `α = −1` for the standard
//! consistent-orientation case covered by this module.
//!
//! Because both sides of the shared boundary evaluate to the same
//! curve, the `P_A[0][·]` / `P_B[·][0]` terms in the two
//! `xtan`'s cancel and the identity reduces to a purely linear
//! equation in **the row-1 near-spoke CPs of both patches plus
//! the two Bernstein coefficients of `β`**.
//!
//! # What we zero
//!
//! Setting every Bernstein coefficient of the identity's LHS to
//! zero gives **8 vector equations** (= 24 scalar equations) per
//! spoke. They are exact: satisfaction at the coefficient level
//! means the polynomial identity holds at *every* `t ∈ [0, 1]`,
//! not just at samples. This is the structural improvement over
//! the sampled solver in [`super::spoke_g1_solver`].
//!
//! # What we solve
//!
//! The system is under-determined (50 unknowns, 24 constraints
//! with default `β` degree 1), so we solve the minimum-norm
//! problem
//!
//! ```text
//!   minimise ‖x − x_seed‖²   subject to   A · x = b,
//! ```
//!
//! where `x_seed` is the CP state from the biseptic-promoted
//! starting point. The result stays close to the seed while
//! satisfying G1 exactly.
//!
//! # Limitations of the single-spoke solver
//!
//! This module handles **one** spoke in isolation. For the
//! N-spoke ring of a three-quad decomposition, shared CPs at the
//! central vertex `C` and at the outer hex corners couple
//! adjacent spokes' constraints — this is the Phase C work. The
//! single-spoke routine here is the Phase B building block and
//! the basis for the per-spoke block of the Phase C global
//! system.

// Phase C (the N-spoke global coupler) is what calls most of
// these symbols in production. Until it lands, the dead-code
// warnings here are an expected consequence of Phase B landing
// as a self-contained, independently-tested unit.
#![allow(dead_code)]
// Index-style loops over 8-wide Bernstein-coefficient or 24-wide
// constraint tables read more clearly than iterator-based form
// because the index is load-bearing in the mathematics (each
// `k` is the Bernstein-basis index of the identity's LHS, each
// row `r` is a scalar constraint row, etc.). The iterator
// transformation `(k, entry)` here would obscure the algebraic
// correspondence.
#![allow(clippy::needless_range_loop)]

use crate::geometry::bernstein_poly::{BernsteinPoly, bernstein_product_scalar_vec3};
use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::gregory::Side;
use crate::geometry::{Point3, Vec3};

/// Bernstein coefficients of `xtan(t)` at the given side of a
/// two-stripe patch. For side `UZero`:
///
/// ```text
///   ∂S/∂u(0, t) = Σ_j 7 · (P[1][j] − P[0][j]) · B_j^7(t)
/// ```
///
/// so Bernstein coefficient `j` is `7 · (P[1][j] − P[0][j])`.
pub fn xtan_bernstein_cps(patch: &RectangularGregoryTwoStripe, side: Side) -> [Vec3; 8] {
    let mut out = [Vec3::zero(); 8];
    for j in 0..8 {
        let (near, far) = match side {
            Side::UZero => (patch.grid[1][j], patch.grid[0][j]),
            Side::UOne => (patch.grid[6][j], patch.grid[7][j]),
            Side::VZero => (patch.grid[j][1], patch.grid[j][0]),
            Side::VOne => (patch.grid[j][6], patch.grid[j][7]),
        };
        out[j] = Vec3::new(
            7.0 * (near.x - far.x),
            7.0 * (near.y - far.y),
            7.0 * (near.z - far.z),
        );
    }
    out
}

/// Boundary control points at the given side (the 8 CPs that
/// define the spoke curve itself).
pub fn boundary_bernstein_cps(patch: &RectangularGregoryTwoStripe, side: Side) -> [Point3; 8] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 8];
    for j in 0..8 {
        out[j] = match side {
            Side::UZero => patch.grid[0][j],
            Side::UOne => patch.grid[7][j],
            Side::VZero => patch.grid[j][0],
            Side::VOne => patch.grid[j][7],
        };
    }
    out
}

/// Boundary tangent `T(t) = C'(t)` as a degree-6 vector Bernstein
/// polynomial, derived from the 8 boundary CPs.
pub fn boundary_tangent_poly(boundary_cps: &[Point3; 8]) -> BernsteinPoly<Vec3> {
    // C(t) = Σ_j P_j · B_j^7(t) (degree 7 vector — viewing Point3
    // as displacement). C'(t) = Σ_{j=0..6} 7·(P_{j+1} − P_j)·B_j^6(t).
    let mut coeffs = Vec::with_capacity(7);
    for j in 0..7 {
        let p0 = boundary_cps[j];
        let p1 = boundary_cps[j + 1];
        coeffs.push(Vec3::new(
            7.0 * (p1.x - p0.x),
            7.0 * (p1.y - p0.y),
            7.0 * (p1.z - p0.z),
        ));
    }
    BernsteinPoly::new(coeffs)
}

/// Coefficient-level G1 identity for a single spoke, cached for
/// constraint assembly. `alpha = −1` is baked in (the standard
/// consistent-orientation case). For reversed-orientation spokes,
/// `b_reversed = true` triggers coefficient reversal on B's xtan.
///
/// With β degree 1, there are 8 vector Bernstein coefficient
/// equations; stacking their `{x, y, z}` scalar equations gives
/// 24 scalar rows.
pub struct SingleSpokePolyG1 {
    /// Which grid entry of A is the j-th free row-1 CP. Helps
    /// the caller map solver unknowns back onto the grid after
    /// the solve.
    pub a_row1_indices: [(usize, usize); 8],
    /// Same for B.
    pub b_row1_indices: [(usize, usize); 8],
    /// 8 Bernstein equations. `equations[k]` corresponds to the
    /// `k`-th Bernstein coefficient of the identity's LHS.
    pub equations: [SpokeEquation; 8],
}

/// One vector-equation row (8 per spoke). The equation form is
///
/// ```text
///   a_coef · P_A[a_row1_index[j_a]]
/// + b_coef · P_B[b_row1_index[j_b]]
/// + β_coef_0 · β_0 · T_vec_0
/// + β_coef_1 · β_1 · T_vec_1
/// = rhs   (known Vec3)
/// ```
///
/// Each `β_coef_i` is a scalar weight (from the Bernstein product
/// formula) and `T_vec_i` are Vec3 values (Bernstein coefficients
/// of `T(t)`). We keep the per-`β`-component product factored
/// this way so it stays unambiguous in the constraint-matrix
/// assembly.
#[derive(Clone, Copy, Debug)]
pub struct SpokeEquation {
    /// Grid index on A whose row-1 CP enters this equation (with
    /// coefficient `+7`).
    pub a_idx: (usize, usize),
    /// Grid index on B whose row-1 CP enters this equation (with
    /// coefficient `α · 7 = −7`).
    pub b_idx: (usize, usize),
    /// `β₀`'s Vec3 coefficient in this Bernstein equation.
    pub beta0_weight: Vec3,
    /// `β₁`'s Vec3 coefficient in this Bernstein equation.
    pub beta1_weight: Vec3,
    /// Known RHS (always zero in exact form — boundary
    /// cancellation — but left explicit for completeness and
    /// eventual degree-mismatch escalation).
    pub rhs: Vec3,
}

/// Build the `SingleSpokePolyG1` for a spoke given both patches'
/// geometries and sides. `b_reversed` controls coefficient-order
/// reversal on B. `α = −1` is baked in.
pub fn build_single_spoke_poly_g1(
    patch_a: &RectangularGregoryTwoStripe,
    side_a: Side,
    side_b: Side,
    b_reversed: bool,
) -> SingleSpokePolyG1 {
    let a_row1_indices = row1_indices_for_side(side_a);
    let b_row1_indices_raw = row1_indices_for_side(side_b);

    let boundary_a = boundary_bernstein_cps(patch_a, side_a);
    let t_poly = boundary_tangent_poly(&boundary_a);

    // β · T with β = B_0^1 = (1 − t) (takes β_0 = 1, β_1 = 0).
    let beta0 = BernsteinPoly::<f64>::new(vec![1.0, 0.0]);
    // β · T with β = B_1^1 = t.
    let beta1 = BernsteinPoly::<f64>::new(vec![0.0, 1.0]);

    let beta0_t = bernstein_product_scalar_vec3(&beta0, &t_poly);
    let beta1_t = bernstein_product_scalar_vec3(&beta1, &t_poly);

    let mut equations = [SpokeEquation {
        a_idx: (0, 0),
        b_idx: (0, 0),
        beta0_weight: Vec3::zero(),
        beta1_weight: Vec3::zero(),
        rhs: Vec3::zero(),
    }; 8];

    for k in 0..8 {
        let a_k = a_row1_indices[k];
        let b_k_raw = if b_reversed { 7 - k } else { k };
        let b_k = b_row1_indices_raw[b_k_raw];

        // RHS = 0 exactly, because P_A[0][k] and P_B[·][0] (the
        // boundary CPs) coincide as points on the shared spoke.
        // We leave the explicit RHS in place because a future
        // Phase D escalation (β degree > 1) may lift LHS degree
        // above 7 and the upper coefficients then need an
        // explicit boundary contribution.
        equations[k] = SpokeEquation {
            a_idx: a_k,
            b_idx: b_k,
            beta0_weight: beta0_t.coeff(k),
            beta1_weight: beta1_t.coeff(k),
            rhs: Vec3::zero(),
        };
    }

    SingleSpokePolyG1 {
        a_row1_indices,
        b_row1_indices: {
            let mut out = [(0usize, 0usize); 8];
            for k in 0..8 {
                let b_k_raw = if b_reversed { 7 - k } else { k };
                out[k] = b_row1_indices_raw[b_k_raw];
            }
            out
        },
        equations,
    }
}

/// Grid-index lookup for the 8 row-1-near-spoke CP positions on
/// the given side. Indexed by Bernstein basis index `j ∈ 0..8`.
///
/// **Caveat**: positions `j = 0` and `j = 7` sit on the TWO other
/// boundaries of the patch (the hex edges adjacent to the spoke's
/// endpoints at `t = 0` / `t = 1`). They are shared CPs of the
/// patch's other boundary curves — moving them breaks G0 with
/// the adjacent-spoke boundary. Callers should use
/// [`free_row1_indices_for_side`] instead, which returns only
/// the genuinely free interior positions `j ∈ 1..7`.
pub fn row1_indices_for_side(side: Side) -> [(usize, usize); 8] {
    let mut out = [(0, 0); 8];
    for j in 0..8 {
        out[j] = match side {
            Side::UZero => (1, j),
            Side::UOne => (6, j),
            Side::VZero => (j, 1),
            Side::VOne => (j, 6),
        };
    }
    out
}

/// Genuinely free row-1 near-spoke CPs: `j ∈ 1..7` (excluding the
/// boundary-adjacent endpoints `j = 0` and `j = 7`). These are
/// the positions the polynomial-identity ring solver is allowed
/// to move without breaking G0 with adjacent spokes.
pub fn free_row1_indices_for_side(side: Side) -> [(usize, usize); 6] {
    let mut out = [(0, 0); 6];
    for (k, j) in (1..7).enumerate() {
        out[k] = match side {
            Side::UZero => (1, j),
            Side::UOne => (6, j),
            Side::VZero => (j, 1),
            Side::VOne => (j, 6),
        };
    }
    out
}

/// Set an entry of a two-stripe grid. Private helper to localise
/// the `.grid[i][j] = …` assignment pattern.
fn grid_set(patch: &mut RectangularGregoryTwoStripe, idx: (usize, usize), p: Point3) {
    patch.grid[idx.0][idx.1] = p;
}

/// Read an entry of a two-stripe grid.
fn grid_get(patch: &RectangularGregoryTwoStripe, idx: (usize, usize)) -> Point3 {
    patch.grid[idx.0][idx.1]
}

/// Outcome of a single-spoke polynomial-identity solve.
#[derive(Clone, Debug)]
pub struct SpokePolyG1Solve {
    /// New values for A's row-1 near-spoke CPs (8 CPs).
    pub a_row1_new: [Point3; 8],
    /// New values for B's row-1 near-spoke CPs (8 CPs).
    pub b_row1_new: [Point3; 8],
    /// Blending polynomial Bernstein coefficients `(β_0, β_1)`.
    pub beta: (f64, f64),
    /// Post-solve residual: `max_k ‖LHS(k)‖₂` across the 8
    /// Bernstein equations. Should be at machine precision when
    /// the system is well-conditioned and has a solution.
    pub residual_max: f64,
}

/// Solve the polynomial-identity G1 constraint for a single
/// spoke in isolation. Minimum-norm solution from the current
/// `(patch_a, patch_b)` CP state. Returns `None` if the KKT
/// linear system is singular.
///
/// **Does not mutate** the input patches — the caller applies
/// [`SpokePolyG1Solve::a_row1_new`] / `b_row1_new` explicitly,
/// which lets the Phase C global solver batch multiple spokes'
/// updates after a combined solve.
pub fn solve_single_spoke_poly_g1(
    patch_a: &RectangularGregoryTwoStripe,
    side_a: Side,
    patch_b: &RectangularGregoryTwoStripe,
    side_b: Side,
    b_reversed: bool,
) -> Option<SpokePolyG1Solve> {
    // Auto-detect α from the initial geometry. When the two xtan
    // vectors point in the same geometric direction (dot > 0),
    // `α = −1` makes the identity `xtan_A − xtan_B = β·T`; when
    // they oppose (dot < 0), `α = +1` so `xtan_A + xtan_B = β·T`.
    // Empirically this is how the K-P form handles mirrored-side
    // junctions without needing to tag the orientation.
    let alpha = {
        let xtan_a_mid = xtan_bernstein_cps(patch_a, side_a)[3];
        let xtan_b_raw = xtan_bernstein_cps(patch_b, side_b);
        let xtan_b_mid = if b_reversed {
            xtan_b_raw[7 - 3]
        } else {
            xtan_b_raw[3]
        };
        let dot =
            xtan_a_mid.x * xtan_b_mid.x + xtan_a_mid.y * xtan_b_mid.y + xtan_a_mid.z * xtan_b_mid.z;
        if dot >= 0.0 { -1.0_f64 } else { 1.0_f64 }
    };
    let system = build_single_spoke_poly_g1(patch_a, side_a, side_b, b_reversed);

    // Unknowns layout:
    //   [0..24):  P_A[row1_A[0..8]].{x, y, z}  (8 CPs × 3 coords)
    //   [24..48): P_B[row1_B[0..8]].{x, y, z}
    //   [48..50): β_0, β_1
    const N_UNK: usize = 50;

    // Seed vector (current CP state + β = 0, 0).
    let mut seed = [0.0_f64; N_UNK];
    for k in 0..8 {
        let a_p = grid_get(patch_a, system.a_row1_indices[k]);
        seed[3 * k] = a_p.x;
        seed[3 * k + 1] = a_p.y;
        seed[3 * k + 2] = a_p.z;
    }
    for k in 0..8 {
        let b_p = grid_get(patch_b, system.b_row1_indices[k]);
        seed[24 + 3 * k] = b_p.x;
        seed[24 + 3 * k + 1] = b_p.y;
        seed[24 + 3 * k + 2] = b_p.z;
    }

    // 24 scalar constraints (8 Bernstein eqs × 3 coords).
    // Identity in inward xtan convention:
    //   7·A_free[k] + α·7·B_free[k] + (β·T)[k] = (1+α)·7·C[k].
    let mut a_mat = vec![vec![0.0_f64; N_UNK]; 24];
    let mut rhs = [0.0_f64; 24];
    let boundary_cps_a = boundary_bernstein_cps(patch_a, side_a);
    let rhs_factor = (1.0 + alpha) * 7.0;
    for k in 0..8 {
        let eq = &system.equations[k];
        let a_offset = 3 * k;
        let b_offset = 24 + 3 * k;
        let c_k = boundary_cps_a[k];
        for (c_idx, (w0, w1, c_c)) in [
            (0, (eq.beta0_weight.x, eq.beta1_weight.x, c_k.x)),
            (1, (eq.beta0_weight.y, eq.beta1_weight.y, c_k.y)),
            (2, (eq.beta0_weight.z, eq.beta1_weight.z, c_k.z)),
        ] {
            let row = 3 * k + c_idx;
            a_mat[row][a_offset + c_idx] = 7.0;
            a_mat[row][b_offset + c_idx] = alpha * 7.0;
            a_mat[row][48] = w0;
            a_mat[row][49] = w1;
            rhs[row] = rhs_factor * c_c;
        }
    }

    // Min-norm toward seed:
    //   minimise ‖x − seed‖²  s.t.  A·x = rhs
    //   => x = seed + Aᵀ·λ,   with A · Aᵀ · λ = rhs − A · seed
    let mut a_seed = [0.0_f64; 24];
    for r in 0..24 {
        for (c, &v) in a_mat[r].iter().enumerate() {
            a_seed[r] += v * seed[c];
        }
    }
    let rhs_adj: Vec<f64> = (0..24).map(|r| rhs[r] - a_seed[r]).collect();

    // Form AAᵀ (24 × 24).
    let mut aat = vec![vec![0.0_f64; 24]; 24];
    let mut trace = 0.0_f64;
    for r in 0..24 {
        for c in r..24 {
            let mut s = 0.0_f64;
            for k in 0..N_UNK {
                s += a_mat[r][k] * a_mat[c][k];
            }
            aat[r][c] = s;
            aat[c][r] = s;
        }
        trace += aat[r][r];
    }
    // Small Tikhonov for numerical safety.
    let lambda = 1e-12 * (trace / 24.0).max(1.0);
    for r in 0..24 {
        aat[r][r] += lambda;
    }
    let lambdas = solve_linear_system(aat, rhs_adj)?;

    // x = seed + Aᵀ · λ.
    let mut x = seed;
    for r in 0..24 {
        for (k, xk) in x.iter_mut().enumerate() {
            *xk += lambdas[r] * a_mat[r][k];
        }
    }

    // Unpack.
    let mut a_row1_new = [Point3::new(0.0, 0.0, 0.0); 8];
    let mut b_row1_new = [Point3::new(0.0, 0.0, 0.0); 8];
    for k in 0..8 {
        a_row1_new[k] = Point3::new(x[3 * k], x[3 * k + 1], x[3 * k + 2]);
        b_row1_new[k] = Point3::new(x[24 + 3 * k], x[24 + 3 * k + 1], x[24 + 3 * k + 2]);
    }
    let beta = (x[48], x[49]);

    // Post-solve residual.
    let mut residual_max = 0.0_f64;
    for k in 0..8 {
        let mut lhs_val = [0.0_f64; 3];
        for c in 0..3 {
            let row = 3 * k + c;
            let mut acc = -rhs[row];
            for (j, &v) in a_mat[row].iter().enumerate() {
                acc += v * x[j];
            }
            lhs_val[c] = acc;
        }
        let mag = (lhs_val[0].powi(2) + lhs_val[1].powi(2) + lhs_val[2].powi(2)).sqrt();
        if mag > residual_max {
            residual_max = mag;
        }
    }

    Some(SpokePolyG1Solve {
        a_row1_new,
        b_row1_new,
        beta,
        residual_max,
    })
}

/// Apply a single-spoke solve result back to the patches in place.
pub fn apply_single_spoke_solve(
    patch_a: &mut RectangularGregoryTwoStripe,
    patch_b: &mut RectangularGregoryTwoStripe,
    side_a: Side,
    side_b: Side,
    b_reversed: bool,
    result: &SpokePolyG1Solve,
) {
    let a_indices = row1_indices_for_side(side_a);
    let b_indices_raw = row1_indices_for_side(side_b);
    for k in 0..8 {
        grid_set(patch_a, a_indices[k], result.a_row1_new[k]);
        let b_k_raw = if b_reversed { 7 - k } else { k };
        grid_set(patch_b, b_indices_raw[b_k_raw], result.b_row1_new[k]);
    }
}

use super::spoke_poly_g1_linalg::solve_linear_system;

#[cfg(test)]
#[path = "spoke_poly_g1_solver_tests.rs"]
mod tests;
