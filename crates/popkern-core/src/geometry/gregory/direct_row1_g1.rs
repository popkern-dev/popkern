//! Direct-construction: the row-1 G1 formula for a single
//! shared spoke.
//!
//! This is the **core math** of Phase G. With the blending
//! polynomials `α(t), β(t)` pre-determined by
//! [`super::direct_blend_coeffs::compute_ring_blend_coeffs`] (via
//! the corner tangent-plane geometry), the Chiyokura-Kimura G1
//! polynomial identity
//!
//! ```text
//!   xtan_A(t) + α(t) · xtan_B(t) + β(t) · T(t) ≡ 0
//! ```
//!
//! becomes a **linear system** in the row-1 near-spoke CPs of
//! both patches. Setting every Bernstein coefficient of the
//! LHS (degree 8 after elevation of the deg-7 terms) to zero
//! yields 9 vector equations per spoke — exact G1 along the
//! whole spoke, not just at sample points.
//!
//! # Sources
//!
//! * Chiyokura, H. & Kimura, F. (1983). "Design of Solids with
//!   Free-form Surfaces". *Computer Graphics* 17(3), §3.2.
//! * Farin, G. (1993). "Curves and Surfaces for CAGD", Ch. 18
//!   §4.
//! * Peters, J. (2000). Two-stripe bidegree-7 extension —
//!   provides extra row-1 DOFs at the plain-Bezier
//!   corner-adjacent positions.
//!
//! All pre-1990 or contemporaneous academic (no vendor
//! patent assignments known).

#![allow(clippy::needless_range_loop)]
//!
//! # Degree bookkeeping
//!
//! Bidegree-7 two-stripe:
//!
//! | Term       | Native degree | Coefficient structure                       |
//! |------------|---------------|---------------------------------------------|
//! | `xtan_A`   | 7             | linear in unknown row-1 CPs on side A       |
//! | `xtan_B`   | 7             | linear in unknown row-1 CPs on side B       |
//! | `T`        | 6             | all known (from spoke boundary CPs)         |
//! | `α·xtan_B` | 8             | linear in unknowns (α is linear scalar)     |
//! | `β·T`      | 7             | all known                                   |
//!
//! Elevating `xtan_A` (7 → 8) and `β·T` (7 → 8) puts
//! everything on the common degree-8 Bernstein basis: **9
//! coefficient equations per spoke**.

use crate::geometry::bernstein_poly::{BernsteinPoly, bernstein_product_scalar_vec3};
use crate::geometry::gregory::direct_blend_coeffs::SpokeBlendCoeffs;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::SpokeBoundaryData;
use crate::geometry::gregory::direct_linalg::{binomial, solve_dense};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::{Point3, Vec3};

/// Inputs to [`solve_single_spoke_row1_g1`]. Both row-1 arrays
/// are supplied in **spoke-t order**: index `k` on each side
/// corresponds to Bernstein-basis index `k` of the spoke's
/// canonical parameter. The caller is responsible for
/// reversing patch B's grid indices before packaging
/// (patch B's VZero side uses `b_reversed = true` in the
/// three-quad convention).
#[derive(Debug, Clone)]
pub struct SingleSpokeRow1Input {
    /// Boundary data for this spoke (curve CPs + corner xtans).
    pub spoke: SpokeBoundaryData,
    /// Blend coefficients from [`SpokeBlendCoeffs`].
    pub blend_coeffs: SpokeBlendCoeffs,
    /// 8 row-1 CPs on patch A near its UOne boundary, in
    /// spoke-t order. `[0]` and `[7]` are pinned boundary CPs;
    /// `[1..=6]` are seed values that the solve refines.
    pub patch_a_row1: [Point3; 8],
    /// 8 row-1 CPs on patch B near its VZero boundary (with
    /// the b_reversed flip already applied), in spoke-t order.
    pub patch_b_row1: [Point3; 8],
}

/// Output of [`solve_single_spoke_row1_g1`]. The two row-1
/// arrays carry only the free indices `1..=6`; the caller
/// preserves the pinned `[0]` / `[7]` entries.
#[derive(Debug, Clone)]
pub struct SingleSpokeRow1Output {
    /// Updated row-1 CPs on side A for `j ∈ 1..=6`.
    pub patch_a_row1_free: [Point3; 6],
    /// Updated row-1 CPs on side B for `k' ∈ 1..=6`
    /// (spoke-t-indexed).
    pub patch_b_row1_free: [Point3; 6],
    /// Max `‖LHS(k)‖₂` across the 9 Bernstein coefficient
    /// equations. Machine precision when the input is
    /// G1-consistent.
    pub residual_max: f64,
}

/// Solve the row-1 G1 direct formula for ONE spoke.
///
/// Uses the Bernstein-polynomial product + degree-elevation
/// machinery from [`crate::geometry::bernstein_poly`] to form
/// the LHS polynomial symbolically (coefficients expressed as
/// linear-in-unknowns vectors), then assembles and solves the
/// min-norm KKT linear system.
///
/// For a SINGLE spoke in isolation, each side's free row-1 CPs
/// are independent of the other spokes' CPs, so there's no
/// global coupling here — that's the job of the ring-level
/// caller in a future commit.
pub fn solve_single_spoke_row1_g1(
    input: &SingleSpokeRow1Input,
    options: &DirectConstructionOptions,
) -> Result<SingleSpokeRow1Output, DirectConstructionError> {
    // 1. Build the constant polynomials (boundary tangent T,
    //    α·xtan_B constant parts, β·T) and split each xtan into
    //    "constant" + "linear-in-unknowns" components.
    let q = &input.spoke.boundary_cps;
    let bc = input.blend_coeffs;

    // xtan_A: 8 Bernstein Vec3 coefficients (degree 7).
    // Constant part (includes the boundary CPs Q_A[j'] always,
    // and the PINNED patch_a_row1 values at j'=0, 7).
    let mut xtan_a_const = [Vec3::zero(); 8];
    for j in 0..8 {
        let q_j = q[j];
        let row1_j = input.patch_a_row1[j];
        xtan_a_const[j] = Vec3::new(
            7.0 * (q_j.x - row1_j.x),
            7.0 * (q_j.y - row1_j.y),
            7.0 * (q_j.z - row1_j.z),
        );
    }
    // Zero-out the UNKNOWN contributions at j' ∈ 1..=6 so we
    // can recover the "const without the unknown row-1 part"
    // later — we treat the unknowns' -7 coefficients separately.
    for j in 1..=6 {
        let q_j = q[j];
        xtan_a_const[j] = Vec3::new(7.0 * q_j.x, 7.0 * q_j.y, 7.0 * q_j.z);
    }

    // Same for xtan_B (sign convention: +7·row1 - 7·Q).
    let mut xtan_b_const = [Vec3::zero(); 8];
    for k in 0..8 {
        let q_k = q[k];
        let row1_k = input.patch_b_row1[k];
        xtan_b_const[k] = Vec3::new(
            7.0 * (row1_k.x - q_k.x),
            7.0 * (row1_k.y - q_k.y),
            7.0 * (row1_k.z - q_k.z),
        );
    }
    for k in 1..=6 {
        let q_k = q[k];
        xtan_b_const[k] = Vec3::new(-7.0 * q_k.x, -7.0 * q_k.y, -7.0 * q_k.z);
    }

    // T(t): degree-6 vector Bernstein poly from spoke boundary
    // CPs. T(t) = 7·Σ_k (Q[k+1] − Q[k])·B_k^6(t).
    let mut t_coeffs = Vec::with_capacity(7);
    for k in 0..7 {
        let a = q[k];
        let b = q[k + 1];
        t_coeffs.push(Vec3::new(
            7.0 * (b.x - a.x),
            7.0 * (b.y - a.y),
            7.0 * (b.z - a.z),
        ));
    }
    let t_poly = BernsteinPoly::<Vec3>::new(t_coeffs);

    // α(t), β(t): linear scalar polys.
    let alpha_poly = BernsteinPoly::<f64>::new(vec![bc.alpha_0, bc.alpha_1]);
    let beta_poly = BernsteinPoly::<f64>::new(vec![bc.beta_0, bc.beta_1]);

    // β·T: degree-7 vector poly (all constant).
    let beta_t = bernstein_product_scalar_vec3(&beta_poly, &t_poly);

    // 2. Form the LHS constant polynomial at degree 8.
    let xtan_a_const_poly = BernsteinPoly::<Vec3>::new(xtan_a_const.to_vec()).degree_elevate(1); // 7 → 8
    let xtan_b_const_poly = BernsteinPoly::<Vec3>::new(xtan_b_const.to_vec());
    let alpha_xtan_b_const = bernstein_product_scalar_vec3(&alpha_poly, &xtan_b_const_poly); // deg 8
    let beta_t_elev = beta_t.degree_elevate(1); // 7 → 8

    let lhs_const = xtan_a_const_poly
        .add_elevated(&alpha_xtan_b_const)
        .add_elevated(&beta_t_elev);
    debug_assert_eq!(lhs_const.degree(), 8);

    // 3. Per-unknown contribution weights at degree 8.
    // Side A unknown U_A[j] contributes -7 at degree-7 xtan_A
    // index j, which degree-elevates to xtan_A at degree 8:
    //   -7·(8-j)/8 at index j, -7·(j+1)/8 at index j+1.
    // Side B unknown U_B[k'] contributes +7 at degree-7 xtan_B
    // index k'; after multiplying by α (linear scalar):
    //   +7·(C(7,k')/C(8,k'))·α_0 at degree-8 index k',
    //   +7·(C(7,k')/C(8,k'+1))·α_1 at degree-8 index k'+1.
    let weights_a: [[f64; 9]; 6] = {
        let mut out = [[0.0_f64; 9]; 6];
        for (idx, j) in (1..=6).enumerate() {
            // U_A[j] contributes -7 at degree-7 index j of xtan_A.
            // Elevate: -7 at degree-8 index j with weight (8-j)/8,
            //          -7 at degree-8 index j+1 with weight (j+1)/8.
            let w_low = (8 - j) as f64 / 8.0;
            let w_high = (j + 1) as f64 / 8.0;
            out[idx][j] = -7.0 * w_low;
            out[idx][j + 1] = -7.0 * w_high;
        }
        out
    };

    let weights_b: [[f64; 9]; 6] = {
        let mut out = [[0.0_f64; 9]; 6];
        for (idx, kp) in (1..=6).enumerate() {
            // U_B[k'] contributes +7 at degree-7 index k' of xtan_B.
            // After α·xtan_B: appears at degree-8 indices k' and k'+1.
            let c_7_kp = binomial(7, kp);
            let c_8_kp = binomial(8, kp);
            let c_8_kpp1 = binomial(8, kp + 1);
            let w_at_kp = 7.0 * (c_7_kp / c_8_kp) * bc.alpha_0;
            let w_at_kpp1 = 7.0 * (c_7_kp / c_8_kpp1) * bc.alpha_1;
            out[idx][kp] = w_at_kp;
            out[idx][kp + 1] = w_at_kpp1;
        }
        out
    };

    // 4. Assemble the linear system. Unknowns flatten to 36
    // scalars (12 Point3 × 3 coords). Rows: 9 Bernstein × 3
    // coords = 27.
    const N_UNK: usize = 36;
    const N_ROWS: usize = 27;
    let mut a_mat = vec![vec![0.0_f64; N_UNK]; N_ROWS];
    let mut rhs = vec![0.0_f64; N_ROWS];
    for k in 0..9 {
        let coeff = lhs_const.coeff(k);
        for c in 0..3 {
            let row = 3 * k + c;
            for idx_a in 0..6 {
                a_mat[row][3 * idx_a + c] = weights_a[idx_a][k];
            }
            for idx_b in 0..6 {
                a_mat[row][18 + 3 * idx_b + c] = weights_b[idx_b][k];
            }
            rhs[row] = match c {
                0 => -coeff.x,
                1 => -coeff.y,
                _ => -coeff.z,
            };
        }
    }

    // Seed vector: the biseptic-promoted input row-1 values.
    let mut seed = [0.0_f64; N_UNK];
    for (idx_a, j) in (1..=6).enumerate() {
        let p = input.patch_a_row1[j];
        seed[3 * idx_a] = p.x;
        seed[3 * idx_a + 1] = p.y;
        seed[3 * idx_a + 2] = p.z;
    }
    for (idx_b, k) in (1..=6).enumerate() {
        let p = input.patch_b_row1[k];
        seed[18 + 3 * idx_b] = p.x;
        seed[18 + 3 * idx_b + 1] = p.y;
        seed[18 + 3 * idx_b + 2] = p.z;
    }

    // 5. Min-norm KKT solve.
    let mut a_seed = [0.0_f64; N_ROWS];
    for r in 0..N_ROWS {
        let mut s = 0.0_f64;
        for (k, v) in a_mat[r].iter().enumerate() {
            s += v * seed[k];
        }
        a_seed[r] = s;
    }
    let rhs_adj: Vec<f64> = (0..N_ROWS).map(|r| rhs[r] - a_seed[r]).collect();

    let mut aat = vec![vec![0.0_f64; N_ROWS]; N_ROWS];
    let mut trace = 0.0_f64;
    for r in 0..N_ROWS {
        for c in r..N_ROWS {
            let mut s = 0.0_f64;
            for k in 0..N_UNK {
                s += a_mat[r][k] * a_mat[c][k];
            }
            aat[r][c] = s;
            aat[c][r] = s;
        }
        trace += aat[r][r];
    }
    let lambda_reg = 1e-12 * (trace / (N_ROWS as f64)).max(1.0);
    for r in 0..N_ROWS {
        aat[r][r] += lambda_reg;
    }
    let Some(lambdas) = solve_dense(aat, rhs_adj) else {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: crate::geometry::gregory::direct_types::CornerKind::Central,
            smallest_sv: 0.0,
            tolerance: options.corner_singularity_eps,
            edge_hint: format!("row-1 KKT singular on edge {}", input.spoke.edge_id),
        });
    };

    let mut x = seed;
    for r in 0..N_ROWS {
        for k in 0..N_UNK {
            x[k] += lambdas[r] * a_mat[r][k];
        }
    }

    // 6. Unpack + residual.
    let mut patch_a_row1_free = [Point3::new(0.0, 0.0, 0.0); 6];
    for idx_a in 0..6 {
        patch_a_row1_free[idx_a] = Point3::new(x[3 * idx_a], x[3 * idx_a + 1], x[3 * idx_a + 2]);
    }
    let mut patch_b_row1_free = [Point3::new(0.0, 0.0, 0.0); 6];
    for idx_b in 0..6 {
        patch_b_row1_free[idx_b] = Point3::new(
            x[18 + 3 * idx_b],
            x[18 + 3 * idx_b + 1],
            x[18 + 3 * idx_b + 2],
        );
    }

    let mut residual_max = 0.0_f64;
    for k in 0..9 {
        let mut lhs = [0.0_f64; 3];
        for c in 0..3 {
            let row = 3 * k + c;
            let mut v = -rhs[row];
            for (j, coef) in a_mat[row].iter().enumerate() {
                v += coef * x[j];
            }
            lhs[c] = v;
        }
        let mag = (lhs[0].powi(2) + lhs[1].powi(2) + lhs[2].powi(2)).sqrt();
        if mag > residual_max {
            residual_max = mag;
        }
    }

    Ok(SingleSpokeRow1Output {
        patch_a_row1_free,
        patch_b_row1_free,
        residual_max,
    })
}

#[cfg(test)]
#[path = "direct_row1_g1_tests.rs"]
mod tests;
