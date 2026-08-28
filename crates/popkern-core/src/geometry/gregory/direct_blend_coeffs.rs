//! Direct-construction: compute the blending-polynomial
//! coefficients `α(t)` and `β(t)` from corner data.
//!
//! The G1 polynomial identity across a shared spoke is
//!
//! ```text
//!   xtan_A(t) + α(t) · xtan_B(t) + β(t) · T(t) ≡ 0 ,
//! ```
//!
//! where `xtan_A`, `xtan_B` are cross-tangent polynomials and
//! `T(t) = C'(t)` is the boundary tangent. A polynomial-identity direct-formula attempt tried to
//! use `α = constant` and failed: on asymmetric geometry the
//! coplanarity ratio varies with `t`.
//!
//! **The key Phase-G insight** (per Chiyokura-Kimura 1983 §3.2
//! and Farin 1993 Ch. 18): `α(t)` and `β(t)` become linear
//! polynomials whose endpoint values are **determined from the
//! corner tangent-plane geometry** — not solver outputs. Once
//! `α`, `β` are pinned, the remaining row-1 CPs satisfy a
//! linear system with an exact closed-form solution.
//!
//! This module computes `(α_0, α_1, β_0, β_1)` per spoke from
//! the corner data in the input. Output is consumed by the
//! row-1 G1 direct-formula module (row-1 G1 direct formulas).
//!
//! # Endpoint formula
//!
//! At `t = 0` (central vertex `C`):
//!
//!   `xtan_A(0) + α_0 · xtan_B(0) + β_0 · T(0) = 0`
//!
//! where all three vectors are known from input geometry
//! (G1-at-`C` ensures coplanarity, so there's a unique
//! `(α_0, β_0)` in the spanning coefficients). Solve the 3-D
//! vector equation as a 2×2 linear system in `(α_0, β_0)` by
//! projecting onto any two independent directions perpendicular
//! to the coplanar-sum null direction.
//!
//! Symmetric at `t = 1`.

use crate::geometry::Vec3;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::gregory::direct_types::CornerKind;

/// Per-spoke blending polynomial coefficients.
///
/// `α(t)` and `β(t)` are linear in Bernstein-basis form:
///   `α(t) = α_0 · B_0^1(t) + α_1 · B_1^1(t) = α_0·(1−t) + α_1·t`
///   `β(t) = β_0 · (1−t) + β_1 · t`
#[derive(Debug, Clone, Copy)]
pub struct SpokeBlendCoeffs {
    /// `α(0)`.
    pub alpha_0: f64,
    /// `α(1)`.
    pub alpha_1: f64,
    /// `β(0)`.
    pub beta_0: f64,
    /// `β(1)`.
    pub beta_1: f64,
    /// Max of the 2×2 endpoint-system conditioning numbers.
    /// Useful for `near_degeneracy_score` telemetry.
    pub condition_number: f64,
}

/// Compute [`SpokeBlendCoeffs`] for every spoke in the ring.
///
/// Output order matches [`DirectBoundaryInput::spokes`].
pub fn compute_ring_blend_coeffs(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec<SpokeBlendCoeffs>, DirectConstructionError> {
    let n = input.n();
    let mut out = Vec::with_capacity(n);
    for (k, spoke) in input.spokes.iter().enumerate() {
        // Boundary tangent at t=0 and t=1, from the spoke's
        // 8-CP Bezier: T(0) = 7·(P[1] − P[0]), T(1) = 7·(P[7] − P[6]).
        let cps = &spoke.boundary_cps;
        let t_at_0 = Vec3::new(
            7.0 * (cps[1].x - cps[0].x),
            7.0 * (cps[1].y - cps[0].y),
            7.0 * (cps[1].z - cps[0].z),
        );
        let t_at_1 = Vec3::new(
            7.0 * (cps[7].x - cps[6].x),
            7.0 * (cps[7].y - cps[6].y),
            7.0 * (cps[7].z - cps[6].z),
        );

        // Endpoint coefficient solve at t=0 (central vertex).
        let (alpha_0, beta_0, cond_0) = solve_endpoint_coeffs(
            spoke.xtan_a_at_c,
            spoke.xtan_b_at_c,
            t_at_0,
            options,
            k,
            CornerKind::Central,
        )?;
        // Endpoint coefficient solve at t=1 (outer vertex).
        let (alpha_1, beta_1, cond_1) = solve_endpoint_coeffs(
            spoke.xtan_a_at_v,
            spoke.xtan_b_at_v,
            t_at_1,
            options,
            k,
            CornerKind::OuterHex,
        )?;

        out.push(SpokeBlendCoeffs {
            alpha_0,
            alpha_1,
            beta_0,
            beta_1,
            condition_number: cond_0.max(cond_1),
        });
    }
    Ok(out)
}

/// Solve the endpoint equation `xtan_A + α·xtan_B + β·T = 0` for
/// `(α, β)` given the three vectors. Uses least-squares on the
/// 3×2 system, returning the conditioning number so the caller
/// can flag near-degenerate corners.
fn solve_endpoint_coeffs(
    xtan_a: Vec3,
    xtan_b: Vec3,
    t_vec: Vec3,
    options: &DirectConstructionOptions,
    spoke_idx: usize,
    corner_kind: CornerKind,
) -> Result<(f64, f64, f64), DirectConstructionError> {
    // Linear system:
    //   [ xtan_B.x  T.x ] [α]   [ -xtan_A.x ]
    //   [ xtan_B.y  T.y ] [β] = [ -xtan_A.y ]
    //   [ xtan_B.z  T.z ]       [ -xtan_A.z ]
    //
    // 3 equations, 2 unknowns. For G1-consistent geometry, the
    // system is rank-2 and exactly solvable. Solve by normal
    // equations:
    //   AᵀA · x = Aᵀ b
    //
    // AᵀA is 2×2 SPD:
    //   [ |xtan_B|²   xtan_B·T ]
    //   [ xtan_B·T    |T|²    ]
    let a11 = xtan_b.dot(&xtan_b);
    let a12 = xtan_b.dot(&t_vec);
    let a22 = t_vec.dot(&t_vec);

    let b1 = -xtan_a.dot(&xtan_b);
    let b2 = -xtan_a.dot(&t_vec);

    let det = a11 * a22 - a12 * a12;
    let scale = (a11 + a22).max(1e-30);
    let cond = (det / (scale * scale)).abs();

    if det.abs() < options.corner_singularity_eps * scale * scale {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: spoke_idx,
            kind: corner_kind,
            smallest_sv: det.abs().sqrt(),
            tolerance: options.corner_singularity_eps * scale,
            edge_hint: format!("spoke {spoke_idx} endpoint xtan/T linearly dependent"),
        });
    }

    let alpha = (a22 * b1 - a12 * b2) / det;
    let beta = (a11 * b2 - a12 * b1) / det;
    Ok((alpha, beta, cond))
}

#[cfg(test)]
#[path = "direct_blend_coeffs_tests.rs"]
mod tests;
