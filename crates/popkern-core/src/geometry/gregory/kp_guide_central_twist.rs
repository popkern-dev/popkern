//! Phase H — Shirman-Séquin twist computation at the central
//! vertex `C` for the bicubic guide surface.
//!
//! Given the C-corner tangent data from
//! [`extract_central_vertex_hermite`](super::extract_central_vertex_hermite)
//! (with tangents already projected onto the common tangent
//! plane at `C`), this module computes the **mixed partial**
//! `∂²G_k/∂u∂v` at `C` for every sub-patch, subject to the
//! Shirman-Séquin N-way G1-along-spoke constraint.
//!
//! # Why proper twist matters
//!
//! Without twist (twist = 0), the bicubic guide's inner CP at
//! `c1` is placed at the parallelogram-rule position — which
//! only guarantees G1 at the spoke endpoint. For G1 to hold
//! along the **entire** shared spoke between adjacent guides,
//! the inner CPs must satisfy a Lagrangian-QP constraint (the
//! classical Shirman-Séquin 1987 N-way enclosure). This module
//! solves that QP for the bicubic-guide case, then derives the
//! mixed partial value from the solved inner CP via the
//! standard Hermite identity.
//!
//! # Algorithm
//!
//! Per sub-patch `k`:
//!
//! 1. Compute the parallelogram-rule inner-CP seed at `c1`:
//!
//!    ```text
//!      grid[2][1]_seed = C + tangent_v_k / 3 − tangent_u_k / 3
//!    ```
//!
//! 2. Project all N seeds onto the 2-D common tangent plane with
//!    normal `n_C` (already established by
//!    [`solve_central_vertex_corner`](super::solve_central_vertex_corner)).
//!
//! 3. Solve the N×N circulant Lagrangian QP for the multipliers
//!    `λ_k` that enforce `(Y_k − Y_{k+1}) · p_k = 0` at every
//!    shared spoke `k`, where `p_k` is the 2-D spoke-
//!    perpendicular direction and `Y_k` the adjusted inner CP.
//!
//! 4. Lift the adjusted inner CPs back to 3-D.
//!
//! 5. Derive the mixed partial at `c1` via the Hermite-bicubic
//!    inversion:
//!
//!    ```text
//!      grid[2][1]_k = C − tu_k / 3 + tv_k / 3 − twist_k / 9
//!      ⟹ twist_k = 9·(C − grid[2][1]_k) − 3·tu_k + 3·tv_k
//!    ```
//!
//! The Lagrangian-QP core is identical to
//! [`solve_central_vertex_corner`](super::solve_central_vertex_corner)
//! — the math is degree-agnostic; only the parallelogram-seed
//! formula and the final inverse derivation differ.

use super::direct_errors::DirectConstructionError;
use super::direct_input::DirectBoundaryInput;
use super::direct_linalg::{orthonormal_tangent_basis, solve_dense};
use super::direct_options::DirectConstructionOptions;
use super::direct_types::CornerKind;
use super::kp_guide_build::GuideCornerData;
use super::kp_guide_central_hermite::extract_central_vertex_hermite;
use crate::geometry::{Point3, Vec3};

/// Return [`GuideCornerData`] for every sub-patch's `c1` corner
/// with a **properly Shirman-Séquin-solved twist** instead of
/// zero.
///
/// On success the returned Hermite data has identical `position`
/// and tangent fields to
/// [`extract_central_vertex_hermite`](super::extract_central_vertex_hermite)
/// (the twist correction doesn't change tangents — it only adds
/// the mixed-partial value). The `twist` field is populated with
/// the G1-along-spoke-correct value.
pub fn extract_central_vertex_hermite_with_twist(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec<GuideCornerData>, DirectConstructionError> {
    // Delegate to the tangent-only extractor for position,
    // tangent_u, tangent_v (common-tangent-plane projection
    // already applied there).
    let base = extract_central_vertex_hermite(input, options)?;
    let n = base.len();
    let c = input.central_vertex;

    // Common tangent plane normal — same as inside the base
    // extractor. Derive locally to avoid re-running the full
    // Shirman-Séquin solve (the common normal is cheap).
    let n_c = common_normal(input, options)?;

    // 1. Parallelogram seeds for each bicubic inner CP at c1.
    let mut raw_parallelogram: Vec<Point3> = Vec::with_capacity(n);
    for bk in &base {
        raw_parallelogram.push(Point3::new(
            c.x + bk.tangent_v.x / 3.0 - bk.tangent_u.x / 3.0,
            c.y + bk.tangent_v.y / 3.0 - bk.tangent_u.y / 3.0,
            c.z + bk.tangent_v.z / 3.0 - bk.tangent_u.z / 3.0,
        ));
    }

    // 2. Per-spoke direction at C (into the tangent plane). For
    //    the bicubic guide, at c1 the "shared spoke" with
    //    sub-patch k+1 is guide_k's u = 1 edge, whose tangent
    //    at v = 0 is `tangent_v_k`. So σ_k (the 3-D spoke
    //    direction between sub-patches k and k+1 as seen from
    //    sub-patch k) is `tangent_v_k`.
    let spoke_directions: Vec<Vec3> = base.iter().map(|b| b.tangent_v).collect();

    // 3. Lagrangian-QP: identical to solve_central_vertex_corner
    //    but taking the bicubic seeds as input.
    let inner_cps = solve_lagrangian_qp(c, n_c, &raw_parallelogram, &spoke_directions, options)?;

    // 4. Derive twist per sub-patch and return full Hermite.
    let mut out = Vec::with_capacity(n);
    for (k, bk) in base.iter().enumerate() {
        let inner = inner_cps[k];
        // twist = 9·(C − inner) − 3·tu + 3·tv
        let twist = Vec3::new(
            9.0 * (c.x - inner.x) - 3.0 * bk.tangent_u.x + 3.0 * bk.tangent_v.x,
            9.0 * (c.y - inner.y) - 3.0 * bk.tangent_u.y + 3.0 * bk.tangent_v.y,
            9.0 * (c.z - inner.z) - 3.0 * bk.tangent_u.z + 3.0 * bk.tangent_v.z,
        );
        out.push(GuideCornerData {
            position: bk.position,
            tangent_u: bk.tangent_u,
            tangent_v: bk.tangent_v,
            twist,
        });
    }
    Ok(out)
}

/// Compute the common tangent-plane normal at `C` as the
/// normalised sum of the input-face tangent-plane normals. Same
/// formula as inside
/// [`solve_central_vertex_corner`](super::solve_central_vertex_corner);
/// duplicated here to avoid a second full QP solve.
fn common_normal(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec3, DirectConstructionError> {
    let n = input.n();
    let mut sum = Vec3::zero();
    for tps in &input.central_vertex_tps {
        sum = sum + tps.normal();
    }
    let len = sum.length();
    let threshold = options.corner_singularity_eps.max(1e-6) * (n as f64);
    if len < threshold {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: len,
            tolerance: threshold,
            edge_hint: "input-face normals at central vertex sum to ~0".into(),
        });
    }
    sum.normalized()
        .ok_or_else(|| DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: len,
            tolerance: options.corner_singularity_eps,
            edge_hint: "averaged normal not normalisable".into(),
        })
}

/// Solve the circulant Lagrangian QP for the adjusted inner CPs
/// given parallelogram seeds and spoke directions.
///
/// This is a degree-agnostic port of the core algorithm inside
/// [`solve_central_vertex_corner`](super::solve_central_vertex_corner).
/// Both the biseptic and bicubic solves consume the same
/// primitive — only the seed computation and the output
/// conversion differ.
fn solve_lagrangian_qp(
    c: Point3,
    n_c: Vec3,
    raw_parallelogram: &[Point3],
    spoke_directions: &[Vec3],
    options: &DirectConstructionOptions,
) -> Result<Vec<Point3>, DirectConstructionError> {
    let n = raw_parallelogram.len();
    let (t1, t2) = orthonormal_tangent_basis(n_c);

    // 2-D projections of the seeds.
    let raw_2d: Vec<[f64; 2]> = raw_parallelogram
        .iter()
        .map(|p| {
            let off = Vec3::new(p.x - c.x, p.y - c.y, p.z - c.z);
            [off.dot(&t1), off.dot(&t2)]
        })
        .collect();

    // 2-D spoke perpendiculars p_k (rotate σ_k 90° CCW).
    let mut p: Vec<[f64; 2]> = Vec::with_capacity(n);
    let mut active: Vec<bool> = Vec::with_capacity(n);
    for sigma in spoke_directions {
        let sx = sigma.dot(&t1);
        let sy = sigma.dot(&t2);
        let mag = (sx * sx + sy * sy).sqrt();
        if mag < options.corner_singularity_eps {
            p.push([0.0, 0.0]);
            active.push(false);
        } else {
            p.push([-sy / mag, sx / mag]);
            active.push(true);
        }
    }
    if !active.iter().any(|&a| a) {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: 0.0,
            tolerance: options.corner_singularity_eps,
            edge_hint: "all spokes project to zero in central tangent plane".into(),
        });
    }

    // N × N circulant tridiagonal λ-system.
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut rhs = vec![0.0_f64; n];
    for k in 0..n {
        if !active[k] {
            matrix[k][k] = 1.0;
            continue;
        }
        let k_prev = (k + n - 1) % n;
        let k_next = (k + 1) % n;
        matrix[k][k] = 2.0;
        if active[k_prev] {
            matrix[k][k_prev] -= p[k_prev][0] * p[k][0] + p[k_prev][1] * p[k][1];
        }
        if active[k_next] {
            matrix[k][k_next] -= p[k_next][0] * p[k][0] + p[k_next][1] * p[k][1];
        }
        let diff = [
            raw_2d[k][0] - raw_2d[k_next][0],
            raw_2d[k][1] - raw_2d[k_next][1],
        ];
        rhs[k] = diff[0] * p[k][0] + diff[1] * p[k][1];
    }

    let Some(lambdas) = solve_dense(matrix, rhs) else {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: 0.0,
            tolerance: options.corner_singularity_eps,
            edge_hint: "Lagrangian QP λ-system singular".into(),
        });
    };

    // Recover adjusted Y_k in 2-D then lift to 3-D.
    let mut inner_cps = Vec::with_capacity(n);
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let mut yx = raw_2d[k][0];
        let mut yy = raw_2d[k][1];
        if active[k] {
            yx -= lambdas[k] * p[k][0];
            yy -= lambdas[k] * p[k][1];
        }
        if active[k_prev] {
            yx += lambdas[k_prev] * p[k_prev][0];
            yy += lambdas[k_prev] * p[k_prev][1];
        }
        inner_cps.push(Point3::new(
            c.x + yx * t1.x + yy * t2.x,
            c.y + yx * t1.y + yy * t2.y,
            c.z + yx * t1.z + yy * t2.z,
        ));
    }
    Ok(inner_cps)
}

#[cfg(test)]
#[path = "kp_guide_central_twist_tests.rs"]
mod tests;
