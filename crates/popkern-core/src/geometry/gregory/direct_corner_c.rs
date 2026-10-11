//! Direct-construction: the central-vertex (`C`) corner.
//!
//! For a ring of `N` two-stripe sub-patches meeting at a shared
//! vertex `C`, the **corner-adjacent inner CP** of each
//! sub-patch (for patch `k`, this is `grid[6][1]`) controls the
//! tangent plane at `C` along with the surrounding boundary CPs.
//! This module emits those `N` inner CPs in closed form so that
//! **all N tangent planes at `C` agree** — the Shirman-Séquin
//! vertex-enclosure condition from the direct-construction
//! literature.
//!
//! # Source
//!
//! Algorithm: Shirman & Séquin (1987). "Local surface
//! interpolation with Bézier patches". *Computer Aided
//! Geometric Design* 4(4), extended to N-fold corners per
//! Chiyokura-Kimura (1983) §3.3. The implementation mirrors
//! the existing version in
//! `crates/popkern-blend/src/blend/gregory_vertex_enclosure.rs`
//! (which operates on the earlier bicubic Gregory rung) — we
//! port the same Lagrangian QP to `popkern-core` so the
//! direct-construction framework doesn't have a reverse
//! dependency on `popkern-blend`.
//!
//! # Mathematical setup (recap)
//!
//! Work in the 2-D tangent plane at `C`, spanned by an
//! orthonormal basis `(t1, t2)` with normal `n_C`. Per
//! sub-patch `k`:
//!
//! * `R_k` — the raw "parallelogram rule" inner-CP position
//!   (projected into the tangent plane) — the minimum-twist
//!   starting point.
//! * `Y_k` — the solver's output, also in the tangent plane.
//!
//! For each shared spoke `k` between sub-patches `k` and
//! `(k+1) mod N`, let `σ_k` be the 2-D spoke direction in the
//! tangent plane and `p_k = perp(σ_k)` the perpendicular. The
//! G1-along-spoke condition at `C` is
//!
//!   `(Y_k − Y_{k+1}) · p_k = 0`.
//!
//! `N` scalar constraints, `2N` DOFs. Under-determined; the
//! solver minimises `Σ ‖Y_k − R_k‖²` subject to the
//! constraints (least-perturbation from the parallelogram
//! seed).
//!
//! The Lagrangian QP reduces to a circulant tridiagonal `N×N`
//! system in the multipliers `λ_k` — Gaussian elimination for
//! `N ≤ 8` is deterministic and fast.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::{DirectBoundaryInput, TangentPlane};
use crate::geometry::gregory::direct_linalg::{orthonormal_tangent_basis, solve_dense};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::gregory::direct_types::CornerKind;
use crate::geometry::{Point3, Vec3};

/// Per-patch output of the central-corner solve: the
/// `grid[6][1]` position for each sub-patch in the ring.
#[derive(Debug, Clone)]
pub struct CentralCornerInnerCps {
    /// `grid[6][1]` per patch, indexed in the same order as
    /// [`DirectBoundaryInput::spokes`].
    pub inner_cps: Vec<Point3>,
    /// `n_C` — outward unit normal of the common tangent plane
    /// at `C`, averaged from the input face tangent planes.
    pub common_normal: Vec3,
}

/// Solve the central-vertex N-way corner: return the
/// `grid[6][1]` position for every sub-patch so that all `N`
/// tangent planes at `C` coincide.
///
/// Takes only the fields of [`DirectBoundaryInput`] that matter
/// here — the central vertex, each spoke's first boundary
/// derivative at `t=0`, and the per-patch input face tangent
/// planes.
///
/// Errors:
/// - [`DirectConstructionError::CornerTangentPlaneSingular`] if
///   the averaged normal is zero-length (all input face normals
///   cancel — geometrically a degenerate vertex).
pub fn solve_central_vertex_corner(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<CentralCornerInnerCps, DirectConstructionError> {
    let n = input.n();
    let c = input.central_vertex;

    // 1. Compute the common tangent-plane normal `n_C` by
    //    averaging the input-face normals. Each `TangentPlane`
    //    exposes its outward normal as `u × v`.
    let mut n_sum = Vec3::new(0.0, 0.0, 0.0);
    for tps in &input.central_vertex_tps {
        let norm = tps.normal();
        n_sum = Vec3::new(n_sum.x + norm.x, n_sum.y + norm.y, n_sum.z + norm.z);
    }
    // Compare the averaged-normal magnitude to the input-count-
    // normalised tolerance. If each face contributes a unit
    // normal, a "typical" sum magnitude is O(N); we reject when
    // the sum is much smaller than one normal's worth, meaning
    // the face normals cancel geometrically.
    let n_c_raw_len = n_sum.length();
    let n_c_threshold = options.corner_singularity_eps.max(1e-6) * (n as f64);
    if n_c_raw_len < n_c_threshold {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: n_c_raw_len,
            tolerance: n_c_threshold,
            edge_hint: "input-face normals at central vertex sum to ~0".into(),
        });
    }
    let Some(n_c) = n_sum.normalized() else {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: n_c_raw_len,
            tolerance: options.corner_singularity_eps,
            edge_hint: "all faces at central vertex".into(),
        });
    };

    // 2. Build the parallelogram-rule seed for each sub-patch's
    //    corner-adjacent inner CP `grid[6][1]`:
    //       R_k = P[7][1] + P[6][0] − P[7][0]   (with P[7][0] = C)
    //
    //    Index mapping — spokes' `boundary_cps` are indexed in
    //    the spoke's canonical parameter `t ∈ [0, 1]` from `C`
    //    outward, so `boundary_cps[j]` is at spoke-t = j/7.
    //
    //    * `P[7][1]` (patch k, u=1, v=1/7): on patch k's UOne
    //      boundary (spoke k, not reversed), at spoke-t = 1/7.
    //      → `spokes[k].boundary_cps[1]`.
    //    * `P[6][0]` (patch k, u=6/7, v=0): on patch k's VZero
    //      boundary (spoke k-1, **reversed**), at patch k's
    //      u=6/7 → spoke-t = 1 − 6/7 = 1/7.
    //      → `spokes[k-1].boundary_cps[1]`.
    //
    //    (The earlier version of this block mistakenly used
    //    `boundary_cps[6]`, which is at spoke-t = 6/7 and maps to
    //    `grid[1][0]` on patch k, not `grid[6][0]`. Caught during
    //    the 2026-04-19 self-review.)
    let mut raw_parallelogram = Vec::with_capacity(n);
    for k in 0..n {
        let p_7_1 = input.spokes[k].boundary_cps[1];
        let k_prev = (k + n - 1) % n;
        let p_6_0 = input.spokes[k_prev].boundary_cps[1];
        // Parallelogram: P[7][1] + P[6][0] − P[7][0] = P[7][1] + (P[6][0] − C).
        let r_k = Point3::new(
            p_7_1.x + (p_6_0.x - c.x),
            p_7_1.y + (p_6_0.y - c.y),
            p_7_1.z + (p_6_0.z - c.z),
        );
        raw_parallelogram.push(r_k);
    }

    // 3. Per-spoke 3-D direction vectors, pointing from `C`
    //    outward along each spoke. Use a tangent at the `t=0`
    //    end — the first Bernstein-derivative coefficient of
    //    the spoke curve:
    //       T(0) = 7 · (P[1] − P[0])   (degree-6 poly, evaluated at t=0 only needs the first coeff)
    //    Take this as the 3-D spoke direction; normalisation
    //    happens inside the solve.
    let mut spoke_directions = Vec::with_capacity(n);
    for k in 0..n {
        let cps = &input.spokes[k].boundary_cps;
        spoke_directions.push(Vec3::new(
            7.0 * (cps[1].x - cps[0].x),
            7.0 * (cps[1].y - cps[0].y),
            7.0 * (cps[1].z - cps[0].z),
        ));
    }

    // 4. Port of the Lagrangian-QP solver from
    //    popkern-blend::blend::gregory_vertex_enclosure.
    let (t1, t2) = orthonormal_tangent_basis(n_c);

    // Project parallelogram seeds into 2-D tangent-plane coords.
    let raw_2d: Vec<[f64; 2]> = raw_parallelogram
        .iter()
        .map(|p| {
            let off = Vec3::new(p.x - c.x, p.y - c.y, p.z - c.z);
            [off.dot(&t1), off.dot(&t2)]
        })
        .collect();

    // 2-D per-spoke perpendiculars `p_k`.
    let mut p: Vec<[f64; 2]> = Vec::with_capacity(n);
    let mut active: Vec<bool> = Vec::with_capacity(n);
    for sigma in &spoke_directions {
        let sx = sigma.dot(&t1);
        let sy = sigma.dot(&t2);
        let mag = (sx * sx + sy * sy).sqrt();
        if mag < options.corner_singularity_eps {
            p.push([0.0, 0.0]);
            active.push(false);
        } else {
            let ux = sx / mag;
            let uy = sy / mag;
            // perp = rotate 90° CCW: (-uy, ux).
            p.push([-uy, ux]);
            active.push(true);
        }
    }
    if !active.iter().any(|&a| a) {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: 0.0,
            tolerance: options.corner_singularity_eps,
            edge_hint: "all spokes project to zero in the central tangent plane".into(),
        });
    }

    // Build the N×N λ system. Inactive spokes pin λ_k = 0.
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

    // Recover Y_k = R_k − λ_k p_k + λ_{k-1} p_{k-1}, lift 2-D → 3-D.
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

    Ok(CentralCornerInnerCps {
        inner_cps,
        common_normal: n_c,
    })
}

/// Projector helper for external callers: drop a point onto a
/// tangent plane centred at `c` with normal `n`.
#[allow(dead_code)]
pub(crate) fn project_to_tangent_plane(p: Point3, c: Point3, n: Vec3) -> Point3 {
    let off = Vec3::new(p.x - c.x, p.y - c.y, p.z - c.z);
    let d = off.dot(&n);
    Point3::new(p.x - d * n.x, p.y - d * n.y, p.z - d * n.z)
}

/// Helper: the common tangent plane at `C` derived from the
/// input's per-patch tangent planes. Returns the averaged
/// normal and an orthonormal basis.
pub fn central_plane_from_input(input: &DirectBoundaryInput) -> Option<TangentPlane> {
    let mut n_sum = Vec3::new(0.0, 0.0, 0.0);
    for tps in &input.central_vertex_tps {
        let norm = tps.normal();
        n_sum = Vec3::new(n_sum.x + norm.x, n_sum.y + norm.y, n_sum.z + norm.z);
    }
    let n_unit = n_sum.normalized()?;
    TangentPlane::from_normal(n_unit)
}

#[cfg(test)]
#[path = "direct_corner_c_tests.rs"]
mod tests;
