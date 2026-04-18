//! §5.3 item #3 — full Shirman-Séquin vertex-enclosure G1 at the
//! central vertex, extended from G1-at-C only to G1 ALONG each
//! internal spoke.
//!
//! For `N` rectangular Gregory sub-patches sharing a central
//! vertex `C`, the corner-3 inner-CP of each sub-patch controls
//! the cross-tangent at `C` along the two spokes adjacent to that
//! sub-patch. A spoke between adjacent sub-patches `k` and
//! `(k+1) mod N` is shared — at every point along it the two
//! sub-patches must agree on the tangent plane for G1 continuity
//! (and the reflection-line / zebra diagnostics in §5.3 item #11
//! to show zero kink).
//!
//! # Mathematical setup
//!
//! Work in the 2-D tangent plane at `C`, spanned by an orthonormal
//! basis `(t1, t2)` with normal `n_C`. For each sub-patch `k`:
//!
//! * `R_k` — the raw parallelogram inner-CP position, projected
//!   into the tangent plane.
//! * `Y_k` — the solver's output, also in the tangent plane.
//!
//! For each shared spoke `k` (between sub-patches `k` and
//! `(k+1) mod N`), let `σ_k` = unit direction from `C` along the
//! spoke in the tangent plane, and `p_k = perp(σ_k)` = unit
//! perpendicular. The G1-along-spoke condition at `C` is that the
//! cross-tangent VECTORS on either side of the spoke (which sit
//! perpendicular to `σ_k` within the tangent plane) have the
//! same component along `p_k`. In terms of the inner-CP offsets:
//! `(Y_k − Y_{k+1}) · p_k = 0`.
//!
//! This is one scalar constraint per spoke — `N` constraints
//! across `N` 2-D unknowns (2N DOF). The system is
//! under-determined, so we minimise `Σ ||Y_k − R_k||²` subject to
//! the constraints (least-perturbation from the parallelogram
//! raw).
//!
//! # Lagrangian QP
//!
//! Lagrange multipliers `λ_k`, one per constraint:
//! `L = (1/2) Σ_k ‖Y_k − R_k‖² + Σ_k λ_k (Y_k − Y_{k+1}) · p_k`.
//!
//! Stationarity:
//! `∂L/∂Y_k = (Y_k − R_k) + λ_k p_k − λ_{k−1} p_{k−1} = 0` gives
//! `Y_k = R_k − λ_k p_k + λ_{k−1} p_{k−1}`.
//!
//! Substituting into the constraint `(Y_k − Y_{k+1}) · p_k = 0`
//! and using `p_k · p_k = 1` gives the **circulant tridiagonal**
//! system in the `λ_k`:
//! `2 λ_k − cos γ_{k−1,k} · λ_{k−1} − cos γ_{k,k+1} · λ_{k+1}
//! = (R_k − R_{k+1}) · p_k`
//! where `cos γ_{i,j} = p_i · p_j`. Solved directly by Gaussian
//! elimination for `N ≤ 8`.
//!
//! # Reduction to the original mean-shift
//!
//! When the spokes are evenly spaced (regular `N`-gon in the
//! tangent plane), the constraint matrix becomes circulant with
//! equal off-diagonals and the λ_k are constant across `k`; Y_k
//! then equals `R_k − mean(R·σ)` along the spoke direction — i.e.
//! the same mean-shift the previous implementation produced.
//! For irregular spoke angles the new solver correctly couples
//! the Y_k, eliminating the spoke-direction residual that the
//! mean-shift left behind.
//!
//! # Rank-deficient fallback
//!
//! The tangent basis is well-defined when `‖n_C‖² ≥ 1e-18`.
//! Individual spokes can also be degenerate (zero projection into
//! tangent plane) — those constraints are dropped. A totally
//! degenerate input (no usable spokes) returns `None`.

use popkern_core::geometry::{Point3, Vec3};

/// Result of the vertex-enclosure solve: one projected inner CP
/// position per sub-patch (used for BOTH `P_u` and `P_v` of the
/// Gregory corner-3 pair in the symmetric-Chiyokura form).
pub(super) type EnclosureSolve = Vec<Point3>;

/// Shirman-Séquin spoke-wise G1 solve for the corner-3 inner CPs
/// at a central vertex shared by `N` sub-patches.
///
/// Inputs:
/// * `central_point` — the shared corner position `C`.
/// * `n_c` — the outward tangent-plane normal at `C`.
/// * `raw_parallelogram` — pre-SS corner-3 inner CPs (length `N`).
/// * `spoke_directions` — 3D direction from `C` along each shared
///   spoke, indexed so `spoke_directions[k]` is the spoke shared
///   between sub-patches `k` and `(k+1) mod N`. Length `N`.
///
/// Returns `Some(shifted)` with one `Point3` per sub-patch that
/// lies in the tangent plane at `C` AND satisfies the spoke-wise
/// G1 constraints. Returns `None` if `n_c` is degenerate
/// (`‖n_c‖² < 1e-18`) or all spokes project to zero in the
/// tangent plane.
pub(super) fn shirman_sequin_vertex_enclosure_solve(
    central_point: Point3,
    n_c: Vec3,
    raw_parallelogram: &[Point3],
    spoke_directions: &[Vec3],
) -> Option<EnclosureSolve> {
    if n_c.length_squared() < 1e-18 {
        return None;
    }
    let n = raw_parallelogram.len();
    if n == 0 || spoke_directions.len() != n {
        return None;
    }
    let (t1, t2) = orthonormal_tangent_basis(n_c);
    let c_pt = central_point;

    // Project inputs to 2D tangent-plane coordinates.
    let raw_2d: Vec<[f64; 2]> = raw_parallelogram
        .iter()
        .map(|p| {
            let off = Vec3::new(p.x - c_pt.x, p.y - c_pt.y, p.z - c_pt.z);
            [off.dot(&t1), off.dot(&t2)]
        })
        .collect();

    // Project spoke directions to 2D + unit-normalise; build
    // per-spoke perpendicular unit vector `p_k`.
    let mut p: Vec<[f64; 2]> = Vec::with_capacity(n);
    let mut active: Vec<bool> = Vec::with_capacity(n);
    for sigma in spoke_directions {
        let sx = sigma.dot(&t1);
        let sy = sigma.dot(&t2);
        let mag = (sx * sx + sy * sy).sqrt();
        if mag < 1e-12 {
            p.push([0.0, 0.0]);
            active.push(false);
        } else {
            let ux = sx / mag;
            let uy = sy / mag;
            // perp = rotate 90° CCW: (-uy, ux)
            p.push([-uy, ux]);
            active.push(true);
        }
    }
    if !active.iter().any(|&a| a) {
        return None;
    }

    // Build the N×N system for λ. For inactive spokes, pin
    // λ_k = 0 by setting the row to the identity.
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut rhs = vec![0.0_f64; n];
    for k in 0..n {
        if !active[k] {
            matrix[k][k] = 1.0;
            rhs[k] = 0.0;
            continue;
        }
        let k_prev = (k + n - 1) % n;
        let k_next = (k + 1) % n;
        // Diagonal = 2 · (p_k · p_k) = 2.
        matrix[k][k] = 2.0;
        if active[k_prev] {
            // Coupling with λ_{k-1} via p_{k-1} · p_k.
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

    let lambdas = solve_linear_system(matrix, rhs)?;

    // Recover Y_k = R_k - λ_k p_k + λ_{k-1} p_{k-1}, then lift to 3D.
    let mut out = Vec::with_capacity(n);
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
        out.push(Point3::new(
            c_pt.x + yx * t1.x + yy * t2.x,
            c_pt.y + yx * t1.y + yy * t2.y,
            c_pt.z + yx * t1.z + yy * t2.z,
        ));
    }
    Some(out)
}

/// Gaussian elimination with partial pivoting for small dense
/// linear systems. Works in-place on `matrix` + `rhs`; returns
/// `None` if the matrix is singular to machine precision.
fn solve_linear_system(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = m.len();
    for i in 0..n {
        // Partial pivot.
        let mut pivot = i;
        let mut best = m[i][i].abs();
        for (r, row) in m.iter().enumerate().skip(i + 1) {
            let v = row[i].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if pivot != i {
            m.swap(i, pivot);
            b.swap(i, pivot);
        }
        // Eliminate below.
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            let pivot_row = m[i].clone();
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                m[r][c] -= factor * pc;
            }
            b[r] -= factor * b[i];
        }
    }
    // Back-substitute.
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for (c, xc) in x.iter().enumerate().skip(i + 1) {
            s -= m[i][c] * xc;
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

/// Build an orthonormal tangent basis `(t1, t2)` in the plane
/// perpendicular to `n`. `t1` is derived from the world-axis
/// least aligned with `n`; `t2 = n × t1`.
pub(super) fn orthonormal_tangent_basis(n: Vec3) -> (Vec3, Vec3) {
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let dot = seed.dot(&n);
    let raw_t1 = Vec3::new(seed.x - dot * n.x, seed.y - dot * n.y, seed.z - dot * n.z);
    let t1 = raw_t1.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let t2 = n.cross(&t1);
    (t1, t2)
}
