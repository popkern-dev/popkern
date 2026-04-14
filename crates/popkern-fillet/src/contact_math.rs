//! Math primitives for the rolling-ball contact solver: signed
//! normal, Weingarten map, 4×4 linear solve, condition estimate,
//! Jacobian/residual builder, Levenberg-Marquardt regularized step.
//!
//! Split from `contact.rs` to keep individual files under the 500-line
//! structural cap. All functions here are `pub(crate)` so they remain
//! callable from the solver entry points and unit-testable from the
//! sibling test files.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};

use crate::contact::BallSide;

/// Compute the unit normal `n̂` at `(u, v)` on `surface`, multiplied by
/// the side sign. Returns `None` if the raw normal is degenerate
/// (length below the numerical floor).
pub(crate) fn signed_unit_normal(
    surface: &dyn SurfaceEval,
    u: f64,
    v: f64,
    side: BallSide,
    numerical_floor: f64,
) -> Option<Vec3> {
    let raw = surface.normal_at(u, v);
    let len = raw.length();
    if len < numerical_floor {
        return None;
    }
    Some(raw * (side.sign() / len))
}

/// Weingarten map: `dn̂/du` at `(u, v)`.
///
/// The shape operator (Weingarten map) in the `{S_u, S_v}` basis is
/// `I⁻¹ · II`, where `I = [[E, F], [F, G]]` is the first fundamental
/// form and `II = [[L, M], [M, N]]` is the second fundamental form.
/// Expanding the inverse and pushing through the multiplication:
///
/// `dn̂/du = -[(G·L − F·M)·S_u + (E·M − F·L)·S_v] / (EG − F²)`
///
/// Coefficients: `E = S_u·S_u`, `F = S_u·S_v`, `G = S_v·S_v`,
/// `L = S_uu·n̂`, `M = S_uv·n̂`, `N = S_vv·n̂`. The second-form
/// coefficients are computed against the *signed* normal `n̂`, so
/// the result already carries the [`BallSide`] sign.
///
/// Returns `Vec3::zero()` if the surface metric is degenerate
/// (`EG − F² < numerical_floor²`).
pub(crate) fn weingarten_dn_du(
    surface: &dyn SurfaceEval,
    u: f64,
    v: f64,
    n_hat: &Vec3,
    numerical_floor: f64,
) -> Vec3 {
    let s_u = surface.partial_u(u, v);
    let s_v = surface.partial_v(u, v);
    let s_uu = surface.second_partial_uu(u, v);
    let s_uv = surface.second_partial_uv(u, v);

    let e = s_u.dot(&s_u);
    let f = s_u.dot(&s_v);
    let g = s_v.dot(&s_v);
    let metric_det = e * g - f * f;
    if metric_det < numerical_floor * numerical_floor {
        return Vec3::zero();
    }

    let l = s_uu.dot(n_hat);
    let m = s_uv.dot(n_hat);

    let coef_u = (g * l - f * m) / metric_det;
    let coef_v = (e * m - f * l) / metric_det;
    (s_u * coef_u + s_v * coef_v) * -1.0
}

/// Weingarten map: `dn̂/dv` at `(u, v)`.
///
/// `dn̂/dv = -[(G·M − F·N)·S_u + (E·N − F·M)·S_v] / (EG − F²)`
///
/// See [`weingarten_dn_du`] for the derivation and coefficient
/// definitions.
pub(crate) fn weingarten_dn_dv(
    surface: &dyn SurfaceEval,
    u: f64,
    v: f64,
    n_hat: &Vec3,
    numerical_floor: f64,
) -> Vec3 {
    let s_u = surface.partial_u(u, v);
    let s_v = surface.partial_v(u, v);
    let s_uv = surface.second_partial_uv(u, v);
    let s_vv = surface.second_partial_vv(u, v);

    let e = s_u.dot(&s_u);
    let f = s_u.dot(&s_v);
    let g = s_v.dot(&s_v);
    let metric_det = e * g - f * f;
    if metric_det < numerical_floor * numerical_floor {
        return Vec3::zero();
    }

    let m = s_uv.dot(n_hat);
    let n_coef = s_vv.dot(n_hat);

    let coef_u = (g * m - f * n_coef) / metric_det;
    let coef_v = (e * n_coef - f * m) / metric_det;
    (s_u * coef_u + s_v * coef_v) * -1.0
}

/// Solve a 4×4 linear system `A · x = b` by Gaussian elimination
/// with partial pivoting. Returns `None` if any pivot falls below
/// the numerical floor (rank-deficient matrix).
#[allow(clippy::needless_range_loop)] // Gaussian elimination is row-and-column indexing by definition.
pub(crate) fn solve_4x4(
    mut a: [[f64; 4]; 4],
    mut b: [f64; 4],
    numerical_floor: f64,
) -> Option<[f64; 4]> {
    for k in 0..4 {
        let mut pivot_row = k;
        let mut pivot_val = a[k][k].abs();
        for r in (k + 1)..4 {
            let v = a[r][k].abs();
            if v > pivot_val {
                pivot_val = v;
                pivot_row = r;
            }
        }
        if pivot_val < numerical_floor {
            return None;
        }
        if pivot_row != k {
            a.swap(k, pivot_row);
            b.swap(k, pivot_row);
        }
        for r in (k + 1)..4 {
            let factor = a[r][k] / a[k][k];
            for c in k..4 {
                a[r][c] -= factor * a[k][c];
            }
            b[r] -= factor * b[k];
        }
    }
    let mut x = [0.0_f64; 4];
    for k in (0..4).rev() {
        let mut sum = b[k];
        for c in (k + 1)..4 {
            sum -= a[k][c] * x[c];
        }
        x[k] = sum / a[k][k];
    }
    Some(x)
}

/// Cheap 1-norm condition number estimate for a 4×4 matrix.
///
/// Returns `||A||_1 · ||A^{-1}||_1` where `||·||_1` is the maximum
/// absolute column sum. Returns `f64::INFINITY` if the matrix is
/// singular.
pub(crate) fn condition_estimate_4x4(a: &[[f64; 4]; 4]) -> f64 {
    fn one_norm(m: &[[f64; 4]; 4]) -> f64 {
        let mut max_col = 0.0_f64;
        for col in 0..4 {
            let col_sum: f64 = m.iter().map(|row| row[col].abs()).sum();
            if col_sum > max_col {
                max_col = col_sum;
            }
        }
        max_col
    }
    let a_norm = one_norm(a);
    if a_norm == 0.0 {
        return f64::INFINITY;
    }
    let mut a_inv = [[0.0_f64; 4]; 4];
    for col in 0..4 {
        let mut e = [0.0; 4];
        e[col] = 1.0;
        match solve_4x4(*a, e, 1e-300) {
            Some(x) => {
                for (r, x_r) in x.iter().enumerate() {
                    a_inv[r][col] = *x_r;
                }
            }
            None => return f64::INFINITY,
        }
    }
    let a_inv_norm = one_norm(&a_inv);
    a_norm * a_inv_norm
}

/// Build the 4×4 Jacobian of the contact + spine-plane system at
/// `(u0, v0, u1, v1)`.
///
/// Rows 0..3 are the contact equality (3D vector equation
/// `C0 - C1 = 0`); row 3 is the spine-plane scalar constraint
/// `(C0 - p_spine) · spine_plane_normal = 0`.
///
/// Returns `(jacobian, residual, c0, c1, n0, n1)` or `None` if
/// either surface produces a degenerate normal.
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub(crate) fn build_jacobian_and_residual(
    surface_0: &dyn SurfaceEval,
    surface_1: &dyn SurfaceEval,
    side_0: BallSide,
    side_1: BallSide,
    radius: f64,
    uv0: (f64, f64),
    uv1: (f64, f64),
    spine_point: Point3,
    spine_plane_normal: Vec3,
    numerical_floor: f64,
) -> Option<([[f64; 4]; 4], [f64; 4], Point3, Point3, Vec3, Vec3)> {
    let n0 = signed_unit_normal(surface_0, uv0.0, uv0.1, side_0, numerical_floor)?;
    let n1 = signed_unit_normal(surface_1, uv1.0, uv1.1, side_1, numerical_floor)?;

    let p0 = surface_0.point_at(uv0.0, uv0.1);
    let p1 = surface_1.point_at(uv1.0, uv1.1);
    let c0 = p0 + n0 * radius;
    let c1 = p1 + n1 * radius;

    let s0_u = surface_0.partial_u(uv0.0, uv0.1);
    let s0_v = surface_0.partial_v(uv0.0, uv0.1);
    let s1_u = surface_1.partial_u(uv1.0, uv1.1);
    let s1_v = surface_1.partial_v(uv1.0, uv1.1);

    let dn0_du = weingarten_dn_du(surface_0, uv0.0, uv0.1, &n0, numerical_floor);
    let dn0_dv = weingarten_dn_dv(surface_0, uv0.0, uv0.1, &n0, numerical_floor);
    let dn1_du = weingarten_dn_du(surface_1, uv1.0, uv1.1, &n1, numerical_floor);
    let dn1_dv = weingarten_dn_dv(surface_1, uv1.0, uv1.1, &n1, numerical_floor);

    // Contact-equality Jacobian columns: ∂(C0 - C1) / ∂x
    let dc_du0 = s0_u + dn0_du * radius;
    let dc_dv0 = s0_v + dn0_dv * radius;
    let dc_du1 = (s1_u + dn1_du * radius) * -1.0;
    let dc_dv1 = (s1_v + dn1_dv * radius) * -1.0;

    let mut j = [[0.0_f64; 4]; 4];
    j[0][0] = dc_du0.x;
    j[1][0] = dc_du0.y;
    j[2][0] = dc_du0.z;
    j[0][1] = dc_dv0.x;
    j[1][1] = dc_dv0.y;
    j[2][1] = dc_dv0.z;
    j[0][2] = dc_du1.x;
    j[1][2] = dc_du1.y;
    j[2][2] = dc_du1.z;
    j[0][3] = dc_dv1.x;
    j[1][3] = dc_dv1.y;
    j[2][3] = dc_dv1.z;

    // Spine-plane row: g(x) = (C0 - p_spine) · spine_plane_normal
    j[3][0] = dc_du0.dot(&spine_plane_normal);
    j[3][1] = dc_dv0.dot(&spine_plane_normal);
    j[3][2] = 0.0;
    j[3][3] = 0.0;

    let f = c0 - c1;
    let g = (c0 - spine_point).dot(&spine_plane_normal);
    let residual = [f.x, f.y, f.z, g];

    Some((j, residual, c0, c1, n0, n1))
}

/// Levenberg-Marquardt regularized step.
///
/// Computes `(JᵀJ + λI)⁻¹·(−Jᵀ·F)` for `λ` starting at `initial_lambda`
/// and growing 10× per failed solve, up to a cap of 1.0. Returns the
/// first successful step or `None` if all damping levels failed.
#[allow(clippy::needless_range_loop)] // Matrix-vector arithmetic is row-and-column indexing.
pub(crate) fn lm_step(
    jacobian: &[[f64; 4]; 4],
    residual: &[f64; 4],
    initial_lambda: f64,
    numerical_floor: f64,
) -> Option<[f64; 4]> {
    let mut jt_j = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut s = 0.0;
            for k in 0..4 {
                s += jacobian[k][i] * jacobian[k][j];
            }
            jt_j[i][j] = s;
        }
    }
    let mut neg_jt_f = [0.0_f64; 4];
    for i in 0..4 {
        let mut s = 0.0;
        for k in 0..4 {
            s += jacobian[k][i] * residual[k];
        }
        neg_jt_f[i] = -s;
    }

    let mut lambda = initial_lambda.max(1.0e-9);
    while lambda <= 1.0 {
        let mut a = jt_j;
        for d in 0..4 {
            a[d][d] += lambda;
        }
        if let Some(x) = solve_4x4(a, neg_jt_f, numerical_floor) {
            return Some(x);
        }
        lambda *= 10.0;
    }
    None
}
