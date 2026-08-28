//! Analytic ring-Jacobian helpers for bi-nonic guide optimization.

use popkern_core::geometry::gregory::{
    RectangularGregoryBiNonic, bernstein_binonic_all, bernstein_binonic_deriv_all,
};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

/// Flat index into the ring-CP vector: patch `p`, cell `(i, j)`,
/// axis `a ∈ {0=x, 1=y, 2=z}`. `pub(crate)` because the
/// analytic-Jacobian helpers are shared with the sibling
/// optimizers (bi-quintic, exact-boundary, chart-orientation
/// polish) that build selector linearizations over ring CPs.
#[inline]
pub(crate) fn ring_idx(p: usize, i: usize, j: usize, a: usize) -> usize {
    p * 300 + (i * 10 + j) * 3 + a
}

/// Flatten a bi-nonic ring's CPs to a row-major `D_ring` vector
/// `[patch_0_grid (row-major), patch_1_grid, ...]` with each CP
/// stored as `(x, y, z)`. Shared with `kp_biquintic_optimizer`
/// for the phase-2 analytic-Jacobian precompute.
pub(crate) fn flatten_ring(ring: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let mut out = Vec::with_capacity(ring.len() * 300);
    for patch in ring {
        for i in 0..10 {
            for j in 0..10 {
                out.push(patch.grid[i][j].x);
                out.push(patch.grid[i][j].y);
                out.push(patch.grid[i][j].z);
            }
        }
    }
    out
}

/// Analytic Jacobian `J[sample, k] = ∂f[sample]/∂guide_CP_k`.
///
/// Both phase 2 (bi-quintic guides) and phase 3 (bi-nonic guides)
/// produce bi-nonic rings; this function takes the ring and the
/// pre-computed `M = ∂ring_CP/∂guide_CP` matrix and produces the
/// residual-vector Jacobian via the chain rule. It does NOT care
/// about the guide degree — all that matters is that `m.len() ==
/// ring.len() * 300` and each row of `m` has exactly `d_guide`
/// columns.
///
/// Cost: `O(samples · D_ring_relevant · D_guide)`. For the bi-nonic
/// output that's `50N · 30 · D_guide` per iteration.
pub(crate) fn analytic_ring_jacobian(
    ring: &[RectangularGregoryBiNonic],
    m: &[Vec<f64>],
    d_guide: usize,
) -> Vec<Vec<f64>> {
    let n = ring.len();
    let total_samples = n * 50;
    let mut jac = vec![vec![0.0; d_guide]; total_samples];

    for spoke in 0..n {
        let a_idx = spoke;
        let b_idx = (spoke + 1) % n;
        let a = &ring[a_idx];
        let b = &ring[b_idx];

        for s in 0..50 {
            let t = 0.05 + 0.90 * (s as f64) / 49.0;
            let t_b = 1.0 - t;

            // Bernstein bases at this sample — constants the
            // chain-rule inner loop needs.
            let b_u = bernstein_binonic_all(t);
            let b_prime = bernstein_binonic_deriv_all(t);
            let b_v = bernstein_binonic_all(t_b);

            // Reference vectors at current iterate.
            let u = a.partial_u(1.0, t);
            let v = b.partial_v(t_b, 0.0);
            let w = a.partial_v(1.0, t);

            let cross_uv = u.cross(&v);
            let triple = cross_uv.x * w.x + cross_uv.y * w.y + cross_uv.z * w.z;
            let sign_t = if triple >= 0.0 { 1.0 } else { -1.0 };

            let nu_sq = u.x * u.x + u.y * u.y + u.z * u.z;
            let nv_sq = v.x * v.x + v.y * v.y + v.z * v.z;
            let nw_sq = w.x * w.x + w.y * w.y + w.z * w.z;
            let denom = nu_sq.sqrt() * nv_sq.sqrt() * nw_sq.sqrt();

            let row = spoke * 50 + s;
            if denom < 1e-300 {
                continue; // matches residual_vector's ill-conditioned zero
            }

            let f_val = triple.abs() / denom;

            // ∂f/∂u = sign_t · (v × w) / denom − f · u / |u|²  (similarly for v, w)
            let v_cross_w = v.cross(&w);
            let w_cross_u = w.cross(&u);
            let df_du = Vec3::new(
                sign_t * v_cross_w.x / denom - f_val * u.x / nu_sq,
                sign_t * v_cross_w.y / denom - f_val * u.y / nu_sq,
                sign_t * v_cross_w.z / denom - f_val * u.z / nu_sq,
            );
            let df_dv = Vec3::new(
                sign_t * w_cross_u.x / denom - f_val * v.x / nv_sq,
                sign_t * w_cross_u.y / denom - f_val * v.y / nv_sq,
                sign_t * w_cross_u.z / denom - f_val * v.z / nv_sq,
            );
            let df_dw = Vec3::new(
                sign_t * cross_uv.x / denom - f_val * w.x / nw_sq,
                sign_t * cross_uv.y / denom - f_val * w.y / nw_sq,
                sign_t * cross_uv.z / denom - f_val * w.z / nw_sq,
            );

            // Row via chain rule.
            // ∂u_α/∂g_k = 9 · Σ_j B_j(t) · (M[a,9,j,α; k] − M[a,8,j,α; k])
            // ∂v_α/∂g_k = 9 · Σ_i B_i(t_b) · (M[b,i,1,α; k] − M[b,i,0,α; k])
            // ∂w_α/∂g_k =      Σ_j B'_j(t) · M[a,9,j,α; k]
            let row_slice = &mut jac[row];
            for k in 0..d_guide {
                let mut du_x = 0.0_f64;
                let mut du_y = 0.0_f64;
                let mut du_z = 0.0_f64;
                let mut dv_x = 0.0_f64;
                let mut dv_y = 0.0_f64;
                let mut dv_z = 0.0_f64;
                let mut dw_x = 0.0_f64;
                let mut dw_y = 0.0_f64;
                let mut dw_z = 0.0_f64;
                for jj in 0..10 {
                    let bj = b_u[jj];
                    let r9x = m[ring_idx(a_idx, 9, jj, 0)][k];
                    let r9y = m[ring_idx(a_idx, 9, jj, 1)][k];
                    let r9z = m[ring_idx(a_idx, 9, jj, 2)][k];
                    let r8x = m[ring_idx(a_idx, 8, jj, 0)][k];
                    let r8y = m[ring_idx(a_idx, 8, jj, 1)][k];
                    let r8z = m[ring_idx(a_idx, 8, jj, 2)][k];
                    du_x += bj * (r9x - r8x);
                    du_y += bj * (r9y - r8y);
                    du_z += bj * (r9z - r8z);

                    let bp = b_prime[jj];
                    dw_x += bp * r9x;
                    dw_y += bp * r9y;
                    dw_z += bp * r9z;
                }
                for ii in 0..10 {
                    let bi = b_v[ii];
                    let rb1x = m[ring_idx(b_idx, ii, 1, 0)][k];
                    let rb1y = m[ring_idx(b_idx, ii, 1, 1)][k];
                    let rb1z = m[ring_idx(b_idx, ii, 1, 2)][k];
                    let rb0x = m[ring_idx(b_idx, ii, 0, 0)][k];
                    let rb0y = m[ring_idx(b_idx, ii, 0, 1)][k];
                    let rb0z = m[ring_idx(b_idx, ii, 0, 2)][k];
                    dv_x += bi * (rb1x - rb0x);
                    dv_y += bi * (rb1y - rb0y);
                    dv_z += bi * (rb1z - rb0z);
                }
                du_x *= 9.0;
                du_y *= 9.0;
                du_z *= 9.0;
                dv_x *= 9.0;
                dv_y *= 9.0;
                dv_z *= 9.0;

                row_slice[k] = df_du.x * du_x
                    + df_du.y * du_y
                    + df_du.z * du_z
                    + df_dv.x * dv_x
                    + df_dv.y * dv_y
                    + df_dv.z * dv_z
                    + df_dw.x * dw_x
                    + df_dw.y * dw_y
                    + df_dw.z * dw_z;
            }
        }
    }

    jac
}

/// One outer-fillet residual sample with the cap normal held fixed:
/// project `p` onto `fillet`, take the fillet normal there, and
/// return `|n_cap × n_f| / (|n_cap| · |n_f|)` — exactly the
/// per-sample quantity `outer_fillet_residual_vector` computes.
/// `None` when the projection fails or the normals are
/// ill-conditioned (the residual pushes 0 there).
fn outer_sample_sine(n_cap: &Vec3, fillet: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    let (u_f, v_f) = project_point_to_uv(fillet, p)?;
    let n_f = fillet
        .partial_u(u_f, v_f)
        .cross(&fillet.partial_v(u_f, v_f));
    let denom = n_cap.length() * n_f.length();
    if denom > 1e-300 {
        Some(n_cap.cross(&n_f).length() / denom)
    } else {
        None
    }
}

/// Analytic Jacobian block for the the outer-fillet-G1 work outer-fillet residual
/// (PK-275): `J[sample, k] = ∂f[sample]/∂guide_CP_k` for the
/// `50 · N` cap-vs-fillet G1 samples appended to the residual vector
/// when `outer_fillets` is active. Rows align with
/// `outer_fillet_residual_vector`'s ordering (patch-major, 50
/// samples each); patches without a fillet entry produce zero rows,
/// matching the residual's zero padding.
///
/// # Derivation
///
/// Per sample, `f = |n_c × n_f| / (|n_c| · |n_f|)` where
/// `n_c = ∂cap/∂u × ∂cap/∂v` at `(0, t)` and `n_f` is the fillet
/// normal at the projection of `p = cap(0, t)`. Three chains reach
/// the guide DoF, all through the constant ring linearization `M`
/// (the ring build is linear in guide CPs):
///
/// - `∂f/∂u` and `∂f/∂v` (cap tangents) are closed-form: with
///   `q = ∂f/∂n_c = (n_f × ĉ)/(|n_c||n_f|) − f · n_c/|n_c|²`,
///   `∂f/∂u = v × q` and `∂f/∂v = q × u` (from `n_c = u × v`).
/// - `∂f/∂p` (footpoint motion moves the fillet normal) is a
///   central finite difference of `outer_sample_sine` in the three
///   spatial axes — a 3-dimensional FD of a cheap smooth map, not
///   the full-DoF FD this block replaces. Fully-analytic would need
///   the fillet's shape operator at the footpoint; the 6 projection
///   calls per sample are negligible against the ring rebuild the
///   FD fallback needed per DoF.
///
/// Boundary-row chain rule (mirror of the inner block's rows 8/9 at
/// `u = 1`, here rows 0/1 at `u = 0`):
///
/// - `∂u_α/∂g_k = 9 · Σ_j B_j(t) · (M[p,1,j,α; k] − M[p,0,j,α; k])`
/// - `∂v_α/∂g_k =     Σ_j B'_j(t) · M[p,0,j,α; k]`
/// - `∂p_α/∂g_k =     Σ_j B_j(t) · M[p,0,j,α; k]`
pub(crate) fn outer_fillet_jacobian(
    ring: &[RectangularGregoryBiNonic],
    outer_fillets: &[SurfaceGeometry],
    m: &[Vec<f64>],
    d_guide: usize,
) -> Vec<Vec<f64>> {
    outer_fillet_jacobian_sampled(ring, outer_fillets, m, d_guide, 50)
}

/// Sample-count-generalized form of [`outer_fillet_jacobian`]: the
/// phase-3 optimizer samples the outer sine at 50 points per patch,
/// the chart-orientation polish at 25 (`outer_fillet_sin` grid).
/// The `t` grid is `0.05 + 0.90 · s/(samples − 1)` in both cases.
pub(crate) fn outer_fillet_jacobian_sampled(
    ring: &[RectangularGregoryBiNonic],
    outer_fillets: &[SurfaceGeometry],
    m: &[Vec<f64>],
    d_guide: usize,
    samples: usize,
) -> Vec<Vec<f64>> {
    let n = ring.len();
    let mut jac = vec![vec![0.0; d_guide]; n * samples];

    for (patch_idx, cap) in ring.iter().enumerate() {
        let Some(fillet) = outer_fillets.get(patch_idx) else {
            continue; // zero rows match the residual's zero padding
        };
        for s in 0..samples {
            let t = 0.05 + 0.90 * (s as f64) / (samples as f64 - 1.0);
            let b_v = bernstein_binonic_all(t);
            let b_prime = bernstein_binonic_deriv_all(t);

            let u = cap.partial_u(0.0, t);
            let v = cap.partial_v(0.0, t);
            let p = cap.point_at(0.0, t);
            let n_c = u.cross(&v);

            let Some((u_f, v_f)) = project_point_to_uv(fillet, &p) else {
                continue; // residual is 0 here; leave the row zero
            };
            let n_f = fillet
                .partial_u(u_f, v_f)
                .cross(&fillet.partial_v(u_f, v_f));

            let b_len = n_c.length();
            let c_len = n_f.length();
            let denom = b_len * c_len;
            if denom < 1e-300 {
                continue; // matches the residual's ill-conditioned zero
            }
            let c = n_c.cross(&n_f);
            let a_len = c.length();
            let f_val = a_len / denom;

            // q = ∂f/∂n_c with n_f held fixed. At exact tangency
            // (|c| → 0) the |·| kink makes ĉ undefined; the residual
            // is at its minimum there and the row's cap-tangent
            // contribution vanishes in the limit, so drop it.
            let q = if a_len > 1e-14 * denom {
                let c_hat = Vec3::new(c.x / a_len, c.y / a_len, c.z / a_len);
                let nf_x_chat = n_f.cross(&c_hat);
                let b_sq = b_len * b_len;
                Vec3::new(
                    nf_x_chat.x / denom - f_val * n_c.x / b_sq,
                    nf_x_chat.y / denom - f_val * n_c.y / b_sq,
                    nf_x_chat.z / denom - f_val * n_c.z / b_sq,
                )
            } else {
                Vec3::zero()
            };
            let df_du = v.cross(&q);
            let df_dv = q.cross(&u);

            // ∂f/∂p by central difference of the projection chain.
            let h = 1e-6 * (1.0 + (p - Point3::new(0.0, 0.0, 0.0)).length());
            let mut df_dp = [0.0_f64; 3];
            for (axis, slot) in df_dp.iter_mut().enumerate() {
                let mut p_plus = p;
                let mut p_minus = p;
                match axis {
                    0 => {
                        p_plus.x += h;
                        p_minus.x -= h;
                    }
                    1 => {
                        p_plus.y += h;
                        p_minus.y -= h;
                    }
                    _ => {
                        p_plus.z += h;
                        p_minus.z -= h;
                    }
                }
                if let (Some(f_plus), Some(f_minus)) = (
                    outer_sample_sine(&n_c, fillet, &p_plus),
                    outer_sample_sine(&n_c, fillet, &p_minus),
                ) {
                    *slot = (f_plus - f_minus) / (2.0 * h);
                }
            }

            let row_slice = &mut jac[patch_idx * samples + s];
            for k in 0..d_guide {
                let mut du_x = 0.0_f64;
                let mut du_y = 0.0_f64;
                let mut du_z = 0.0_f64;
                let mut dv_x = 0.0_f64;
                let mut dv_y = 0.0_f64;
                let mut dv_z = 0.0_f64;
                let mut dp_x = 0.0_f64;
                let mut dp_y = 0.0_f64;
                let mut dp_z = 0.0_f64;
                for jj in 0..10 {
                    let bj = b_v[jj];
                    let bp = b_prime[jj];
                    let r0x = m[ring_idx(patch_idx, 0, jj, 0)][k];
                    let r0y = m[ring_idx(patch_idx, 0, jj, 1)][k];
                    let r0z = m[ring_idx(patch_idx, 0, jj, 2)][k];
                    let r1x = m[ring_idx(patch_idx, 1, jj, 0)][k];
                    let r1y = m[ring_idx(patch_idx, 1, jj, 1)][k];
                    let r1z = m[ring_idx(patch_idx, 1, jj, 2)][k];
                    du_x += bj * (r1x - r0x);
                    du_y += bj * (r1y - r0y);
                    du_z += bj * (r1z - r0z);
                    dv_x += bp * r0x;
                    dv_y += bp * r0y;
                    dv_z += bp * r0z;
                    dp_x += bj * r0x;
                    dp_y += bj * r0y;
                    dp_z += bj * r0z;
                }
                du_x *= 9.0;
                du_y *= 9.0;
                du_z *= 9.0;

                row_slice[k] = df_du.x * du_x
                    + df_du.y * du_y
                    + df_du.z * du_z
                    + df_dv.x * dv_x
                    + df_dv.y * dv_y
                    + df_dv.z * dv_z
                    + df_dp[0] * dp_x
                    + df_dp[1] * dp_y
                    + df_dp[2] * dp_z;
            }
        }
    }

    jac
}
