//! Levenberg-Marquardt solver for the K-P ring-residual
//! nonlinear least-squares problem.
//!
//! # Why L-M here
//!
//! The ring G1 residual is structurally a sum of squared per-
//! sample coplanarities: `r(x) = Σ_i f_i(x)²`, so it's a natural
//! NLS problem. The existing FD gradient-
//! descent (`kp_*_optimizer.rs`) steps in `-∇r`; LM steps in
//! `-(JᵀJ + λI)⁻¹ Jᵀf`, incorporating curvature from the
//! Gauss-Newton approximation. On cube_extreme the GD stack
//! bottoms out with ~60% of steps rejected by the fold-safety
//! guard; LM's adaptive damping is designed to handle exactly
//! that case (λ small → GN, λ large → SD, auto-switch based on
//! step quality).
//!
//! # Scope
//!
//! Generic enough for both the bi-quintic (phase 2, 48N DoF)
//! and bi-nonic (phase 3, 192N DoF) cases. The caller provides:
//!
//! - Initial parameter vector `x0` (flat Vec<f64>).
//! - A residual function `x → (f, fold_safe)` returning the
//!   per-sample residual vector and a fold-safety flag.
//! - Hyperparameters: `fd_eps`, `lambda_init`, `lambda_max`,
//!   `max_iter`, `target_residual`.
//!
//! The FD Jacobian is computed by perturbing each DoF by
//! `fd_eps` and storing per-sample differences (the same
//! ring-build count as the existing GD, just with more data
//! retained per build).

/// One L-M optimization pass.
pub(super) struct LmConfig {
    pub fd_eps: f64,
    pub lambda_init: f64,
    pub lambda_max: f64,
    pub max_iter: usize,
    pub target_residual: f64,
}

pub(super) struct LmResult {
    pub x: Vec<f64>,
    #[allow(dead_code)]
    pub final_residual: f64,
    #[allow(dead_code)]
    pub iters: usize,
    #[allow(dead_code)]
    pub accepted: usize,
}

/// Callback trait. Splits residual evaluation from fold-safety
/// checking because the latter can be expensive (e.g. surface-
/// surface self-intersection on bidegree-9 Gregory patches takes
/// ~100 ms per call). The L-M driver calls:
///
/// - `residual(x)` on every FD perturbation (576 calls per
///   Jacobian iteration for the bi-nonic phase 3).
/// - `fold_safe(x)` only on *candidate* step vectors — typically
///   1-10 calls per iteration, inside the damping loop.
pub(super) trait LmEval {
    /// Fast path: residuals only. Returns `None` if the build
    /// fails (e.g. ring assembly errored).
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>>;

    /// Slow path: has the fold-safety property? Only called on
    /// candidate step vectors.
    fn fold_safe(&mut self, x: &[f64]) -> bool;

    /// Optional: analytic Jacobian. Returns an `m × n` matrix
    /// where `m = residual(x).len()` and `n = x.len()`. Return
    /// `None` to use FD fallback. Implementations that can
    /// exploit linear structure in the residual chain (as in the
    /// K-P ring builder, where ring CPs are linear in guide CPs)
    /// get order-of-magnitude speedup over FD.
    fn jacobian(&mut self, _x: &[f64]) -> Option<Vec<Vec<f64>>> {
        None
    }
}

/// Run the L-M loop. Returns the best (lowest-residual) fold-
/// safe iterate found during the search.
pub(super) fn run_lm<E: LmEval>(
    mut eval: E,
    x0: Vec<f64>,
    config: &LmConfig,
    trace_tag: &str,
) -> LmResult {
    let n_dof = x0.len();
    let mut x = x0;

    // Initial residual. If the starting point isn't fold-safe,
    // we return it unchanged — the caller is responsible for
    // providing a feasible seed.
    let Some(mut f) = eval.residual(&x) else {
        return LmResult {
            x,
            final_residual: f64::INFINITY,
            iters: 0,
            accepted: 0,
        };
    };
    let mut residual_sq = sum_sq(&f);
    if !eval.fold_safe(&x) {
        return LmResult {
            x,
            final_residual: residual_sq,
            iters: 0,
            accepted: 0,
        };
    }

    let mut best_x = x.clone();
    let mut best_residual = residual_sq;
    let mut lambda = config.lambda_init;
    let mut iters = 0usize;
    let mut accepted = 0usize;

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!("[lm {trace_tag}] init: residual={residual_sq:.3e} D={n_dof}");
    }

    for _ in 0..config.max_iter {
        iters += 1;
        if residual_sq < config.target_residual {
            break;
        }

        // Analytic Jacobian if the trait provides one, else FD.
        let m = f.len();
        let jac: Vec<Vec<f64>> = match eval.jacobian(&x) {
            Some(j) if j.len() == m && j.iter().all(|row| row.len() == n_dof) => j,
            _ => {
                let mut jac: Vec<Vec<f64>> = vec![vec![0.0; n_dof]; m];
                for j in 0..n_dof {
                    let saved = x[j];
                    x[j] = saved + config.fd_eps;
                    if let Some(f_pert) = eval.residual(&x)
                        && f_pert.len() == m
                    {
                        for (i, row) in jac.iter_mut().enumerate() {
                            row[j] = (f_pert[i] - f[i]) / config.fd_eps;
                        }
                    }
                    x[j] = saved;
                }
                jac
            }
        };

        // g = Jᵀ f  (size n_dof)
        let mut g = vec![0.0; n_dof];
        for (i, fi) in f.iter().enumerate() {
            for (j, gj) in g.iter_mut().enumerate() {
                *gj += jac[i][j] * fi;
            }
        }

        // H = JᵀJ  (n_dof × n_dof), lower triangle only — we
        // mirror on solve.
        let mut h = vec![vec![0.0_f64; n_dof]; n_dof];
        for jac_row in jac.iter() {
            for a in 0..n_dof {
                let v = jac_row[a];
                if v == 0.0 {
                    continue;
                }
                for (b, hab) in h[a].iter_mut().enumerate().take(a + 1) {
                    *hab += v * jac_row[b];
                }
            }
        }
        // Mirror to upper triangle.
        #[allow(clippy::needless_range_loop)]
        for a in 0..n_dof {
            for b in (a + 1)..n_dof {
                h[a][b] = h[b][a];
            }
        }

        let grad_norm = g.iter().map(|v| v * v).sum::<f64>().sqrt();
        if grad_norm < 1e-12 {
            break;
        }

        // Inner damping loop: grow λ until we find an accepted
        // step or hit the cap.
        let mut accepted_this_iter = false;
        while lambda <= config.lambda_max {
            let mut a = h.clone();
            #[allow(clippy::needless_range_loop)]
            for d in 0..n_dof {
                a[d][d] += lambda;
            }
            let Some(delta) = solve_dense(a, neg(&g)) else {
                lambda *= 10.0;
                continue;
            };

            let mut x_cand = x.clone();
            for (xi, di) in x_cand.iter_mut().zip(delta.iter()) {
                *xi += *di;
            }
            let Some(f_cand) = eval.residual(&x_cand) else {
                lambda *= 10.0;
                continue;
            };
            let r_cand = sum_sq(&f_cand);
            // Only check fold-safety on improving candidates —
            // avoids the expensive SI check on every rejected
            // damping attempt.
            if r_cand < residual_sq && eval.fold_safe(&x_cand) {
                x = x_cand;
                f = f_cand;
                residual_sq = r_cand;
                if r_cand < best_residual {
                    best_residual = r_cand;
                    best_x = x.clone();
                }
                lambda = (lambda / 3.0).max(1e-12);
                accepted += 1;
                accepted_this_iter = true;
                break;
            } else {
                lambda *= 10.0;
            }
        }
        if !accepted_this_iter {
            // All damping levels exhausted — converged or stuck.
            break;
        }
    }

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!(
            "[lm {trace_tag}] done: residual={best_residual:.3e} iters={iters} accepted={accepted}"
        );
    }

    LmResult {
        x: best_x,
        final_residual: best_residual,
        iters,
        accepted,
    }
}

fn sum_sq(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

fn neg(v: &[f64]) -> Vec<f64> {
    v.iter().map(|x| -x).collect()
}

/// In-place Gaussian elimination with partial pivoting for a
/// dense square SPD-plus-damping system `A · x = b`. Mirrors
/// `popkern-core::geometry::gregory::direct_linalg::solve_dense`
/// (that helper is `pub(crate)` to its crate, so we keep a local
/// copy here rather than widen its visibility across the fillet
/// boundary).
fn solve_dense(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = m.len();
    for i in 0..n {
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
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            let pivot_row = m[i].clone();
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                m[r][c] -= factor * pc;
            }
            b[r] -= factor * b[i];
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    /// 2-D Rosenbrock-like residual: f = [x-1, y-x²]. Minimum
    /// at (1, 1) with zero residual. Standard NLS sanity test.
    struct Rosenbrock;
    impl LmEval for Rosenbrock {
        fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
            let f0 = x[0] - 1.0;
            let f1 = x[1] - x[0] * x[0];
            Some(vec![f0, f1])
        }
        fn fold_safe(&mut self, _x: &[f64]) -> bool {
            true
        }
    }

    #[test]
    fn lm_solves_rosenbrock_from_origin() {
        let result = run_lm(
            Rosenbrock,
            vec![0.0, 0.0],
            &LmConfig {
                fd_eps: 1e-6,
                lambda_init: 1e-3,
                lambda_max: 1e6,
                max_iter: 50,
                target_residual: 1e-12,
            },
            "rosenbrock",
        );
        assert!(
            result.final_residual < 1e-10,
            "expected convergence, got residual={:.3e}",
            result.final_residual
        );
        assert!((result.x[0] - 1.0).abs() < 1e-5);
        assert!((result.x[1] - 1.0).abs() < 1e-5);
    }

    /// Diagonal linear LSQ: f_i = a_i · x_i - b_i. Converges
    /// in one Gauss-Newton step when λ is small.
    struct DiagonalLsq {
        a: Vec<f64>,
        b: Vec<f64>,
    }
    impl LmEval for DiagonalLsq {
        fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
            let f: Vec<f64> = self
                .a
                .iter()
                .zip(self.b.iter())
                .zip(x.iter())
                .map(|((a, b), x)| a * x - b)
                .collect();
            Some(f)
        }
        fn fold_safe(&mut self, _x: &[f64]) -> bool {
            true
        }
    }

    #[test]
    fn lm_solves_linear_lsq_in_few_iters() {
        let result = run_lm(
            DiagonalLsq {
                a: vec![1.0, 2.0, 3.0],
                b: vec![1.0, 4.0, 9.0],
            },
            vec![0.0, 0.0, 0.0],
            &LmConfig {
                fd_eps: 1e-6,
                lambda_init: 1e-6,
                lambda_max: 1e6,
                max_iter: 20,
                target_residual: 1e-16,
            },
            "diag-lsq",
        );
        assert!(result.final_residual < 1e-12);
        assert!((result.x[0] - 1.0).abs() < 1e-6);
        assert!((result.x[1] - 2.0).abs() < 1e-6);
        assert!((result.x[2] - 3.0).abs() < 1e-6);
    }
}
