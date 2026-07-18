//! KP20 eq (21) normalization fit (PK-134 paper-fidelity archaeology).
//!
//! Appendix A prints the central cubic expansion τ (eq 22) in a
//! normalized chart: the boundary legs' interior BB-coefficients sit
//! at the affine third-points (eq 21) and the whole cubic layer is a
//! one-parameter family in the tabulated µ0 (eq 23). The kernel's
//! paper τ SOLVE (WP2) produces the same construction in its own
//! sector-scale chart, so the printed table must equal the solved
//! cubic jet at `o` composed with a DOMAIN 3-jet ψ — the legitimate
//! reformulation freedom: ψ(o) = o, commuting with the sector's
//! diagonal symmetry (swap a ↔ b in domain and range).
//!
//! Method (PK-134 resolution): µ is FIXED at the printed table and a
//! 9-parameter σ-equivariant domain 3-jet ψ is fitted by Gauss–Newton
//! against the 20 scalar jet equations — 11 equations of
//! over-determination. A residual at solver scale therefore PROVES
//! the printed table is the solved τ up to the reformulation gauge:
//! measured 1e-10..1e-9 at every printed valence n ∈ {3, 5..=10}.
//!
//! Falsified along the way (recorded so nobody re-derives):
//! - fitting µ freely alongside ψ matches the family at 2.9e-10 but
//!   µ0 lands at ≈ 8.9–9.6 — the eq (22) template family is closed
//!   under domain 3-jets, so free µ is GAUGE, not data;
//! - restricting ψ to the quadratic class with the eq (23) couplings
//!   active strands Gauss–Newton in a wrong basin (residual 0.46).

use super::super::kp20_bi6_tau::Kp20TauMap;
use super::super::kp20_bi6_tau_solve_assemble::Kp20PlanarTauPaperSolveError;
use super::super::kp20_paper_tau::{Kp20PaperTauReport, solve_kp20_paper_tau};

/// Monomial order for total-degree-3 2D jets: coefficients of
/// `a^i b^j` for the listed `(i, j)`.
const MONOMIALS: [(usize, usize); 10] = [
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (3, 0),
];

/// A 2-component total-degree-3 polynomial jet in monomial form.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Jet2 {
    /// `m[k][axis]` = coefficient of monomial `MONOMIALS[k]`.
    m: [[f64; 2]; 10],
}

impl Jet2 {
    fn zero() -> Self {
        Self { m: [[0.0; 2]; 10] }
    }
}

fn mono_index(i: usize, j: usize) -> Option<usize> {
    MONOMIALS.iter().position(|&(a, b)| a == i && b == j)
}

/// Multiply two scalar degree-≤3 polynomials, truncating above total
/// degree 3.
fn mul_trunc3(p: &[f64; 10], q: &[f64; 10]) -> [f64; 10] {
    let mut out = [0.0_f64; 10];
    for (kp, &(ip, jp)) in MONOMIALS.iter().enumerate() {
        if p[kp] == 0.0 {
            continue;
        }
        for (kq, &(iq, jq)) in MONOMIALS.iter().enumerate() {
            if q[kq] == 0.0 {
                continue;
            }
            if let Some(k) = mono_index(ip + iq, jp + jq) {
                out[k] += p[kp] * q[kq];
            }
        }
    }
    out
}

/// Compose a jet `s` with the domain map `(psi_a, psi_b)` (scalar
/// polynomials with zero constant term), truncated to total degree 3.
fn compose(s: &Jet2, psi_a: &[f64; 10], psi_b: &[f64; 10]) -> Jet2 {
    // Powers of psi_a and psi_b up to 3, truncated.
    let one = {
        let mut p = [0.0_f64; 10];
        p[0] = 1.0;
        p
    };
    let mut pow_a = [one, *psi_a, mul_trunc3(psi_a, psi_a), [0.0; 10]];
    pow_a[3] = mul_trunc3(&pow_a[2], psi_a);
    let mut pow_b = [one, *psi_b, mul_trunc3(psi_b, psi_b), [0.0; 10]];
    pow_b[3] = mul_trunc3(&pow_b[2], psi_b);

    let mut out = Jet2::zero();
    for (k, &(i, j)) in MONOMIALS.iter().enumerate() {
        let basis = mul_trunc3(&pow_a[i], &pow_b[j]);
        for axis in 0..2 {
            let c = s.m[k][axis];
            if c == 0.0 {
                continue;
            }
            for t in 0..10 {
                out.m[t][axis] += c * basis[t];
            }
        }
    }
    out
}

/// ψ parameters: linear `[x, y]`, quadratic `[q1, q2, q3]`, cubic
/// `[c1, c2, c3, c4]`, all σ-equivariant (swap a↔b in domain AND
/// range):
///   ψ_a = x·a + y·b + q1·a² + q2·ab + q3·b² + c1·a³ + c2·a²b + c3·ab² + c4·b³
///   ψ_b = y·a + x·b + q3·a² + q2·ab + q1·b² + c4·a³ + c3·a²b + c2·ab² + c1·b³
fn psi_from_params(theta: &[f64; 9]) -> ([f64; 10], [f64; 10]) {
    let [x, y, q1, q2, q3, c1, c2, c3, c4] = *theta;
    let mut a = [0.0_f64; 10];
    let mut b = [0.0_f64; 10];
    let set = |arr: &mut [f64; 10], i: usize, j: usize, v: f64| {
        if let Some(k) = mono_index(i, j) {
            arr[k] = v;
        }
    };
    set(&mut a, 1, 0, x);
    set(&mut a, 0, 1, y);
    set(&mut a, 2, 0, q1);
    set(&mut a, 1, 1, q2);
    set(&mut a, 0, 2, q3);
    set(&mut a, 3, 0, c1);
    set(&mut a, 2, 1, c2);
    set(&mut a, 1, 2, c3);
    set(&mut a, 0, 3, c4);
    set(&mut b, 1, 0, y);
    set(&mut b, 0, 1, x);
    set(&mut b, 2, 0, q3);
    set(&mut b, 1, 1, q2);
    set(&mut b, 0, 2, q1);
    set(&mut b, 3, 0, c4);
    set(&mut b, 2, 1, c3);
    set(&mut b, 1, 2, c2);
    set(&mut b, 0, 3, c1);
    (a, b)
}

/// Monomial jet of a [`Kp20TauMap`] at the origin corner of its
/// domain triangle (its BB corner `(0,0)`).
fn tau_map_jet(map: &Kp20TauMap) -> Jet2 {
    let jet = map.jet3(0.0, 0.0);
    let mut out = Jet2::zero();
    let fact = |k: usize| -> f64 { (1..=k).map(|v| v as f64).product::<f64>().max(1.0) };
    for (k, &(i, j)) in MONOMIALS.iter().enumerate() {
        for axis in 0..2 {
            out.m[k][axis] = jet.d[i][j][axis] / (fact(i) * fact(j));
        }
    }
    out
}

/// Monomial jet of the SOLVED τ at `o`, in sector scale, with the
/// `(a, b)` index assignment optionally transposed (the solve chart
/// may transpose the printed convention — same ambiguity the
/// `jet_delta` measure resolves by trying both).
fn solved_jet(report: &Kp20PaperTauReport, swap: bool) -> Jet2 {
    use super::super::kp20_bi6_cap::Kp20PieceId;
    let partials = super::bi6_corner_partials(report.tau.piece(Kp20PieceId::NearOriginLeft));
    let mut out = Jet2::zero();
    let fact = |k: usize| -> f64 { (1..=k).map(|v| v as f64).product::<f64>().max(1.0) };
    for (k, &(i, j)) in MONOMIALS.iter().enumerate() {
        let (pi, pj) = if swap { (j, i) } else { (i, j) };
        for axis in 0..2 {
            // Piece P1 spans [0,½]²: sector-scale derivative =
            // piece derivative × 2^(order).
            out.m[k][axis] =
                partials[pi][pj][axis] * 2.0_f64.powi((i + j) as i32) / (fact(i) * fact(j));
        }
    }
    out
}

/// Result of the joint (ψ, µ) fit for one valence.
#[derive(Debug, Clone, Copy)]
pub struct Kp20Eq21MuFit {
    /// Valence.
    pub n: usize,
    /// Sup-norm of the final 20-equation residual.
    pub residual_inf: f64,
    /// Fitted `[µ0, µ1, µ2, µ3]`.
    pub mu: [f64; 4],
    /// Whether the transposed index assignment won.
    pub swapped: bool,
}

/// Residual vector: composed solved jet minus the eq (22) template
/// with the given µ.
fn eq23_mu(n: usize, mu0: f64) -> [f64; 4] {
    let c = (core::f64::consts::TAU / n as f64).cos();
    [
        mu0,
        1.0 / 3.0,
        2.0 / 3.0 - 61.0 * c / 135.0 + 4.0 * c * mu0 / 9.0,
        29.0 / 270.0 + 2.0 * mu0 / 9.0,
    ]
}

fn residual(n: usize, solved: &Jet2, theta: &[f64; 9], mu: &[f64; 4]) -> Result<[f64; 20], String> {
    let (pa, pb) = psi_from_params(theta);
    let composed = compose(solved, &pa, &pb);
    let template = Kp20TauMap::for_sector_paper_with_mu(n, 0, mu[0], mu[1], mu[2], mu[3])?;
    let target = tau_map_jet(&template);
    let mut r = [0.0_f64; 20];
    for k in 0..10 {
        for axis in 0..2 {
            r[2 * k + axis] = composed.m[k][axis] - target.m[k][axis];
        }
    }
    Ok(r)
}

/// Fit (ψ, µ) for one valence by Gauss–Newton with finite-difference
/// Jacobians. Tries both index assignments and keeps the better fit.
pub fn kp20_paper_tau_eq21_mu_fit(n: usize) -> Result<Kp20Eq21MuFit, Kp20PlanarTauPaperSolveError> {
    let report = solve_kp20_paper_tau(n)?;
    let printed = super::super::kp20_bi6_math::Kp20TauParams::for_valence(n)
        .map_err(Kp20PlanarTauPaperSolveError::TensorBorder)?;
    let mut best: Option<Kp20Eq21MuFit> = None;
    for swap in [false, true] {
        let solved = solved_jet(&report, swap);
        // Seed: identity-scale ψ from the solved linear diagonal, µ
        // from the printed table.
        let Some(lin_idx) = mono_index(1, 0) else {
            continue;
        };
        let lin = solved.m[lin_idx];
        let scale = lin[1].abs().max(lin[0].abs()).max(1e-12);
        let mut theta = [1.0 / scale, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mu = eq23_mu(n, printed.mu0);
        let mut r = match residual(n, &solved, &theta, &mu) {
            Ok(r) => r,
            Err(_) => continue,
        };
        for _iter in 0..60 {
            // Finite-difference Jacobian over the ACTIVE unknowns
            // (quadratic ψ: 5 + 4 µ; cubic ψ frozen at zero).
            let mut jac = [[0.0_f64; 9]; 20];
            let h = 1.0e-7;
            for p in 0..9 {
                let mut theta_p = theta;
                theta_p[p] += h;
                let rp = match residual(n, &solved, &theta_p, &mu) {
                    Ok(r) => r,
                    Err(_) => break,
                };
                for (row, jrow) in jac.iter_mut().enumerate() {
                    jrow[p] = (rp[row] - r[row]) / h;
                }
            }
            // Normal equations (9×9) with light Tikhonov.
            let mut ata = [[0.0_f64; 9]; 9];
            let mut atb = [0.0_f64; 9];
            for i in 0..9 {
                for j in 0..9 {
                    for jac_row in jac.iter() {
                        ata[i][j] += jac_row[i] * jac_row[j];
                    }
                }
                for (row, jac_row) in jac.iter().enumerate() {
                    atb[i] -= jac_row[i] * r[row];
                }
                ata[i][i] += 1.0e-14;
            }
            let step = match solve_9(ata, atb) {
                Some(s) => s,
                None => break,
            };
            for p in 0..9 {
                theta[p] += step[p];
            }
            r = match residual(n, &solved, &theta, &mu) {
                Ok(r) => r,
                Err(_) => break,
            };
            let step_inf = step.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            if step_inf < 1.0e-14 {
                break;
            }
        }
        let residual_inf = r.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
        let fit = Kp20Eq21MuFit {
            n,
            residual_inf,
            mu,
            swapped: swap,
        };
        if best.is_none_or(|b| residual_inf < b.residual_inf) {
            best = Some(fit);
        }
    }
    best.ok_or(Kp20PlanarTauPaperSolveError::TensorBorder(
        "eq21 fit: no candidate".into(),
    ))
}

/// Gaussian elimination with partial pivoting for the 9×9 normal
/// equations.
fn solve_9(mut a: [[f64; 9]; 9], mut b: [f64; 9]) -> Option<[f64; 9]> {
    for col in 0..9 {
        let mut piv = col;
        for row in (col + 1)..9 {
            if a[row][col].abs() > a[piv][col].abs() {
                piv = row;
            }
        }
        if a[piv][col].abs() < 1.0e-18 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for row in (col + 1)..9 {
            let f = a[row][col] / a[col][col];
            for k in col..9 {
                a[row][k] -= f * a[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = [0.0_f64; 9];
    for col in (0..9).rev() {
        let mut acc = b[col];
        for k in (col + 1)..9 {
            acc -= a[col][k] * x[k];
        }
        x[col] = acc / a[col][col];
    }
    Some(x)
}

#[cfg(test)]
#[path = "kp20_paper_tau_eq21_fit_tests.rs"]
mod tests;
