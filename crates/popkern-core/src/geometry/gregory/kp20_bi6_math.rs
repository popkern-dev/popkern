//! KP20 bi-6 macro-patch standard formulas.
//!
//! This module intentionally lives beside, not inside, the KP15
//! implementation. KP20 replaces the single-bi-6-per-sector
//! construction with a 2x2 bi-6 macro-patch per sector, but many
//! formulas can be verified independently before the full cap
//! assembler exists.

/// KP20 equation (23) reparameterisation coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20TauParams {
    /// Valence.
    pub n: usize,
    /// `c = cos(2π/n)`.
    pub c: f64,
    /// KP20 table value `κ_n`.
    pub kappa: f64,
    /// `µ_0 = 1 + κ_n / 105`.
    pub mu0: f64,
    /// `µ_1 = 1/3`.
    pub mu1: f64,
    /// `µ_2 = 2/3 - 61c/135 + 4cµ_0/9`.
    pub mu2: f64,
    /// `µ_3 = 29/270 + 2µ_0/9`.
    pub mu3: f64,
}

impl Kp20TauParams {
    /// Build KP20 equation (23) parameters for the paper-tabulated
    /// valences `n ∈ {3, 5, 6, 7, 8, 9, 10}`.
    pub fn for_valence(n: usize) -> Result<Self, String> {
        let kappa = kappa_for_valence(n)
            .ok_or_else(|| format!("KP20 tau table is defined for n = 3, 5..=10; got {n}"))?;
        let c = (core::f64::consts::TAU / n as f64).cos();
        let mu0 = 1.0 + kappa / 105.0;
        let mu1 = 1.0 / 3.0;
        let mu2 = 2.0 / 3.0 - 61.0 * c / 135.0 + 4.0 * c * mu0 / 9.0;
        let mu3 = 29.0 / 270.0 + 2.0 * mu0 / 9.0;
        Ok(Self {
            n,
            c,
            kappa,
            mu0,
            mu1,
            mu2,
            mu3,
        })
    }
}

/// KP20 equation (23) table.
pub fn kappa_for_valence(n: usize) -> Option<f64> {
    match n {
        3 => Some(185.0),
        5 => Some(1033.0),
        6 => Some(352.0),
        7 => Some(-145.0),
        8 => Some(-469.0),
        9 => Some(-687.0),
        10 => Some(-840.0),
        _ => None,
    }
}

/// KP18/KP20 central limit weights for `c_5` and `c_6`.
pub fn kp20_limit_weights(n: usize) -> Result<(f64, f64), String> {
    if n == 3 {
        Ok((5.0 / 96.0, 1.0 / 6.0))
    } else if n > 4 {
        let gamma5 = 1.0 / (n as f64 * (n as f64 + 5.0));
        Ok((gamma5, 4.0 * gamma5))
    } else {
        Err(format!(
            "KP20 limit weights are for extraordinary n = 3 or n > 4; got {n}"
        ))
    }
}

/// KP20 equation (10) `b(u)` BB coefficients for the
/// sector-separating reparameterisation.
pub fn kp20_b_coeffs_near_o(c: f64) -> [f64; 3] {
    [2.0 * c, 2.0 * c, c]
}

/// KP20 equation (10) `b(u)` BB coefficients for the piece closer
/// to the outer midpoint `m`.
pub fn kp20_b_coeffs_near_m(c: f64) -> [f64; 3] {
    [c, 0.0, 0.0]
}

/// KP20 equation (13) strip weights from the first `b` coefficient.
pub fn kp20_g2_strip_weights(b0: f64) -> (f64, f64, f64) {
    let w1 = -b0 / 100.0;
    let w2 = -w1;
    let w3 = b0 / 30.0;
    (w1, w2, w3)
}

/// One-dimensional Hermite-to-degree-6 BB weights for orders 0..3
/// at the lower endpoint.
pub const HERMITE_TO_BI6_W: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0 / 6.0, 0.0, 0.0],
    [1.0, 1.0 / 3.0, 1.0 / 30.0, 0.0],
    [1.0, 1.0 / 2.0, 1.0 / 10.0, 1.0 / 120.0],
];

/// Convert a 4x4 scalar Hermite jet `H[p][q] = ∂^p_u ∂^q_v f(0,0)`
/// to the corresponding 4x4 corner block of a bi-6 BB grid.
pub fn hermite_to_bi6_corner_block(h: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let mut value = 0.0;
            for (p, row) in h.iter().enumerate().take(i + 1) {
                for (q, &hpq) in row.iter().enumerate().take(j + 1) {
                    value += HERMITE_TO_BI6_W[i][p] * HERMITE_TO_BI6_W[j][q] * hpq;
                }
            }
            out[i][j] = value;
        }
    }
    out
}

#[cfg(test)]
#[path = "kp20_bi6_math_tests.rs"]
mod tests;
