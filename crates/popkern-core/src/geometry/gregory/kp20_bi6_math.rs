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
    /// `µ_0 = 1 + κ_n / 10⁵`.
    ///
    /// The divisor is ten to the fifth (100000), per the printed
    /// exponent in KP20 eq (23) page 10 (verified at 300 dpi). A
    /// pdftotext extraction flattens the superscript to `105`, which
    /// makes µ₀ ≈ 10.84 for n = 5 instead of the correct ≈ 1.01033
    /// and turns τ¹ from a near-identity perturbation into a ~10×
    /// distortion. See
    /// `docs/architecture/kp20-paper-faithful-recovery-plan.md` (RC1).
    pub mu0: f64,
    /// `µ_1 = 1/3`.
    pub mu1: f64,
    /// `µ_2 = 2/3 - 61c/135 + 4cµ_0/9`.
    pub mu2: f64,
    /// `µ_3 = 29/270 + 2µ_0/9`.
    pub mu3: f64,
}

/// KP18/KP20 Catmull-Clark characteristic-ring scalars.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20CharacteristicScalars {
    /// Valence.
    pub n: usize,
    /// `c = cos(2π/n)`.
    pub c: f64,
    /// Catmull-Clark subdominant eigenvalue
    /// `(c + 5 + sqrt((c + 1)(c + 9))) / 16`.
    pub lambda: f64,
    /// `a_n = (1 + c) / 4`.
    pub a: f64,
    /// `b_n = (1 + c) / 8`.
    pub b: f64,
    /// One KP18/Stam normalisation of the six-entry characteristic
    /// eigenvector at Fourier frequency one.
    pub eigenvector: [f64; 6],
}

impl Kp20TauParams {
    /// Build KP20 equation (23) parameters for the paper-tabulated
    /// valences `n ∈ {3, 5, 6, 7, 8, 9, 10}`.
    pub fn for_valence(n: usize) -> Result<Self, String> {
        let kappa = kappa_for_valence(n)
            .ok_or_else(|| format!("KP20 tau table is defined for n = 3, 5..=10; got {n}"))?;
        let c = (core::f64::consts::TAU / n as f64).cos();
        let mu0 = 1.0 + kappa / 1.0e5;
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

impl Kp20CharacteristicScalars {
    /// Build KP18/KP20 characteristic-ring scalars for `n >= 3`.
    pub fn for_valence(n: usize) -> Result<Self, String> {
        if n < 3 {
            return Err(format!(
                "KP20 characteristic-ring scalars require n >= 3; got {n}"
            ));
        }
        let c = (core::f64::consts::TAU / n as f64).cos();
        let lambda = (c + 5.0 + ((c + 1.0) * (c + 9.0)).sqrt()) / 16.0;
        let a = (1.0 + c) / 4.0;
        let b = (1.0 + c) / 8.0;
        let eigenvector = kp20_stam_characteristic_eigenvector(lambda, a, b);
        Ok(Self {
            n,
            c,
            lambda,
            a,
            b,
            eigenvector,
        })
    }
}

/// KP18/PR08 Catmull-Clark subdominant eigenvalue `λ_n`.
pub fn kp20_cc_subdominant_eigenvalue(n: usize) -> Result<f64, String> {
    Ok(Kp20CharacteristicScalars::for_valence(n)?.lambda)
}

/// KP18/Stam six-entry characteristic eigenvector normalisation.
pub fn kp20_stam_characteristic_eigenvector(lambda: f64, a: f64, b: f64) -> [f64; 6] {
    let edge = (a - lambda) * (12.0 * lambda * lambda - 8.0 * lambda + 1.0);
    let spoke = 4.0 * lambda * (3.0 * lambda - 1.0);
    let face = 16.0 * lambda * lambda * b;
    [edge, spoke, face, spoke, edge, face]
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

/// KP20 free scalar count before splitting each sector into a 2x2
/// macro patch.
///
/// The count is `7n + 10`: seven scalar degrees of freedom per
/// sector plus the central ten-dimensional cap block. KP20 keeps this
/// count as the baseline cap space before the midpoint refinement adds
/// interior freedom.
pub fn kp20_s0_free_count(n: usize) -> Result<usize, String> {
    if n < 3 {
        return Err(format!("KP20 free-count formulas require n >= 3; got {n}"));
    }
    Ok(7 * n + 10)
}

/// KP20 free scalar count after the 2x2 macro-patch refinement.
///
/// The refined count is `16n + 10`, so the refinement contributes
/// nine additional scalar degrees of freedom per sector while
/// preserving the same central ten-dimensional block.
pub fn kp20_refined_s0_free_count(n: usize) -> Result<usize, String> {
    if n < 3 {
        return Err(format!("KP20 free-count formulas require n >= 3; got {n}"));
    }
    Ok(16 * n + 10)
}

/// Additional scalar freedom introduced by KP20's midpoint
/// refinement, equal to `(16n + 10) - (7n + 10) = 9n`.
pub fn kp20_refinement_added_free_count(n: usize) -> Result<usize, String> {
    Ok(kp20_refined_s0_free_count(n)? - kp20_s0_free_count(n)?)
}

/// KP18 total guide scalar freedom, including the ten cubic-expansion
/// coefficients.
pub fn kp18_guide_total_free_count(n: usize) -> Result<usize, String> {
    if n < 3 {
        return Err(format!(
            "KP18 guide free-count formulas require n >= 3; got {n}"
        ));
    }
    Ok(12 * n + 10)
}

/// KP18 guide scalar freedom solved by the outer-coefficient LSQ after
/// the cubic expansion is fixed.
pub fn kp18_guide_lsq_free_count(n: usize) -> Result<usize, String> {
    if n < 3 {
        return Err(format!(
            "KP18 guide free-count formulas require n >= 3; got {n}"
        ));
    }
    Ok(12 * n + 9)
}

/// KP18 guide boundary-coupling free scalar count.
pub fn kp18_guide_boundary_free_count(n: usize) -> Result<usize, String> {
    if n < 3 {
        return Err(format!(
            "KP18 guide free-count formulas require n >= 3; got {n}"
        ));
    }
    Ok(6 * n)
}

/// KP18 guide interior-only free scalar count.
pub fn kp18_guide_interior_free_count(n: usize) -> Result<usize, String> {
    if n < 3 {
        return Err(format!(
            "KP18 guide free-count formulas require n >= 3; got {n}"
        ));
    }
    Ok(6 * n)
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

/// KP18 Fig. 5 scalar C² join rule for the coefficient immediately
/// left of a Bezier-curve join.
pub fn kp18_c2_join_left(a_d_minus_2: f64, a_d: f64, b_2: f64) -> f64 {
    0.25 * a_d_minus_2 + a_d - 0.25 * b_2
}

/// KP18 Fig. 5 mirrored scalar C² join rule for the coefficient
/// immediately right of a Bezier-curve join.
pub fn kp18_c2_join_right(a_d_minus_2: f64, b_0: f64, b_2: f64) -> f64 {
    0.25 * b_2 + b_0 - 0.25 * a_d_minus_2
}

/// One-dimensional Hermite-to-degree-5 BB weights for orders 0..2
/// at the lower endpoint.
pub const HERMITE_TO_BI5_W: [[f64; 3]; 3] = [
    [1.0, 0.0, 0.0],
    [1.0, 1.0 / 5.0, 0.0],
    [1.0, 2.0 / 5.0, 1.0 / 20.0],
];

/// Convert a 3x3 scalar Hermite jet `H[p][q] = ∂^p_u ∂^q_v f(0,0)`
/// to the corresponding 3x3 corner block of a bi-5 BB grid.
pub fn hermite_to_bi5_corner_block(h: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut value = 0.0;
            for (p, row) in h.iter().enumerate().take(i + 1) {
                for (q, &hpq) in row.iter().enumerate().take(j + 1) {
                    value += HERMITE_TO_BI5_W[i][p] * HERMITE_TO_BI5_W[j][q] * hpq;
                }
            }
            out[i][j] = value;
        }
    }
    out
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
