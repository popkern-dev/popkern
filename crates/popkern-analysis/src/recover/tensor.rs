//! Symmetric stress tensors in Voigt order (xx, yy, zz, yz, xz, xy):
//! von Mises, principal stresses and the complementary energy density.

use popkern_math::PortableMath;

/// A symmetric stress or strain tensor, Voigt order
/// (xx, yy, zz, yz, xz, xy), tensor (not engineering) shears.
pub type Voigt = [f64; 6];

fn sq(x: f64) -> f64 {
    x * x
}

/// Von Mises equivalent stress.
pub fn von_mises(s: &Voigt) -> f64 {
    let normal = sq(s[0] - s[1]) + sq(s[1] - s[2]) + sq(s[2] - s[0]);
    (0.5 * normal + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])).sqrt()
}

/// Principal stresses, descending, by the closed-form trigonometric
/// solution of the characteristic cubic.
pub fn principal(s: &Voigt) -> [f64; 3] {
    let (xx, yy, zz, yz, xz, xy) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    let p1 = xy * xy + xz * xz + yz * yz;
    let q = (xx + yy + zz) / 3.0;
    let p2 = sq(xx - q) + sq(yy - q) + sq(zz - q) + 2.0 * p1;
    if p2 <= 0.0 {
        return [q; 3];
    }
    let p = (p2 / 6.0).sqrt();
    let b = [
        (xx - q) / p,
        (yy - q) / p,
        (zz - q) / p,
        yz / p,
        xz / p,
        xy / p,
    ];
    let det = b[0] * (b[1] * b[2] - b[3] * b[3]) - b[5] * (b[5] * b[2] - b[3] * b[4])
        + b[4] * (b[5] * b[3] - b[1] * b[4]);
    let r = (0.5 * det).clamp(-1.0, 1.0);
    let phi = r.pacos() / 3.0;
    let l1 = q + 2.0 * p * phi.pcos();
    let l3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).pcos();
    [l1, 3.0 * q - l1 - l3, l3]
}

/// Complementary energy density σ : C⁻¹ : σ for an isotropic material,
/// (1/E)[(1 + ν) σ:σ − ν (tr σ)²].
pub fn energy_density(s: &Voigt, youngs: f64, poisson: f64) -> f64 {
    let ss =
        s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
    let tr = s[0] + s[1] + s[2];
    ((1.0 + poisson) * ss - poisson * tr * tr) / youngs
}

/// Stress from a displacement gradient H[i][j] = ∂u_i/∂x_j.
pub fn stress_from_gradient(h: &[[f64; 3]; 3], lambda: f64, mu: f64) -> Voigt {
    let tr = h[0][0] + h[1][1] + h[2][2];
    [
        lambda * tr + 2.0 * mu * h[0][0],
        lambda * tr + 2.0 * mu * h[1][1],
        lambda * tr + 2.0 * mu * h[2][2],
        mu * (h[1][2] + h[2][1]),
        mu * (h[0][2] + h[2][0]),
        mu * (h[0][1] + h[1][0]),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn principal_stresses_of_known_tensors() {
        let d = principal(&[3.0, -1.0, 2.0, 0.0, 0.0, 0.0]);
        assert!(
            (d[0] - 3.0).abs() < 1e-14 && (d[1] - 2.0).abs() < 1e-14 && (d[2] + 1.0).abs() < 1e-14
        );
        // Pure shear τ in xy: principal ±τ and 0.
        let s = principal(&[0.0, 0.0, 0.0, 0.0, 0.0, 2.5]);
        assert!((s[0] - 2.5).abs() < 1e-14 && s[1].abs() < 1e-14 && (s[2] + 2.5).abs() < 1e-14);
        assert!((von_mises(&[0.0, 0.0, 0.0, 0.0, 0.0, 2.5]) - 2.5 * 3f64.sqrt()).abs() < 1e-14);
        // Uniaxial: σ_vM = |σ|.
        assert!((von_mises(&[7.0, 0.0, 0.0, 0.0, 0.0, 0.0]) - 7.0).abs() < 1e-14);
        // A general tensor: invariants are preserved.
        let t = [1.0, 2.0, -0.5, 0.3, -0.7, 0.9];
        let p = principal(&t);
        assert!((p.iter().sum::<f64>() - (1.0 + 2.0 - 0.5)).abs() < 1e-13);
        assert!(p[0] >= p[1] && p[1] >= p[2]);
    }
}
