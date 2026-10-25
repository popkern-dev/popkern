//! The stress field inside one cell: σ_h(ξ) = D(B(ξ)u + G̃(ξ)α), a
//! polynomial of degree ≤ 1 per local axis, so it can be evaluated
//! anywhere in the cell and integrated exactly against the moments.

use crate::element::poly::{Poly27, eval, shape_derivative};
use crate::element::vector::{NODAL_DOFS, VectorCell};
use crate::moments::moment_index;

use super::tensor::Voigt;

/// The six stress components of a cell as polynomials in local
/// coordinates.
pub type StressPoly = [Poly27; 6];

/// Shape-function derivative polynomials, built once.
#[derive(Debug, Clone)]
pub struct DerivativeTables {
    dn: [[Poly27; 3]; 8],
}

impl Default for DerivativeTables {
    fn default() -> Self {
        Self::new()
    }
}

impl DerivativeTables {
    /// Builds the tables.
    pub fn new() -> DerivativeTables {
        DerivativeTables {
            dn: std::array::from_fn(|a| std::array::from_fn(|k| shape_derivative(a, k))),
        }
    }
}

/// Fragment centroid in local coordinates, from the moments, clamped
/// to the cell (a near-empty fragment's m₁₀₀/m₀₀₀ is round-off and can
/// land far outside it).
pub fn centroid(m: &[f64; 27]) -> [f64; 3] {
    if m[0] <= 0.0 {
        return [0.0; 3];
    }
    [
        m[moment_index(1, 0, 0)] / m[0],
        m[moment_index(0, 1, 0)] / m[0],
        m[moment_index(0, 0, 1)] / m[0],
    ]
    .map(|c| {
        if c.is_finite() {
            c.clamp(-1.0, 1.0)
        } else {
            0.0
        }
    })
}

/// σ_h as polynomials for nodal displacements `u` of a cell of size
/// `h` with moments `m`, using the cell operator's recovery (when EAS
/// was used) for the enhanced dofs.
pub fn stress_poly(
    t: &DerivativeTables,
    op: &VectorCell,
    m: &[f64; 27],
    u: &[f64; NODAL_DOFS],
    h: f64,
    lambda: f64,
    mu: f64,
) -> StressPoly {
    let scale = 2.0 / h;
    // H[i][j] = ∂u_i/∂x_j as polynomials.
    let mut grad: [[Poly27; 3]; 3] = [[[0.0; 27]; 3]; 3];
    for a in 0..8 {
        for i in 0..3 {
            let ua = u[3 * a + i];
            if ua == 0.0 {
                continue;
            }
            for (j, slot) in grad[i].iter_mut().enumerate() {
                for (x, y) in slot.iter_mut().zip(&t.dn[a][j]) {
                    *x += scale * ua * y;
                }
            }
        }
    }
    if let Some(recovery) = &op.recovery {
        let mean = centroid(m);
        for (row, rec) in recovery.iter().enumerate() {
            let alpha: f64 = rec.iter().zip(u).map(|(r, x)| r * x).sum();
            let (bubble, i) = (row / 3, row % 3);
            // ∂/∂x_m of α(1 − ξ_m²), relative to the fragment mean:
            // −2α (ξ_m − ξ̄_m)·(2/h).
            let slot = &mut grad[i][bubble];
            let mut e = [0; 3];
            e[bubble] = 1;
            slot[moment_index(e[0], e[1], e[2])] += -2.0 * scale * alpha;
            slot[0] += 2.0 * scale * alpha * mean[bubble];
        }
    }
    let combine = |k: usize| -> Poly27 {
        std::array::from_fn(|c| {
            let h = |i: usize, j: usize| grad[i][j][c];
            let tr = h(0, 0) + h(1, 1) + h(2, 2);
            match k {
                0 => lambda * tr + 2.0 * mu * h(0, 0),
                1 => lambda * tr + 2.0 * mu * h(1, 1),
                2 => lambda * tr + 2.0 * mu * h(2, 2),
                3 => mu * (h(1, 2) + h(2, 1)),
                4 => mu * (h(0, 2) + h(2, 0)),
                _ => mu * (h(0, 1) + h(1, 0)),
            }
        })
    };
    std::array::from_fn(combine)
}

/// σ_h at local point `xi`.
pub fn stress_at(p: &StressPoly, xi: [f64; 3]) -> Voigt {
    std::array::from_fn(|k| eval(&p[k], xi))
}

/// Adds a thermal strain's part to σ_h (`stress_poly`'s output): the
/// enhanced dofs' affine part α₀ (from the cell's material `parts`) and
/// −(3λ + 2μ) α θ I with the (λ, μ) the stress is read in
/// (`element::thermal_strain`).
#[allow(clippy::too_many_arguments)]
pub fn add_thermal(
    poly: &mut StressPoly,
    op: &VectorCell,
    m: &[f64; 27],
    parts: &[(&[f64; 27], f64, f64)],
    h: f64,
    (lambda, mu): (f64, f64),
    alpha: f64,
    theta: &[f64; 8],
) {
    use crate::element::thermal_strain::{cell_thermal_parts, theta_poly};
    let t = cell_thermal_parts(op, m, parts, h, alpha, theta);
    let scale = 2.0 / h;
    let mean = crate::element::vector::fragment_mean(m);
    // The α₀ bubbles' gradient: ∂u_i/∂x_m = α₀ · (−2)(ξ_m − ξ̄_m)·(2/h).
    let mut grad: [[Poly27; 3]; 3] = [[[0.0; 27]; 3]; 3];
    for (r, a0) in t.alpha0.iter().enumerate() {
        let (bubble, i) = (r / 3, r % 3);
        let mut e = [0; 3];
        e[bubble] = 1;
        let slot = &mut grad[i][bubble];
        slot[moment_index(e[0], e[1], e[2])] += -2.0 * scale * a0;
        slot[0] += 2.0 * scale * a0 * mean[bubble];
    }
    let tp = theta_poly(theta);
    let c = (3.0 * lambda + 2.0 * mu) * alpha;
    for (k, p) in poly.iter_mut().enumerate() {
        for (q, x) in p.iter_mut().enumerate() {
            let g = |i: usize, j: usize| grad[i][j][q];
            let tr = g(0, 0) + g(1, 1) + g(2, 2);
            *x += match k {
                0 => lambda * tr + 2.0 * mu * g(0, 0) - c * tp[q],
                1 => lambda * tr + 2.0 * mu * g(1, 1) - c * tp[q],
                2 => lambda * tr + 2.0 * mu * g(2, 2) - c * tp[q],
                3 => mu * (g(1, 2) + g(2, 1)),
                4 => mu * (g(0, 2) + g(2, 0)),
                _ => mu * (g(0, 1) + g(1, 0)),
            };
        }
    }
}
