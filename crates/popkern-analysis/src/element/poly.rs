//! Polynomials in cell-local (ξ, η, ζ) with every exponent ≤ 2, stored
//! as 27 coefficients in moment order (`a + 3b + 9c`), so that
//! ∫_frag p dV = Σ_k p[k] · m_k for the fragment moments m.

use std::sync::OnceLock;

use crate::grid::CELL_CORNERS;
use crate::moments::moment_index;

/// A polynomial with exponents ≤ 2 per axis.
pub type Poly27 = [f64; 27];

/// The product of two polynomials whose product stays within degree 2
/// per axis. Terms beyond are a caller error; they are dropped and
/// flagged by the debug assertion.
pub fn mul(p: &Poly27, q: &Poly27) -> Poly27 {
    let mut out = [0.0; 27];
    for c1 in 0..3 {
        for b1 in 0..3 {
            for a1 in 0..3 {
                let x = p[moment_index(a1, b1, c1)];
                if x == 0.0 {
                    continue;
                }
                for c2 in 0..3 {
                    for b2 in 0..3 {
                        for a2 in 0..3 {
                            let y = q[moment_index(a2, b2, c2)];
                            if y == 0.0 {
                                continue;
                            }
                            let (a, b, c) = (a1 + a2, b1 + b2, c1 + c2);
                            debug_assert!(a <= 2 && b <= 2 && c <= 2, "degree overflow");
                            if a <= 2 && b <= 2 && c <= 2 {
                                out[moment_index(a, b, c)] += x * y;
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

/// Σ_k p[k] m[k]: the integral of `p` against moments `m`.
pub fn integrate(p: &Poly27, m: &[f64; 27]) -> f64 {
    p.iter().zip(m).map(|(x, y)| x * y).sum()
}

/// Sign of corner `a` on each axis: −1 for offset 0, +1 for offset 1.
pub fn corner_signs(a: usize) -> [f64; 3] {
    CELL_CORNERS[a].map(|o| if o == 0 { -1.0 } else { 1.0 })
}

/// (1 + s t) as a polynomial in the variable of `axis`.
fn linear_factor(axis: usize, s: f64) -> Poly27 {
    let mut p = [0.0; 27];
    p[0] = 1.0;
    let mut e = [0; 3];
    e[axis] = 1;
    p[moment_index(e[0], e[1], e[2])] = s;
    p
}

/// Each corner's shape function and its three derivatives, built once:
/// every cell, vertex and sample reads them.
fn table() -> &'static [[Poly27; 4]; 8] {
    static TABLE: OnceLock<[[Poly27; 4]; 8]> = OnceLock::new();
    TABLE.get_or_init(|| {
        std::array::from_fn(|a| {
            [
                shape_product(a),
                derivative_product(a, 0),
                derivative_product(a, 1),
                derivative_product(a, 2),
            ]
        })
    })
}

/// The trilinear shape function of corner `a` (0..8, zero otherwise):
/// N_a = ⅛ (1 + ξ_a ξ)(1 + η_a η)(1 + ζ_a ζ).
pub fn shape(a: usize) -> Poly27 {
    table().get(a).map_or([0.0; 27], |t| t[0])
}

/// ∂N_a/∂(local coordinate `axis`) (zero outside 0..8 and 0..3).
pub fn shape_derivative(a: usize, axis: usize) -> Poly27 {
    match (table().get(a), axis) {
        (Some(t), 0..3) => t[1 + axis],
        _ => [0.0; 27],
    }
}

/// [`shape`] as the product of its three factors.
fn shape_product(a: usize) -> Poly27 {
    let s = corner_signs(a);
    let xy = mul(&linear_factor(0, s[0]), &linear_factor(1, s[1]));
    let mut n = mul(&xy, &linear_factor(2, s[2]));
    for c in &mut n {
        *c *= 0.125;
    }
    n
}

/// [`shape_derivative`] as the product of its three factors.
fn derivative_product(a: usize, axis: usize) -> Poly27 {
    let s = corner_signs(a);
    let mut factors = [
        linear_factor(0, s[0]),
        linear_factor(1, s[1]),
        linear_factor(2, s[2]),
    ];
    let mut d = [0.0; 27];
    d[0] = s[axis];
    factors[axis] = d;
    let mut out = mul(&mul(&factors[0], &factors[1]), &factors[2]);
    for c in &mut out {
        *c *= 0.125;
    }
    out
}

/// Evaluates `p` at local point `x`.
pub fn eval(p: &Poly27, x: [f64; 3]) -> f64 {
    let pw = |t: f64| [1.0, t, t * t];
    let (px, py, pz) = (pw(x[0]), pw(x[1]), pw(x[2]));
    let mut sum = 0.0;
    for c in 0..3 {
        for b in 0..3 {
            for a in 0..3 {
                sum += p[moment_index(a, b, c)] * px[a] * py[b] * pz[c];
            }
        }
    }
    sum
}

#[cfg(test)]
#[path = "poly_tests.rs"]
mod tests;
