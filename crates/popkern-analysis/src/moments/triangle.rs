//! Exact integrals of monomials over a triangle.
//!
//! Each coordinate is a linear form in the barycentrics λ₀, λ₁, λ₂, so
//! a monomial ξ^α η^β ζ^γ is a homogeneous polynomial of degree
//! d = α + β + γ in λ, and
//! ∫_T λ₀^p λ₁^q λ₂^r dA = 2A · p! q! r! / (d + 2)!.
//! The result is exact up to round-off.

/// Highest power of ξ needed (G_a for a ≤ 2 has degree 3).
pub const MAX_XI: usize = 3;
/// Highest power of η and ζ needed.
pub const MAX_YZ: usize = 2;

const MAX_DEGREE: usize = MAX_XI + 2 * MAX_YZ;

/// n! for n up to MAX_DEGREE + 2.
const FACTORIAL: [f64; MAX_DEGREE + 3] = {
    let mut f = [1.0; MAX_DEGREE + 3];
    let mut i = 1;
    while i < f.len() {
        f[i] = f[i - 1] * i as f64;
        i += 1;
    }
    f
};

/// Position of the coefficient of λ₀^p λ₁^q λ₂^(deg−p−q): rows of
/// increasing p, row p holding q = 0..=deg−p.
fn slot(deg: usize, p: usize, q: usize) -> usize {
    p * (deg + 1) - p * p.saturating_sub(1) / 2 + q
}

/// A homogeneous polynomial of degree `deg` in λ₀, λ₁, λ₂.
#[derive(Debug, Clone)]
struct Homogeneous {
    deg: usize,
    c: Vec<f64>,
}

impl Homogeneous {
    fn one() -> Homogeneous {
        Homogeneous {
            deg: 0,
            c: vec![1.0],
        }
    }

    /// x₀λ₀ + x₁λ₁ + x₂λ₂.
    fn linear(x: [f64; 3]) -> Homogeneous {
        let mut c = vec![0.0; 3];
        c[slot(1, 1, 0)] = x[0];
        c[slot(1, 0, 1)] = x[1];
        c[slot(1, 0, 0)] = x[2];
        Homogeneous { deg: 1, c }
    }

    fn times(&self, other: &Homogeneous) -> Homogeneous {
        let deg = self.deg + other.deg;
        let mut c = vec![0.0; (deg + 1) * (deg + 2) / 2];
        for p1 in 0..=self.deg {
            for q1 in 0..=(self.deg - p1) {
                let a = self.c[slot(self.deg, p1, q1)];
                for p2 in 0..=other.deg {
                    for q2 in 0..=(other.deg - p2) {
                        c[slot(deg, p1 + p2, q1 + q2)] += a * other.c[slot(other.deg, p2, q2)];
                    }
                }
            }
        }
        Homogeneous { deg, c }
    }

    /// ∫ over a triangle of area `area`.
    fn integrate(&self, area: f64) -> f64 {
        let d = self.deg;
        let mut sum = 0.0;
        for (p, fp) in FACTORIAL.iter().enumerate().take(d + 1) {
            for (q, fq) in FACTORIAL.iter().enumerate().take(d - p + 1) {
                sum += self.c[slot(d, p, q)] * fp * fq * FACTORIAL[d - p - q];
            }
        }
        2.0 * area * sum / FACTORIAL[d + 2]
    }
}

/// ∫ ξ^α η^β ζ^γ dA over the triangle with local-coordinate vertices
/// `v` and physical area `area`, for α ≤ 3, β ≤ 2, γ ≤ 2, indexed
/// `[α][β][γ]`.
pub fn monomial_integrals(v: [[f64; 3]; 3], area: f64) -> [[[f64; 3]; 3]; 4] {
    let powers = |k: usize, n: usize| {
        let l = Homogeneous::linear([v[0][k], v[1][k], v[2][k]]);
        let mut out = vec![Homogeneous::one()];
        for i in 1..=n {
            let next = out[i - 1].times(&l);
            out.push(next);
        }
        out
    };
    let xi = powers(0, MAX_XI);
    let eta = powers(1, MAX_YZ);
    let zeta = powers(2, MAX_YZ);
    let mut out = [[[0.0; 3]; 3]; 4];
    for (a, xa) in xi.iter().enumerate() {
        for (b, yb) in eta.iter().enumerate() {
            let xy = xa.times(yb);
            for (c, zc) in zeta.iter().enumerate() {
                out[a][b][c] = xy.times(zc).integrate(area);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slots_are_dense_and_distinct() {
        for deg in 0..=MAX_DEGREE {
            let mut seen = Vec::new();
            for p in 0..=deg {
                for q in 0..=(deg - p) {
                    seen.push(slot(deg, p, q));
                }
            }
            let n = (deg + 1) * (deg + 2) / 2;
            seen.sort_unstable();
            assert_eq!(seen, (0..n).collect::<Vec<_>>(), "deg {deg}");
        }
    }

    #[test]
    fn known_triangle_integrals() {
        // Unit right triangle in the ξη plane: ∫1 = 1/2, ∫ξ = 1/6,
        // ∫ξ² = 1/12, ∫ξη = 1/24, ∫ξ³ = 1/20, ∫ξ³η² = 3!2!/7! = 1/420.
        let t = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let m = monomial_integrals(t, 0.5);
        let close = |a: f64, b: f64| (a - b).abs() < 1e-16;
        assert!(close(m[0][0][0], 0.5));
        assert!(close(m[1][0][0], 1.0 / 6.0));
        assert!(close(m[2][0][0], 1.0 / 12.0));
        assert!(close(m[1][1][0], 1.0 / 24.0));
        assert!(close(m[3][0][0], 1.0 / 20.0));
        assert!(close(m[3][2][0], 1.0 / 420.0));
        assert!(close(m[0][0][1], 0.0));
    }
}
