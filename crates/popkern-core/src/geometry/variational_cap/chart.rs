//! The characteristic map φ_n (Loop & Schaefer 2008, §3.1).
//!
//! φ_n maps the unit square bilinearly onto the kite with corners
//! (0,0), (1,0), (1, tan(π/n)) and (cos(2π/n), sin(2π/n)): the sector
//! of a regular n-valent vertex. In Bernstein form,
//! `φ_x = b¹(u)ᵀ [[0, c],[1, 1]] b¹(v)`, `φ_y = b¹(u)ᵀ [[0, s],[0, τ]]
//! b¹(v)` with c = cos(2π/n), s = sin(2π/n), τ = tan(π/n), i.e.
//! `φ_x = u + c·v·(1 − u)` and `φ_y = v·(s·(1 − u) + τ·u)`.
//!
//! Its correspondence maps `φ_n⁻¹ ∘ r_n⁻¹ ∘ φ_n` between neighbouring
//! sectors satisfy the cocycle condition at the centre, which is what
//! makes the spoke conditions of [`super::spoke`] consistent there.

/// φ_n for one valence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct CharacteristicMap {
    /// cos(2π/n): the coefficient of the spoke conditions.
    pub c: f64,
    s: f64,
    tau: f64,
}

impl CharacteristicMap {
    /// φ_n for valence `n` (n ≥ 3).
    pub fn new(n: usize) -> Self {
        let angle = 2.0 * std::f64::consts::PI / n as f64;
        CharacteristicMap {
            c: angle.cos(),
            s: angle.sin(),
            tau: (std::f64::consts::PI / n as f64).tan(),
        }
    }

    /// Jacobian columns `(φ_u, φ_v)` at `(u, v)`.
    pub fn jacobian(&self, u: f64, v: f64) -> [[f64; 2]; 2] {
        let xu = 1.0 - self.c * v;
        let yu = v * (self.tau - self.s);
        let xv = self.c * (1.0 - u);
        let yv = self.s * (1.0 - u) + self.tau * u;
        [[xu, yu], [xv, yv]]
    }

    /// The metric `JᵀJ` at `(u, v)` as `[g_uu, g_uv, g_vv]`.
    pub fn metric(&self, u: f64, v: f64) -> [f64; 3] {
        let [a, b] = self.jacobian(u, v);
        [
            a[0] * a[0] + a[1] * a[1],
            a[0] * b[0] + a[1] * b[1],
            b[0] * b[0] + b[1] * b[1],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn corners_land_on_the_kite() {
        for n in 3..=24 {
            let phi = CharacteristicMap::new(n);
            // φ(u, v) from the closed form.
            let eval = |u: f64, v: f64| {
                (
                    u + phi.c * v * (1.0 - u),
                    v * (phi.s * (1.0 - u) + phi.tau * u),
                )
            };
            let (x, y) = eval(0.0, 1.0);
            assert!((x - phi.c).abs() < 1e-15 && (y - phi.s).abs() < 1e-15);
            let (x, y) = eval(1.0, 1.0);
            assert!((x - 1.0).abs() < 1e-15 && (y - phi.tau).abs() < 1e-15);
            // The metric at the centre measures the sector angle 2π/n.
            let [guu, guv, gvv] = phi.metric(0.0, 0.0);
            let cos_angle = guv / (guu * gvv).sqrt();
            let expected = (2.0 * std::f64::consts::PI / n as f64).cos();
            assert!((cos_angle - expected).abs() < 1e-14, "n = {n}");
        }
    }
}
