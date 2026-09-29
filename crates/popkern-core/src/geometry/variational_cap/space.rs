//! The sector spline space: n tensor-product B-spline sectors of one
//! degree, sharing one clamped knot vector in both directions.
//!
//! Coefficients are indexed per coordinate as one flat vector:
//! `sector·m² + i·m + j`, with `i` along u and `j` along v: the
//! [`NurbsSurface`] storage (`i_u · num_v + i_v`, v varying fastest). Every
//! surface quantity the fit uses is a linear functional of that vector,
//! represented sparsely as `(index, weight)` pairs.

use crate::geometry::nurbs_ops::{basis::basis_funs_and_two_derivs, find_span};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

/// A sparse linear functional over one coordinate's coefficients.
pub(super) type Functional = Vec<(usize, f64)>;

/// Basis functions with their first two derivatives at one parameter.
#[derive(Debug, Clone)]
pub(super) struct Basis1 {
    /// Global index of the first non-zero function.
    pub first: usize,
    /// `[value, d/dt, d²/dt²]` per non-zero function.
    pub funcs: Vec<[f64; 3]>,
}

/// The sector space for one valence, degree and knot vector.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SectorSpace {
    /// Valence (number of sectors).
    pub n: usize,
    /// Degree in both directions.
    pub degree: usize,
    /// Clamped knot vector on [0, 1], shared by u and v.
    pub knots: Vec<f64>,
    /// Control points per direction.
    pub m: usize,
}

impl SectorSpace {
    /// The space with `knots` (clamped, degree + 1 end multiplicity).
    pub fn new(n: usize, degree: usize, knots: Vec<f64>) -> Self {
        let m = knots.len() - degree - 1;
        SectorSpace {
            n,
            degree,
            knots,
            m,
        }
    }

    /// The single-span (Bézier) space.
    pub fn bezier(n: usize, degree: usize) -> Self {
        let mut knots = vec![0.0; degree + 1];
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        Self::new(n, degree, knots)
    }

    /// Coefficients per coordinate over all sectors.
    pub fn total(&self) -> usize {
        self.n * self.m * self.m
    }

    /// Flat index of control point `(i, j)` of `sector`.
    pub fn index(&self, sector: usize, i: usize, j: usize) -> usize {
        sector * self.m * self.m + i * self.m + j
    }

    /// The distinct knot spans `[a, b)` of `K`.
    pub fn spans(&self) -> Vec<(f64, f64)> {
        self.knots
            .windows(2)
            .filter(|w| w[1] > w[0])
            .map(|w| (w[0], w[1]))
            .collect()
    }

    /// Basis values and derivatives at `t` in [0, 1].
    pub fn basis(&self, t: f64) -> Basis1 {
        let span = find_span(self.degree, &self.knots, self.m, t);
        let (b0, b1, b2) = basis_funs_and_two_derivs(self.degree, &self.knots, span, t);
        Basis1 {
            first: span - self.degree,
            funcs: (0..=self.degree).map(|r| [b0[r], b1[r], b2[r]]).collect(),
        }
    }

    /// The functional `∂^{du+dv} S_sector / ∂u^du ∂v^dv` at `(u, v)`
    /// (`du, dv <= 2`), from precomputed bases.
    pub fn functional_from(
        &self,
        sector: usize,
        bu: &Basis1,
        bv: &Basis1,
        du: usize,
        dv: usize,
    ) -> Functional {
        let mut out = Vec::with_capacity(bu.funcs.len() * bv.funcs.len());
        for (rj, fv) in bv.funcs.iter().enumerate() {
            for (ri, fu) in bu.funcs.iter().enumerate() {
                let w = fu[du] * fv[dv];
                if w != 0.0 {
                    out.push((self.index(sector, bu.first + ri, bv.first + rj), w));
                }
            }
        }
        out
    }

    /// The functional `∂^{du+dv} S_sector / ∂u^du ∂v^dv` at `(u, v)`.
    pub fn functional(&self, sector: usize, u: f64, v: f64, du: usize, dv: usize) -> Functional {
        let bu = self.basis(u);
        let bv = self.basis(v);
        self.functional_from(sector, &bu, &bv, du, dv)
    }

    /// Wrap `sector`'s coefficients as a surface.
    pub fn sector_surface(&self, coords: &[Vec<f64>; 3], sector: usize) -> NurbsSurface {
        let base = sector * self.m * self.m;
        let control_points = (0..self.m * self.m)
            .map(|k| {
                Point3::new(
                    coords[0][base + k],
                    coords[1][base + k],
                    coords[2][base + k],
                )
            })
            .collect();
        NurbsSurface {
            degree_u: self.degree,
            degree_v: self.degree,
            knots_u: self.knots.clone(),
            knots_v: self.knots.clone(),
            control_points,
            num_u: self.m,
            num_v: self.m,
            weights: vec![1.0; self.m * self.m],
        }
    }
}

/// Apply a functional to one coordinate's coefficients.
pub(super) fn apply(f: &Functional, x: &[f64]) -> f64 {
    f.iter().map(|&(i, w)| w * x[i]).sum()
}

/// Apply a functional to all three coordinates.
pub(super) fn apply3(f: &Functional, coords: &[Vec<f64>; 3]) -> Vec3 {
    Vec3::new(
        apply(f, &coords[0]),
        apply(f, &coords[1]),
        apply(f, &coords[2]),
    )
}

/// Gauss–Legendre nodes and weights on [0, 1] (`count` <= 12).
pub(crate) fn gauss_legendre(count: usize) -> Vec<(f64, f64)> {
    // Newton on the Legendre polynomial from Chebyshev-like seeds.
    let n = count.max(1);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let pn = if n == 1 { x } else { p1 };
            let pm = if n == 1 { 1.0 } else { p0 };
            dp = n as f64 * (x * pn - pm) / (x * x - 1.0);
            let dx = pn / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        let w = 2.0 / ((1.0 - x * x) * dp * dp);
        out.push((0.5 * (1.0 - x), 0.5 * w));
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceEval;

    #[test]
    fn gauss_legendre_integrates_polynomials_exactly() {
        for n in 1..=12 {
            let nodes = gauss_legendre(n);
            for p in 0..(2 * n) {
                let integral: f64 = nodes.iter().map(|&(x, w)| w * x.powi(p as i32)).sum();
                assert!(
                    (integral - 1.0 / (p as f64 + 1.0)).abs() < 1e-13,
                    "n {n} p {p}"
                );
            }
        }
    }

    #[test]
    fn functionals_match_the_surface_evaluator() {
        let mut space = SectorSpace::bezier(3, 5);
        space = SectorSpace::new(space.n, space.degree, {
            let mut k = space.knots.clone();
            k.insert(6, 0.4);
            k
        });
        let total = space.total();
        let coords: [Vec<f64>; 3] = [
            (0..total).map(|k| (k as f64 * 0.37).sin()).collect(),
            (0..total).map(|k| (k as f64 * 0.11).cos()).collect(),
            (0..total).map(|k| (k as f64 * 0.23).sin() * 0.5).collect(),
        ];
        let surf = space.sector_surface(&coords, 1);
        for &(u, v) in &[(0.1, 0.7), (0.4, 0.4), (0.9, 0.05), (1.0, 1.0)] {
            let p = apply3(&space.functional(1, u, v, 0, 0), &coords);
            let q = surf.point_at(u, v);
            let d = (p - Vec3::new(q.x, q.y, q.z)).length();
            assert!(d < 1e-12, "({u}, {v}): {d:.3e} {p:?} {q:?}");
            let su = apply3(&space.functional(1, u, v, 1, 0), &coords);
            assert!((su - surf.partial_u(u, v)).length() < 1e-9);
            let sv = apply3(&space.functional(1, u, v, 0, 1), &coords);
            assert!((sv - surf.partial_v(u, v)).length() < 1e-9);
        }
    }
}
