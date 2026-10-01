//! The sector spline space: n tensor-product B-spline sectors of one
//! degree, sharing one clamped knot vector in both directions.
//!
//! Coefficients are indexed per coordinate as one flat vector:
//! `sector·m² + i·m + j`, with `i` along u and `j` along v: the
//! [`NurbsSurface`] storage (`i_u · num_v + i_v`, v varying fastest). Every
//! surface quantity the fit uses is a linear functional of that vector,
//! represented sparsely as `(index, weight)` pairs.

use super::corner_cells::{
    CELLS, CORNER_CELL_DEGREE, CornerOwner, cell_slot, claim_weight, complement, dual_jet_patch,
    product,
};
use crate::geometry::SurfaceGeometry;
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::nurbs_ops::{basis::basis_funs_and_two_derivs, find_span};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

/// A sparse linear functional over one coordinate's coefficients.
pub(super) type Functional = Vec<(usize, f64)>;

/// Basis functions with their first two derivatives at one parameter.
#[derive(Debug, Clone)]
pub(super) struct Basis1 {
    /// The parameter.
    pub t: f64,
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
    /// Dual-jet corner cells at every sector's rim corner (PK-595, see
    /// `corner_cells`): four v-claims per sector after its `m²` plain
    /// coefficients.
    pub corner_cells: bool,
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
            corner_cells: false,
        }
    }

    /// The same space with dual-jet corner cells: degree 9, one span (a
    /// sector is then a dual-jet patch, see `corner_cells`); other
    /// spaces are returned unchanged.
    pub fn with_corner_cells(mut self) -> Self {
        self.corner_cells = self.spans().len() == 1 && self.degree == CORNER_CELL_DEGREE;
        self
    }

    /// Coefficients per coordinate in one sector.
    pub fn per(&self) -> usize {
        self.m * self.m + if self.corner_cells { CELLS } else { 0 }
    }

    /// Flat index of dual-cell `slot`'s v-claim in `sector`.
    pub fn claim_index(&self, sector: usize, slot: usize) -> usize {
        sector * self.per() + self.m * self.m + slot
    }

    /// The single-span (Bézier) space.
    pub fn bezier(n: usize, degree: usize) -> Self {
        let mut knots = vec![0.0; degree + 1];
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        Self::new(n, degree, knots)
    }

    /// Coefficients per coordinate over all sectors.
    pub fn total(&self) -> usize {
        self.n * self.per()
    }

    /// Flat index of control point `(i, j)` of `sector` (the u-claim at
    /// a dual cell).
    pub fn index(&self, sector: usize, i: usize, j: usize) -> usize {
        sector * self.per() + i * self.m + j
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
            t,
            first: span - self.degree,
            funcs: (0..=self.degree).map(|r| [b0[r], b1[r], b2[r]]).collect(),
        }
    }

    /// The basis products of one sector at `(bu.t, bv.t)`: local index
    /// and `[f, ∂u, ∂v, ∂uu, ∂uv, ∂vv]`. A dual cell contributes its
    /// u-claim at the plain slot and its v-claim after the `m²` plain
    /// slots, each with the claim weight folded in.
    pub fn local_supports(
        &self,
        bu: &Basis1,
        bv: &Basis1,
        owner: CornerOwner,
    ) -> Vec<(usize, [f64; 6])> {
        self.supports_in_order(bu, bv, owner, false)
    }

    /// [`Self::local_supports`] in u-major (`v_major == false`, the
    /// energy's order) or v-major order (the functionals' order): each
    /// caller keeps its historical summation order, so a space without
    /// corner cells evaluates bit for bit as before PK-595.
    fn supports_in_order(
        &self,
        bu: &Basis1,
        bv: &Basis1,
        owner: CornerOwner,
        v_major: bool,
    ) -> Vec<(usize, [f64; 6])> {
        let m = self.m;
        let (nu, nv) = (bu.funcs.len(), bv.funcs.len());
        let mut out = Vec::with_capacity(nu * nv + CELLS);
        let weight = self.corner_cells.then(|| claim_weight(bu.t, bv.t, owner));
        let (outer_len, inner_len) = if v_major { (nv, nu) } else { (nu, nv) };
        for outer in 0..outer_len {
            for inner in 0..inner_len {
                let (ri, rj) = if v_major {
                    (inner, outer)
                } else {
                    (outer, inner)
                };
                let (fu, fv) = (bu.funcs[ri], bv.funcs[rj]);
                let (i, j) = (bu.first + ri, bv.first + rj);
                let b = [
                    fu[0] * fv[0],
                    fu[1] * fv[0],
                    fu[0] * fv[1],
                    fu[2] * fv[0],
                    fu[1] * fv[1],
                    fu[0] * fv[2],
                ];
                match (weight, cell_slot(m, i, j)) {
                    (Some(a), Some(slot)) => {
                        out.push((i * m + j, product(b, a)));
                        out.push((m * m + slot, product(b, complement(a))));
                    }
                    _ => out.push((i * m + j, b)),
                }
            }
        }
        out
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
        self.functional_owned(sector, bu, bv, du, dv, CornerOwner::Interior)
    }

    /// [`Self::functional_from`] for a sample owned by a rim (decides
    /// the dual cells' claim at the corner itself).
    pub fn functional_owned(
        &self,
        sector: usize,
        bu: &Basis1,
        bv: &Basis1,
        du: usize,
        dv: usize,
        owner: CornerOwner,
    ) -> Functional {
        let component = match (du, dv) {
            (0, 0) => 0,
            (1, 0) => 1,
            (0, 1) => 2,
            (2, 0) => 3,
            (1, 1) => 4,
            _ => 5,
        };
        let base = sector * self.per();
        self.supports_in_order(bu, bv, owner, true)
            .into_iter()
            .filter_map(|(l, d)| (d[component] != 0.0).then_some((base + l, d[component])))
            .collect()
    }

    /// The functional `∂^{du+dv} S_sector / ∂u^du ∂v^dv` at `(u, v)`.
    #[cfg(test)]
    pub fn functional(&self, sector: usize, u: f64, v: f64, du: usize, dv: usize) -> Functional {
        self.functional_at(sector, u, v, du, dv, CornerOwner::Interior)
    }

    /// [`Self::functional`] for a sample a rim owns.
    pub fn functional_at(
        &self,
        sector: usize,
        u: f64,
        v: f64,
        du: usize,
        dv: usize,
        owner: CornerOwner,
    ) -> Functional {
        let bu = self.basis(u);
        let bv = self.basis(v);
        self.functional_owned(sector, &bu, &bv, du, dv, owner)
    }

    /// The face surface of `sector`: the dual-jet patch with corner
    /// cells, else the sector itself.
    pub fn face_surface(&self, coords: &[Vec<f64>; 3], sector: usize) -> SurfaceGeometry {
        if !self.corner_cells {
            return SurfaceGeometry::Nurbs(self.sector_surface(coords, sector));
        }
        let at = |k: usize| Point3::new(coords[0][k], coords[1][k], coords[2][k]);
        SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(Box::new(
            dual_jet_patch(&|i, j| at(self.index(sector, i, j)), &|slot| {
                at(self.claim_index(sector, slot))
            }),
        )))
    }

    /// Wrap `sector`'s plain coefficients (the u-claims at dual cells) as
    /// a surface. With corner cells it carries the sector's boundary
    /// curves exactly (the claims are interior cells), not its interior.
    pub fn sector_surface(&self, coords: &[Vec<f64>; 3], sector: usize) -> NurbsSurface {
        let base = sector * self.per();
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
