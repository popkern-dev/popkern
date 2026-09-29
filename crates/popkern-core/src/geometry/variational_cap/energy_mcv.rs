//! The mean-curvature-variation integrand's pieces: the frozen
//! reference fields, `H = N·Δ_g x` as a functional, and its gradient by
//! central differences (see `energy.rs` for the energy itself).

use super::energy::supports;
use super::space::SectorSpace;
use crate::geometry::Vec3;

/// The frozen reference fields at one parameter point: unit normal,
/// metric `[g_uu, g_uv, g_vv]`, and the Laplace–Beltrami coefficients
/// `Δf = a_uu f_uu + 2 a_uv f_uv + a_vv f_vv + b_u f_u + b_v f_v`
/// (`a = g⁻¹`, `b_k = −g^{ij} Γ^k_ij`).
#[derive(Debug, Clone, Copy)]
pub(super) struct Reference {
    pub normal: Vec3,
    pub metric: [f64; 3],
    pub lap2: [f64; 3],
    pub lap1: [f64; 2],
}

impl Reference {
    /// Reference fields from a surface's first and second partials.
    pub fn from_partials(su: Vec3, sv: Vec3, suu: Vec3, suv: Vec3, svv: Vec3) -> Option<Self> {
        let normal = su.cross(&sv).normalized()?;
        let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let det = e * g - f * f;
        if det.is_nan() || det <= 1e-24 * e * g {
            return None;
        }
        let inv = [g / det, -f / det, e / det];
        // Γ^k_ij = g^{kl} <x_ij, x_l>.
        let gamma = |xij: Vec3| {
            let (p, q) = (xij.dot(&su), xij.dot(&sv));
            [inv[0] * p + inv[1] * q, inv[1] * p + inv[2] * q]
        };
        let (guu, guv, gvv) = (gamma(suu), gamma(suv), gamma(svv));
        let lap1 = [
            -(inv[0] * guu[0] + 2.0 * inv[1] * guv[0] + inv[2] * gvv[0]),
            -(inv[0] * guu[1] + 2.0 * inv[1] * guv[1] + inv[2] * gvv[1]),
        ];
        Some(Reference {
            normal,
            metric: [e, f, g],
            lap2: inv,
            lap1,
        })
    }
}

/// A sparse functional over one sector's `(coord, local)` indices.
pub(super) type Sparse = Vec<(usize, f64)>;

/// `H = N·Δ_g x` at `(u, v)` as a functional, with the reference fields
/// taken at the same point.
fn mean_curvature(
    space: &SectorSpace,
    sector: usize,
    u: f64,
    v: f64,
    reference: &dyn Fn(usize, f64, f64) -> Option<Reference>,
) -> Option<Sparse> {
    let r = reference(sector, u, v)?;
    let m = space.m;
    let per = m * m;
    let bu = space.basis(u);
    let bv = space.basis(v);
    let nc = [r.normal.x, r.normal.y, r.normal.z];
    let mut out = Vec::new();
    for (l, d) in supports(m, &bu, &bv) {
        let lap = r.lap2[0] * d[3]
            + 2.0 * r.lap2[1] * d[4]
            + r.lap2[2] * d[5]
            + r.lap1[0] * d[1]
            + r.lap1[1] * d[2];
        for (c, n) in nc.iter().enumerate() {
            out.push((c * per + l, n * lap));
        }
    }
    Some(out)
}

/// `(∂_u H, ∂_v H)` by central differences of the (linear) functional
/// `H`, each side with its own frozen reference. The step is 1e-3 of
/// the span, inside the span for every Gauss node used.
pub(super) fn mean_curvature_gradient(
    space: &SectorSpace,
    sector: usize,
    u: f64,
    v: f64,
    h: f64,
    reference: &dyn Fn(usize, f64, f64) -> Option<Reference>,
) -> Option<[Sparse; 2]> {
    let at = |a: f64, b: f64| mean_curvature(space, sector, a, b, reference);
    let (up, um) = (at(u + h, v)?, at(u - h, v)?);
    let (vp, vm) = (at(u, v + h)?, at(u, v - h)?);
    let s = 0.5 / h;
    Some([
        combine(&[(s, &up), (-s, &um)]),
        combine(&[(s, &vp), (-s, &vm)]),
    ])
}

/// `Σ c_i f_i`, merged by index.
pub(super) fn combine(terms: &[(f64, &Sparse)]) -> Sparse {
    let mut all: Sparse = terms
        .iter()
        .flat_map(|(c, f)| f.iter().map(move |&(i, w)| (i, c * w)))
        .collect();
    all.sort_by_key(|&(i, _)| i);
    let mut merged: Sparse = Vec::with_capacity(all.len());
    for (i, w) in all {
        match merged.last_mut() {
            Some(last) if last.0 == i => last.1 += w,
            _ => merged.push((i, w)),
        }
    }
    merged
}

/// `block += w · row rowᵀ`.
pub(super) fn outer(block: &mut [f64], size: usize, row: &Sparse, w: f64) {
    for &(i, a) in row {
        for &(j, b) in row {
            block[i * size + j] += w * a * b;
        }
    }
}
