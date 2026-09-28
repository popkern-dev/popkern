//! Gauss–Newton rows for the fitted rim terms, over the null-space
//! coordinates `y = [y_x; y_y; y_z]` (each of length `k`).
//!
//! With `w` the cap's cross-boundary derivative at a rim sample (`S_u`
//! on a trailing half, `S_v` on a leading half), `w'` the second one,
//! `N` and `II` the oriented neighbour normal and form, and `w_k` the
//! previous iterate's `w`:
//!
//! ```text
//! G1:  N·w / |w_k|                                   = 0
//! G2:  (N·w' − 2 II(w_k, w) + II(w_k, w_k)) / |w_k|² = 0
//! ```
//!
//! G2 is `N·w' − II(w, w)` linearised about `w_k`: with G0 exact and G1
//! met, the along-rim and mixed curvature terms agree automatically, so
//! the transverse term is the only independent G2 condition. Each row is
//! divided by its target, so the least-squares term reads Σ (r/target)².

use faer::Mat;

use super::measure::{RimSample, cap_jet, rim_uv};
use super::nullspace::HardSystem;
use super::space::{Functional, SectorSpace};
use crate::geometry::Vec3;

/// A functional reduced to the null-space coordinates of one
/// coordinate: `f·x = a·y + c`.
pub(super) struct Reduced {
    pub a: Vec<f64>,
    pub c: f64,
}

/// Reduce `f` for one coordinate with particular free values `xp` and
/// fixed values `fixed`.
pub(super) fn reduce(f: &Functional, hs: &HardSystem, xp: &[f64], fixed: &[f64]) -> Reduced {
    let k = hs.dim();
    let mut a = vec![0.0; k];
    let mut c = 0.0;
    for &(i, w) in f {
        match hs.slots[i] {
            (true, r) => c += w * fixed[r],
            (false, r) => {
                c += w * xp[r];
                for (j, slot) in a.iter_mut().enumerate() {
                    *slot += w * hs.z.read(r, j);
                }
            }
        }
    }
    Reduced { a, c }
}

/// Everything the rows need about the hard system for all coordinates.
pub(super) struct Reduction<'a> {
    pub hs: &'a HardSystem,
    pub xp: &'a [Vec<f64>; 3],
    pub fixed: &'a [Vec<f64>; 3],
}

/// Rows `(coefficients over 3k, rhs)` for every sample.
pub(super) fn rim_rows(
    space: &SectorSpace,
    red: &Reduction<'_>,
    coords: &[Vec<f64>; 3],
    samples: &[RimSample],
    g1_target: f64,
    g2_target: Option<f64>,
) -> (Mat<f64>, Vec<f64>) {
    let k = red.hs.dim();
    let per = if g2_target.is_some() { 2 } else { 1 };
    let mut rows = Mat::<f64>::zeros(samples.len() * per, 3 * k);
    let mut rhs = vec![0.0; samples.len() * per];
    for (si, sample) in samples.iter().enumerate() {
        let ((u, v), cross_u) = rim_uv(sample);
        let (d1, d2) = if cross_u {
            ((1, 0), (2, 0))
        } else {
            ((0, 1), (0, 2))
        };
        let fw = space.functional(sample.sector, u, v, d1.0, d1.1);
        let fw2 = space.functional(sample.sector, u, v, d2.0, d2.1);
        let jet = cap_jet(space, coords, sample.sector, u, v);
        let wk = if cross_u { jet.su } else { jet.sv };
        let len = wk.length().max(1e-12);
        let n = sample.jet.normal;
        let comps = |v: Vec3| [v.x, v.y, v.z];
        // G1.
        let r = si * per;
        let scale1 = 1.0 / (len * g1_target);
        let mut constant = 0.0;
        for (c, nc) in comps(n).iter().enumerate() {
            let red_w = reduce(&fw, red.hs, &red.xp[c], &red.fixed[c]);
            for j in 0..k {
                rows.write(r, c * k + j, scale1 * nc * red_w.a[j]);
            }
            constant += nc * red_w.c;
        }
        rhs[r] = -scale1 * constant;
        // G2.
        if let Some(t2) = g2_target {
            let r = r + 1;
            let qwk = sample.jet.apply(wk);
            let scale2 = 1.0 / (len * len * t2);
            let mut constant = sample.jet.ii(wk, wk);
            for (c, (nc, qc)) in comps(n).iter().zip(comps(qwk)).enumerate() {
                let red_w = reduce(&fw, red.hs, &red.xp[c], &red.fixed[c]);
                let red_w2 = reduce(&fw2, red.hs, &red.xp[c], &red.fixed[c]);
                for j in 0..k {
                    let coef = nc * red_w2.a[j] - 2.0 * qc * red_w.a[j];
                    rows.write(r, c * k + j, scale2 * coef);
                }
                constant += nc * red_w2.c - 2.0 * qc * red_w.c;
            }
            rhs[r] = -scale2 * constant;
        }
    }
    (rows, rhs)
}
