//! Gauss–Newton rows for the fitted rim terms, over the null-space
//! coordinates `y = [y_x; y_y; y_z]` (each of length `k`).
//!
//! With `w` the cap's cross-boundary derivative at a rim sample (`S_u`
//! on a trailing half, `S_v` on a leading half), `w'` the second one,
//! `N` and `II` the oriented neighbour normal and form, and `w_k` the
//! previous iterate's `w`:
//!
//! ```text
//! G1:        N·w / |w_k|                                   = 0
//! G2 cross:  (N·w' − 2 II(w_k, w) + II(w_k, w_k)) / |w_k|² = 0
//! G2 mixed:  (N·∂w/∂t − II(T, w)) / (|w_k| |T|)            = 0
//! ```
//!
//! with `T` the (fixed) rim tangent. The cross row is `N·w' − II(w, w)`
//! linearised about `w_k`. The mixed row is exact and linear. With G0
//! and G1 exact along the whole rim the mixed term would follow from
//! differentiating G1 along the rim, but a fitted G1 barely constrains
//! the twist: measured on the pentagon flagship, the corner twist
//! missed the strip by 0.53 (normalised mixed curvature) while G1 read
//! 5e-4. The along-rim term is fixed by the rim (Meusnier). Each row is
//! divided by its target, so the least-squares term reads Σ (r/target)².

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::measure::{RimSample, cap_jet_owned, rim_owner, rim_uv};
use super::nullspace::HardSystem;
use super::space::{Functional, SectorSpace};
use crate::geometry::Vec3;

/// Extra weight on the mixed row at a setback corner (the rim end
/// `x = 1`, where two sides meet). Each side's G2 data fixes the corner
/// twist the same way (the two neighbours' second forms agree there),
/// but the neighbours' curvature slopes conflict at third order for any
/// polynomial cap (see the design doc, "setback corner theorem"). At
/// unit weight the fit spent the corner twist to trade against that
/// conflict and kept a constant corner error under every refinement
/// (measured on the pentagon flagship: mixed curvature −0.078 against
/// 0 at every level). Pinning the twist moves the conflict into the
/// cross term beside the corner, where it decays linearly with the
/// distance to the corner and so shrinks under corner refinement.
const CORNER_TWIST_WEIGHT: f64 = 1.0e3;

/// A functional reduced to the null-space coordinates of one
/// coordinate: `f·x = a·y + c`.
pub(super) struct Reduced {
    /// Null-space coefficients (the same for every coordinate).
    pub a: Vec<f64>,
    /// Constant term per coordinate.
    pub c: [f64; 3],
}

/// Reduce `f` onto the null space; the constant uses each coordinate's
/// particular free values and fixed values.
pub(super) fn reduce(f: &Functional, red: &Reduction<'_>) -> Reduced {
    let hs = red.hs;
    let k = hs.dim();
    let mut a = vec![0.0; k];
    let mut c = [0.0; 3];
    for &(i, w) in f {
        match hs.slots[i] {
            (true, r) => {
                for (cc, slot) in c.iter_mut().enumerate() {
                    *slot += w * red.fixed[cc][r];
                }
            }
            (false, r) => {
                for (cc, slot) in c.iter_mut().enumerate() {
                    *slot += w * red.xp[cc][r];
                }
                hs.add_row(r, w, &mut a);
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

/// How the transversality prior sizes its target across the rim.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CrossDerivativeTarget {
    /// The length of the start's cross derivative within the
    /// neighbour's tangent plane (at least a quarter of its full
    /// length): the start's own speed across the rim, where the start
    /// meets the neighbour near tangentially.
    #[default]
    InPlaneLength,
    /// The start's full length, rotated into the tangent plane. Where the
    /// start crosses the rim steeply to the neighbour's tangent plane its
    /// in-plane part collapses, and a short target pins a short cross
    /// derivative: the cap turns into the neighbour inside a thin rim
    /// layer (PK-625, measured on the octagon h3 wide apex: 44 / r_min
    /// at the rim against 0.3 inside, the corner-patched cap reading 37;
    /// with the start's length, 2.33).
    StartLength,
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
    let per = if g2_target.is_some() { 3 } else { 1 };
    let mut rows = Mat::<f64>::zeros(samples.len() * per, 3 * k);
    let mut rhs = vec![0.0; samples.len() * per];
    for (si, sample) in samples.iter().enumerate() {
        let ((u, v), cross_u) = rim_uv(sample);
        let (d1, d2) = if cross_u {
            ((1, 0), (2, 0))
        } else {
            ((0, 1), (0, 2))
        };
        let owner = rim_owner(sample);
        let fw = space.functional_at(sample.sector, u, v, d1.0, d1.1, owner);
        let fw2 = space.functional_at(sample.sector, u, v, d2.0, d2.1, owner);
        let jet = cap_jet_owned(space, coords, sample.sector, u, v, owner);
        let wk = if cross_u { jet.su } else { jet.sv };
        let len = wk.length().max(1e-12);
        let n = sample.jet.normal;
        let comps = |v: Vec3| [v.x, v.y, v.z];
        // G1.
        let r = si * per;
        let scale1 = 1.0 / (len * g1_target);
        let mut constant = 0.0;
        let red_w = reduce(&fw, red);
        for (c, nc) in comps(n).iter().enumerate() {
            for j in 0..k {
                rows.write(r, c * k + j, scale1 * nc * red_w.a[j]);
            }
            constant += nc * red_w.c[c];
        }
        rhs[r] = -scale1 * constant;
        // G2.
        if let Some(t2) = g2_target {
            let r = r + 1;
            let qwk = sample.jet.apply(wk);
            let scale2 = 1.0 / (len * len * t2);
            let mut constant = sample.jet.ii(wk, wk);
            let red_w2 = reduce(&fw2, red);
            for (c, (nc, qc)) in comps(n).iter().zip(comps(qwk)).enumerate() {
                for j in 0..k {
                    let coef = nc * red_w2.a[j] - 2.0 * qc * red_w.a[j];
                    rows.write(r, c * k + j, scale2 * coef);
                }
                constant += nc * red_w2.c[c] - 2.0 * qc * red_w.c[c];
            }
            rhs[r] = -scale2 * constant;
            // Mixed: N·S_uv − II(T, w), T the rim tangent.
            let r = r + 1;
            let tangent = if cross_u { jet.sv } else { jet.su };
            let fuv = space.functional_at(sample.sector, u, v, 1, 1, owner);
            let qt = sample.jet.apply(tangent);
            // With corner cells each rim owns its corner twist, so the
            // pin that moved the theorem's conflict off the corner is
            // not needed (PK-595).
            let corner = if sample.x == 1.0 && !space.corner_cells {
                CORNER_TWIST_WEIGHT
            } else {
                1.0
            };
            let scale3 = corner / (len * tangent.length().max(1e-12) * t2);
            let mut constant = 0.0;
            let red_uv = reduce(&fuv, red);
            for (c, (nc, qc)) in comps(n).iter().zip(comps(qt)).enumerate() {
                for j in 0..k {
                    rows.write(r, c * k + j, scale3 * (nc * red_uv.a[j] - qc * red_w.a[j]));
                }
                constant += nc * red_uv.c[c] - qc * red_w.c[c];
            }
            rhs[r] = -scale3 * constant;
        }
    }
    (rows, rhs)
}

/// Transversality prior rows: each rim sample's cross derivative `w`
/// pulled towards `w_t`, per coordinate, divided by `|w_t|`.
///
/// The G1 rows are homogeneous in `w`, so nothing else fixes how long
/// the cross derivative is: measured on the cube corner 0.05/0.2/0.3,
/// the G1 step shrank it at a rim midpoint until the Jacobian nearly
/// vanished (regularity 0.35 -> 1.9e-2 -> 1.5e-3). The targets are the
/// start surface's cross derivatives projected into each neighbour's
/// tangent plane, so the prior never competes with G1.
pub(super) fn transversal_rows(
    space: &SectorSpace,
    red: &Reduction<'_>,
    samples: &[RimSample],
    targets: &[Vec3],
    weight: f64,
) -> (Mat<f64>, Vec<f64>) {
    let k = red.hs.dim();
    let mut rows = Mat::<f64>::zeros(samples.len() * 3, 3 * k);
    let mut rhs = vec![0.0; samples.len() * 3];
    for (si, (sample, wt)) in samples.iter().zip(targets).enumerate() {
        let ((u, v), cross_u) = rim_uv(sample);
        let (du, dv) = if cross_u { (1, 0) } else { (0, 1) };
        let fw = space.functional_at(sample.sector, u, v, du, dv, rim_owner(sample));
        let red_w = reduce(&fw, red);
        let scale = weight / wt.length().max(1e-12);
        for (c, t) in [wt.x, wt.y, wt.z].into_iter().enumerate() {
            let r = 3 * si + c;
            for j in 0..k {
                rows.write(r, c * k + j, scale * red_w.a[j]);
            }
            rhs[r] = scale * (t - red_w.c[c]);
        }
    }
    (rows, rhs)
}

/// The prior's targets: across each rim, into the neighbour, in its
/// tangent plane.
///
/// With the sides counter-clockwise about the outward normal N, the cap
/// lies to the left of the loop tangent T, so the cross derivative (u or
/// v increasing towards the rim) must point along `e_out = T × N`. The
/// linear G1 row `N·w = 0` cannot tell `w` from `−w`: on the cube corner
/// 0.05/0.2/0.3 the G1 step produced a cap folded back over its strip
/// along a third of one rim (cap normal against the strip's, dot −0.73).
/// The target's component along the rim is the start's (the parameter
/// flow); its length across the rim follows `rule`.
pub(super) fn transversal_targets(
    space: &SectorSpace,
    coords: &[Vec<f64>; 3],
    samples: &[RimSample],
    rule: CrossDerivativeTarget,
) -> Vec<Vec3> {
    samples
        .iter()
        .map(|sample| {
            let ((u, v), cross_u) = rim_uv(sample);
            let j = cap_jet_owned(space, coords, sample.sector, u, v, rim_owner(sample));
            // Trailing half (u = 1): the rim runs with the loop along +v.
            // Leading half (v = 1): the rim runs against it along +u.
            let (w, t_loop) = if cross_u {
                (j.su, j.sv)
            } else {
                (j.sv, j.su * -1.0)
            };
            let n = sample.jet.normal;
            let Some(t) = t_loop.normalized() else {
                return w - n * n.dot(&w);
            };
            let e_out = t.cross(&n);
            let along = t * t.dot(&w);
            let across = match rule {
                CrossDerivativeTarget::InPlaneLength => {
                    (w - n * n.dot(&w) - along).length().max(0.25 * w.length())
                }
                CrossDerivativeTarget::StartLength => (w.dot(&w) - along.dot(&along))
                    .max(0.0)
                    .sqrt()
                    .max(0.25 * w.length()),
            };
            e_out * across + along
        })
        .collect()
}
