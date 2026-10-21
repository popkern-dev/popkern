//! Bent correspondences for the triangle chordal bound (PK-908).
//!
//! The affine correspondence pairs `S(y)` with `L(y)` at the same
//! parameter, so a surface that curves inside its own tangent plane (a
//! torus's u-circle at the top of its tube) is charged the chord's
//! in-plane sag although it passes close to the chord. A bent
//! correspondence pairs `S(φ(y))` with `L(y)`, where `φ` is the
//! quadratic Bézier map of the parameter triangle with its corners
//! fixed and each edge `k` (vertices `k, k + 1`) bowed by `δ_k` at its
//! midpoint:
//!
//! `φ(y) = y + Σ_k 4 λ_k λ_{k+1} δ_k`
//!
//! (`λ` the barycentric coordinates). A bent edge is given by its
//! control point `c_k = m_k + 2 δ_k` (`m_k` the edge midpoint), a float,
//! so the map is exact for the caller's predicates and `δ_k` is
//! enclosed from it; a straight edge has `δ_k = 0` exactly (a rounded
//! midpoint control would bow it by an ulp). On edge `k` only its own
//! bubble survives, so two triangles bending a shared edge with the
//! same control agree there, and a straight edge is the identity. Waldron's bound applies to `S∘φ − L` with
//!
//! `∂_i∂_j (S∘φ) = Σ_kl S_kl ∂_iφ_k ∂_jφ_l + Σ_k S_k ∂_i∂_jφ_k`,
//!
//! from second-order jets of `S` over the enclosure of `φ`. Any `δ`
//! gives a valid bound for this triangle's pairing; whether the pairings
//! of a whole mesh cover the face is the caller's (see the mesh crate's
//! `chordal_bend`).

use super::super::interval::Iv;
use super::super::interval_vec::IvVec3;
use super::super::surface_jet::SurfaceJet;
use super::chord_triangle::Uv;

/// The bent map of one parameter triangle, enclosed.
#[derive(Debug, Clone)]
pub(super) struct Bend {
    w: [Uv; 3],
    /// `∂λ_i / ∂(u, v)`.
    grad: [[Iv; 2]; 3],
    /// Midpoint displacement `δ_k` of edge `k`.
    d: [[Iv; 2]; 3],
    det: Iv,
}

impl Bend {
    /// The map whose edge `k` has the control point `controls[k]`
    /// (`None`: straight); `None` for a degenerate triangle (no
    /// barycentric coordinates).
    pub(super) fn new(w: [Uv; 3], controls: [Option<Uv>; 3]) -> Option<Bend> {
        let p = |x: f64| Iv::point(x);
        let (e1, e2) = (
            [p(w[1].0) - p(w[0].0), p(w[1].1) - p(w[0].1)],
            [p(w[2].0) - p(w[0].0), p(w[2].1) - p(w[0].1)],
        );
        let det = e1[0] * e2[1] - e2[0] * e1[1];
        if det.contains_zero() || !det.is_finite() {
            return None;
        }
        // λ1 = (q − w0) × e2 / det, λ2 = e1 × (q − w0) / det.
        let g1 = [e2[1] / det, -e2[0] / det];
        let g2 = [-e1[1] / det, e1[0] / det];
        let g0 = [-(g1[0] + g2[0]), -(g1[1] + g2[1])];
        // δ_k = c_k / 2 − (w_k + w_{k+1}) / 4.
        let mut d = [[Iv::point(0.0); 2]; 3];
        for (k, dk) in d.iter_mut().enumerate() {
            let Some(c) = controls[k] else { continue };
            let (a, b) = (w[k], w[(k + 1) % 3]);
            *dk = [
                p(c.0).scale(0.5) - (p(a.0) + p(b.0)).scale(0.25),
                p(c.1).scale(0.5) - (p(a.1) + p(b.1)).scale(0.25),
            ];
        }
        Some(Bend {
            w,
            grad: [g0, g1, g2],
            d,
            det,
        })
    }

    /// Barycentric coordinates of `q`, enclosed.
    fn lambda(&self, q: Uv) -> [Iv; 3] {
        let p = |x: f64| Iv::point(x);
        let w = &self.w;
        let r = [p(q.0) - p(w[0].0), p(q.1) - p(w[0].1)];
        let (e1, e2) = (
            [p(w[1].0) - p(w[0].0), p(w[1].1) - p(w[0].1)],
            [p(w[2].0) - p(w[0].0), p(w[2].1) - p(w[0].1)],
        );
        let l1 = (r[0] * e2[1] - e2[0] * r[1]) / self.det;
        let l2 = (e1[0] * r[1] - r[0] * e1[1]) / self.det;
        [p(1.0) - l1 - l2, l1, l2]
    }

    /// Ranges of the barycentric coordinates over piece `t` (each is
    /// affine, so its range is the hull of its vertex values).
    fn lambda_over(&self, t: &[Uv]) -> [Iv; 3] {
        let mut out: Option<[Iv; 3]> = None;
        for &q in t {
            let l = self.lambda(q);
            out = Some(match out {
                None => l,
                Some(acc) => [acc[0].hull(l[0]), acc[1].hull(l[1]), acc[2].hull(l[2])],
            });
        }
        out.unwrap_or([Iv::point(0.0); 3])
    }

    /// `φ` over the points `t` (their hull, with the bubbles over the
    /// hull's barycentric ranges): a box holding `φ(conv t)`.
    pub(super) fn map_box(&self, t: &[Uv]) -> (Iv, Iv) {
        let l = self.lambda_over(t);
        let hull = |k: usize| {
            let xs = t.iter().map(|q| if k == 0 { q.0 } else { q.1 });
            Iv::new(
                xs.clone().fold(f64::INFINITY, f64::min),
                xs.fold(f64::NEG_INFINITY, f64::max),
            )
        };
        let mut out = [hull(0), hull(1)];
        for k in 0..3 {
            let bubble = (l[k] * l[(k + 1) % 3]).scale(4.0);
            for (o, d) in out.iter_mut().zip(self.d[k]) {
                *o = *o + bubble * d;
            }
        }
        (out[0], out[1])
    }

    /// `∂φ_k / ∂y_j` over piece `t`, as `[k][j]`.
    fn jacobian(&self, t: &[Uv]) -> [[Iv; 2]; 2] {
        let l = self.lambda_over(t);
        let mut j = [
            [Iv::point(1.0), Iv::point(0.0)],
            [Iv::point(0.0), Iv::point(1.0)],
        ];
        for e in 0..3 {
            let (a, b) = (e, (e + 1) % 3);
            for jj in 0..2 {
                let dbubble = (l[a] * self.grad[b][jj] + l[b] * self.grad[a][jj]).scale(4.0);
                for (k, row) in j.iter_mut().enumerate() {
                    row[jj] = row[jj] + dbubble * self.d[e][k];
                }
            }
        }
        j
    }

    /// `∂_i∂_j φ_k` (constant), as `[k][i][j]`.
    fn hessian(&self) -> [[[Iv; 2]; 2]; 2] {
        let mut h = [[[Iv::point(0.0); 2]; 2]; 2];
        for e in 0..3 {
            let (a, b) = (e, (e + 1) % 3);
            for i in 0..2 {
                for jj in 0..2 {
                    let dd = (self.grad[a][i] * self.grad[b][jj]
                        + self.grad[a][jj] * self.grad[b][i])
                        .scale(4.0);
                    for (k, hk) in h.iter_mut().enumerate() {
                        hk[i][jj] = hk[i][jj] + dd * self.d[e][k];
                    }
                }
            }
        }
        h
    }

    /// Bounds `(|G_11|, |G_12|, |G_22|)` of the second partials of
    /// `S∘φ` over piece `t`, from a jet of `S` (to order 2) over a box
    /// holding `φ(t)`.
    pub(super) fn second_partials(&self, t: &[Uv], jet: &SurfaceJet) -> [f64; 3] {
        let j = self.jacobian(t);
        let h = self.hessian();
        let (s_u, s_v) = (jet.partial(1, 0), jet.partial(0, 1));
        let (s_uu, s_uv, s_vv) = (jet.partial(2, 0), jet.partial(1, 1), jet.partial(0, 2));
        let g = |a: usize, b: usize| -> IvVec3 {
            s_uu.scale(j[0][a] * j[0][b])
                + s_uv.scale(j[0][a] * j[1][b] + j[1][a] * j[0][b])
                + s_vv.scale(j[1][a] * j[1][b])
                + s_u.scale(h[0][a][b])
                + s_v.scale(h[1][a][b])
        };
        [g(0, 0), g(0, 1), g(1, 1)].map(|x| x.norm().hi)
    }
}

#[cfg(test)]
#[path = "chord_bend_tests.rs"]
mod tests;
