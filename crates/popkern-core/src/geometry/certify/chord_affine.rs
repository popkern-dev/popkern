//! The affine map of a parameter triangle onto its emitted triangle,
//! the pairing the triangle chordal bound compares the surface with
//! (see `chord_triangle`).

use super::super::interval::Iv;
use super::super::interval_vec::IvVec3;
use super::super::surface_jet::ParamMap;
use super::chord_triangle::Uv;
use crate::geometry::Point3;

/// The affine map of the parameter triangle onto the emitted one.
pub(super) struct Affine {
    w: [Uv; 3],
    p: [Point3; 3],
    det: f64,
}

impl Affine {
    pub(super) fn new(w: [Uv; 3], p: [Point3; 3]) -> Self {
        let det = (w[1].0 - w[0].0) * (w[2].1 - w[0].1) - (w[2].0 - w[0].0) * (w[1].1 - w[0].1);
        Affine { w, p, det }
    }

    /// `L` as `l0 + u·lu + v·lv` in the parameters shifted by `shift`
    /// (a spline chart's), enclosed; `None` for a degenerate triangle.
    pub(super) fn param_map(&self, shift: Uv) -> Option<ParamMap> {
        let p = |x: f64| Iv::point(x);
        let w = &self.w;
        let (e1, e2) = (
            [p(w[1].0) - p(w[0].0), p(w[1].1) - p(w[0].1)],
            [p(w[2].0) - p(w[0].0), p(w[2].1) - p(w[0].1)],
        );
        let det = e1[0] * e2[1] - e2[0] * e1[1];
        if det.contains_zero() || !det.is_finite() {
            return None;
        }
        // ∇λ1 = (e2.v, −e2.u) / det, ∇λ2 = (−e1.v, e1.u) / det.
        let (g1, g2) = ([e2[1] / det, -e2[0] / det], [-e1[1] / det, e1[0] / det]);
        let p0 = IvVec3::from_point(self.p[0]);
        let (d1, d2) = (
            IvVec3::from_point(self.p[1]) - p0,
            IvVec3::from_point(self.p[2]) - p0,
        );
        let lu = d1.scale(g1[0]) + d2.scale(g2[0]);
        let lv = d1.scale(g1[1]) + d2.scale(g2[1]);
        // L(c) = p0 + lu (c_u − s_u − w0_u) + lv (c_v − s_v − w0_v).
        let l0 = p0 - lu.scale(p(w[0].0) + p(shift.0)) - lv.scale(p(w[0].1) + p(shift.1));
        Some(ParamMap {
            u0: p(0.0),
            v0: p(0.0),
            hu: p(1.0),
            hv: p(1.0),
            p: l0,
            q: lu,
            r: lv,
            t: IvVec3::new(p(0.0), p(0.0), p(0.0)),
        })
    }

    /// `L(q)`; a degenerate parameter triangle maps every point to the
    /// nearest vertex's image (its pieces are then single points or
    /// segments whose discrepancies carry the bound).
    pub(super) fn at(&self, q: Uv) -> Point3 {
        let w = &self.w;
        if self.det.abs() <= f64::MIN_POSITIVE {
            let d = |a: Uv| (a.0 - q.0).powi(2) + (a.1 - q.1).powi(2);
            let k = (0..3)
                .min_by(|&i, &j| d(w[i]).total_cmp(&d(w[j])))
                .unwrap_or(0);
            return self.p[k];
        }
        let l1 =
            ((q.0 - w[0].0) * (w[2].1 - w[0].1) - (w[2].0 - w[0].0) * (q.1 - w[0].1)) / self.det;
        let l2 =
            ((w[1].0 - w[0].0) * (q.1 - w[0].1) - (q.0 - w[0].0) * (w[1].1 - w[0].1)) / self.det;
        let l0 = 1.0 - l1 - l2;
        let p = &self.p;
        Point3::new(
            l0 * p[0].x + l1 * p[1].x + l2 * p[2].x,
            l0 * p[0].y + l1 * p[1].y + l2 * p[2].y,
            l0 * p[0].z + l1 * p[1].z + l2 * p[2].z,
        )
    }
}
