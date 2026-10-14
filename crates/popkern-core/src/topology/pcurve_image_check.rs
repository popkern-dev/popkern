//! Image-off-edge measures for kept pcurves (split from
//! `pcurve_recompute_helpers` under the file-size cap): does an
//! existing pcurve's surface image still trace its edge?

use crate::topology::BRep;

impl BRep {
    /// Every surface image of `pcurve` (three interior stations) lies
    /// within `gate` of `heh`'s edge curve, as measured by
    /// [`BRep::pcurve_image_off_edge`]. Decided first by the samples
    /// that search would visit at each station's own fraction (and its
    /// reverse): those are among its 257 samples, so a bound under the
    /// gate is the search's verdict without the search (it cost ~90 µs
    /// per half-edge, every kept pcurve of every recomputed face —
    /// PK-802: 22 ms on a 64-slot disc's cap face).
    pub(super) fn pcurve_image_within(
        &self,
        heh: crate::topology::HalfEdgeHandle,
        pcurve: &crate::geometry::pcurve::Pcurve,
        surface: &crate::geometry::surface::SurfaceGeometry,
        gate: f64,
    ) -> bool {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::surface::SurfaceEval;
        let Some(curve) = self
            .half_edges
            .get(heh)
            .and_then(|h| self.edges.get(h.edge))
            .and_then(|e| e.curve.as_ref())
        else {
            return 0.0 <= gate;
        };
        let (c0, c1) = (curve.t_start, curve.t_end);
        let (p0, p1) = pcurve.parameter_range();
        let bounded = [0.25, 0.5, 0.75].iter().all(|f| {
            let q = pcurve.point_at(p0 + (p1 - p0) * f);
            let image = surface.point_at(q.x, q.y);
            [*f, 1.0 - f].iter().any(|k| {
                curve
                    .geometry
                    .point_at(c0 + (c1 - c0) * k)
                    .distance_to(&image)
                    <= gate
            })
        });
        bounded || self.pcurve_image_off_edge(heh, pcurve, surface) <= gate
    }

    /// Worst distance from the surface images of `pcurve` (three
    /// interior stations) to `heh`'s edge curve: the nearest of 256
    /// edge samples, refined by golden section on its neighbourhood.
    /// 0 when the edge has no curve.
    fn pcurve_image_off_edge(
        &self,
        heh: crate::topology::HalfEdgeHandle,
        pcurve: &crate::geometry::pcurve::Pcurve,
        surface: &crate::geometry::surface::SurfaceGeometry,
    ) -> f64 {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::surface::SurfaceEval;
        const SAMPLES: usize = 256;
        let Some(curve) = self
            .half_edges
            .get(heh)
            .and_then(|h| self.edges.get(h.edge))
            .and_then(|e| e.curve.as_ref())
        else {
            return 0.0;
        };
        let (c0, c1) = (curve.t_start, curve.t_end);
        let at = |k: f64| c0 + (c1 - c0) * k;
        let (p0, p1) = pcurve.parameter_range();
        [0.25, 0.5, 0.75]
            .iter()
            .map(|f| {
                let q = pcurve.point_at(p0 + (p1 - p0) * f);
                let image = surface.point_at(q.x, q.y);
                let dist = |k: f64| curve.geometry.point_at(at(k)).distance_to(&image);
                // First minimal sample (each distance evaluated once).
                let mut best = (0usize, f64::INFINITY);
                for i in 0..=SAMPLES {
                    let d = dist(i as f64 / SAMPLES as f64);
                    if d.total_cmp(&best.1).is_lt() {
                        best = (i, d);
                    }
                }
                let best = best.0;
                let (mut lo, mut hi) = (
                    (best as f64 - 1.0).max(0.0) / SAMPLES as f64,
                    (best as f64 + 1.0).min(SAMPLES as f64) / SAMPLES as f64,
                );
                for _ in 0..40 {
                    let m1 = lo + (hi - lo) * 0.382;
                    let m2 = lo + (hi - lo) * 0.618;
                    if dist(m1) <= dist(m2) {
                        hi = m2;
                    } else {
                        lo = m1;
                    }
                }
                dist(0.5 * (lo + hi))
            })
            .fold(0.0, f64::max)
    }
}
