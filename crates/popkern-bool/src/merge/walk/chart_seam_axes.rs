//! Seam axes of a face on a doubly-periodic chart (PK-459 stage 2;
//! docs/architecture/torus-cut-domain-chart-walk.md §3.2).
//!
//! A closed profile revolved a full circle is ONE face wrapping both
//! periods of its torus chart, bounded by its seam square: a u-seam
//! pair (constant u) and a v-seam pair (constant v). Such a face is
//! walked on its cut domain D = [u₀, u₀ + Pu] × [v₀, v₀ + Pv], the
//! seams as D's borders. Recognition reads the face's STORED seam
//! half-edges (both sides of the edge on this face) and their own
//! pcurves, never sampled boundary gaps: a u-seam's full-period v run
//! sampled at a few stations leaves gaps that read as a v-band
//! (measured 2026-10-22: the bare torus face took the band path and
//! served loops winding −1 round u).

use popkern_core::topology::{BRep, FaceHandle};

/// A chart axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ChartAxis {
    /// The u parameter.
    U,
    /// The v parameter.
    V,
}

/// The face's seam coordinates per axis: `u` is the azimuth of a
/// seam constant in u (the u-seam), `v` the height of a seam constant
/// in v (the v-seam). Each is folded into `[0, period)`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct SeamAxes {
    /// The u-seam's u, when the face carries one.
    pub u: Option<f64>,
    /// The v-seam's v, when the face carries one.
    pub v: Option<f64>,
}

impl SeamAxes {
    /// The cut domain's low corner `(u₀, v₀)` when the face carries a
    /// seam on BOTH axes.
    pub(crate) fn corner(&self) -> Option<(f64, f64)> {
        Some((self.u?, self.v?))
    }
}

/// Which axis a chart run is constant on, within `eps` of the
/// period-scale: `U` when every sample shares one u, `V` likewise.
/// `None` for runs that move on both axes (or neither).
pub(crate) fn constant_axis(uv: &[(f64, f64)], eps_u: f64, eps_v: f64) -> Option<ChartAxis> {
    let spread = |f: &dyn Fn(&(f64, f64)) -> f64| {
        uv.iter()
            .map(f)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (lo.min(x), hi.max(x))
            })
    };
    let (u_lo, u_hi) = spread(&|q| q.0);
    let (v_lo, v_hi) = spread(&|q| q.1);
    let (du, dv) = (u_hi - u_lo, v_hi - v_lo);
    match (du <= eps_u, dv <= eps_v) {
        (true, false) => Some(ChartAxis::U),
        (false, true) => Some(ChartAxis::V),
        _ => None,
    }
}

/// The seam axes of `face`, whose surface has chart periods
/// `(pu, pv)`. Reads every stored half-edge whose twin lies on the
/// same face (a seam) and classifies it by its own pcurve's constant
/// axis. Faces without pcurves on their seams report no axis there.
pub(crate) fn face_seam_axes(brep: &BRep, face: FaceHandle, periods: (f64, f64)) -> SeamAxes {
    let (pu, pv) = periods;
    let mut axes = SeamAxes::default();
    let Some(fd) = brep.faces.get(face) else {
        return axes;
    };
    let (eps_u, eps_v) = (pu * 1.0e-6, pv * 1.0e-6);
    for &lh in std::iter::once(&fd.outer_loop).chain(fd.inner_loops.iter()) {
        let Some(lp) = brep.loops.get(lh) else {
            continue;
        };
        let start = lp.half_edge;
        let mut he = start;
        // A loop never holds more half-edges than the body does.
        for _ in 0..brep.half_edges.len() {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            let twin_on_face = brep
                .half_edges
                .get(hd.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .is_some_and(|tl| tl.face == face);
            if twin_on_face && let Some(pc) = hd.pcurve.as_ref() {
                let (t0, t1) = pc.parameter_range();
                let uv: Vec<(f64, f64)> = (0..=8)
                    .map(|i| pc.point_at(t0 + (t1 - t0) * i as f64 / 8.0))
                    .filter(|q| q.x.is_finite() && q.y.is_finite())
                    .map(|q| (q.x, q.y))
                    .collect();
                if uv.len() >= 2 {
                    match constant_axis(&uv, eps_u, eps_v) {
                        Some(ChartAxis::U) => axes.u = axes.u.or(Some(uv[0].0.rem_euclid(pu))),
                        Some(ChartAxis::V) => axes.v = axes.v.or(Some(uv[0].1.rem_euclid(pv))),
                        None => {}
                    }
                }
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    axes
}
