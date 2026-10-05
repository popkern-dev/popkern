//! Parameter windows on a CLOSED record source for the far-tier
//! endpoint re-trim (PK-673).
//!
//! A clip window on a closed record may run THROUGH the source's
//! parameter seam (`clip_curve_wrapped`): the kept arc of a coaxial
//! circle crossing a multi-turn thread band straddles the seam, so
//! its two ends project onto the source at, e.g., t = 0.761 and
//! t = 0.009 while the piece runs 0.761 → 1.009. Read raw, the
//! re-trim's window [t*, t_other] ran the wrong way round the loop
//! (decreasing parameters, the complement arc), the refit declined
//! silently, and both ends stayed 6–7e-3 from their exact junctions
//! — the crest arc of the default M3 nut's top mouth chamfer met
//! neither flank record and the cellular body went non-manifold.
//!
//! On a closed source every parameter is a class modulo the period;
//! the piece's ends are unwrapped so its traversal passes through its
//! own midpoint, junction parameters are lifted next to the end they
//! correct, and evaluation folds back into the source domain. Open
//! sources (no period) pass through every helper unchanged.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv_near;
use popkern_core::tolerance::Tolerance;

/// Chart samples per re-authored pcurve of a seam re-trim (the far
/// tier's re-trim density).
const RETRIM_CHART_SAMPLES: usize = 17;

/// The parameter period of a record source: `Some` only on a closed
/// source with a non-degenerate domain.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SourceWindow {
    s0: f64,
    period: Option<f64>,
}

impl SourceWindow {
    /// The window of `source`'s 3D curve.
    pub(crate) fn of(source: &CurveOnSurface) -> Self {
        let (s0, s1) = source.curve_3d.parameter_range();
        let span = s1 - s0;
        let period = (source.is_closed && span.is_finite() && span > 0.0).then_some(span);
        Self { s0, period }
    }

    /// Fold an unwrapped parameter back into the source domain.
    pub(crate) fn wrap(&self, t: f64) -> f64 {
        match self.period {
            Some(p) => self.s0 + (t - self.s0).rem_euclid(p),
            None => t,
        }
    }

    /// The lift of `t` nearest `anchor`.
    pub(crate) fn lift_near(&self, t: f64, anchor: f64) -> f64 {
        match self.period {
            Some(p) => t + ((anchor - t) / p).round() * p,
            None => t,
        }
    }

    /// Unwrap a piece's end parameters (`t_this` kept) so the
    /// traversal `t_this → t_other` passes through the lift of the
    /// piece's midpoint parameter `t_mid` and spans less than one
    /// period. A piece that does not run through the seam keeps its
    /// raw parameters exactly; an ambiguous configuration (both ends
    /// one class) is returned raw.
    pub(crate) fn unwrap_piece(&self, t_this: f64, t_other: f64, t_mid: f64) -> (f64, f64) {
        let Some(p) = self.period else {
            return (t_this, t_other);
        };
        let between = |a: f64, b: f64| {
            let (lo, hi) = (a.min(b), a.max(b));
            let m = self.lift_near(t_mid, 0.5 * (lo + hi));
            hi - lo < p && m > lo && m < hi
        };
        if between(t_this, t_other) {
            return (t_this, t_other);
        }
        for k in [-1.0, 1.0] {
            let lifted = t_other + k * p;
            if between(t_this, lifted) {
                return (t_this, lifted);
            }
        }
        (t_this, t_other)
    }

    /// Evaluate `curve` at an unwrapped parameter.
    pub(crate) fn point_at(&self, curve: &NurbsCurve, t: f64) -> Point3 {
        curve.point_at(self.wrap(t))
    }

    /// Make a run of chart samples continuous across the source seam:
    /// on a closed source each sample is folded by whole chart periods
    /// toward its predecessor (a sample that already continues its
    /// predecessor is untouched). No-op on open sources.
    pub(crate) fn unwrap_chart_run(&self, run: &mut [Point3], surface: &SurfaceGeometry) {
        if self.period.is_none() {
            return;
        }
        let (pu, pv) = surface.parameter_periodicity();
        for i in 1..run.len() {
            let prev = run[i - 1];
            let q = &mut run[i];
            q.x = fold_toward(q.x, prev.x, pu);
            q.y = fold_toward(q.y, prev.y, pv);
        }
    }
}

impl SourceWindow {
    /// Exact sub-curve of the closed source `curve` over the
    /// unwrapped window `[lo, hi]` when that window runs THROUGH the
    /// source seam: the head to the domain end and the tail from the
    /// domain start, spliced (C0 at the seam point, the
    /// `clip_curve_wrapped` construction). The result's domain is
    /// the window itself (shifted by whole periods). `None` for open
    /// sources and for windows that do not cross the seam.
    pub(crate) fn seam_crossing_subcurve(
        &self,
        curve: &NurbsCurve,
        lo: f64,
        hi: f64,
    ) -> Option<NurbsCurve> {
        let p = self.period?;
        if !(hi > lo && hi - lo < p) {
            return None;
        }
        let lo_w = self.wrap(lo);
        let hi_w = lo_w + (hi - lo);
        let s1 = self.s0 + p;
        if hi_w <= s1 {
            return None;
        }
        let head = crate::clip::clip_curve_to_crossings(curve, lo_w, s1)?;
        let tail = crate::clip::clip_curve_to_crossings(curve, self.s0, hi_w - p)?;
        let join = curve.point_at(self.s0).distance_to(&curve.point_at(s1)) + 1e-9;
        head.concatenate_clamped(&tail, join).ok()
    }

    /// Re-trim a piece EXACTLY when its corrected window `[n0, n1]`
    /// (unwrapped, traversal order) runs through the closed source's
    /// seam: the 3D curve is the source's own sub-curve
    /// ([`Self::seam_crossing_subcurve`]) with its ends pinned onto
    /// `j` (this end) and `other` (the far end), and each chart
    /// pcurve is re-authored by projecting the 3D samples onto its
    /// surface, continuous across the chart seam. Re-sampling and
    /// re-fitting across the seam instead bent the 3D curve 1.3e-5
    /// off the source exactly AT the seam (the closed source's
    /// parameter speed jumps there), and the seam-crossing junction
    /// node no longer split the record (PK-673). Declines (`false`,
    /// piece untouched) unless the window crosses the seam, runs in
    /// source order, and both pinned ends already lie on the source
    /// within the linear tolerance — then the end moves are below
    /// tolerance and cannot bend the curve.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn exact_seam_retrim(
        &self,
        piece: &mut CurveOnSurface,
        source: &CurveOnSurface,
        (n0, n1): (f64, f64),
        end: bool,
        j: Point3,
        other: Point3,
        (surf_a, surf_b): (&SurfaceGeometry, &SurfaceGeometry),
        tolerance: &Tolerance,
    ) -> bool {
        if n1 <= n0 {
            return false;
        }
        let Some(mut curve) = self.seam_crossing_subcurve(&source.curve_3d, n0, n1) else {
            return false;
        };
        let (start, finish) = if end { (other, j) } else { (j, other) };
        let (d0, d1) = curve.parameter_range();
        if curve.point_at(d0).distance_to(&start) > tolerance.linear()
            || curve.point_at(d1).distance_to(&finish) > tolerance.linear()
        {
            return false;
        }
        if let Some(cp) = curve.control_points.first_mut() {
            *cp = start;
        }
        if let Some(cp) = curve.control_points.last_mut() {
            *cp = finish;
        }
        let n = RETRIM_CHART_SAMPLES;
        let params: Vec<f64> = (0..n)
            .map(|k| d0 + (d1 - d0) * k as f64 / (n - 1) as f64)
            .collect();
        let mut charts = Vec::with_capacity(2);
        for (src_uv, surf) in [(&source.curve_uv1, surf_a), (&source.curve_uv2, surf_b)] {
            let mut run = Vec::with_capacity(n);
            for &t in &params {
                let hint = self.point_at(src_uv, t);
                let Some((u, v)) =
                    project_point_to_uv_near(surf, &curve.point_at(t), hint.x, hint.y)
                else {
                    return false;
                };
                run.push(Point3::new(u, v, 0.0));
            }
            // Anchor the run on the source chart's branch, then keep
            // it continuous across the chart seam.
            let (pu, pv) = surf.parameter_periodicity();
            let anchor = self.point_at(src_uv, params[0]);
            run[0].x = fold_toward(run[0].x, anchor.x, pu);
            run[0].y = fold_toward(run[0].y, anchor.y, pv);
            for i in 1..run.len() {
                let prev = run[i - 1];
                run[i].x = fold_toward(run[i].x, prev.x, pu);
                run[i].y = fold_toward(run[i].y, prev.y, pv);
            }
            let Ok(fitted) = popkern_core::geometry::fitting::fit_curve(&run, &params, 3, n) else {
                return false;
            };
            charts.push(fitted);
        }
        let (Some(uv2), Some(uv1)) = (charts.pop(), charts.pop()) else {
            return false;
        };
        piece.curve_3d = curve;
        piece.curve_uv1 = uv1;
        piece.curve_uv2 = uv2;
        true
    }
}

/// Fold `x` by whole `period`s to the branch nearest `anchor`.
pub(crate) fn fold_toward(x: f64, anchor: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p.is_finite() && p > 0.0 => x + ((anchor - x) / p).round() * p,
        _ => x,
    }
}

#[cfg(test)]
#[path = "closed_source_window_tests.rs"]
mod tests;
