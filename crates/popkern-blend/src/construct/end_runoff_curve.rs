//! Runoff-curve construction for fillet end runoff (PK-518).
//!
//! The runoff edge is the strip ∩ end-plane intersection between the
//! two trim-crossing anchors. It was a degree-1 polyline through 33
//! exact crossings: its chords stay in the end plane (an affine
//! combination of in-plane points) but cut inside the curved strip by
//! the chord sagitta, so the seam measured a positional gap of
//! L²/(8R) — 2.9e-5 on the pentagonal pyramid, 1.0e-4 on the square
//! one, 8.0e-4 on the triangular one at r = 0.16 — with nothing
//! recording it, which floored every such fillet's continuity
//! evidence at G0.
//!
//! The curve is now a cubic interpolant of exact crossings whose
//! accuracy is MEASURED, not assumed: after each fit the exact
//! crossing at every node midpoint is compared with the curve there,
//! and the node count doubles until the worst deviation is within the
//! target (a tenth of the model's linear tolerance) or the density
//! cap is reached. The interpolant lies in the end plane exactly (its
//! control points are affine combinations of in-plane nodes), so the
//! strip side carries the whole residual, which is returned for the
//! caller to record as the edge's tolerance when it exceeds the
//! model's linear tolerance.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, SurfaceGeometry, Vec3};

/// Node counts tried, coarse to fine (2^k + 1 so every refinement
/// keeps the previous nodes).
const NODE_COUNTS: [usize; 4] = [33, 65, 129, 257];

/// Exact strip ∩ end-plane crossings between two anchors, as a
/// function of the station fraction `f ∈ [0, 1]`.
pub(super) struct RunoffCrossing<'a> {
    /// The strip surface.
    pub surface: &'a SurfaceGeometry,
    /// A point on the end plane.
    pub plane_origin: Point3,
    /// The end plane's normal (any length).
    pub plane_normal: Vec3,
    /// Anchor parameters `(u, v)` at `f = 0` and `f = 1`.
    pub anchors: ((f64, f64), (f64, f64)),
    /// The strip face's parameter box `((u0, u1), (v0, v1))`.
    pub parameter_box: ((f64, f64), (f64, f64)),
}

impl RunoffCrossing<'_> {
    fn signed(&self, p: Point3) -> f64 {
        (p - self.plane_origin).dot(&self.plane_normal) / self.plane_normal.length()
    }

    /// The exact crossing at station fraction `f`.
    ///
    /// Analytic cylinders: the curve is graph-like over the ANGLE (at
    /// a symmetric corner both anchors share the axial value and an
    /// axial march degenerates to the anchor chord — PK-097), so `u`
    /// is marched and the plane equation solved for `v` in closed
    /// form. Other strips march the parameter along which the anchors
    /// are genuinely separated (box-normalized) and bisect the other
    /// from the strip's mid-surface toward the overshooting border
    /// (PK-098 hardening); a station with no crossing takes the border
    /// point projected to the plane.
    pub(super) fn at(&self, f: f64) -> Point3 {
        let ((ua, va), (ub, vb)) = self.anchors;
        let ((bu0, bu1), (bv0, bv1)) = self.parameter_box;
        let n = self.plane_normal;
        if let SurfaceGeometry::Cylinder { axis, .. } = self.surface {
            let denom = axis.dot(&n);
            if denom.abs() > 1.0e-9 {
                let ring = self.surface.point_at(ua + (ub - ua) * f, 0.0);
                return ring + *axis * ((self.plane_origin - ring).dot(&n) / denom);
            }
        }
        let du_n = ((ua - ub) / (bu1 - bu0).max(1.0e-12)).abs();
        let dv_n = ((va - vb) / (bv1 - bv0).max(1.0e-12)).abs();
        let march_v = dv_n >= du_n;
        let (sa, sb, (b0, b1), (ca, cb)) = if march_v {
            (va, vb, (bu0, bu1), (ua, ub))
        } else {
            (ua, ub, (bv0, bv1), (va, vb))
        };
        let eval = |station: f64, cross: f64| -> Point3 {
            if march_v {
                self.surface.point_at(cross, station)
            } else {
                self.surface.point_at(station, cross)
            }
        };
        let station = sa + (sb - sa) * f;
        let mid_box = 0.5 * (b0 + b1);
        let border = if (ca + cb) * 0.5 > mid_box { b1 } else { b0 };
        let (mut lo, mut hi) = (mid_box, border);
        let d_lo = self.signed(eval(station, lo));
        let d_hi = self.signed(eval(station, hi));
        if d_lo.abs() < 1.0e-12 {
            return eval(station, lo);
        }
        if d_lo.signum() == d_hi.signum() {
            let p = eval(station, hi);
            let d = self.signed(p);
            let unit = n * (1.0 / n.length());
            return p - unit * d;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.signed(eval(station, mid)).signum() == d_lo.signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        eval(station, 0.5 * (lo + hi))
    }
}

/// Build the runoff curve from `p_prev` to `p_next` (the topology
/// anchors, which the curve's ends are snapped to) to within
/// `target` of the exact crossings where the density cap allows.
/// Returns the curve (parameter range `[0, 1]`) and its measured
/// worst midpoint deviation from the exact crossings.
pub(super) fn build_runoff_curve(
    crossing: &RunoffCrossing<'_>,
    p_prev: Point3,
    p_next: Point3,
    target: f64,
) -> Result<(NurbsCurve, f64), String> {
    let mut best: Option<(NurbsCurve, f64)> = None;
    for &count in &NODE_COUNTS {
        let params: Vec<f64> = (0..count).map(|i| i as f64 / (count - 1) as f64).collect();
        let mut points: Vec<Point3> = params.iter().map(|&f| crossing.at(f)).collect();
        points[0] = p_prev;
        points[count - 1] = p_next;
        let curve = popkern_core::geometry::fitting::interpolate_curve(&points, &params, 3)
            .map_err(|e| format!("runoff curve interpolation: {e}"))?;
        let deviation = params
            .windows(2)
            .map(|pair| {
                let f = 0.5 * (pair[0] + pair[1]);
                use popkern_core::geometry::curve::CurveEval;
                curve.point_at(f).distance_to(&crossing.at(f))
            })
            .fold(0.0_f64, f64::max);
        let done = deviation <= target;
        if best.as_ref().is_none_or(|(_, d)| deviation < *d) {
            best = Some((curve, deviation));
        }
        if done {
            break;
        }
    }
    best.ok_or_else(|| "runoff curve: no node count attempted".to_string())
}
