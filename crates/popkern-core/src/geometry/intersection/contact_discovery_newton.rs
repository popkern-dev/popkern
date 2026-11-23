//! Damped Newton on the gap's gradient (split from `contact_discovery.rs`
//! for the file-size limit).

use super::ParamWindow;
use super::{Gap, NEWTON_STEPS, full_periods};
use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_math::PortableMath;

/// Damped Newton on `∇g = 0` from `seed`, derivatives by central
/// differences of step `h` (parameter units), kept inside `w`; the point
/// reached is polished with the gradient's fourth-order stencil.
pub(super) fn newton_critical(
    s1: &SurfaceGeometry,
    gap: &Gap<'_>,
    seed: (f64, f64),
    h: f64,
    w: ParamWindow,
) -> Option<((f64, f64), f64)> {
    // Each parameter's difference step in proportion to its window: one
    // step for an angle and a length (a cylinder at scale 1000: 0.2 in
    // both) differenced the angle over a fifth of a radian, and Newton
    // stalled 0.1 from the touch with its gap six tolerances.
    let (wu, wv) = ((w.0.1 - w.0.0).abs(), (w.1.1 - w.1.0).abs());
    let wide = wu.max(wv);
    let (hu, hv) = if wide > 0.0 && wide.is_finite() {
        (h * (wu / wide).max(1e-3), h * (wv / wide).max(1e-3))
    } else {
        (h, h)
    };
    let walk = Walk { s1, gap, w };
    let mut uv = walk.run(seed, hu, hv, NEWTON_STEPS, false)?;
    // A central difference's error grows with its step squared: the
    // window's step (a sixth of a length unit at scale 1000) left the
    // point 4e-5 from a cone's touch on a cylinder, where the normals
    // already part by more than the tolerance allows. The gradient's
    // fourth-order stencil at the same step removes it; a shorter step
    // would instead raise the rounding noise, and on a flat touch (the
    // gap quartic) Newton stops where the gradient sinks under it.
    if let Some(polished) = walk.run(uv, hu, hv, NEWTON_STEPS / 4, true) {
        uv = polished;
    }
    Some((uv, walk.g(uv.0, uv.1)?))
}

/// One Newton walk's surroundings.
struct Walk<'a> {
    s1: &'a SurfaceGeometry,
    gap: &'a Gap<'a>,
    w: ParamWindow,
}

impl Walk<'_> {
    fn g(&self, u: f64, v: f64) -> Option<f64> {
        self.gap.at(self.s1, (u, v)).map(|(x, _)| x)
    }

    /// The gradient at `(u, v)` by central differences of steps `hu`,
    /// `hv`, of the fourth order when `fine`.
    fn gradient(&self, u: f64, v: f64, hu: f64, hv: f64, fine: bool) -> Option<(f64, f64)> {
        let g = |u: f64, v: f64| self.g(u, v);
        let d = |f: &dyn Fn(f64) -> Option<f64>, h: f64| -> Option<f64> {
            let (p, m) = (f(h)?, f(-h)?);
            if !fine {
                return Some((p - m) / (2.0 * h));
            }
            let (p2, m2) = (f(2.0 * h)?, f(-2.0 * h)?);
            Some((8.0 * (p - m) - (p2 - m2)) / (12.0 * h))
        };
        Some((d(&|t| g(u + t, v), hu)?, d(&|t| g(u, v + t), hv)?))
    }

    /// Levenberg–Marquardt from `start` with steps `hu`, `hv`, the
    /// gradient of the fourth order when `fine`.
    fn run(
        &self,
        start: (f64, f64),
        hu: f64,
        hv: f64,
        steps: usize,
        fine: bool,
    ) -> Option<(f64, f64)> {
        let (s1, w) = (self.s1, self.w);
        let g = |u: f64, v: f64| self.g(u, v);
        // A step leaving the window wraps across a whole period, else
        // stops at the edge (a contact at a pole lies on it).
        let (wrap_u, wrap_v) = full_periods(s1, w);
        let (pu, pv) = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(s1);
        let settle = |x: f64, (lo, hi): (f64, f64), wrap: bool, period: Option<f64>| {
            let (lo, hi) = (lo.min(hi), lo.max(hi));
            match period.filter(|_| wrap) {
                Some(p) => lo + crate::geometry::chart::ops::wrap_into(x - lo, 0.0, p),
                None => x.clamp(lo, hi),
            }
        };
        let (mut u, mut v) = start;
        let mut mu = 1e-6;
        for _ in 0..steps {
            let c = g(u, v)?;
            let (e, wv, n, s) = (g(u + hu, v)?, g(u - hu, v)?, g(u, v + hv)?, g(u, v - hv)?);
            let (ne, nw, se, sw) = (
                g(u + hu, v + hv)?,
                g(u - hu, v + hv)?,
                g(u + hu, v - hv)?,
                g(u - hu, v - hv)?,
            );
            let (gu, gv) = self.gradient(u, v, hu, hv, fine)?;
            let huu = (e - 2.0 * c + wv) / (hu * hu);
            let hvv = (n - 2.0 * c + s) / (hv * hv);
            let huv = (ne - nw - se + sw) / (4.0 * hu * hv);
            // The chart's speeds: lengths on the surface per unit of
            // each parameter.
            let speed = |a: Point3, b: Point3, step: f64| a.distance_to(&b) / (2.0 * step);
            let su = speed(s1.point_at(u + hu, v), s1.point_at(u - hu, v), hu);
            let sv = speed(s1.point_at(u, v + hv), s1.point_at(u, v - hv), hv);
            let floor = 1e-3 * su.max(sv).max(f64::MIN_POSITIVE);
            let (su, sv) = (su.max(floor), sv.max(floor));
            // Levenberg–Marquardt on the gradient: (H + μ|H|) Δ = −∇g,
            // the damping per unit length too. Taken in parameter units
            // (|H| the larger of an angle's curvature per radian² and a
            // length's per length², at scale 0.001 a millionfold apart),
            // a millionth of it cancelled a saddle's curvature along the
            // angle: Newton stepped three times past a tangent crossing,
            // and then the wrong way.
            let scale = (huu.abs() / (su * su))
                .max(hvv.abs() / (sv * sv))
                .max(1e-300);
            let (a, b, d) = (huu + mu * scale * su * su, huv, hvv + mu * scale * sv * sv);
            let det = a * d - b * b;
            if det.abs() <= 1e-300 {
                mu *= 10.0;
                continue;
            }
            let du = -(d * gu - b * gv) / det;
            let dv = -(-b * gu + a * gv) / det;
            let (nu, nv) = (
                settle(u + du, w.0, wrap_u, pu),
                settle(v + dv, w.1, wrap_v, pv),
            );
            // Accept a step that does not grow the gradient, measured
            // per unit length on the surface: in parameter units an
            // angle's component at radius 1000 outweighed a length's a
            // thousandfold, and every Newton step along a shallow valley
            // read as worse (the damping then crawled, stalling with the
            // gap 2.5).
            let grad_now = (gu / su).phypot(gv / sv);
            let next = self.gradient(nu, nv, hu, hv, fine)?;
            let grad_next = (next.0 / su).phypot(next.1 / sv);
            if grad_next <= grad_now {
                (u, v) = (nu, nv);
                mu = (mu * 0.3).max(1e-12);
                if du.abs() <= 1e-6 * hu.max(1e-12) && dv.abs() <= 1e-6 * hv.max(1e-12) {
                    break;
                }
            } else {
                mu *= 10.0;
                if mu > 1e8 {
                    break;
                }
            }
        }
        Some((u, v))
    }
}
