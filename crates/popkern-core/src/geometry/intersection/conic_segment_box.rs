//! Exact axis-aligned boxes of conic curve segments (PK-785).
//!
//! The general surface-curve solver culls (curve segment, surface cell)
//! pairs by box overlap. A conic's world coordinates are sinusoids of
//! its angle, `c + U cos φ + V sin φ`, so a segment's box is exact from
//! its end angles and the stationary angles between them.
//!
//! The surface side keeps its padded lattice box: for quadric cells
//! the sag/3 tail bound covers (measured: exact quadric cell boxes left
//! the PK-781/785 sweep unchanged once the curve side was exact, and
//! were removed as unmeasured).

use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// `c·cos t + d·sin t + e·t + f`: a sinusoid plus a linear part.
#[derive(Clone, Copy)]
struct Wave {
    c: f64,
    d: f64,
    e: f64,
    f: f64,
}

impl Wave {
    fn at(&self, t: f64) -> f64 {
        self.c * t.pcos() + self.d * t.psin() + self.e * t + self.f
    }

    /// `(min, max)` over `[t0, t1]`: the endpoints and, for the
    /// sinusoid part, every stationary angle inside the window.
    fn range(&self, (t0, t1): (f64, f64)) -> (f64, f64) {
        let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        let mut min = self.at(lo).min(self.at(hi));
        let mut max = self.at(lo).max(self.at(hi));
        if self.c != 0.0 || self.d != 0.0 {
            // d/dt (c cos t + d sin t) = 0 at t = atan2(d, c) + kπ.
            let base = self.d.patan2(self.c);
            let pi = std::f64::consts::PI;
            let mut k = ((lo - base) / pi).ceil();
            while base + k * pi <= hi {
                let v = self.at(base + k * pi);
                min = min.min(v);
                max = max.max(v);
                k += 1.0;
            }
        }
        (min, max)
    }
}

fn component(v: Vec3, k: usize) -> f64 {
    [v.x, v.y, v.z][k]
}

fn point_component(p: Point3, k: usize) -> f64 {
    [p.x, p.y, p.z][k]
}

/// The exact box of one Bézier segment of a conic (circle, arc or
/// ellipse), padded for roundoff; `None` for other curves. The segment's
/// angular range comes from its own end and mid points, which lie on
/// the conic, so the box is the conic's over exactly that range. A
/// 16-sample box missed an arc's apex by up to R(1 − cos(Δ/2)) — 2.4e-3
/// on a radius-5 arc — while a rod of radius 0.02 placed 0.999 r past
/// the arc is pierced only 2e-5 deep: every pair was culled (PK-785).
pub(crate) fn exact_conic_segment_box(
    curve: &crate::geometry::curve::CurveGeometry,
    segment: &crate::geometry::nurbs_curve::NurbsCurve,
    t_range: (f64, f64),
) -> Option<(Point3, Point3)> {
    use crate::geometry::curve::{CurveEval, CurveGeometry};
    let (center, u, v) = match curve {
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        }
        | CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            let (e1, e2) = crate::geometry::analytic_surfaces::orthonormal_basis(*normal);
            (*center, e1 * *radius, e2 * *radius)
        }
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => (*center, *u_basis, *v_basis),
        _ => return None,
    };
    // p − c = U cos φ + V sin φ, solved through the 2×2 Gram system.
    let (uu, uv, vv) = (u.dot(&u), u.dot(&v), v.dot(&v));
    let det = uu * vv - uv * uv;
    if det.abs() <= 1e-24 * uu.max(vv).powi(2) {
        return None;
    }
    let angle = |t: f64| {
        let w = segment.point_at(t) - center;
        let (wu, wv) = (w.dot(&u), w.dot(&v));
        let cos = (vv * wu - uv * wv) / det;
        let sin = (uu * wv - uv * wu) / det;
        sin.patan2(cos)
    };
    let pi = std::f64::consts::PI;
    let near = |phi: f64, to: f64| phi + ((to - phi) / (2.0 * pi)).round() * 2.0 * pi;
    let a = angle(t_range.0);
    let m = near(angle(0.5 * (t_range.0 + t_range.1)), a);
    let b = near(angle(t_range.1), m);
    if !((a <= m && m <= b) || (b <= m && m <= a)) {
        return None;
    }
    let mut lo = [0.0; 3];
    let mut hi = [0.0; 3];
    for k in 0..3 {
        let wave = Wave {
            c: component(u, k),
            d: component(v, k),
            e: 0.0,
            f: point_component(center, k),
        };
        let (x, y) = wave.range((a, b));
        let pad = 1e-12 * (1.0 + x.abs().max(y.abs()));
        lo[k] = x - pad;
        hi[k] = y + pad;
    }
    Some((
        Point3::new(lo[0], lo[1], lo[2]),
        Point3::new(hi[0], hi[1], hi[2]),
    ))
}

#[cfg(test)]
#[path = "conic_segment_box_tests.rs"]
mod tests;
