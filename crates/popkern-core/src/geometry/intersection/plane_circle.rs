//! Exact plane × circle (or arc) intersection (PK-758).
//!
//! A circular curve's parameter is its angle, so the signed distance to
//! a plane is `C + A cos t + B sin t`: a closed form with both roots,
//! a single tangent root, or an identical zero when the circle lies in
//! the plane. The general NURBS solver seeded one Newton start per
//! Bézier segment and found only one of two crossings when both fell
//! in the same segment (measured: a ring's notch-rim arc crossing a
//! boss's end cap at z = 148.6 and 176.4, seeded at 162.5 between
//! them). Roots are lifted by whole turns into the requested window.

use popkern_math::PortableMath;

use super::SurfaceCurveIntersection;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// The plane, by origin, unit normal and in-plane axes.
pub(crate) struct PlaneFrame<'a> {
    pub origin: &'a Point3,
    pub normal: &'a Vec3,
    pub u_axis: &'a Vec3,
    pub v_axis: &'a Vec3,
}

/// Roots of `curve` (a circle or arc, angle-parameterised) on the plane
/// within `t_range` and the plane's `(u, v)` window; `None` when
/// `curve` is not circular.
pub(crate) fn intersect_plane_circle(
    plane: PlaneFrame<'_>,
    curve: &CurveGeometry,
    t_range: (f64, f64),
    (u_range, v_range): ((f64, f64), (f64, f64)),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceCurveIntersection>> {
    let center = match curve {
        CurveGeometry::Circle { center, .. } | CurveGeometry::Arc { center, .. } => *center,
        _ => return None,
    };
    let n = plane.normal.normalized()?;
    // Signed distance d(t) = C + A cos t + B sin t, read off the curve's
    // own evaluation (its basis is the curve's, not re-derived here).
    let c = n.dot(&(center - *plane.origin));
    let a = n.dot(&(curve.point_at(0.0) - center));
    let b = n.dot(&(curve.point_at(std::f64::consts::FRAC_PI_2) - center));
    let amp = a.phypot(b);
    let lin = tolerance.linear();
    if amp <= lin {
        // Parallel planes: in the plane (an overlap) or clear of it.
        return Some(if c.abs() <= lin {
            vec![SurfaceCurveIntersection::Overlap { t_range }]
        } else {
            Vec::new()
        });
    }
    let ratio = -c / amp;
    if ratio.abs() > 1.0 + lin / amp {
        return Some(Vec::new());
    }
    let base = b.patan2(a);
    // Tangency is a distance band, as for line × cylinder (PK-789): a
    // circle whose extreme lies past the plane by less than the linear
    // tolerance touches it once, at that extreme. Its two secant roots
    // sit 2·sqrt(2·r·depth) apart there (9.5e-4 for r = 4.5 at depth
    // 1e-7), far beyond the tolerance along the curve although every
    // point between them is within it of the plane (PK-887).
    // Open at a whole tolerance, judged past the computation's own
    // roundoff: a cap typed one tolerance inside (58 − 57.999999 =
    // 9.99999997e-7) is a secant, as line × cylinder reads it.
    let scale = super::coordinate_scale(&[center, *plane.origin]) + amp;
    let tangent = amp - c.abs() < tolerance.tangency_band(scale);
    let half = if tangent {
        if ratio > 0.0 {
            0.0
        } else {
            std::f64::consts::PI
        }
    } else {
        ratio.clamp(-1.0, 1.0).pacos()
    };
    let mut roots = vec![base + half];
    if !tangent {
        roots.push(base - half);
    }
    let tau = std::f64::consts::TAU;
    let eps = tolerance.parametric_base();
    let mut out = Vec::new();
    for root in roots {
        // Every whole-turn lift inside the window.
        let mut t = root + ((t_range.0 - root) / tau).ceil() * tau;
        if t - tau >= t_range.0 - eps {
            t -= tau;
        }
        while t <= t_range.1 + eps {
            let point = curve.point_at(t);
            let d = point - *plane.origin;
            let (u, v) = (d.dot(plane.u_axis), d.dot(plane.v_axis));
            if u >= u_range.0 - eps
                && u <= u_range.1 + eps
                && v >= v_range.0 - eps
                && v <= v_range.1 + eps
            {
                out.push(SurfaceCurveIntersection::Point { point, t, u, v });
            }
            t += tau;
        }
    }
    Some(out)
}

#[cfg(test)]
#[path = "plane_circle_tests.rs"]
mod tests;
