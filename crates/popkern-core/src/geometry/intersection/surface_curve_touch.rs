//! Tangential contact of a curve with a surface: one root, at the touch.
//!
//! Where a curve touches a surface (a double root of the distance), it
//! stays within the linear tolerance of the surface over a band of
//! width about `2·√(2·r·tol)` around the touch — 6e-4 mm for a circle
//! of radius 40 against a cylinder at the default tolerance. Newton
//! seeded from neighbouring grid cells converges anywhere in that band:
//! one contact arrived as two "crossings" far enough apart to survive
//! de-duplication, or as one crossing 3e-4 off the touch. Measured
//! (PK-732): the equal-radius tee's branch end circle touches the main
//! cylinder at the two singular points; the circle was split twice at
//! each and the walk stranded the intersection records ending between
//! the splits.
//!
//! Every root moves to the curve's closest approach to the surface
//! within its band — the stretch around it where the curve stays within
//! the tolerance. A transversal crossing's band is `tol / sin θ` wide
//! and its closest approach is the crossing itself, so it stays put; a
//! touch's roots all land on the touch and collapse to one. On a closed
//! curve the band runs through the parameter seam: a touch AT the seam
//! (a closed edge's own vertex) lands on the window's end.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv_near;
use crate::tolerance::Tolerance;

/// One root: curve parameter, surface (u, v), point.
pub(super) type Root = (f64, f64, f64, Point3);

/// The distance from the curve at `t` to the surface (projected near
/// `uv`), with the projection.
fn distance(
    curve: &CurveGeometry,
    surface: &SurfaceGeometry,
    t: f64,
    uv: (f64, f64),
) -> Option<(f64, (f64, f64), Point3)> {
    let p = curve.point_at(t);
    let (u, v) = project_point_to_uv_near(surface, &p, uv.0, uv.1)?;
    Some((surface.point_at(u, v).distance_to(&p), (u, v), p))
}

/// The curve's parameter window, closed or not.
#[derive(Clone, Copy)]
struct Window {
    lo: f64,
    hi: f64,
    closed: bool,
}

impl Window {
    /// `t` brought into the window (through the seam when closed).
    fn wrap(self, t: f64) -> f64 {
        if !self.closed {
            return t.clamp(self.lo, self.hi);
        }
        let span = self.hi - self.lo;
        if t > self.hi {
            t - span
        } else if t < self.lo {
            t + span
        } else {
            t
        }
    }
}

/// The edge of `root`'s band in direction `dir` (±1): the farthest
/// parameter, by doubling steps then bisection, where the curve still
/// lies within `band` of the surface. At most a quarter of the window.
fn band_edge(
    curve: &CurveGeometry,
    surface: &SurfaceGeometry,
    root: &Root,
    window: Window,
    band: f64,
    dir: f64,
) -> f64 {
    let uv = (root.1, root.2);
    let within =
        |t: f64| distance(curve, surface, window.wrap(t), uv).is_some_and(|(d, _, _)| d <= band);
    let span = window.hi - window.lo;
    let reach = if window.closed {
        0.25 * span
    } else if dir > 0.0 {
        window.hi - root.0
    } else {
        root.0 - window.lo
    };
    let mut inside = 0.0;
    let mut step = span * 1e-9;
    let mut outside = None;
    while step <= reach {
        if within(root.0 + dir * step) {
            inside = step;
            step *= 2.0;
        } else {
            outside = Some(step);
            break;
        }
    }
    let Some(mut out) = outside else {
        return root.0 + dir * inside;
    };
    for _ in 0..40 {
        let mid = 0.5 * (inside + out);
        if within(root.0 + dir * mid) {
            inside = mid;
        } else {
            out = mid;
        }
    }
    root.0 + dir * inside
}

/// The closest approach of the curve to the surface over `[a, b]`
/// (golden section on the distance; `a`, `b` may run past the window
/// on a closed curve). The root itself when nothing is closer.
fn closest_approach(
    curve: &CurveGeometry,
    surface: &SurfaceGeometry,
    root: &Root,
    window: Window,
    (a, b): (f64, f64),
) -> Root {
    let uv = (root.1, root.2);
    let d = |t: f64| distance(curve, surface, window.wrap(t), uv).map_or(f64::INFINITY, |r| r.0);
    let g = 0.5 * (5.0_f64.sqrt() - 1.0);
    let (mut lo, mut hi) = (a, b);
    let mut x1 = hi - g * (hi - lo);
    let mut x2 = lo + g * (hi - lo);
    let (mut f1, mut f2) = (d(x1), d(x2));
    for _ in 0..80 {
        if f1 <= f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - g * (hi - lo);
            f1 = d(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + g * (hi - lo);
            f2 = d(x2);
        }
    }
    let t = window.wrap(0.5 * (lo + hi));
    match distance(curve, surface, t, uv) {
        Some((gap, (u, v), p)) if gap <= d(root.0) => (t, u, v, p),
        _ => *root,
    }
}

/// Move each root to its band's closest approach and collapse roots
/// that land on one contact. `t_range` is the curve's window; it is
/// closed when its ends meet.
pub(super) fn consolidate_touches(
    curve: &CurveGeometry,
    surface: &SurfaceGeometry,
    roots: Vec<Root>,
    t_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<Root> {
    let band = tolerance.linear();
    let window = Window {
        lo: t_range.0,
        hi: t_range.1,
        closed: curve
            .point_at(t_range.0)
            .distance_to(&curve.point_at(t_range.1))
            <= band,
    };
    let mut out: Vec<Root> = Vec::with_capacity(roots.len());
    for root in roots {
        let left = band_edge(curve, surface, &root, window, band, -1.0);
        let right = band_edge(curve, surface, &root, window, band, 1.0);
        let touch = closest_approach(curve, surface, &root, window, (left, right));
        // One contact: a root landing on one already kept (in space —
        // a seam touch lands at either window end).
        if !out
            .iter()
            .any(|kept| kept.3.distance_to(&touch.3) <= 10.0 * band)
        {
            out.push(touch);
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

#[cfg(test)]
#[path = "surface_curve_touch_tests.rs"]
mod tests;
