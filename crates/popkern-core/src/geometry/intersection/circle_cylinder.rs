//! Exact circle (or arc) × cylinder intersection with a tolerance-band
//! tangency (PK-887).
//!
//! The general NURBS solver read a circle lying in a plane tangent to
//! a cylinder, and dipping past it by less than the linear tolerance,
//! as two secant crossings 2·sqrt(2·r·depth) apart, while the plane ×
//! cylinder solver reads that plane as tangent (no curve) and the line
//! and plane × circle solvers read their contacts as one touch at the
//! extreme. The boolean then split one edge at two nodes 1.1e-2 apart
//! that nothing joined (measured: a web's ring arc against a boss 5e-7
//! proud of the web face).
//!
//! Along an angle-parameterised circle `P(t) = C + a cos t + b sin t`
//! the squared distance to the axis is a degree-2 trigonometric
//! polynomial, so its derivative has at most four zeros per turn. The
//! window is split at those critical points into monotone pieces, each
//! holding at most one root; an excursion across the wall between two
//! roots shallower than the tolerance is one contact at its extreme,
//! and an approach from inside within the tolerance touches once. An
//! approach from outside misses, as line × cylinder decides it.

use popkern_math::PortableMath;

use super::SurfaceCurveIntersection;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Derivative samples per turn when bracketing the critical points.
const SAMPLES_PER_TURN: f64 = 256.0;
/// Bisection steps on a bracket (2^-60 of a sample spacing).
const BISECTIONS: usize = 60;

/// The cylinder, by origin, axis and radius.
pub(crate) struct CylinderFrame<'a> {
    pub origin: &'a Point3,
    pub axis: &'a Vec3,
    pub radius: f64,
}

/// Roots of `curve` (a circle or arc, angle-parameterised) on the
/// cylinder within `t_range` and the cylinder's `(u, v)` window;
/// `None` when `curve` is not circular or the frame is degenerate.
pub(crate) fn intersect_circle_cylinder(
    cylinder: CylinderFrame<'_>,
    curve: &CurveGeometry,
    t_range: (f64, f64),
    (u_range, v_range): ((f64, f64), (f64, f64)),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceCurveIntersection>> {
    let center = match curve {
        CurveGeometry::Circle { center, .. } | CurveGeometry::Arc { center, .. } => *center,
        _ => return None,
    };
    let axis = cylinder.axis.normalized()?;
    let perp = |w: Vec3| w - axis * w.dot(&axis);
    let c = perp(center - *cylinder.origin);
    let a = perp(curve.point_at(0.0) - center);
    let b = perp(curve.point_at(std::f64::consts::FRAC_PI_2) - center);
    let (lo, hi) = (t_range.0.min(t_range.1), t_range.0.max(t_range.1));
    if !(lo.is_finite() && hi.is_finite()) {
        return None;
    }
    let q = |t: f64| c + a * t.pcos() + b * t.psin();
    let d = |t: f64| q(t).length() - cylinder.radius;
    // A full turn has no ends: start it where the circle is farthest
    // from the wall, so an in-band dip is never split by the window
    // seam into two one-sided pieces. Roots map back into the window.
    let tau = std::f64::consts::TAU;
    let full_turn = hi - lo >= tau * (1.0 - 1e-12);
    let (t0, t1) = if full_turn {
        let far = (0..64)
            .map(|k| lo + tau * k as f64 / 64.0)
            .max_by(|x, y| d(*x).abs().total_cmp(&d(*y).abs()))
            .unwrap_or(lo);
        (far, far + tau)
    } else {
        (lo, hi)
    };
    let home = |t: f64| {
        if full_turn {
            lo + (t - lo).rem_euclid(tau)
        } else {
            t
        }
    };
    // d/dt of |q|²/2: zero exactly where the distance is extremal.
    let slope = |t: f64| q(t).dot(&(b * t.pcos() - a * t.psin()));

    // Breakpoints: the window ends and every critical point inside.
    let steps = ((t1 - t0) / std::f64::consts::TAU * SAMPLES_PER_TURN)
        .ceil()
        .max(8.0) as usize;
    let at = |k: usize| t0 + (t1 - t0) * k as f64 / steps as f64;
    let mut breaks = vec![t0];
    for k in 0..steps {
        let (lo, hi) = (at(k), at(k + 1));
        if let Some(t) = bracket(&slope, lo, hi) {
            breaks.push(t);
        }
    }
    breaks.push(t1);

    let lin = tolerance.linear();
    // A circle on the wall (a coaxial rim) is an overlap: its critical
    // points and roots are roundoff, so the general solver owns it.
    if (0..=steps).all(|k| d(at(k)).abs() < lin) {
        return None;
    }
    let scale = super::coordinate_scale(&[center, *cylinder.origin]) + a.length() + b.length();
    let band = tolerance.tangency_band(scale + cylinder.radius);
    let noise = lin - band;
    // One root per monotone piece (None when the piece keeps its side).
    let roots: Vec<Option<f64>> = breaks.windows(2).map(|w| bracket(&d, w[0], w[1])).collect();
    let mut keep: Vec<f64> = Vec::new();
    let mut taken = vec![false; roots.len()];
    for i in 1..breaks.len() - 1 {
        let e = breaks[i];
        let depth = d(e);
        let shallow = depth.abs() < band;
        let crossed = roots[i - 1].is_some() && roots[i].is_some();
        if shallow && crossed && !taken[i - 1] && !taken[i] {
            // An excursion across the wall shallower than the tolerance.
            keep.push(e);
            taken[i - 1] = true;
            taken[i] = true;
        } else if roots[i - 1].is_none()
            && roots[i].is_none()
            && depth <= noise
            && depth > -band
            && depth >= d(breaks[i - 1])
            && depth >= d(breaks[i + 1])
        {
            // Approaching from inside within the tolerance: one touch.
            keep.push(e);
        }
    }
    keep.extend(
        roots
            .iter()
            .zip(&taken)
            .filter(|(_, taken)| !**taken)
            .filter_map(|(r, _)| *r),
    );
    keep.sort_by(f64::total_cmp);

    let (perp1, perp2) = orthonormal_basis(axis);
    let eps = tolerance.parametric_base();
    let mut out = Vec::new();
    for t in keep.into_iter().map(home) {
        let point = curve.point_at(t);
        let rel = point - *cylinder.origin;
        let v = rel.dot(&axis);
        let mut u = rel.dot(&perp2).patan2(rel.dot(&perp1));
        u += ((u_range.0 - u) / tau).ceil() * tau;
        if u - tau >= u_range.0 - eps {
            u -= tau;
        }
        let seen = out.iter().any(|h| match h {
            SurfaceCurveIntersection::Point { t: s, .. } => (s - t).abs() <= eps,
            SurfaceCurveIntersection::Overlap { .. } => false,
        });
        if seen || u > u_range.1 + tolerance.angular() || v < v_range.0 - lin || v > v_range.1 + lin
        {
            continue;
        }
        out.push(SurfaceCurveIntersection::Point { point, t, u, v });
    }
    Some(out)
}

/// The zero of `f` in `[lo, hi]` by bisection, when `f` changes sign
/// there (a zero exactly at `hi` belongs to the next bracket).
fn bracket(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64) -> Option<f64> {
    let (mut lo, mut hi) = (lo, hi);
    let (mut flo, fhi) = (f(lo), f(hi));
    if flo == 0.0 {
        return Some(lo);
    }
    if flo.signum() == fhi.signum() || fhi == 0.0 {
        return None;
    }
    for _ in 0..BISECTIONS {
        let mid = 0.5 * (lo + hi);
        let fm = f(mid);
        if fm == 0.0 {
            return Some(mid);
        }
        if fm.signum() == flo.signum() {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
#[path = "circle_cylinder_tests.rs"]
mod tests;
