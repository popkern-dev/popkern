//! Certified chordal deviation (PK-505, design
//! `docs/architecture/certified-chordal-bound.md`): rigorous upper
//! bounds on how far an emitted chord or triangle lies from the curve
//! piece or surface patch it stands for.
//!
//! - A boundary segment `p0 p1` for the curve over `[t0, t1]`: the
//!   rational Bernstein hull (NURBS, positive weights) or a closed form
//!   (analytic kinds), plus the endpoints' distance from the curve.
//! - A triangle with parameter vertices `w_i` and emitted vertices
//!   `p_i`: `(M/2) R²` (Waldron) with `M` the certified second-order
//!   bound over the triangle and `R` its smallest enclosing circle,
//!   plus the vertices' distance from the surface; bisected to tighten.
//!
//! Every bound covers the exact geometry: interval enclosures for the
//! geometry, and float arithmetic on emitted coordinates inflated by
//! [`ROUNDING`].

use super::CertifyDecline;
use super::bernstein::{self, Coef};
use super::interval::Iv;
use super::interval_vec::IvVec3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, PortableMath};

#[path = "chord_bend.rs"]
mod chord_bend;
#[path = "chord_rectangle.rs"]
mod chord_rectangle;
#[path = "chord_triangle.rs"]
mod chord_triangle;
pub use chord_rectangle::{quad_pair_offset, rectangle_chord_bound};
pub use chord_triangle::{C0Breaks, bent_triangle_chord_bound, triangle_chord_bound};

/// Relative inflation of float-computed distances (a few ulps of the
/// coordinates, generously: emitted coordinates are combined by at
/// most a dozen operations).
pub const ROUNDING: f64 = 1e-12;

/// Bisection depth of a curve piece's hull (2^depth leaves at most).
const CURVE_DEPTH: usize = 10;

/// One flattened segment: its parameter interval and certified bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlatSegment {
    /// Start parameter.
    pub t0: f64,
    /// End parameter.
    pub t1: f64,
    /// Certified upper bound of the curve's distance from the chord
    /// `C(t0) C(t1)` (both ways).
    pub bound: f64,
}

/// Why a curve could not be flattened within the tolerance.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum FlattenError {
    /// The curve has no certified bound.
    #[error("no certified bound: {0}")]
    Decline(#[from] CertifyDecline),
    /// A segment shorter than the depth limit still missed the
    /// tolerance (its bound is the one reached).
    #[error("segment [{t0}, {t1}] stays at bound {bound} after the depth limit")]
    DepthLimit {
        /// Start parameter.
        t0: f64,
        /// End parameter.
        t1: f64,
        /// The best bound reached.
        bound: f64,
    },
}

/// The upper end of `x`, inflated for float rounding.
fn up(x: f64, scale: f64) -> f64 {
    x + ROUNDING * (x.abs() + scale)
}

fn dist_to_segment(q: Point3, a: Point3, b: Point3) -> f64 {
    let ab = b - a;
    let len2 = ab.dot(&ab);
    let t = if len2 > 0.0 {
        ((q - a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    q.distance_to(&(a + ab * t))
}

/// Largest distance from the box `c` to the segment `a b`: distance to
/// a convex set is convex, so it peaks at a corner.
fn box_to_segment(c: IvVec3, a: Point3, b: Point3) -> f64 {
    let mut worst = 0.0_f64;
    for k in 0..8 {
        let pick = |iv: Iv, bit: usize| if k & bit == 0 { iv.lo } else { iv.hi };
        let q = Point3::new(pick(c.x, 1), pick(c.y, 2), pick(c.z, 4));
        worst = worst.max(dist_to_segment(q, a, b));
    }
    worst
}

/// Distance upper bound between a point and an enclosure.
pub(crate) fn point_to_box(p: Point3, c: IvVec3) -> f64 {
    let d = c - IvVec3::from_point(p);
    d.norm().hi
}

/// A certified upper bound of the distance (both ways) between the
/// curve over `[t0, t1]` and the segment `p0 p1`, tightened by
/// bisection until it is at most `target` or the depth runs out:
/// [`curve_piece_bound`] plus the endpoints' distance from the curve
/// ([`curve_point_distance`]).
///
/// # Errors
///
/// [`CertifyDecline`] when the curve has no rigorous enclosure (a
/// NURBS weight enclosure reaching zero, non-finite data).
pub fn curve_chord_bound(
    curve: &CurveGeometry,
    t0: f64,
    t1: f64,
    p0: Point3,
    p1: Point3,
    target: f64,
) -> Result<f64, CertifyDecline> {
    let ends = curve_point_distance(curve, t0, p0)?.max(curve_point_distance(curve, t1, p1)?);
    let body = curve_piece_bound(curve, t0, t1, p0, p1, target)?;
    let scale = magnitude(p0).max(magnitude(p1));
    Ok(up(body + ends, scale))
}

/// Upper bound of the distance from `p` to the curve's point at `t`.
///
/// # Errors
///
/// As [`curve_chord_bound`].
pub fn curve_point_distance(
    curve: &CurveGeometry,
    t: f64,
    p: Point3,
) -> Result<f64, CertifyDecline> {
    Ok(point_to_box(p, curve_point(curve, t)?))
}

/// One way: an upper bound of the distance from every point of the
/// curve over `[t0, t1]` to the segment `p0 p1`. The other way (every
/// segment point near the curve) follows when the piece's ends are
/// near `p0` and `p1`: the piece's projection onto the segment covers
/// it, so add [`curve_point_distance`] of the ends.
///
/// # Errors
///
/// As [`curve_chord_bound`].
pub fn curve_piece_bound(
    curve: &CurveGeometry,
    t0: f64,
    t1: f64,
    p0: Point3,
    p1: Point3,
    target: f64,
) -> Result<f64, CertifyDecline> {
    let (a, b) = (t0.min(t1), t0.max(t1));
    let scale = magnitude(p0).max(magnitude(p1));
    let body = match curve {
        CurveGeometry::Nurbs(n) => nurbs_hull_bound(n, a, b, p0, p1, target)?,
        _ => analytic_piece(curve, a, b, p0, p1, target, CURVE_DEPTH)?,
    };
    Ok(up(body, scale))
}

/// Deviation of an analytic piece over `[a, b]` from its own chord
/// `C(a) C(b)`: the sagitta for circles, `h²/8 · sup |C''|` for the
/// other conics.
fn own_chord(curve: &CurveGeometry, a: f64, b: f64) -> f64 {
    match curve {
        CurveGeometry::Circle { radius, .. } | CurveGeometry::Arc { radius, .. } => {
            let half = Iv::point(b - a).scale(0.5);
            (Iv::point(*radius) * (Iv::point(1.0) - Iv::cos(half))).hi
        }
        CurveGeometry::Ellipse {
            u_basis, v_basis, ..
        } => {
            // |C''| = |cos t u + sin t v| ≤ sqrt(|u|² + |v|²).
            let m = (IvVec3::from_vec(*u_basis).norm_sqr() + IvVec3::from_vec(*v_basis).norm_sqr())
                .sqrt()
                .hi;
            second_derivative_bound(m, b - a)
        }
        CurveGeometry::Parabola { quadratic, .. } => {
            second_derivative_bound(2.0 * IvVec3::from_vec(*quadratic).norm().hi, b - a)
        }
        CurveGeometry::Hyperbola {
            semi_major,
            semi_minor,
            ..
        } => {
            // |C''| = |cosh t a + sinh t b| ≤ cosh m |a| + sinh m |b|.
            let m = a.abs().max(b.abs());
            let (ch, sh) = (up(m.pcosh(), 1.0), up(m.psinh(), 1.0));
            let bound = ch * IvVec3::from_vec(*semi_major).norm().hi
                + sh * IvVec3::from_vec(*semi_minor).norm().hi;
            second_derivative_bound(bound, b - a)
        }
        // A line piece is its own chord.
        CurveGeometry::Line { .. } | CurveGeometry::Nurbs(_) => 0.0,
    }
}

/// An analytic piece over `[a, b]` against the segment: its own-chord
/// deviation plus its chord ends' distance from the segment, bisected
/// while above `target` (each half is within its own bound of the
/// segment, so the worse half bounds the piece).
fn analytic_piece(
    curve: &CurveGeometry,
    a: f64,
    b: f64,
    p0: Point3,
    p1: Point3,
    target: f64,
    depth: usize,
) -> Result<f64, CertifyDecline> {
    let end = |t: f64| -> Result<f64, CertifyDecline> {
        Ok(box_to_segment(curve_point(curve, t)?, p0, p1))
    };
    let own = own_chord(curve, a, b) + end(a)?.max(end(b)?);
    if own <= target || depth == 0 || b - a <= f64::EPSILON * (1.0 + a.abs()) {
        return Ok(own);
    }
    let m = 0.5 * (a + b);
    let left = analytic_piece(curve, a, m, p0, p1, target, depth - 1)?;
    if left >= own {
        return Ok(own);
    }
    let right = analytic_piece(curve, m, b, p0, p1, target, depth - 1)?;
    Ok(own.min(left.max(right)))
}

/// The 1-D interpolation bound `h²/8 · sup |C''|`.
fn second_derivative_bound(m: f64, h: f64) -> f64 {
    (Iv::point(h).sqr() * Iv::point(m) / Iv::point(8.0)).hi
}

/// An enclosure of the curve's point at `t`.
fn curve_point(curve: &CurveGeometry, t: f64) -> Result<IvVec3, CertifyDecline> {
    match curve {
        CurveGeometry::Nurbs(n) => {
            let (span, net) = nurbs_span_net(n, t)?;
            let (lo, hi) = (n.knots[span], n.knots[span + 1]);
            project(bernstein::de_casteljau_point(
                &net,
                bernstein::local_param(t, lo, hi),
            ))
        }
        other => {
            // Analytic evaluation in floats: a few ulps, enclosed by a
            // relative pad (the bound adds this distance, so a pad
            // covers it rigorously at that scale).
            let p = crate::geometry::curve::CurveEval::point_at(other, t);
            let pad = ROUNDING * (1.0 + magnitude(p));
            let iv = |x: f64| Iv::new(x - pad, x + pad);
            let out = IvVec3::new(iv(p.x), iv(p.y), iv(p.z));
            if out.is_finite() {
                Ok(out)
            } else {
                Err(CertifyDecline::NonFinite)
            }
        }
    }
}

fn homogeneous(n: &NurbsCurve, k: usize) -> Result<Coef<4>, CertifyDecline> {
    let (p, w) = (n.control_points[k], n.weights[k]);
    if !(w > 0.0 && w.is_finite()) {
        return Err(CertifyDecline::NonPositiveWeight);
    }
    let w = Iv::point(w);
    Ok([
        Iv::point(p.x) * w,
        Iv::point(p.y) * w,
        Iv::point(p.z) * w,
        w,
    ])
}

/// The knot span holding `t` (the last non-empty one at the end) and
/// its Bernstein net.
fn nurbs_span_net(n: &NurbsCurve, t: f64) -> Result<(usize, Vec<Coef<4>>), CertifyDecline> {
    let spans = nurbs_spans(n)?;
    let span = spans
        .iter()
        .copied()
        .find(|&s| t < n.knots[s + 1])
        .or_else(|| spans.last().copied())
        .ok_or(CertifyDecline::MalformedNurbs)?;
    Ok((span, span_net(n, span)?))
}

fn nurbs_spans(n: &NurbsCurve) -> Result<Vec<usize>, CertifyDecline> {
    let p = n.degree;
    if n.weights.len() != n.control_points.len() || n.knots.len() != n.control_points.len() + p + 1
    {
        return Err(CertifyDecline::MalformedNurbs);
    }
    Ok((p..n.control_points.len())
        .filter(|&s| n.knots[s] < n.knots[s + 1])
        .collect())
}

fn span_net(n: &NurbsCurve, span: usize) -> Result<Vec<Coef<4>>, CertifyDecline> {
    let p = n.degree;
    let ctrl: Vec<Coef<4>> = (span - p..=span)
        .map(|k| homogeneous(n, k))
        .collect::<Result<_, _>>()?;
    Ok(bernstein::bspline_span_to_bernstein(
        p, &n.knots, span, &ctrl,
    ))
}

/// The projected point of a homogeneous coefficient.
fn project(c: Coef<4>) -> Result<IvVec3, CertifyDecline> {
    if c[3].lo <= 0.0 {
        return Err(CertifyDecline::DenominatorContainsZero);
    }
    let out = IvVec3::new(c[0] / c[3], c[1] / c[3], c[2] / c[3]);
    if out.is_finite() {
        Ok(out)
    } else {
        Err(CertifyDecline::NonFinite)
    }
}

/// The hull bound of the NURBS curve over `[a, b]` against `p0 p1`:
/// per knot span, the projected control points of the restricted net
/// (positive weights: the curve lies in their convex hull), bisected
/// while the bound exceeds `target`.
fn nurbs_hull_bound(
    n: &NurbsCurve,
    a: f64,
    b: f64,
    p0: Point3,
    p1: Point3,
    target: f64,
) -> Result<f64, CertifyDecline> {
    let mut worst = 0.0_f64;
    for span in nurbs_spans(n)? {
        let (lo, hi) = (n.knots[span], n.knots[span + 1]);
        let (x0, x1) = (a.max(lo), b.min(hi));
        if x0 >= x1 && !(a == b && (lo..=hi).contains(&a)) {
            continue;
        }
        let net = span_net(n, span)?;
        worst = worst.max(piece_bound(
            &net,
            (lo, hi),
            x0,
            x1,
            p0,
            p1,
            target,
            CURVE_DEPTH,
        )?);
    }
    Ok(worst)
}

/// The hull of the net restricted to `[x0, x1]` against the segment,
/// bisected while above `target` (the smaller of a parent's bound and
/// its children's worst is kept: both are bounds).
#[allow(clippy::too_many_arguments)]
fn piece_bound(
    net: &[Coef<4>],
    span: (f64, f64),
    x0: f64,
    x1: f64,
    p0: Point3,
    p1: Point3,
    target: f64,
    depth: usize,
) -> Result<f64, CertifyDecline> {
    let local = |x: f64| bernstein::local_param(x, span.0, span.1);
    let restricted = bernstein::restrict(net, local(x0), local(x1));
    let mut own = 0.0_f64;
    for c in restricted {
        own = own.max(box_to_segment(project(c)?, p0, p1));
    }
    if own <= target || depth == 0 || x1 - x0 <= f64::EPSILON * (1.0 + x0.abs()) {
        return Ok(own);
    }
    let m = 0.5 * (x0 + x1);
    let left = piece_bound(net, span, x0, m, p0, p1, target, depth - 1)?;
    if left >= own {
        return Ok(own);
    }
    let right = piece_bound(net, span, m, x1, p0, p1, target, depth - 1)?;
    Ok(own.min(left.max(right)))
}

/// Flatten the curve over `[t0, t1]` into segments whose certified
/// bound is at most `tol`, bisecting the parameter interval (at most
/// `2^depth` segments).
///
/// # Errors
///
/// [`FlattenError::Decline`] without a certified bound;
/// [`FlattenError::DepthLimit`] when a segment at the depth limit
/// still misses `tol`.
pub fn certified_flatten(
    curve: &CurveGeometry,
    t0: f64,
    t1: f64,
    tol: f64,
    depth: usize,
) -> Result<Vec<FlatSegment>, FlattenError> {
    let point = |t: f64| -> Result<Point3, CertifyDecline> {
        let c = curve_point(curve, t)?;
        Ok(Point3::new(c.x.mid(), c.y.mid(), c.z.mid()))
    };
    let mut out = Vec::new();
    let mut stack = vec![(t0, t1, depth)];
    while let Some((a, b, d)) = stack.pop() {
        let bound = curve_chord_bound(curve, a, b, point(a)?, point(b)?, tol)?;
        if bound <= tol {
            out.push(FlatSegment {
                t0: a,
                t1: b,
                bound,
            });
        } else if d == 0 {
            return Err(FlattenError::DepthLimit {
                t0: a,
                t1: b,
                bound,
            });
        } else {
            let m = 0.5 * (a + b);
            stack.push((m, b, d - 1));
            stack.push((a, m, d - 1));
        }
    }
    Ok(out)
}

/// The distance of `p` from the origin (the scale of its rounding).
pub(crate) fn magnitude(p: Point3) -> f64 {
    p.x.abs().max(p.y.abs()).max(p.z.abs())
}

#[cfg(test)]
#[path = "chord_tests.rs"]
mod tests;
