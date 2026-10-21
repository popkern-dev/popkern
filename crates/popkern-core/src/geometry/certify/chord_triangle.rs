//! The triangle half of the certified chordal bound (see `chord`).
//!
//! For parameter vertices `w_i` and the affine map `L` with
//! `L(w_i) = p_i`, every `x` in the triangle satisfies
//! `|L(x) − S(x)| ≤ (M/2) R² + max_i |p_i − S(w_i)|`, where `M`
//! bounds `|D²_d S|` over the triangle for unit directions `d` and `R`
//! is the radius of its smallest enclosing circle (Waldron). The same
//! form holds on any sub-triangle against the same `L`, with the
//! discrepancies `|L(q_j) − S(q_j)|` at its own vertices, so the
//! triangle is cut along C⁰ parameter lines (where Taylor's remainder
//! does not hold) and bisected into four while the bound exceeds the
//! target.

use super::super::CertifyDecline;
use super::super::interval::Iv;
use super::super::surface_jet::SurfaceJetSource;
use super::chord_bend::Bend;
use super::{ROUNDING, magnitude, point_to_box};
use crate::geometry::Point3;

/// Bisection depth of a triangle (4^depth sub-triangles at most).
const TRIANGLE_DEPTH: usize = 3;

/// Bisection depth of a bent triangle: bent bounds are taken only for
/// triangles the others leave over the tolerance, and land near it (the
/// double-boss fillet torus read 5.007e-3 at depth 3, against a 5e-3
/// chord height and a 4.82e-3 sampled pairing).
const BENT_DEPTH: usize = 4;

/// Bisection depth toward a zero of the weight function (a Gregory
/// corner), where only the position hull encloses: each level halves
/// the corner piece, whose bound shrinks with its size.
pub(super) const CORNER_DEPTH: usize = 20;

/// Pieces near weight zeros one bound may visit (each corner subtree
/// splits every sibling, so depth alone does not bound the work).
/// Measured on pentagon_tower's Gregory caps (PK-904): the bound is the
/// same from 4 to 4096 pieces (2.3e-2, set by the pieces beside the
/// activated corner, where the jets' quotient diverges as the box
/// shrinks), while the time grows from 0.1 s to 25 s per triangle;
/// at 2 pieces 13 of the 20 caps' bounds loosen (PK-905).
pub(super) const CORNER_WORK: usize = 4;

/// Where the surface is not C¹, and how its chart wraps: the triangle
/// is cut along these lines before Taylor bounds apply, and pieces
/// outside a spline's domain are shifted by its period.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct C0Breaks {
    /// `u` values of C⁰ (or worse) parameter lines.
    pub u: Vec<f64>,
    /// `v` values of C⁰ (or worse) parameter lines.
    pub v: Vec<f64>,
    /// `(start, period)` of a periodic `u` domain a piece is shifted
    /// into before its jets are taken (spline sources only).
    pub wrap_u: Option<(f64, f64)>,
    /// `(start, period)` of a periodic `v` domain.
    pub wrap_v: Option<(f64, f64)>,
}

pub(super) type Uv = (f64, f64);

/// A certified upper bound of the distance between the triangle
/// `p0 p1 p2` and the surface patch over its parameter triangle `w`
/// (each point of either within the bound of the other through the
/// affine correspondence), tightened toward `target`.
///
/// # Errors
///
/// [`CertifyDecline`] when the surface has no enclosure over a piece
/// (non-finite jets, a spline box outside its knots).
pub fn triangle_chord_bound(
    source: &mut SurfaceJetSource,
    breaks: &C0Breaks,
    w: [Uv; 3],
    p: [Point3; 3],
    target: f64,
) -> Result<f64, CertifyDecline> {
    let affine = Affine::new(w, p);
    let mut worst = 0.0_f64;
    let mut work = CORNER_WORK;
    for piece in cut(vec![w], breaks) {
        let ctx = Ctx {
            breaks,
            affine: &affine,
            bend: None,
        };
        worst = worst.max(top_bound(source, &ctx, piece, target, &mut work)?);
    }
    Ok(worst + ROUNDING * (worst + scale_of(&p)))
}

/// The bound of [`triangle_chord_bound`] under the bent correspondence
/// `y ↦ (S(φ(y)), L(y))`, `φ` the quadratic map whose edge `k`
/// (vertices `k`, `k + 1`) has the control point `controls[k]` (`None`
/// for a straight edge; see `chord_bend`). `None` when
/// the bend does not apply: a degenerate triangle, or one a C⁰ line of
/// the surface crosses (in `y` or under `φ`).
///
/// # Errors
///
/// As [`triangle_chord_bound`].
pub fn bent_triangle_chord_bound(
    source: &mut SurfaceJetSource,
    breaks: &C0Breaks,
    w: [Uv; 3],
    p: [Point3; 3],
    controls: [Option<Uv>; 3],
    target: f64,
) -> Result<Option<f64>, CertifyDecline> {
    let Some(bend) = Bend::new(w, controls) else {
        return Ok(None);
    };
    let (bu, bv) = bend.map_box(&w);
    let crosses = |lines: &[f64], x: Iv| lines.iter().any(|&k| x.lo < k && k < x.hi);
    if cut(vec![w], breaks).len() != 1 || crosses(&breaks.u, bu) || crosses(&breaks.v, bv) {
        return Ok(None);
    }
    let affine = Affine::new(w, p);
    let ctx = Ctx {
        breaks,
        affine: &affine,
        bend: Some(&bend),
    };
    let mut work = CORNER_WORK;
    let worst = top_bound(source, &ctx, w, target, &mut work)?;
    Ok(Some(worst + ROUNDING * (worst + scale_of(&p))))
}

fn scale_of(p: &[Point3; 3]) -> f64 {
    p.iter().map(|q| magnitude(*q)).fold(0.0, f64::max)
}

/// What every piece of one triangle's bound shares.
struct Ctx<'a> {
    breaks: &'a C0Breaks,
    affine: &'a Affine,
    /// The bent map, or `None` for the affine correspondence.
    bend: Option<&'a Bend>,
}

/// A top-level piece: its vertex discrepancies, then its bound.
fn top_bound(
    source: &mut SurfaceJetSource,
    ctx: &Ctx<'_>,
    piece: [Uv; 3],
    target: f64,
    work: &mut usize,
) -> Result<f64, CertifyDecline> {
    let mut disc = [0.0; 3];
    for (d, q) in disc.iter_mut().zip(piece) {
        *d = discrepancy(source, ctx, q)?;
    }
    piece_bound(
        source,
        ctx,
        (piece, disc),
        target,
        (
            if ctx.bend.is_some() {
                BENT_DEPTH
            } else {
                TRIANGLE_DEPTH
            },
            CORNER_DEPTH,
            false,
        ),
        work,
    )
}

/// The affine map of the parameter triangle onto the emitted one.
struct Affine {
    w: [Uv; 3],
    p: [Point3; 3],
    det: f64,
}

impl Affine {
    fn new(w: [Uv; 3], p: [Point3; 3]) -> Self {
        let det = (w[1].0 - w[0].0) * (w[2].1 - w[0].1) - (w[2].0 - w[0].0) * (w[1].1 - w[0].1);
        Affine { w, p, det }
    }

    /// `L(q)`; a degenerate parameter triangle maps every point to the
    /// nearest vertex's image (its pieces are then single points or
    /// segments whose discrepancies carry the bound).
    fn at(&self, q: Uv) -> Point3 {
        let w = &self.w;
        if self.det.abs() <= f64::MIN_POSITIVE {
            let d = |a: Uv| (a.0 - q.0).powi(2) + (a.1 - q.1).powi(2);
            let k = (0..3)
                .min_by(|&i, &j| d(w[i]).total_cmp(&d(w[j])))
                .unwrap_or(0);
            return self.p[k];
        }
        let l1 =
            ((q.0 - w[0].0) * (w[2].1 - w[0].1) - (w[2].0 - w[0].0) * (q.1 - w[0].1)) / self.det;
        let l2 =
            ((w[1].0 - w[0].0) * (q.1 - w[0].1) - (q.0 - w[0].0) * (w[1].1 - w[0].1)) / self.det;
        let l0 = 1.0 - l1 - l2;
        let p = &self.p;
        Point3::new(
            l0 * p[0].x + l1 * p[1].x + l2 * p[2].x,
            l0 * p[0].y + l1 * p[1].y + l2 * p[2].y,
            l0 * p[0].z + l1 * p[1].z + l2 * p[2].z,
        )
    }
}

/// The triangles cut along every break line strictly inside them.
fn cut(mut pieces: Vec<[Uv; 3]>, breaks: &C0Breaks) -> Vec<[Uv; 3]> {
    let lines = breaks
        .u
        .iter()
        .map(|&x| (0, x))
        .chain(breaks.v.iter().map(|&x| (1, x)));
    for (axis, x) in lines {
        let mut next = Vec::with_capacity(pieces.len());
        for t in pieces {
            next.extend(split(t, axis, x));
        }
        pieces = next;
    }
    pieces
}

pub(super) fn coord(q: Uv, axis: usize) -> f64 {
    if axis == 0 { q.0 } else { q.1 }
}

/// Triangle `t` split by the line `coord(·, axis) = x` into triangles
/// on each side (itself when the line misses its interior).
pub(super) fn split(t: [Uv; 3], axis: usize, x: f64) -> Vec<[Uv; 3]> {
    let c: Vec<f64> = t.iter().map(|&q| coord(q, axis) - x).collect();
    if c.iter().all(|&d| d >= 0.0) || c.iter().all(|&d| d <= 0.0) {
        return vec![t];
    }
    let mut sides: [Vec<Uv>; 2] = [Vec::new(), Vec::new()];
    for i in 0..3 {
        let (a, b) = (t[i], t[(i + 1) % 3]);
        let (ca, cb) = (c[i], c[(i + 1) % 3]);
        if ca <= 0.0 {
            sides[0].push(a);
        }
        if ca >= 0.0 {
            sides[1].push(a);
        }
        if (ca < 0.0 && cb > 0.0) || (ca > 0.0 && cb < 0.0) {
            let s = ca / (ca - cb);
            let mut q = (a.0 + s * (b.0 - a.0), a.1 + s * (b.1 - a.1));
            // On the line exactly: the pieces then meet along it.
            if axis == 0 {
                q.0 = x;
            } else {
                q.1 = x;
            }
            sides[0].push(q);
            sides[1].push(q);
        }
    }
    let mut out = Vec::new();
    for poly in sides {
        for k in 1..poly.len().saturating_sub(1) {
            out.push([poly[0], poly[k], poly[k + 1]]);
        }
    }
    out
}

/// Squared radius of the smallest circle holding the triangle.
fn enclosing_radius_sqr(t: &[Uv; 3]) -> f64 {
    let d2 = |a: Uv, b: Uv| (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2);
    let sides = [d2(t[1], t[2]), d2(t[0], t[2]), d2(t[0], t[1])];
    let long = sides.iter().copied().fold(0.0, f64::max);
    let others: f64 = sides.iter().sum::<f64>() - long;
    if long >= others {
        // Right or obtuse: the longest side's midpoint.
        return long / 4.0;
    }
    // Acute: the circumradius, R² = a²b²c² / (16 · area²).
    let area2 = ((t[1].0 - t[0].0) * (t[2].1 - t[0].1) - (t[2].0 - t[0].0) * (t[1].1 - t[0].1))
        .powi(2)
        / 4.0;
    if area2 <= 0.0 {
        return long / 4.0;
    }
    (sides[0] * sides[1] * sides[2] / (16.0 * area2)).max(long / 4.0)
}

/// Shift into a periodic spline domain (all vertices together).
pub(super) fn wrap<const N: usize>(
    t: [Uv; N],
    breaks: &C0Breaks,
    source: &SurfaceJetSource,
) -> [Uv; N] {
    if !matches!(source, SurfaceJetSource::Spline(_)) {
        return t;
    }
    let shift = |axis: usize, w: Option<(f64, f64)>| -> f64 {
        let Some((start, period)) = w.filter(|w| w.1 > 0.0) else {
            return 0.0;
        };
        let lo = t
            .iter()
            .map(|&q| coord(q, axis))
            .fold(f64::INFINITY, f64::min);
        -((lo - start) / period).floor() * period
    };
    let (su, sv) = (shift(0, breaks.wrap_u), shift(1, breaks.wrap_v));
    t.map(|q| (q.0 + su, q.1 + sv))
}

/// The surface parameters a piece's points pair with, as a box in the
/// spline's chart: the hull of `t` (affine correspondence) or of `φ(t)`
/// (bent), shifted into a periodic domain.
fn chart_box(source: &SurfaceJetSource, ctx: &Ctx<'_>, t: &[Uv]) -> (Iv, Iv) {
    let (u, v) = match ctx.bend {
        Some(bend) => bend.map_box(t),
        None => {
            let span = |k: usize| {
                let xs = t.iter().map(|q| coord(*q, k));
                Iv::new(
                    xs.clone().fold(f64::INFINITY, f64::min),
                    xs.fold(f64::NEG_INFINITY, f64::max),
                )
            };
            (span(0), span(1))
        }
    };
    let [lo, hi] = wrap([(u.lo, v.lo), (u.hi, v.hi)], ctx.breaks, source);
    (Iv::new(lo.0, hi.0), Iv::new(lo.1, hi.1))
}

/// `|L(q) − S(φ(q))|` (`φ` the identity unbent), bounded above;
/// infinite at a zero of the weight function, where the quotient
/// declines (a piece with such a vertex has no Taylor bound and takes
/// the position hull).
fn discrepancy(source: &mut SurfaceJetSource, ctx: &Ctx<'_>, q: Uv) -> Result<f64, CertifyDecline> {
    let (u, v) = chart_box(source, ctx, &[q]);
    match source.jet_piecewise_to(u, v, 0) {
        Ok((jet, _)) => Ok(point_to_box(ctx.affine.at(q), jet.partial(0, 0))),
        Err(CertifyDecline::DenominatorContainsZero) => Ok(f64::INFINITY),
        Err(e) => Err(e),
    }
}

/// The bound over piece `t`, whose vertices' discrepancies are `disc`
/// (a piece's children share its corners, so each split encloses only
/// its three midpoints).
fn piece_bound(
    source: &mut SurfaceJetSource,
    ctx: &Ctx<'_>,
    (t, disc): ([Uv; 3], [f64; 3]),
    target: f64,
    (depth, corner, near): (usize, usize, bool),
    work: &mut usize,
) -> Result<f64, CertifyDecline> {
    let taylor = match own_bound(source, ctx, t, disc) {
        Ok(own) if own.is_finite() => Some(own),
        Ok(_) | Err(CertifyDecline::DenominatorContainsZero) => None,
        Err(e) => return Err(e),
    };
    // Near a zero of the weight function (a Gregory corner) the jets
    // decline, or enclose loosely beside it: the position hull bounds
    // those pieces too, and they recurse on the corner budget.
    let near = near || taylor.is_none();
    let own = if near {
        let (u, v) = chart_box(source, ctx, &t);
        let hull = source.position_hull(u, v)?;
        let own = t
            .iter()
            .map(|q| point_to_box(ctx.affine.at(*q), hull))
            .fold(0.0, f64::max);
        taylor.map_or(own, |x| x.min(own))
    } else {
        taylor.unwrap_or(f64::INFINITY)
    };
    let budget = if near {
        // Each piece near the corner spends one unit of the call's
        // work: out of it, the piece keeps its own (valid) bound.
        *work = work.saturating_sub(1);
        if *work == 0 { 0 } else { corner }
    } else {
        depth
    };
    if own <= target || budget == 0 {
        return Ok(own);
    }
    let mid = |i: usize, j: usize| (0.5 * (t[i].0 + t[j].0), 0.5 * (t[i].1 + t[j].1));
    let (m01, m12, m20) = (mid(0, 1), mid(1, 2), mid(2, 0));
    let (d01, d12, d20) = (
        discrepancy(source, ctx, m01)?,
        discrepancy(source, ctx, m12)?,
        discrepancy(source, ctx, m20)?,
    );
    let next = (depth.saturating_sub(1), corner.saturating_sub(1), near);
    let mut children = 0.0_f64;
    for sub in [
        ([t[0], m01, m20], [disc[0], d01, d20]),
        ([m01, t[1], m12], [d01, disc[1], d12]),
        ([m20, m12, t[2]], [d20, d12, disc[2]]),
        ([m01, m12, m20], [d01, d12, d20]),
    ] {
        children = children.max(piece_bound(source, ctx, sub, target, next, work)?);
        if children >= own {
            return Ok(own);
        }
    }
    Ok(children)
}

/// The piece's own Taylor bound (Waldron's, smallest of three metric
/// forms) plus its largest vertex discrepancy.
fn own_bound(
    source: &mut SurfaceJetSource,
    ctx: &Ctx<'_>,
    t: [Uv; 3],
    disc: [f64; 3],
) -> Result<f64, CertifyDecline> {
    let disc = disc.into_iter().fold(0.0_f64, f64::max);
    let (u, v) = chart_box(source, ctx, &t);
    let (jet, _) = source.jet_piecewise_to(u, v, 2)?;
    // Second partials of S (unbent), or of S∘φ by the chain rule.
    let [a, b, c] = match ctx.bend {
        Some(bend) => bend.second_partials(&t, &jet),
        None => [(2, 0), (1, 1), (0, 2)].map(|(i, j)| jet.partial(i, j).norm().hi),
    };
    // Bounds of |D²_d S| = |d1² S_uu + 2 d1 d2 S_uv + d2² S_vv|, all
    // valid, the smallest kept: isotropic, λ_max([[a, b], [b, c]])
    // times |d|²; and diagonal, (a + εb) d1² + (c + b/ε) d2² for any
    // ε > 0 (as 2|d1 d2| ≤ ε d1² + d2²/ε), the metric of the
    // parameters scaled by sqrt(a + εb) and sqrt(c + b/ε). ε = 1, and
    // ε = the triangle's v/u extent ratio, which balances the cross
    // term on long strips (a cylinder's or a far cone's: the isotropic
    // form charges the axial length with the ring curvature).
    let (half_sum, half_diff) = (Iv::point(a + c).scale(0.5), Iv::point(a - c).scale(0.5));
    let m = (half_sum + (half_diff.sqr() + Iv::point(b).sqr()).sqrt()).hi;
    let isotropic = (Iv::point(m) * Iv::point(enclosing_radius_sqr(&t)))
        .scale(0.5)
        .hi;
    let extent = |k: usize| {
        let xs = t.iter().map(|q| coord(*q, k));
        xs.clone().fold(f64::NEG_INFINITY, f64::max) - xs.fold(f64::INFINITY, f64::min)
    };
    let ratio = (extent(1) / extent(0)).clamp(1e-6, 1e6);
    let diagonal = |eps: f64| {
        let (eps, b) = (Iv::point(eps), Iv::point(b));
        let su = (Iv::point(a) + eps * b).sqrt().hi;
        let sv = (Iv::point(c) + b / eps).sqrt().hi;
        let scaled = t.map(|q| (q.0 * su, q.1 * sv));
        Iv::point(enclosing_radius_sqr(&scaled)).scale(0.5).hi
    };
    let quad = if ratio.is_finite() {
        diagonal(ratio)
    } else {
        f64::INFINITY
    };
    Ok(isotropic.min(diagonal(1.0)).min(quad) * (1.0 + 1e-12) + disc)
}
