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
use super::chord_affine::Affine;
use super::chord_bend::Bend;
use super::chord_waldron::waldron;
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

/// The corner budget of a bent bound, whose corner pieces the nets bound
/// to second order (PK-909): pentagon_tower's caps read 8.9e-3 at 4
/// pieces, 4.8e-3 at 64, the same at 1024. Bent bounds run only on the
/// triangles the others leave over the tolerance.
const BENT_CORNER_WORK: usize = 64;

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
    let Some(mut bend) = Bend::new(w, controls) else {
        return Ok(None);
    };
    // A spline's domain (where its position hull applies) clips the
    // image enclosures; not across a periodic chart's wrap.
    if let (Some(domain), None, None) = (source.domain(), breaks.wrap_u, breaks.wrap_v) {
        bend = bend.clipped_to(&controls, domain);
    }
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
    let mut work = BENT_CORNER_WORK;
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

/// Shift into a periodic spline domain (all vertices together).
pub(super) fn wrap<const N: usize>(
    t: [Uv; N],
    breaks: &C0Breaks,
    source: &SurfaceJetSource,
) -> [Uv; N] {
    let (su, sv) = wrap_shift(&t, breaks, source);
    t.map(|q| (q.0 + su, q.1 + sv))
}

/// The period multiples [`wrap`] adds to `t`.
fn wrap_shift(t: &[Uv], breaks: &C0Breaks, source: &SurfaceJetSource) -> Uv {
    if !matches!(source, SurfaceJetSource::Spline(_)) {
        return (0.0, 0.0);
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
    (shift(0, breaks.wrap_u), shift(1, breaks.wrap_v))
}

/// The surface parameters a piece's points pair with, as a box in the
/// spline's chart: the hull of `t` (affine correspondence) or of `φ(t)`
/// (bent), shifted into a periodic domain.
fn chart_box(source: &SurfaceJetSource, ctx: &Ctx<'_>, t: &[Uv]) -> (Iv, Iv) {
    shifted_chart_box(source, ctx, t).0
}

/// [`chart_box`] and the shift it applied.
fn shifted_chart_box(source: &SurfaceJetSource, ctx: &Ctx<'_>, t: &[Uv]) -> ((Iv, Iv), Uv) {
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
    let shift = wrap_shift(&[(u.lo, v.lo), (u.hi, v.hi)], ctx.breaks, source);
    let (lo, hi) = (
        (u.lo + shift.0, v.lo + shift.1),
        (u.hi + shift.0, v.hi + shift.1),
    );
    ((Iv::new(lo.0, hi.0), Iv::new(lo.1, hi.1)), shift)
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
    let mut own = taylor.unwrap_or(f64::INFINITY);
    if near || own > target {
        let ((u, v), shift) = shifted_chart_box(source, ctx, &t);
        if near {
            let hull = source.position_hull(u, v)?;
            own = own.min(
                t.iter()
                    .map(|q| point_to_box(ctx.affine.at(*q), hull))
                    .fold(0.0, f64::max),
            );
        }
        // A spline's nets bound S − L themselves, second order where the
        // position hull is first and unaffected by a weight zero
        // (PK-909); bent pieces through an affine map in the surface's
        // parameters (`Bend::difference`). An optional tightening: where
        // it declines, the others stand.
        let netted = match ctx.bend {
            None => ctx
                .affine
                .param_map(shift)
                .and_then(|l| source.difference_hull(u, v, &l).ok()),
            Some(bend) => bend.difference(source, &t, t.map(|q| ctx.affine.at(q)), (u, v), shift),
        };
        if let Some(d) = netted {
            own = own.min(d);
        }
    }
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
    Ok(waldron(&t, [a, b, c]) + disc)
}
