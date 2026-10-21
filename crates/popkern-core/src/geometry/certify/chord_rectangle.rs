//! The parameter-rectangle half of the certified chordal bound (see
//! `chord`): an emitted quad, or a pole fan's triangle, whose parameter
//! domain is a rectangle mapped bilinearly onto its corners.
//!
//! With `I` the bilinear map onto the emitted corners and `I_R S` the
//! bilinear interpolant of `S` at a sub-rectangle `R`'s corners,
//! `|I − S| ≤ |I − I_R S| + |I_R S − S|`. The first is a bilinear
//! interpolant of the corner discrepancies (at most their largest);
//! the second, as `I_R = I_u I_v` with each a convex combination, is at
//! most `h_u²/8 · sup |S_uu| + h_v²/8 · sup |S_vv|` over `R`.

use super::super::CertifyDecline;
use super::super::interval::Iv;
use super::super::surface_jet::SurfaceJetSource;
use super::chord_triangle::{CORNER_DEPTH, CORNER_WORK, Uv, wrap};
use super::{C0Breaks, ROUNDING, magnitude, point_to_box};
use crate::geometry::Point3;

/// Bisection depth of a rectangle (4^depth sub-rectangles at most).
const RECTANGLE_DEPTH: usize = 3;

/// A certified upper bound of the distance between the surface over
/// `u × v` and its bilinear image onto `corners` (`[p(u0, v0),
/// p(u1, v0), p(u0, v1), p(u1, v1)]`), tightened toward `target`.
///
/// # Errors
///
/// [`CertifyDecline`] when the surface has no enclosure over a piece.
pub fn rectangle_chord_bound(
    source: &mut SurfaceJetSource,
    breaks: &C0Breaks,
    u: (f64, f64),
    v: (f64, f64),
    corners: [Point3; 4],
    target: f64,
) -> Result<f64, CertifyDecline> {
    let map = Bilinear { u, v, corners };
    let scale = corners.iter().map(|q| magnitude(*q)).fold(0.0, f64::max);
    let cuts = |lo: f64, hi: f64, lines: &[f64]| -> Vec<f64> {
        let mut xs = vec![lo];
        let mut inner: Vec<f64> = lines
            .iter()
            .copied()
            .filter(|&x| x > lo && x < hi)
            .collect();
        inner.sort_by(f64::total_cmp);
        xs.extend(inner);
        xs.push(hi);
        xs
    };
    let us = cuts(u.0.min(u.1), u.0.max(u.1), &breaks.u);
    let vs = cuts(v.0.min(v.1), v.0.max(v.1), &breaks.v);
    let mut worst = 0.0_f64;
    let mut work = CORNER_WORK;
    for uu in us.windows(2) {
        for vv in vs.windows(2) {
            let piece = ((uu[0], uu[1]), (vv[0], vv[1]));
            worst = worst.max(piece_bound(
                source,
                breaks,
                &map,
                piece,
                target,
                (RECTANGLE_DEPTH, CORNER_DEPTH, false),
                &mut work,
            )?);
        }
    }
    Ok(worst + ROUNDING * (worst + scale))
}

struct Bilinear {
    u: (f64, f64),
    v: (f64, f64),
    corners: [Point3; 4],
}

impl Bilinear {
    fn at(&self, q: Uv) -> Point3 {
        let r = if self.u.1 != self.u.0 {
            (q.0 - self.u.0) / (self.u.1 - self.u.0)
        } else {
            0.0
        };
        let s = if self.v.1 != self.v.0 {
            (q.1 - self.v.0) / (self.v.1 - self.v.0)
        } else {
            0.0
        };
        let c = &self.corners;
        let w = [(1.0 - r) * (1.0 - s), r * (1.0 - s), (1.0 - r) * s, r * s];
        Point3::new(
            w[0] * c[0].x + w[1] * c[1].x + w[2] * c[2].x + w[3] * c[3].x,
            w[0] * c[0].y + w[1] * c[1].y + w[2] * c[2].y + w[3] * c[3].y,
            w[0] * c[0].z + w[1] * c[1].z + w[2] * c[2].z + w[3] * c[3].z,
        )
    }
}

fn piece_bound(
    source: &mut SurfaceJetSource,
    breaks: &C0Breaks,
    map: &Bilinear,
    (u, v): ((f64, f64), (f64, f64)),
    target: f64,
    (depth, corner, near): (usize, usize, bool),
    work: &mut usize,
) -> Result<f64, CertifyDecline> {
    let corners = [(u.0, v.0), (u.1, v.0), (u.0, v.1), (u.1, v.1)];
    let taylor = match own_bound(source, breaks, map, (u, v)) {
        Ok(own) => Some(own),
        Err(CertifyDecline::DenominatorContainsZero) => None,
        Err(e) => return Err(e),
    };
    // Near a zero of the weight function (a Gregory corner): the
    // position hull bounds the piece too (see `chord_triangle`).
    let near = near || taylor.is_none();
    let own = if near {
        let chart = wrap(corners, breaks, source);
        let (du, dv) = (chart[0].0 - corners[0].0, chart[0].1 - corners[0].1);
        let hull =
            source.position_hull(Iv::new(u.0 + du, u.1 + du), Iv::new(v.0 + dv, v.1 + dv))?;
        let own = corners
            .iter()
            .map(|q| point_to_box(map.at(*q), hull))
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
    let (um, vm) = (0.5 * (u.0 + u.1), 0.5 * (v.0 + v.1));
    let next = (depth.saturating_sub(1), corner.saturating_sub(1), near);
    let mut children = 0.0_f64;
    for sub in [
        ((u.0, um), (v.0, vm)),
        ((um, u.1), (v.0, vm)),
        ((u.0, um), (vm, v.1)),
        ((um, u.1), (vm, v.1)),
    ] {
        children = children.max(piece_bound(source, breaks, map, sub, target, next, work)?);
        if children >= own {
            return Ok(own);
        }
    }
    Ok(children)
}

/// The piece's own bilinear bound plus its corner discrepancies.
fn own_bound(
    source: &mut SurfaceJetSource,
    breaks: &C0Breaks,
    map: &Bilinear,
    (u, v): ((f64, f64), (f64, f64)),
) -> Result<f64, CertifyDecline> {
    let corners = [(u.0, v.0), (u.1, v.0), (u.0, v.1), (u.1, v.1)];
    let chart = wrap(corners, breaks, source);
    let mut disc = 0.0_f64;
    for (q, c) in corners.iter().zip(&chart) {
        let jet = source
            .jet_piecewise_to(Iv::point(c.0), Iv::point(c.1), 0)?
            .0;
        disc = disc.max(point_to_box(map.at(*q), jet.partial(0, 0)));
    }
    let (du, dv) = (chart[0].0 - corners[0].0, chart[0].1 - corners[0].1);
    let (jet, _) =
        source.jet_piecewise_to(Iv::new(u.0 + du, u.1 + du), Iv::new(v.0 + dv, v.1 + dv), 2)?;
    let (a, c) = (jet.partial(2, 0).norm().hi, jet.partial(0, 2).norm().hi);
    let eighth = |h: f64, m: f64| (Iv::point(h).sqr() * Iv::point(m) / Iv::point(8.0)).hi;
    Ok((eighth(u.1 - u.0, a) + eighth(v.1 - v.0, c)) * (1.0 + 1e-12) + disc)
}

/// How far an emitted quad, split into two triangles along a diagonal,
/// is from the bilinear patch through its corners (`[p(u0, v0),
/// p(u1, v0), p(u0, v1), p(u1, v1)]`, split along `p(u0, v0) p(u1, v1)`
/// when `main_diagonal`, else along the other).
///
/// Projected along the normal `n` of either triangle, the quad must be
/// strictly convex and the other triangle must face the same way.
/// Then the bilinear map is a bijection onto the projected quad, and
/// over each projected point both the patch and the triangle pair sit
/// at heights between the corners' (`0` and the off corner's `h`).
/// Pairing points that share a projection maps the rectangle onto the
/// triangles, and no pair is farther apart than `|h|`. This avoids the
/// in-plane shear an affine correspondence charges: a planar trapezoid
/// (two rulings of a cone) is offset 0. `None` when neither triangle
/// gives a valid projection.
pub fn quad_pair_offset(corners: [Point3; 4], main_diagonal: bool) -> Option<f64> {
    let q = [corners[0], corners[1], corners[3], corners[2]];
    // Each triangle of the split and the quad corner it leaves out, as
    // indices into the cyclic `q`.
    let splits: [([usize; 3], usize); 2] = if main_diagonal {
        [([0, 1, 2], 3), ([0, 2, 3], 1)]
    } else {
        [([1, 2, 3], 0), ([1, 3, 0], 2)]
    };
    let scale = q.iter().map(|p| magnitude(*p)).fold(0.0, f64::max);
    let normal = |t: [usize; 3]| (q[t[1]] - q[t[0]]).cross(&(q[t[2]] - q[t[0]]));
    let mut best: Option<f64> = None;
    for (k, (tri, off)) in splits.iter().enumerate() {
        let Some(n) = normal(*tri).normalized() else {
            continue;
        };
        let (other, _) = splits[1 - k];
        if normal(other).dot(&n) <= 0.0 {
            continue;
        }
        let turns: Vec<f64> = (0..4)
            .map(|i| {
                let (a, b, c) = (q[i], q[(i + 1) % 4], q[(i + 2) % 4]);
                (b - a).cross(&(c - b)).dot(&n)
            })
            .collect();
        let edge = (0..4)
            .map(|i| (q[(i + 1) % 4] - q[i]).length())
            .fold(0.0, f64::max);
        let margin = 1e-9 * edge * edge;
        if !(turns.iter().all(|t| *t > margin) || turns.iter().all(|t| *t < -margin)) {
            continue;
        }
        let h = (q[*off] - q[tri[0]]).dot(&n).abs() + ROUNDING * (edge + scale);
        best = Some(best.map_or(h, |b| b.min(h)));
    }
    best
}
