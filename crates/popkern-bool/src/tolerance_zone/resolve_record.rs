//! Distances along records and edges for the vertex table's
//! incompatibility rules, and the record-meeting exemptions (split from
//! `resolve.rs` under the 500-line rule).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;

use super::candidates::{CandidateVertex, Incidence};
use super::resolve::ResolutionContext;

/// Arc length of `curve` between `t0` and `t1`, by 32 chords.
pub(super) fn arc_between(point_at: impl Fn(f64) -> Point3, t0: f64, t1: f64) -> f64 {
    const CHORDS: usize = 32;
    let mut prev = point_at(t0);
    let mut sum = 0.0;
    for k in 1..=CHORDS {
        let p = point_at(t0 + (t1 - t0) * k as f64 / CHORDS as f64);
        sum += prev.distance_to(&p);
        prev = p;
    }
    sum
}

/// Arc length between two parameters of a curve over `(lo, hi)`; on a
/// closed curve the shorter way round.
pub(super) fn separation(
    point_at: impl Fn(f64) -> Point3 + Copy,
    (lo, hi): (f64, f64),
    t: (f64, f64),
    closing: f64,
) -> f64 {
    let (t0, t1) = (t.0.min(t.1), t.0.max(t.1));
    let inner = arc_between(point_at, t0, t1);
    // Closed when the ends meet within `closing`, the candidates'
    // radii: a closed record rotated to start at a vertex is refitted,
    // and its ends then differ by the fit (1e-8 on a midsphere circle),
    // which read as open put its two ends a whole circle apart.
    // Within the radii it closes unless it doubles back on itself: a
    // record curling back to within a tolerance of its own start (its
    // end tangent reversed against its start tangent) is a hairpin, two
    // ends apart along it; a loop closes, smoothly or through a
    // crossing (a lobe of Viviani's figure eight leaves the double point
    // at 45° and returns at −45°).
    let gap = point_at(lo).distance_to(&point_at(hi));
    let doubles_back = || {
        let h = 1e-6 * (hi - lo);
        let (s0, s1) = (point_at(lo), point_at(lo + h));
        let (e0, e1) = (point_at(hi - h), point_at(hi));
        let (ds, de) = (s1 - s0, e1 - e0);
        ds.dot(&de) < -0.9 * ds.length() * de.length()
    };
    let closed = gap <= 1e-12 * (1.0 + inner) || (gap <= closing && !doubles_back());
    if closed {
        let outer = arc_between(point_at, lo, t0) + arc_between(point_at, t1, hi);
        inner.min(outer)
    } else {
        inner
    }
}

/// Whether some parameter of `p` on `record` and some of `q` on it lie
/// within their zones of each other along it.
pub(super) fn near_on_record(
    ctx: &ResolutionContext<'_>,
    p: &CandidateVertex,
    q: &CandidateVertex,
    record: usize,
) -> bool {
    let Some(r) = ctx.records.get(record) else {
        return false;
    };
    let curve = &r.curve.curve_3d;
    let params = |c: &CandidateVertex| -> Vec<f64> {
        c.incidence
            .iter()
            .filter_map(|i| match i {
                Incidence::Record { record: x, s, .. } if *x == record => Some(*s),
                _ => None,
            })
            .collect()
    };
    let (ps, qs) = (params(p), params(q));
    ps.iter().any(|&a| {
        qs.iter().any(|&b| {
            separation(
                |s| curve.point_at(s),
                curve.parameter_range(),
                (a, b),
                p.radius + q.radius,
            ) <= p.zone + q.zone
        })
    })
}
