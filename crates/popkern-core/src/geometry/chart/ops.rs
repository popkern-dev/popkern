//! The operations every consumer uses (design §5): lift, canonical
//! branch, space distance, and the event rule's test for an artifact
//! near a point.

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::descriptor::{Axis, ChartDescriptor};

/// A lifted parameter polyline: continuous across cuts, and its winding
/// about each closed direction (whole periods from first to last point).
#[derive(Debug, Clone, PartialEq)]
pub struct Lifted {
    /// The points, each on the period branch continuous with the last.
    pub points: Vec<(f64, f64)>,
    /// Whole periods turned along u and v.
    pub winding: (i64, i64),
}

/// The closed directions' periods of a descriptor.
fn periods(d: &ChartDescriptor) -> (Option<f64>, Option<f64>) {
    (d.u.map(|c| c.period), d.v.map(|c| c.period))
}

/// Lift `points` across the chart's cuts: each point moves by whole
/// periods to the branch nearest its predecessor in space, read through
/// the parameters (the branch nearest in parameter is nearest in space
/// for a period cut; a pole is crossed by the caller's regular chart).
pub fn lift(d: &ChartDescriptor, points: &[(f64, f64)]) -> Lifted {
    let (pu, pv) = periods(d);
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(points.len());
    for &(u, v) in points {
        let (u, v) = match out.last() {
            Some(&(lu, lv)) => (nearest(u, lu, pu), nearest(v, lv, pv)),
            None => (u, v),
        };
        out.push((u, v));
    }
    let winding = match (out.first(), out.last()) {
        (Some(a), Some(b)) => (whole(b.0 - a.0, pu), whole(b.1 - a.1, pv)),
        _ => (0, 0),
    };
    Lifted {
        points: out,
        winding,
    }
}

fn nearest(t: f64, to: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) => t + ((to - t) / p).round() * p,
        None => t,
    }
}

fn whole(delta: f64, period: Option<f64>) -> i64 {
    period.map_or(0, |p| (delta / p).round() as i64)
}

/// The representative of `uv` in the chart's domain, and the whole periods
/// it moved by.
pub fn canonical(d: &ChartDescriptor, uv: (f64, f64)) -> ((f64, f64), (i64, i64)) {
    let ((u0, u1), (v0, v1)) = d.domain;
    let one = |t: f64, lo: f64, hi: f64, axis: Axis| -> (f64, i64) {
        match d.period(axis) {
            Some(c) if lo.is_finite() && !(lo..=hi).contains(&t) => {
                let k = ((t - lo) / c.period).floor();
                (t - k * c.period, -(k as i64))
            }
            _ => (t, 0),
        }
    };
    let (u, ku) = one(uv.0, u0, u1, Axis::U);
    let (v, kv) = one(uv.1, v0, v1, Axis::V);
    ((u, v), (ku, kv))
}

/// The length in space of a chart step (tolerances are lengths).
pub fn space_distance(s: &SurfaceGeometry, a: (f64, f64), b: (f64, f64)) -> f64 {
    s.point_at(a.0, a.1).distance_to(&s.point_at(b.0, b.1))
}

/// An artifact near a point (the event rule, design §5): a singular row's
/// point or a cut within `radius` in space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Artifact {
    /// A singular row's point, and its distance.
    Singular(Point3, f64),
    /// A cut along `axis`, and the point's distance from it.
    Cut(Axis, f64),
}

/// The nearest artifact within `radius` of the point at `uv`, if any.
pub fn near_artifact(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    uv: (f64, f64),
    radius: f64,
) -> Option<Artifact> {
    let p = s.point_at(uv.0, uv.1);
    let mut best: Option<Artifact> = None;
    let mut keep = |a: Artifact, dist: f64| {
        let current = match best {
            Some(Artifact::Singular(_, d)) | Some(Artifact::Cut(_, d)) => d,
            None => f64::INFINITY,
        };
        if dist <= radius && dist < current {
            best = Some(a);
        }
    };
    for row in &d.singular_rows {
        let dist = row.point.distance_to(&p);
        keep(Artifact::Singular(row.point, dist), dist);
    }
    let ((u0, u1), (v0, v1)) = d.domain;
    if d.u.is_some() && u0.is_finite() {
        let dist = s
            .point_at(u0, uv.1)
            .distance_to(&p)
            .min(s.point_at(u1, uv.1).distance_to(&p));
        keep(Artifact::Cut(Axis::U, dist), dist);
    }
    if d.v.is_some() && v0.is_finite() {
        let dist = s
            .point_at(uv.0, v0)
            .distance_to(&p)
            .min(s.point_at(uv.0, v1).distance_to(&p));
        keep(Artifact::Cut(Axis::V, dist), dist);
    }
    best
}

/// A chart singularity's distortion radius (design §4): the distance from
/// its point within which the chart's speed across the row exceeds ten
/// times its speed along it. Found by doubling a step off the row along
/// the other parameter until the ratio falls to ten, then bisecting.
pub fn distortion_radius(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    row: &super::descriptor::SingularRow,
) -> f64 {
    let ((u0, u1), (v0, v1)) = d.domain;
    let (lo, hi, a, b) = match row.axis {
        Axis::V => (v0, v1, u0, u1),
        Axis::U => (u0, u1, v0, v1),
    };
    let span = if lo.is_finite() && hi.is_finite() {
        hi - lo
    } else {
        1.0
    };
    let inward = if (row.at - lo).abs() <= (row.at - hi).abs() {
        1.0
    } else {
        -1.0
    };
    let t = if a.is_finite() && b.is_finite() {
        0.5 * (a + b)
    } else {
        0.0
    };
    let at = |step: f64| match row.axis {
        Axis::V => (t, row.at + inward * step),
        Axis::U => (row.at + inward * step, t),
    };
    let ratio = |step: f64| {
        let (u, v) = at(step);
        let (su, sv) = s.first_partials(u, v);
        let (along, across) = match row.axis {
            Axis::V => (su.length(), sv.length()),
            Axis::U => (sv.length(), su.length()),
        };
        across / along.max(1e-300)
    };
    let (mut inside, mut outside) = (1e-12 * span, 1e-12 * span);
    while ratio(outside) > 10.0 && outside < 0.5 * span {
        inside = outside;
        outside *= 2.0;
    }
    for _ in 0..60 {
        let mid = 0.5 * (inside + outside);
        if ratio(mid) > 10.0 {
            inside = mid
        } else {
            outside = mid
        }
    }
    let (u, v) = at(0.5 * (inside + outside));
    s.point_at(u, v).distance_to(&row.point)
}

/// The singular row whose point lies within `reach` of `p` (in space), if
/// any, and whether its row collapses the u coordinate (a row `v = at`
/// maps every u to its point) rather than v. The one question a consumer
/// asks about a pole or an apex: in space, never by partial ratios, which
/// read a point a fixed fraction of the radius off the pole as the pole
/// whatever the scale.
pub fn singular_point_near(
    d: &ChartDescriptor,
    p: &Point3,
    reach: f64,
) -> Option<(super::descriptor::SingularRow, bool)> {
    d.singular_rows
        .iter()
        .filter(|r| r.point.distance_to(p) <= reach)
        .min_by(|a, b| a.point.distance_to(p).total_cmp(&b.point.distance_to(p)))
        .map(|r| (*r, r.axis == Axis::V))
}

/// `uv` moved by whole periods to the branch nearest `to` (a lift to a
/// reference point: a triangle's corners onto one branch, a chord end
/// onto its trim's).
pub fn aligned(d: &ChartDescriptor, uv: (f64, f64), to: (f64, f64)) -> (f64, f64) {
    let (pu, pv) = periods(d);
    (nearest(uv.0, to.0, pu), nearest(uv.1, to.1, pv))
}

/// The whole-period shifts `(k_u·T_u, k_v·T_v)` that bring `uv` nearest
/// `to` (the lift amount itself, exact, where `aligned` returns the moved
/// point).
pub fn shift_toward(d: &ChartDescriptor, uv: (f64, f64), to: (f64, f64)) -> (f64, f64) {
    let (pu, pv) = periods(d);
    let shift = |t: f64, to: f64, p: Option<f64>| p.map_or(0.0, |p| p * ((to - t) / p).round());
    (shift(uv.0, to.0, pu), shift(uv.1, to.1, pv))
}

/// A mesh point within this share of the body's coordinate scale of a
/// singular row's point is on it to rounding.
pub const SINGULAR_EXACT: f64 = 1e-12;

/// A face vertex stands for a singular row only within this share of the
/// row's distortion radius (or of its triangle's other corners' distance
/// from the row point): beyond it the chart resolves the vertex's
/// collapsed coordinate, the mesher tiles with it, and reading it as the
/// row double-covers the fan. 1e-3 of the distortion radius is 1e-4 of a
/// sphere's radius. A face vertex within its tolerance but past this is a
/// vertex beside the row: a tiny rounded cube at a tolerance of 0.29 of
/// its corner balls' radius read its ball corners as poles and its volume
/// 1.6e-5 high.
pub const SINGULAR_RESOLVED: f64 = 1e-3;

/// Per singular row of `d`, the distance within which a point is resolved
/// at the row rather than beside it: `share` of the row's distortion
/// radius (beyond it the chart resolves the point's collapsed coordinate,
/// so the point carries it).
pub fn row_resolutions(d: &ChartDescriptor, s: &SurfaceGeometry, share: f64) -> Vec<f64> {
    d.singular_rows
        .iter()
        .map(|row| share * distortion_radius(d, s, row))
        .collect()
}

/// The singular row a mesh corner `points[k]` stands for, if any (design
/// §5, the event rule's measuring form). On the row's point to `exact`
/// (rounding), any corner stands for it. A corner carrying a vertex
/// tolerance stands for it within that tolerance where it is also
/// resolved at the row: within `resolutions[r]` of it, or within `share`
/// of the other corners' distance from it (a triangle small beside the
/// row). Returns the row and whether it collapses u.
pub fn corner_on_singular_row(
    d: &ChartDescriptor,
    resolutions: &[f64],
    points: [Point3; 3],
    k: usize,
    (vertex_tolerance, exact, share): (Option<f64>, f64, f64),
) -> Option<(super::descriptor::SingularRow, bool)> {
    d.singular_rows.iter().enumerate().find_map(|(r, row)| {
        let dist = row.point.distance_to(&points[k]);
        if dist <= exact {
            return Some((*row, row.axis == Axis::V));
        }
        let tolerance = vertex_tolerance?;
        let far = (1..3)
            .map(|j| row.point.distance_to(&points[(k + j) % 3]))
            .fold(f64::INFINITY, f64::min);
        let resolved = resolutions.get(r).copied().unwrap_or(0.0);
        (dist <= tolerance && (dist <= resolved || dist <= share * far))
            .then_some((*row, row.axis == Axis::V))
    })
}
