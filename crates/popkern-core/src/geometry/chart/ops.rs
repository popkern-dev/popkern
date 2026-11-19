//! The operations every consumer uses (design §5): lift, canonical
//! branch, space distance, and the event rule's test for an artifact
//! near a point.

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::descriptor::{Axis, ChartDescriptor};

pub use super::interpolate::{centroid, evaluate, interpolate};

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
/// its point within which the chart's speed along the row has collapsed
/// under a tenth of its speed away from the row (at the middle of the
/// domain across it). On a sphere it is where `sin v` falls under a tenth,
/// as the across/along ratio of ten reads it; but a long thin strip is
/// anisotropic everywhere, and the ratio read a vanishing fillet strip's
/// whole length as distorted, sending every interpolation on it through
/// the pole projection (ten times the tee saddle's time). The collapse is
/// what a chart singularity is.
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
    let along = |step: f64| {
        let (u, v) = at(step);
        let (su, sv) = s.first_partials(u, v);
        match row.axis {
            Axis::V => su.length(),
            Axis::U => sv.length(),
        }
    };
    let reference = along(0.5 * span).max(1e-300);
    let collapsed = |step: f64| along(step) < 0.1 * reference;
    let (mut inside, mut outside) = (1e-12 * span, 1e-12 * span);
    while collapsed(outside) && outside < 0.5 * span {
        inside = outside;
        outside *= 2.0;
    }
    for _ in 0..60 {
        let mid = 0.5 * (inside + outside);
        if collapsed(mid) {
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
pub fn row_resolutions(d: &ChartDescriptor, share: f64) -> Vec<f64> {
    d.singular_rows
        .iter()
        .map(|row| share * row.distortion)
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

/// The surface's unit normal at face parameters `uv` (its point `p`), in
/// the face chart's sense (`S_u × S_v`, as its limit where that
/// vanishes): at a chart singularity (a sphere's pole, a collapsed row on
/// a smooth cap) the row's normal, which every approach shares; `None` at
/// a geometric singularity (a cone's apex: the limit depends on the
/// approach, so no facet can track one normal there) or where the chart
/// is degenerate with no singular row to read. A point is at a singular
/// row within `exact` of its point in space (rounding: a point a hair off
/// a pole reads its own `S_u × S_v` accurately). Non-periodic coordinates
/// outside the domain are clamped into it first (a projection can land a
/// hair past an edge).
pub fn chart_normal(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    uv: (f64, f64),
    p: Point3,
    exact: f64,
) -> Option<crate::geometry::Vec3> {
    // At a chart singularity and within 1e-5 of its collapse radius (where
    // the speed along the row is a millionth of its speed away from it,
    // and `S_u × S_v` is a cancellation, noise in direction) the row's
    // normal, which the true one approaches: raw normals there sent a
    // vanishing fillet strip's display refinement to its growth budget.
    for row in &d.singular_rows {
        let reach = match row.kind {
            super::SingularKind::Chart => exact.max(1e-5 * row.distortion),
            _ => exact,
        };
        if row.point.distance_to(&p) <= reach {
            return match super::regular::pole_projection(d, s, row) {
                Some(super::RegularChart::PoleProjection { e1, e2, .. }) => Some(e1.cross(&e2)),
                _ => None,
            };
        }
    }
    let ((u0, u1), (v0, v1)) = d.domain;
    let (c, _) = canonical(d, uv);
    // A closure takes its representative, an evaluation period is used as
    // given (as `evaluate`), and a coordinate with no period is clamped.
    let place = |x: f64, cx: f64, lo: f64, hi: f64, axis: Axis| match d.period(axis) {
        Some(p) if p.kind == super::PeriodKind::GeometricClosure => cx,
        Some(_) => x,
        None if lo.is_finite() && hi.is_finite() => x.clamp(lo, hi),
        None => x,
    };
    let (u, v) = (
        place(uv.0, c.0, u0, u1, Axis::U),
        place(uv.1, c.1, v0, v1, Axis::V),
    );
    let (su, sv) = s.first_partials(u, v);
    let n = su.cross(&sv);
    let l = n.length();
    if l.is_finite() && l > 1e-300 * su.length().max(sv.length()).max(1e-300) {
        return Some(n * (1.0 / l));
    }
    // Degenerate where the clamp landed: a point a hair past a pole
    // clamps onto the row itself, whose normal is the row's.
    let q = s.point_at(u, v);
    let scale = q.x.abs().max(q.y.abs()).max(q.z.abs()).max(1.0);
    d.singular_rows
        .iter()
        .find(|row| row.point.distance_to(&q) <= exact.max(SINGULAR_EXACT * scale))
        .and_then(|row| match super::regular::pole_projection(d, s, row) {
            Some(super::RegularChart::PoleProjection { e1, e2, .. }) => Some(e1.cross(&e2)),
            _ => None,
        })
}

/// The regular chart for a neighbourhood reaching `radius` in space from
/// `p0` (design §4): the pole projection about the first chart
/// singularity the neighbourhood holds, else the face chart (which is
/// right beside a singularity as long as the neighbourhood does not wrap
/// round it; reading the distortion radius sent a fillet strip's every
/// relaxed star through the projection). A consumer moving points within the neighbourhood works in
/// it (a vertex star relaxed beside a pole, where face-chart moves are
/// distorted without bound).
pub fn neighbourhood_chart(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    p0: &Point3,
    radius: f64,
) -> super::RegularChart {
    d.singular_rows
        .iter()
        .filter(|row| row.kind == super::SingularKind::Chart)
        .find(|row| row.point.distance_to(p0) <= radius)
        .and_then(|row| super::regular::pole_projection(d, s, row))
        .unwrap_or(super::RegularChart::Face)
}

/// The shortest distance between `a` and `b` along an axis of period
/// `period` (no farther than half a period): the one place a consumer
/// compares two periodic parameters.
pub fn periodic_gap(a: f64, b: f64, period: f64) -> f64 {
    let d = popkern_math::periodic::offset_mod(a, b, period);
    d.min(period - d)
}

/// The signed offset from `b` to `a` on the branch nearest `b`, in
/// `(-period / 2, period / 2]`: the step an unwrap adds.
pub fn nearest_offset(a: f64, b: f64, period: f64) -> f64 {
    let d = popkern_math::periodic::offset_mod(a, b, period);
    if d > 0.5 * period { d - period } else { d }
}

/// `x` brought into `[lo, lo + period)` by whole periods (a curve's
/// angle into its parameter window).
pub fn wrap_into(x: f64, lo: f64, period: f64) -> f64 {
    popkern_math::periodic::wrap_from(x, lo, period)
}

/// The offset from `b` to `a` reduced into `[0, period)` (forward along
/// the axis).
pub fn forward_offset(a: f64, b: f64, period: f64) -> f64 {
    popkern_math::periodic::offset_mod(a, b, period)
}

/// `x` reduced into `[0, period)`.
pub fn reduce(x: f64, period: f64) -> f64 {
    popkern_math::periodic::modulo(x, period)
}
