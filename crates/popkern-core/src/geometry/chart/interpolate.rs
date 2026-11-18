//! Interpolation and evaluation in a regular chart (design §5): a point
//! between surface points, a centroid, and evaluation on any period
//! branch.

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::descriptor::{Axis, ChartDescriptor};
use super::ops::{SINGULAR_EXACT, aligned, canonical};

/// The point a fraction `t` of the way from `a` to `b` (each face
/// parameters and their point) along the surface, interpolated in a
/// regular chart (design §5): in the pole projection when either end lies
/// within a chart singularity's distortion radius (a face-chart midpoint
/// of two points either side of a pole lands on the far side of the
/// sphere), else in the face chart with `b` lifted to `a`'s branch. The
/// parameters come back on the branch of the end farther from the row (a
/// mesh's triangles share one lifted branch; the singular point's own
/// collapsed coordinate is arbitrary). A geometric singularity (an apex)
/// is not re-charted: interpolation there runs in the face chart, which
/// the consumer treats as a topological event.
///
/// # Errors
///
/// The pole projection's inverse failing to converge.
pub fn interpolate(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    a: ((f64, f64), Point3),
    b: ((f64, f64), Point3),
    t: f64,
) -> Result<((f64, f64), Point3), super::ChartError> {
    // The face chart first: its point is the surface's and is right
    // wherever it lands near the chord. Only where it strays (two points
    // either side of a pole) or an end sits on the singular point does the
    // pole projection answer; reading it for every point within the
    // distortion radius ran a fillet strip's every split through the
    // projection's inverse.
    let bl = aligned(d, b.0, a.0);
    let face_uv = (a.0.0 * (1.0 - t) + bl.0 * t, a.0.1 * (1.0 - t) + bl.1 * t);
    let face_p = evaluate(d, s, face_uv);
    let chord = a.1 + (b.1 - a.1) * t;
    let strays = face_p.distance_to(&chord) > 0.25 * a.1.distance_to(&b.1);
    let exact = |p: &Point3| SINGULAR_EXACT * p.x.abs().max(p.y.abs()).max(p.z.abs()).max(1.0);
    let on_row = d.singular_rows.iter().any(|row| {
        row.point.distance_to(&a.1) <= exact(&a.1) || row.point.distance_to(&b.1) <= exact(&b.1)
    });
    if !strays && !on_row {
        return Ok((face_uv, face_p));
    }
    let near = d.singular_rows.iter().find(|row| {
        row.kind == super::SingularKind::Chart
            && (row.point.distance_to(&a.1) <= row.distortion
                || row.point.distance_to(&b.1) <= row.distortion)
    });
    if let Some(row) = near
        && let Some(chart) = super::regular::pole_projection(d, s, row)
        && let super::RegularChart::PoleProjection { centre, e1, e2, .. } = &chart
    {
        let project = |p: &Point3| ((*p - *centre).dot(e1), (*p - *centre).dot(e2));
        let (pa, pb) = (project(&a.1), project(&b.1));
        let xy = (pa.0 + (pb.0 - pa.0) * t, pa.1 + (pb.1 - pa.1) * t);
        // Seeded from the end farther from the row (its parameters are
        // meaningful; the singular point's collapsed one is not).
        let seed = if row.point.distance_to(&a.1) >= row.point.distance_to(&b.1) {
            a.0
        } else {
            b.0
        };
        let uv = chart.to_surface_near(s, xy, Some(canonical(d, seed).0))?;
        let reference = if row.point.distance_to(&a.1) >= row.point.distance_to(&b.1) {
            a.0
        } else {
            b.0
        };
        let uv = aligned(d, uv, reference);
        return Ok((uv, evaluate(d, s, uv)));
    }
    // `a (1 − t) + b t`: at t = ½ exactly the plain midpoint `½ (a + b)`
    // (halving is exact), so a face-chart split is unchanged by routing
    // it here.
    let bl = aligned(d, b.0, a.0);
    let uv = (a.0.0 * (1.0 - t) + bl.0 * t, a.0.1 * (1.0 - t) + bl.1 * t);
    Ok((uv, evaluate(d, s, uv)))
}

/// The surface point at chart parameters `uv` on any period branch:
/// evaluated at the canonical representative, since a geometric closure
/// (a clamped closed NURBS) is not periodic in evaluation and extrapolates
/// past its domain. Lifted parameters are chart coordinates; evaluation
/// always comes through here.
pub fn evaluate(d: &ChartDescriptor, s: &SurfaceGeometry, uv: (f64, f64)) -> Point3 {
    // Only a geometric closure needs the representative; an analytic or
    // periodic-form period evaluates periodically, and is evaluated as
    // given (bit-identical to the surface's own evaluation there).
    let (c, _) = canonical(d, uv);
    let closure = |axis: Axis| {
        d.period(axis)
            .is_some_and(|p| p.kind == super::PeriodKind::GeometricClosure)
    };
    let u = if closure(Axis::U) { c.0 } else { uv.0 };
    let v = if closure(Axis::V) { c.1 } else { uv.1 };
    s.point_at(u, v)
}

/// The centroid of `points` (face parameters and their points) on the
/// surface, in a regular chart: the pole projection when any lies within
/// a chart singularity's distortion radius, else the face chart with each
/// lifted to the first's branch (the mean `Σ / n` there, in the given
/// order). Parameters come back on the branch of the point farthest from
/// the row (pole projection) or the first's (face chart).
///
/// # Errors
///
/// The pole projection's inverse failing to converge; an empty slice.
pub fn centroid(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    points: &[((f64, f64), Point3)],
) -> Result<((f64, f64), Point3), super::ChartError> {
    let Some(first) = points.first() else {
        return Err(super::ChartError::InverseDiverged {
            x: 0.0,
            y: 0.0,
            residual: f64::INFINITY,
        });
    };
    let n = points.len() as f64;
    // The face chart first (as `interpolate`): right wherever its point
    // lands near the points' mean in space.
    {
        let (mut su, mut sv) = (0.0, 0.0);
        let (mut mx, mut my, mut mz) = (0.0, 0.0, 0.0);
        for q in points {
            let a = aligned(d, q.0, first.0);
            su += a.0;
            sv += a.1;
            mx += q.1.x;
            my += q.1.y;
            mz += q.1.z;
        }
        let uv = (su * (1.0 / n), sv * (1.0 / n));
        let p = evaluate(d, s, uv);
        let mean = Point3::new(mx / n, my / n, mz / n);
        let spread = points
            .iter()
            .map(|q| q.1.distance_to(&mean))
            .fold(0.0, f64::max);
        let exact = |p: &Point3| SINGULAR_EXACT * p.x.abs().max(p.y.abs()).max(p.z.abs()).max(1.0);
        let on_row = points.iter().any(|q| {
            d.singular_rows
                .iter()
                .any(|row| row.point.distance_to(&q.1) <= exact(&q.1))
        });
        if !on_row && p.distance_to(&mean) <= 0.5 * spread {
            return Ok((uv, p));
        }
    }
    let near = d.singular_rows.iter().find(|row| {
        row.kind == super::SingularKind::Chart
            && points
                .iter()
                .any(|q| row.point.distance_to(&q.1) <= row.distortion)
    });
    if let Some(row) = near
        && let Some(chart) = super::regular::pole_projection(d, s, row)
        && let super::RegularChart::PoleProjection { centre, e1, e2, .. } = &chart
    {
        let (mut x, mut y) = (0.0, 0.0);
        for q in points {
            let r = q.1 - *centre;
            x += r.dot(e1);
            y += r.dot(e2);
        }
        let far = points
            .iter()
            .max_by(|a, b| {
                row.point
                    .distance_to(&a.1)
                    .total_cmp(&row.point.distance_to(&b.1))
            })
            .map_or(first.0, |q| q.0);
        let uv = chart.to_surface_near(s, (x / n, y / n), Some(canonical(d, far).0))?;
        let reference = points
            .iter()
            .max_by(|a, b| {
                row.point
                    .distance_to(&a.1)
                    .total_cmp(&row.point.distance_to(&b.1))
            })
            .map_or(first.0, |q| q.0);
        let uv = aligned(d, uv, reference);
        return Ok((uv, evaluate(d, s, uv)));
    }
    let (mut su, mut sv) = (0.0, 0.0);
    for q in points {
        let a = aligned(d, q.0, first.0);
        su += a.0;
        sv += a.1;
    }
    let uv = (su * (1.0 / n), sv * (1.0 / n));
    Ok((uv, evaluate(d, s, uv)))
}
