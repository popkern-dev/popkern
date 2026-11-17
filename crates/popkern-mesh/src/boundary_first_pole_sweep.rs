//! A curve's chart step beside a pole, read off its own points
//! (PK-1178).
//!
//! The face rings unwrap each station step of an edge as the nearest
//! period image, faithful only while the step turns under half a
//! period; the sampler splits steps turning past a quarter
//! (`EdgeNormalProbe::wide_uv_step`) reading each face's pcurve at the
//! stations' curve fractions. A pcurve lying on its curve need not
//! share its parameter: a ball cut by a plane whose circle passes
//! 1.8e-5 beside its pole, the seam turned 0.4, splits the circle on
//! the seam beside the pole, and the arc leaving there turns -1.96 rad
//! within 3e-5 of its range where its pcurve turns over 2.4e-4. The
//! step from the arc's end to its next station, under the span floor,
//! turned 1.96 rad unsplit, the ring took it as a pole crossing run
//! along the pole row the face flag's way (+4.32 rad), where the row
//! also closes round the pole the face holds, and the sphere face fell
//! to the fallback grid (its area 0.84 over, the volume 5.3e-4 low).
//!
//! Here the step is judged from the curve's points, each inverted onto
//! each incident sphere or cone, the charts with a singular row a short
//! step can turn round (closed-form inversions, so the sampler's
//! passes stay cheap); a point where the chart is singular (the pole
//! itself) has no longitude and judges nothing.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, Edge};

/// Ratio of `|S_u|` to `|S_v|` under which a chart point is its pole.
const POLE_RATIO: f64 = 1e-9;
/// Smallest curve-window fraction a step is split to below the span
/// floor: a curve beside a pole turns its longitude within a stretch as
/// short as its distance from the pole.
const SWEEP_FLOOR: f64 = 1e-12;

/// Whether the sampler splits the step between curve fractions `a` and
/// `b` (`at` the curve's point at a fraction) under its span floor.
pub(super) fn splits_below_floor(
    brep: &BRep,
    edge: &Edge,
    (a, b): (f64, f64),
    at: &dyn Fn(f64) -> Point3,
) -> bool {
    b - a > SWEEP_FLOOR && turns_past_quarter(brep, edge, &at(a), &at(b))
}

/// Whether the curve's step from `pa` to `pb` turns a periodic
/// coordinate of an incident face's chart past a quarter period.
fn turns_past_quarter(brep: &BRep, edge: &Edge, pa: &Point3, pb: &Point3) -> bool {
    let mate = brep
        .half_edges
        .get(edge.half_edge)
        .and_then(|he| he.mate())
        .filter(|m| *m != edge.half_edge);
    [Some(edge.half_edge), mate]
        .into_iter()
        .flatten()
        .any(|heh| {
            let Some(surface) = brep
                .half_edges
                .get(heh)
                .and_then(|he| brep.loops.get(he.loop_))
                .and_then(|l| brep.faces.get(l.face))
                .and_then(|f| f.surface.as_ref())
                .map(|sof| &sof.geometry)
            else {
                return false;
            };
            // The charts with a singular row a short step can turn round:
            // a sphere's poles, a cone's apex.
            if !matches!(
                surface,
                SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. }
            ) {
                return false;
            }
            let (period_u, period_v) =
                popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(surface);
            // A chart point, `None` at a pole (no longitude there).
            let chart = |p: &Point3| {
                let (u, v) = surface.closest_point(p)?;
                let (su, sv) = (
                    surface.partial_u(u, v).length(),
                    surface.partial_v(u, v).length(),
                );
                (su > POLE_RATIO * sv && sv > POLE_RATIO * su).then_some((u, v))
            };
            let (Some(a), Some(b)) = (chart(pa), chart(pb)) else {
                return false;
            };
            [(b.0 - a.0, period_u), (b.1 - a.1, period_v)]
                .iter()
                .any(|&(d, p)| {
                    p.is_some_and(|p| p > 0.0 && (d - (d / p).round() * p).abs() > 0.25 * p)
                })
        })
}
