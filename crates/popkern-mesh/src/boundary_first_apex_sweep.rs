//! A curve's chart step beside a cone's apex, read off its own points
//! (PK-1178).
//!
//! The face rings unwrap each station step of an edge as the nearest
//! period image, faithful only while the step turns under half a
//! period; the sampler splits steps turning past a quarter
//! (`EdgeNormalProbe::wide_uv_step`) reading each face's pcurve at the
//! stations' curve fractions. A pcurve lying on its curve need not
//! share its parameter: beside a singular point a curve turns its
//! longitude within a stretch as short as its distance from it, and a
//! step under the span floor turned past half a period unsplit (found on
//! a ball whose cut circle passed 1.8e-5 beside its pole, a face the
//! atlas now meshes in charts with no singular row; a cone's apex stays
//! on the face chart).
//!
//! Here the step is judged from the curve's points, each inverted onto
//! each incident cone (closed-form inversions, so the sampler's passes
//! stay cheap); a point at the apex itself has no longitude and judges
//! nothing.

use popkern_core::geometry::Point3;
use popkern_core::geometry::chart::{ChartDescriptor, ops};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, Edge};

/// Smallest curve-window fraction a step is split to below the span
/// floor: a curve beside an apex turns its longitude within a stretch as
/// short as its distance from the apex.
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
            // The face chart with a singular row a short step can turn
            // round: a cone's apex (a sphere's poles are the atlas's).
            if !matches!(surface, SurfaceGeometry::Cone { .. }) {
                return false;
            }
            let (period_u, period_v) = ChartDescriptor::evaluation_periods_of(surface);
            let descriptor = ChartDescriptor::of(surface, &brep.tolerance);
            // A chart point, `None` at the apex, within the linear
            // tolerance of it (one place: no longitude there).
            let chart = |p: &Point3| {
                if ops::singular_point_near(&descriptor, p, brep.tolerance.linear()).is_some() {
                    return None;
                }
                surface.closest_point(p)
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
