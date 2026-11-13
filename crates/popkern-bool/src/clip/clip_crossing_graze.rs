//! A record grazing its face's boundary within the tolerance crosses
//! nothing there (PK-1181).
//!
//! Where a record runs within the linear reach of a boundary edge and
//! leaves it again on the side it came from, the intersector reports
//! the two ends of the stretch within its tolerance as two crossings.
//! At the tolerance's resolution the record and the edge are one place
//! along that stretch, and the record never leaves the face (or never
//! enters it): the pair is a graze, not an exit and a re-entry.
//! Measured: two balls overlapping above a box's top, their rims half a
//! tolerance apart, the circle they overlap in passing 4.3e-7 over the
//! first ball's rim at the top; the clip of that circle against the
//! first ball's cap read the stretch's ends, 1.4e-3 apart, as crossings,
//! kept the stretch between them and dropped the circle, and the union
//! served both caps whole (the lens counted twice).
//!
//! A pair is a graze when the record stays within the reach of the edge
//! in space between the two reports and lies on the same side of the
//! edge beyond both, where it is more than two reaches off it (nearer,
//! a side is noise). The side is read in the chart against the edge's
//! pcurve at the record's foot on it, near the reports, where the
//! record and the boundary share a branch. A pair leaving to opposite
//! sides is a crossing (`resolve_contact_runs` makes it one), and a
//! pair whose sides abstain keeps its reports.
//!
//! Only a closed record's lone pair is judged: there the clip keeps the
//! hull between two crossings, the stretch, and drops the record. An
//! open record's graze, or one among other crossings, keeps its reports
//! for the arms that read them (a cap a tolerance short of a ring's
//! extreme touches the ring's rim where another record ends, and its
//! graze reports split it there: dropped, the complex declined).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::HalfEdgeHandle;

use super::clip_uv::UvCrossing;

/// Stations between a pair's reports, each within the reach.
const STATIONS: usize = 8;
/// Doublings of the step beyond a report looking for a decisive side.
const SIDE_STEPS: usize = 40;

/// The `crossings` (sorted by `t_curve`) of a closed record without its
/// lone pair when that pair is a graze.
pub(super) fn drop_grazes(
    (surface, boundary): (Option<&SurfaceGeometry>, &[(HalfEdgeHandle, NurbsCurve)]),
    curve_uv: &NurbsCurve,
    crossings: Vec<UvCrossing>,
    tolerance: &Tolerance,
) -> Vec<UvCrossing> {
    let Some(surface) = surface else {
        return crossings;
    };
    let (ct0, ct1) = curve_uv.parameter_range();
    // A closed record's two crossings: elsewhere the clip's tuned arms
    // read a graze's reports (a record touching a rim where another
    // record meets it ends there, and its reports split it).
    let closed = {
        let (p, q) = (curve_uv.point_at(ct0), curve_uv.point_at(ct1));
        surface
            .point_at(p.x, p.y)
            .distance_to(&surface.point_at(q.x, q.y))
            <= tolerance.linear()
    };
    if crossings.len() != 2 || !closed {
        return crossings;
    }
    let (lo, hi) = (ct0.min(ct1), ct0.max(ct1));
    let span = (hi - lo).max(1e-12);
    let mut out = Vec::with_capacity(crossings.len());
    let mut i = 0;
    while i < crossings.len() {
        let graze = i + 1 < crossings.len()
            && crossings[i + 1].boundary_he == crossings[i].boundary_he
            && boundary
                .iter()
                .find(|(h, _)| *h == crossings[i].boundary_he)
                .is_some_and(|(_, pc)| {
                    let before = if i > 0 { crossings[i - 1].t_curve } else { lo };
                    let after = crossings.get(i + 2).map_or(hi, |c| c.t_curve);
                    let edge = Edge {
                        surface,
                        pcurve: pc,
                        curve: curve_uv,
                    };
                    let (a, b) = (crossings[i].t_curve, crossings[i + 1].t_curve);
                    // A spread pair: reports at one parameter were judged
                    // by the duplicate collapse.
                    b - a >= span * 1e-5
                        && edge.within_between(a, b, tolerance)
                        && edge.same_side_beyond((before, a), (b, after), tolerance)
                });
        if graze {
            i += 2;
        } else {
            out.push(crossings[i].clone());
            i += 1;
        }
    }
    out
}

/// A boundary edge's pcurve and a record's, on their surface.
struct Edge<'a> {
    surface: &'a SurfaceGeometry,
    pcurve: &'a NurbsCurve,
    curve: &'a NurbsCurve,
}

impl Edge<'_> {
    /// The record's signed side of the edge in the chart at `t`, and
    /// its distance off the edge in space there.
    fn side(&self, t: f64) -> Option<(f64, f64)> {
        let q = self.curve.point_at(t);
        let tb = crate::curve_projection::project_onto_curve_dense(self.pcurve, &q)?;
        let (b, tan) = (self.pcurve.point_at(tb), self.pcurve.tangent_at(tb));
        let s = tan.x * (q.y - b.y) - tan.y * (q.x - b.x);
        let d = self
            .surface
            .point_at(q.x, q.y)
            .distance_to(&self.surface.point_at(b.x, b.y));
        Some((s, d))
    }

    /// Whether the record stays within the reach of the edge from `a`
    /// to `b`.
    fn within_between(&self, a: f64, b: f64, tolerance: &Tolerance) -> bool {
        (1..STATIONS).all(|k| {
            let t = a + (b - a) * k as f64 / STATIONS as f64;
            self.side(t)
                .is_some_and(|(_, d)| tolerance.within_linear(d))
        })
    }

    /// The record's side of the edge stepping from `from` toward `limit`
    /// (no farther than halfway), at the first point more than two
    /// reaches off the edge.
    fn decisive_side(
        &self,
        from: f64,
        limit: f64,
        step: f64,
        tolerance: &Tolerance,
    ) -> Option<f64> {
        let band = 2.0 * tolerance.linear_reach();
        let room = 0.5 * (limit - from);
        let mut h = step.copysign(room);
        for _ in 0..SIDE_STEPS {
            if h.abs() > room.abs() {
                return None;
            }
            let (s, d) = self.side(from + h)?;
            if d > band {
                return (s != 0.0).then_some(s.signum());
            }
            h *= 2.0;
        }
        None
    }

    /// Whether the record lies on one side of the edge beyond both ends
    /// of the pair, decisively.
    fn same_side_beyond(
        &self,
        (before, a): (f64, f64),
        (b, after): (f64, f64),
        tolerance: &Tolerance,
    ) -> bool {
        let step = (b - a).max(1e-12);
        match (
            self.decisive_side(a, before, step, tolerance),
            self.decisive_side(b, after, step, tolerance),
        ) {
            (Some(x), Some(y)) => x == y,
            _ => false,
        }
    }
}

#[cfg(test)]
#[path = "clip_crossing_graze_tests.rs"]
mod tests;
