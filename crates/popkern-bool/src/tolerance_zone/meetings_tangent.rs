//! A meeting of two records of one face pair is a tangent crossing
//! (tangency design §4.8, W3c): two transversal intersection curves of
//! one pair of smooth faces cross only where the faces are tangent. The
//! marcher's records stray near such a point (the intersection's
//! Jacobian is singular there): records of a cylinder and a sphere
//! crossing at one exactly met twice, 2.05e-5 apart, two vertex groups
//! for one place, and the cellular result left a 2.5e-6 shell gap. The
//! meeting is refined to the exact contact (`contact_discovery::
//! contact_from`) where it lies within the zone the faces stay within the
//! reach of each other (`√(2·reach / min|λ|)`, the relative curvature's
//! eigenvalues), and stands there, its radius covering how far each
//! record passes from it: every meeting of the crossing then lands on one
//! point and resolves as one vertex.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::contact_classify::ContactScale;
use popkern_core::geometry::intersection::contact_discovery::contact_from;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

use super::candidates::CandidateVertex;
use crate::intersect::IntersectionRecord;

/// The exact tangent crossing a meeting of records `i` and `j` (one face
/// pair) at `near` stands for, with each record's parameter there and the
/// farthest either passes from it; `None` when the records are of
/// different pairs or no contact lies within the zone.
pub(super) fn tangent_crossing(
    (a, b): (&BRep, &BRep),
    records: &[IntersectionRecord],
    (i, j): (usize, usize),
    (near, at): (Point3, [f64; 2]),
    linear: f64,
) -> Option<(Point3, [f64; 2], f64)> {
    let (ri, rj) = (&records[i], &records[j]);
    if (ri.face_a, ri.face_b) != (rj.face_a, rj.face_b) {
        return None;
    }
    let (sa, wa) = surface_window(a, ri.face_a)?;
    let (sb, wb) = surface_window(b, ri.face_b)?;
    let seed = popkern_core::geometry::uv_projection::project_point_to_uv(sa, &near)?;
    let extent = wa.0.1 - wa.0.0;
    let step = 1e-6 * extent.abs().max(wa.1.1 - wa.1.0).max(1e-9);
    let size = sa
        .point_at(wa.0.0, wa.1.0)
        .distance_to(&sa.point_at(wa.0.1, wa.1.1));
    let scale = ContactScale {
        linear,
        extent: size.max(linear),
    };
    let contact = contact_from((sa, wa), (sb, wb), seed, step, scale)?;
    let reach = popkern_core::tolerance::Tolerance::default()
        .with_linear(linear)
        .linear_reach();
    let lambda = contact
        .classification
        .eigenvalues
        .map(|[l0, l1]| l0.abs().min(l1.abs()))
        .unwrap_or(0.0);
    let zone = if lambda > 0.0 {
        (2.0 * reach / lambda).sqrt()
    } else {
        f64::INFINITY
    };
    if contact.point.distance_to(&near) > zone {
        return None;
    }
    let mut params = [0.0; 2];
    let mut pass = 0.0_f64;
    for (k, r) in [ri, rj].into_iter().enumerate() {
        // Near the meeting's own parameter on each branch: a record
        // meeting itself passes the crossing twice.
        let c = &r.curve.curve_3d;
        let (lo, hi) = c.parameter_range();
        let span = 0.05 * (hi - lo);
        let s = popkern_core::geometry::curve::closest_point_in_window(
            c,
            &contact.point,
            (at[k] - span).max(lo),
            (at[k] + span).min(hi),
        );
        params[k] = s;
        pass = pass.max(c.point_at(s).distance_to(&contact.point));
    }
    Some((contact.point, params, pass))
}

/// A surface and a parameter window on it.
type SurfaceWindow<'a> = (&'a SurfaceGeometry, ((f64, f64), (f64, f64)));

/// A face's surface and its trim's parameter window.
fn surface_window(brep: &BRep, face: FaceHandle) -> Option<SurfaceWindow<'_>> {
    let s = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    Some((s, brep.face_trim_uv_window(face)?))
}

/// `cands` with the meetings standing at one tangent crossing (each of
/// `crossings`, within the linear reach) merged into one candidate
/// carrying every branch's incidence: a record passing the crossing
/// twice (a figure-eight) holds both passes at its one vertex, as a
/// record meeting itself does, where two meetings each holding one pass
/// stayed apart (`SeparatedAlongRecord`).
pub(super) fn merge_at_crossings(
    cands: Vec<CandidateVertex>,
    crossings: &[Point3],
    linear: f64,
) -> Vec<CandidateVertex> {
    let reach = popkern_core::tolerance::Tolerance::default()
        .with_linear(linear)
        .linear_reach();
    let at_crossing = |c: &CandidateVertex| {
        crossings
            .iter()
            .position(|p| p.distance_to(&c.position) <= reach)
    };
    let mut out: Vec<CandidateVertex> = Vec::new();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for c in cands {
        let Some(k) = at_crossing(&c) else {
            out.push(c);
            continue;
        };
        match merged
            .iter()
            .find(|(key, _)| crossings[*key].distance_to(&crossings[k]) <= reach)
        {
            Some(&(_, slot)) => {
                let into = &mut out[slot];
                for inc in c.incidence {
                    if !into.incidence.contains(&inc) {
                        into.incidence.push(inc);
                    }
                }
                into.radius = into.radius.max(c.radius);
                into.zone = into.zone.max(c.zone);
                if c.key < into.key {
                    into.key = c.key;
                }
            }
            None => {
                merged.push((k, out.len()));
                out.push(c);
            }
        }
    }
    out
}

/// The records and their meetings (the `arrangement` trace topic).
pub(super) fn trace_meetings(
    records: &[IntersectionRecord],
    found: &[super::meetings::RecordMeeting],
) {
    for (i, r) in records.iter().enumerate() {
        let c = &r.curve.curve_3d;
        let (lo, hi) = c.parameter_range();
        eprintln!(
            "  [arrangement] record {i}: degree {}, {} control points, closed {}, from {:?} to {:?}",
            c.degree,
            c.control_points.len(),
            r.curve.is_closed,
            c.point_at(lo),
            c.point_at(hi)
        );
    }
    for m in found {
        eprintln!(
            "  [arrangement] records {:?} at {:?}: gap {:.2e}, zone {:.2e}, sin {:.2e}",
            m.records, m.point, m.gap, m.zone, m.sin_angle
        );
    }
}

#[cfg(test)]
#[path = "meetings_tangent_tests.rs"]
mod tests;
