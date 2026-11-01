//! The face-side rule (plan §4.3a): two candidates on one face, on
//! opposite sides of a boundary edge of that face, are points of
//! different fragments-to-be and stay distinct.
//!
//! A candidate within its own reach of the edge lies on it (the
//! interior-incidence pass gives it the edge), so the rule bites only
//! where a group grows wider than its members' radii (up to the growth
//! cap) and would span the edge: each member clear of the edge, the
//! edge passing between them.
//!
//! Measured (bool, examples, mesh and blend, table on): no group spans
//! a boundary edge, so the rule changes no served boolean; it keeps the
//! growth fixture's two sides apart.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge, closest_point_in_window};
use popkern_core::topology::{BRep, FaceHandle};

use super::candidates::{CandidateVertex, Incidence};
use super::region::Bounds;

/// The edge curves bounding `face`, every loop.
fn boundary(brep: &BRep, face: FaceHandle) -> Vec<&CurveOnEdge> {
    let mut out = Vec::new();
    let Some(f) = brep.faces.get(face) else {
        return out;
    };
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(c) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()) {
                out.push(c);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}

/// Whether `p` and `q` both lie on a face of one operand and a boundary
/// edge of that face passes between them, each clear of it by more
/// than its radius.
pub(crate) fn opposite_sides(
    operands: [&BRep; 2],
    p: &CandidateVertex,
    q: &CandidateVertex,
) -> bool {
    let gap = p.position.distance_to(&q.position);
    // Each must be clear of the edge by more than its radius and within
    // the gap of it: no edge separates points closer than either radius.
    if gap <= p.radius.max(q.radius) {
        return false;
    }
    p.incidence.iter().any(|ip| {
        let &Incidence::Face { operand, face, .. } = ip else {
            return false;
        };
        let on_both = q
            .incidence
            .iter()
            .any(|iq| matches!(iq, Incidence::Face { operand: o, face: f, .. } if *o == operand && *f == face));
        on_both
            && boundary(operands[operand], face)
                .into_iter()
                .any(|c| separates(c, p, q, gap))
    })
}

/// Whether edge curve `c` passes between `p` and `q` (`gap` apart),
/// each clear of it by more than its radius.
fn separates(c: &CurveOnEdge, p: &CandidateVertex, q: &CandidateVertex, gap: f64) -> bool {
    let mut hull = vec![c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end)];
    crate::face_box::conservative_curve_extent(c, &mut hull);
    let near = Bounds::of(&hull, gap).is_some_and(|b| b.holds(&p.position) && b.holds(&q.position));
    if !near {
        return false;
    }
    let foot = |x: &Point3| {
        c.geometry
            .point_at(closest_point_in_window(&c.geometry, x, c.t_start, c.t_end))
    };
    let (fp, fq) = (foot(&p.position), foot(&q.position));
    let (dp, dq) = (p.position.distance_to(&fp), q.position.distance_to(&fq));
    dp > p.radius
        && dq > q.radius
        && dp <= gap
        && dq <= gap
        && (p.position - fp).dot(&(q.position - fq)) < 0.0
}
