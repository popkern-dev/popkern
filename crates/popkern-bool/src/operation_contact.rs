//! Bodies that only touch (PK-756).
//!
//! Two solids meeting on coincident faces whose normals oppose, and
//! otherwise at most at points, have disjoint interiors: near each
//! contact face the bodies lie on opposite sides of it, no surface of
//! one crosses the other (every record is a point or runs along a
//! contact face's rim, where the contact ends), and containment would
//! need the contact faces cooriented.
//! So their intersection is empty and A minus B is A.
//!
//! The full pipeline got this wrong on the cycloid module's output pin
//! standing on its flange: the pin's top circle met the flange face as
//! a zero-length record, splitting a zero-height sliver off the pin
//! wall, sampled ON the flange plane and voted Inside; the cellular
//! assembly declined (a valence-1 vertex) and the graph-merge fallback
//! served the two coincident end discs as a 301 mm³ "overlap" that the
//! interference check reported as a collision.
//!
//! The union of touching bodies merges their contact faces, so it keeps
//! the full pipeline.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::coincidence::types::CoincidenceOrientation;
use crate::error::BoolOperation;
use crate::intersect::IntersectionResult;

/// Samples per record when measuring its length.
const LENGTH_SAMPLES: usize = 8;

/// The result when `a` and `b` only touch (see the module docs):
/// `Some` for an intersection (empty) or a subtraction (`a`), `None`
/// otherwise — a record crossing anywhere but a contact face's rim, any
/// cooriented coincidence, no coincidence at all, or a union.
pub(crate) fn touching_only(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    intersection: &IntersectionResult,
    tolerance: &Tolerance,
) -> Option<BRep> {
    if matches!(op, BoolOperation::Union) {
        return None;
    }
    let overlaps = &intersection.coincidence.face_overlaps;
    if overlaps.is_empty()
        || overlaps
            .iter()
            .any(|o| o.orientation != CoincidenceOrientation::Opposite)
    {
        return None;
    }
    // Every record must be a point or run along the rim of a contact
    // face (the pin wall meets the flange plane only along its own top
    // edge): such a record is where the contact ends, not a crossing.
    let rims: Vec<&popkern_core::geometry::curve::CurveOnEdge> = overlaps
        .iter()
        .flat_map(|o| [(a, o.face_a), (b, o.face_b)])
        .flat_map(|(brep, face)| face_edges(brep, face))
        .collect();
    let near = tolerance.linear() * 10.0;
    let on_rim = |p: &popkern_core::geometry::Point3| {
        rims.iter().any(|c| {
            let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
            c.geometry
                .closest_point(p)
                .map(|t| t.clamp(lo, hi))
                .is_some_and(|t| c.geometry.point_at(t).distance_to(p) <= near)
        })
    };
    let contact_only = intersection.records.iter().all(|r| {
        let c = &r.curve.curve_3d;
        let (t0, t1) = c.parameter_range();
        let samples: Vec<_> = (0..=LENGTH_SAMPLES)
            .map(|i| c.point_at(t0 + (t1 - t0) * i as f64 / LENGTH_SAMPLES as f64))
            .collect();
        let length: f64 = samples.windows(2).map(|w| w[0].distance_to(&w[1])).sum();
        length <= near || samples.iter().all(on_rim)
    });
    if !contact_only {
        return None;
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [exec] -> touching only ({} opposite contact face(s), {} point record(s))",
            overlaps.len(),
            intersection.records.len()
        );
    }
    match op {
        BoolOperation::Intersection => Some(BRep::new(*tolerance)),
        BoolOperation::Subtraction => Some(a.clone()),
        BoolOperation::Union => None,
    }
}

/// The edge curves bounding `face` (every loop).
fn face_edges(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Vec<&popkern_core::geometry::curve::CurveOnEdge> {
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
