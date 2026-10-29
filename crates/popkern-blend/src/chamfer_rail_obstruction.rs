//! Chamfer strips obstructed by other topology of their faces (PK-1003).
//!
//! **The class.** A planar chamfer of a straight edge removes, from
//! each of its two planar faces, the strip between the edge and the
//! rail (the trim line at the setback). When another feature's
//! boundary reaches into that strip (a boss standing on the top face
//! closer to the edge than the setback, as a hole in the face's
//! boundary or as a notch in it), the per-edge trims ignore it: the
//! bevel runs under the boss and the served body is validator-clean
//! but WRONG (measured: a boss 0.05 behind the edge, d = 0.06, the
//! chamfer ADDED 0.0158 of volume; d = 0.1 added 0.005).
//!
//! **The gate.** Before construction, every vertex of either face
//! (every loop) must lie outside the strip, rail included: not at an
//! in-face distance above `tol` and up to the setback plus `tol` on the
//! rail's side while projecting strictly inside the edge's span (a rail
//! through a vertex splits the face degenerately and served the same
//! wrong body). The edge's own ends, and its neighbours' far ends at
//! the strip's ends, are on its boundary and pass. An obstructed edge
//! fails with [`BlendFailure::InconsistentGeometry`] naming the face
//! and the vertex, and the request declines whole: a clean decline
//! where the construction would otherwise lie. Bevel growth (PK-1001)
//! then retries smaller and its sweep guard declines a bevel that
//! would reach the feature.
//!
//! **Scope.** Planar faces and straight rails only (the classes the
//! measured failure lives in); curved faces and rails are not gated
//! here. **Tolerance:** the caller's linear tolerance × 10.

use crate::construct::EdgeGeometry;
use crate::error::BlendFailure;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// The first vertex of `face` inside the strip between the edge
/// `a → b` and the straight `rail`, if any.
fn strip_intruder(
    brep: &BRep,
    face: FaceHandle,
    a: Point3,
    b: Point3,
    rail: &NurbsCurve,
    tol: f64,
) -> Option<Point3> {
    let f = brep.faces.get(face)?;
    if !matches!(f.surface.as_ref()?.geometry, SurfaceGeometry::Plane { .. }) {
        return None;
    }
    let length = a.distance_to(&b);
    let along = (b - a).normalized()?;
    let perpendicular = |p: Point3| {
        let w = p - a;
        w - along * w.dot(&along)
    };
    // A straight rail parallel to the edge: its offset from the edge.
    let (first, last) = (*rail.control_points.first()?, *rail.control_points.last()?);
    let offset = perpendicular(first);
    if rail
        .control_points
        .iter()
        .any(|p| (perpendicular(*p) - offset).length() > tol)
        || (perpendicular(last) - offset).length() > tol
    {
        return None;
    }
    let setback = offset.length();
    let side = offset.normalized()?;
    for loop_ in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(loop_)?.half_edge;
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let h = brep.half_edges.get(he)?;
            let p = brep.vertices.get(h.origin)?.position;
            let t = (p - a).dot(&along);
            let w = perpendicular(p);
            let across = w.dot(&side);
            if t > tol
                && t < length - tol
                && across > tol
                && across <= setback + tol
                && (w - side * across).length() <= tol
            {
                return Some(p);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    None
}

/// Why `geometry`'s strips are obstructed, as a failure for its edge.
pub(crate) fn obstructed_strip(
    brep: &BRep,
    geometry: &EdgeGeometry,
    tolerance: &Tolerance,
) -> Option<BlendFailure> {
    let tol = (tolerance.linear() * 10.0).max(1e-9);
    let a = brep.vertices.get(geometry.analysis.vertex_start)?.position;
    let b = brep.vertices.get(geometry.analysis.vertex_end)?.position;
    for (face, rail) in [
        (geometry.analysis.face_1, &geometry.trim_3d_1),
        (geometry.analysis.face_2, &geometry.trim_3d_2),
    ] {
        if let Some(p) = strip_intruder(brep, face, a, b, rail, tol) {
            return Some(BlendFailure::InconsistentGeometry {
                reason: format!(
                    "chamfer strip on {face:?} is obstructed: another feature's vertex lies \
                     between the edge and its rail"
                ),
                location: Some(p),
            });
        }
    }
    None
}
