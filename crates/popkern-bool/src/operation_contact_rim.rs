//! Whether a contact face's rim really ends the contact (PK-850).
//!
//! `touching_only` reads a record that runs along a contact face's rim
//! as the place where the contact ends. That holds only when the faces
//! beyond the rim turn AWAY from the other body: a pin standing on a
//! flange meets the flange plane along its own top edge, and the pin
//! wall rises away from the flange. A peg through the contact plane
//! makes the same record (the plane's hole rim is the peg's wall
//! crossing it), but the wall beyond the rim dips INTO the other body:
//! a block under a box with a peg 0.2 into it read as touching, and the
//! subtraction served the block untouched (the peg's 0.63 of overlap
//! lost), the intersection served empty.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use crate::error::PointClassification;

/// How far into a rim's neighbour face the probe steps, as a fraction
/// of the rim edge's chord length: small enough to stay beside the
/// rim, large enough to clear the tolerance band.
const STEP: f64 = 1e-3;

/// Whether some face beside a rim of `face` (one of `x`'s contact faces)
/// reaches into `other`: a point just past the rim on the neighbour
/// face classifies strictly inside `other`.
pub(crate) fn rim_neighbour_enters(
    x: &BRep,
    face: FaceHandle,
    other: &BRep,
    tolerance: &Tolerance,
) -> bool {
    let Some(f) = x.faces.get(face) else {
        return false;
    };
    let boxes = crate::intersect::compute_face_aabbs(other);
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = x.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=x.half_edges.len() {
            let Some(h) = x.half_edges.get(he) else {
                break;
            };
            if let Some(p) = beyond_rim(x, he, tolerance)
                && matches!(
                    crate::classify::point_in_solid_final(&p, other, &boxes, tolerance),
                    Ok(PointClassification::Inside)
                )
            {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!("  [exec] contact rim neighbour enters the other body at {p:?}");
                }
                return true;
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    false
}

/// A point on the face across half-edge `he`'s edge, a small step from
/// the edge's midpoint, strictly inside that face's trim.
fn beyond_rim(
    x: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    tolerance: &Tolerance,
) -> Option<Point3> {
    let h = x.half_edges.get(he)?;
    let neighbour = x.loops.get(x.half_edges.get(h.twin)?.loop_)?.face;
    let surface = &x.faces.get(neighbour)?.surface.as_ref()?.geometry;
    let c = x.edges.get(h.edge)?.curve.as_ref()?;
    let tm = 0.5 * (c.t_start + c.t_end);
    let m = c.geometry.point_at(tm);
    let chord = c
        .geometry
        .point_at(c.t_start)
        .distance_to(&c.geometry.point_at(c.t_end));
    let length = if chord > tolerance.linear() {
        chord
    } else {
        // A closed rim: its quarter-span chord stands in.
        let tq = c.t_start + 0.25 * (c.t_end - c.t_start);
        c.geometry.point_at(tq).distance_to(&m)
    };
    let step = (STEP * length).max(100.0 * tolerance.linear());
    let tangent = c.geometry.tangent_at(tm).normalized()?;
    let (u, v) = project_point_to_uv(surface, &m)?;
    let normal = surface.normal_at(u, v).normalized()?;
    let across = normal.cross(&tangent).normalized()?;
    for side in [1.0, -1.0] {
        let guess = m + across * (side * step);
        let Some((pu, pv)) = project_point_to_uv(surface, &guess) else {
            continue;
        };
        if x.face_uv_trim_containment(neighbour, surface, pu, pv)
            == Some(FaceUvContainment::Interior)
        {
            return Some(surface.point_at(pu, pv));
        }
    }
    None
}
