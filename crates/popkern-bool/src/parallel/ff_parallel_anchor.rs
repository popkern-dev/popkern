//! A closed record touching its pair's boundary is anchored at the touch
//! before the pair's faces clip it (W3b, tangency design §4.8).
//!
//! A closed record meeting a face's boundary only where it touches it (the
//! ball-ball circle grazing a ball's rim on a box top) divides the face
//! through the touch; kept as a closed loop whose touch is an interior
//! point, the face split read it as a loop inside the face. The pair's
//! edge-face hits on edges bounding two faces (`FacePairTask::anchors`:
//! exact, only where the operands meet, never a seam, which a record
//! crosses without leaving its face) say where: a
//! closed record with exactly one hit on it, off its closure, touches the
//! boundary there (a record crossing the boundary meets it at two hits or
//! more), and lying in both faces off the touch. The record is rotated to start at that hit, a closed record
//! still, so both faces clip it by their closure rules and the vertex
//! table closes it at the hit's vertex.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::tolerance::Tolerance;

/// `curve` rotated to start at the one hit lying on it, or `None` when it
/// is open, holds no hit or several, or the hit is its closure.
pub(super) fn anchored(
    curve: &CurveOnSurface,
    seeds: &[Point3],
    ((a, face_a), (b, face_b)): (
        (
            &popkern_core::topology::BRep,
            popkern_core::topology::FaceHandle,
        ),
        (
            &popkern_core::topology::BRep,
            popkern_core::topology::FaceHandle,
        ),
    ),
    tolerance: &Tolerance,
) -> Option<CurveOnSurface> {
    if !curve.is_closed || seeds.is_empty() {
        return None;
    }
    let c3 = &curve.curve_3d;
    let (t_min, t_max) = c3.parameter_range();
    let reach = tolerance.linear_reach();
    let mut on: Vec<(f64, Point3)> = Vec::new();
    for p in seeds {
        let t = closest_point_in_window(c3, p, t_min, t_max);
        if c3.point_at(t).distance_to(p) <= reach
            && !on.iter().any(|(_, q)| q.distance_to(p) <= reach)
        {
            on.push((t, *p));
        }
    }
    let [(t, p)] = on[..] else {
        return None;
    };
    if c3.point_at(t_min).distance_to(&p) <= reach {
        return None;
    }
    let span = t_max - t_min;
    // A touch divides a face the record lies in: its point half a turn
    // from the touch inside both faces (a circle meeting a face only at a
    // touch on its boundary from outside, a plane's circle on a ball
    // resting on the box's edge, is outside, and anchored the clip's
    // closure rules would keep it).
    let far = t + 0.5 * span;
    let far = if far > t_max { far - span } else { far };
    let inside = |brep: &popkern_core::topology::BRep,
                  face,
                  uv: &popkern_core::geometry::nurbs_curve::NurbsCurve| {
        crate::split_inner_loop::point_inside_face_uv(brep, face, uv.point_at(far), tolerance)
    };
    if !(inside(a, face_a, &curve.curve_uv1) && inside(b, face_b, &curve.curve_uv2)) {
        return None;
    }
    let rotate = |c: &popkern_core::geometry::nurbs_curve::NurbsCurve| {
        crate::clip::clip_curve_wrapped(c, t, t + span)
    };
    Some(CurveOnSurface {
        curve_3d: rotate(&curve.curve_3d)?,
        curve_uv1: rotate(&curve.curve_uv1)?,
        curve_uv2: rotate(&curve.curve_uv2)?,
        ..curve.clone()
    })
}
