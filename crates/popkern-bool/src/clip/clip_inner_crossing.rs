//! Whether a record with no outer-boundary crossing crosses one of the
//! face's holes (PK-1024).
//!
//! The crossing finder is deliberately outer-loop-only: the split
//! stage's hole machinery partitions records that run through holes,
//! and expects them at full span. A record crossing no outer edge then
//! reaches the zero-crossing containment vote, which samples five
//! stations and keeps the record on a majority inside the face region.
//! That vote is for records wholly on one side; a record crossing a
//! hole is not, and when most of it lies in the hole (a pin's cap
//! circle 2.7e-3 into a window wall: an arc of 0.19 inside the face,
//! 10 around) every station lands in the hole and the record dropped,
//! stranding the lens (graph-merge and the complex both found it half
//! registered). Crossing a hole's rim is an exact curve intersection,
//! not a sampling question.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::intersect_curves;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Whether `curve_uv` crosses (or touches) an inner loop of `face` in
/// its chart.
pub(super) fn crosses_inner_loop(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    let surface = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    let curve = CurveGeometry::Nurbs(curve_uv.clone());
    let curve_range = curve_uv.parameter_range();
    crate::clip::face_uv_inner_boundaries(brep, face)
        .iter()
        .flatten()
        .any(|(_, pc)| {
            let (a, b) = pc.parameter_range();
            if !(a.is_finite() && b.is_finite()) {
                return false;
            }
            let rim = match surface {
                Some(s) => pc.to_nurbs_on_surface(s, 1e-3 * tolerance.linear()).0,
                None => pc.to_nurbs(),
            };
            let rim_range = rim.parameter_range();
            intersect_curves(
                &CurveGeometry::Nurbs(rim),
                rim_range,
                &curve,
                curve_range,
                tolerance,
            )
            .is_ok_and(|r| !r.intersections.is_empty())
        })
}

#[cfg(test)]
#[path = "clip_inner_crossing_tests.rs"]
mod tests;
