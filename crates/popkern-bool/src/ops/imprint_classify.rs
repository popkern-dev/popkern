//! Where each face of an imprinted target lies relative to the tool
//! (design 5.1, 6): inside, outside, or on one of the tool's faces, and
//! then whether the two face the same way (a flush, coplanar overlap) or
//! opposite ways (a contact patch).
//!
//! Decided per face at a point strictly inside it: on the tool's
//! boundary geometrically first (the ray vote has no on-boundary
//! verdict of its own), then by point classification. A face the
//! imprint did not split wholly on one side has one class by
//! construction: every imprinted face is bounded by the tool's curves.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::PointClassification;

/// Where an imprinted face lies relative to the tool.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FaceClass {
    /// Inside the tool.
    Inside,
    /// Outside the tool.
    Outside,
    /// On a tool face that faces the same way (a flush, coplanar
    /// overlap).
    OnAligned,
    /// On a tool face that faces the other way (a contact patch).
    OnOpposed,
}

/// The class of every face of `target` relative to `tool`; a face with
/// no point to probe or no verdict is left out.
pub(super) fn classify_faces(
    target: &BRep,
    tool: &BRep,
    tolerance: &Tolerance,
) -> Vec<(FaceHandle, FaceClass)> {
    let boxes = crate::intersect::compute_face_aabbs(tool);
    target
        .faces
        .keys()
        .filter_map(|face| {
            let p = crate::operation::contact::face_interior_point(target, face)?;
            Some((face, class_at(target, face, tool, &boxes, &p, tolerance)?))
        })
        .collect()
}

fn class_at(
    target: &BRep,
    face: FaceHandle,
    tool: &BRep,
    boxes: &[(
        FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    p: &Point3,
    tolerance: &Tolerance,
) -> Option<FaceClass> {
    if crate::postprocess::zero_thickness::on_body_boundary(p, tool, boxes, tolerance) {
        let mine = outward_normal(target, face, p)?;
        let theirs = tool
            .faces
            .keys()
            .find_map(|f| {
                on_face(tool, f, p, tolerance.linear()).then(|| outward_normal(tool, f, p))
            })
            .flatten()?;
        return Some(if mine.dot(&theirs) >= 0.0 {
            FaceClass::OnAligned
        } else {
            FaceClass::OnOpposed
        });
    }
    match crate::classify::point_in_solid_final(p, tool, boxes, tolerance) {
        Ok(PointClassification::Inside) => Some(FaceClass::Inside),
        Ok(PointClassification::Outside) => Some(FaceClass::Outside),
        _ => None,
    }
}

/// Whether `p` lies on `face` of `body`: within ten tolerances of its
/// surface and inside its trim.
fn on_face(body: &BRep, face: FaceHandle, p: &Point3, linear: f64) -> bool {
    let Some(surf) = body.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some((u, v)) = surf.geometry.closest_point(p) else {
        return false;
    };
    surf.geometry.point_at(u, v).distance_to(p) <= linear * 10.0
        && crate::classify_uv_loops::point_in_face_region_exact(body, face, u, v) == Some(true)
}

/// The face's outward normal (its orientation flag applied) at `p`.
fn outward_normal(
    body: &BRep,
    face: FaceHandle,
    p: &Point3,
) -> Option<popkern_core::geometry::Vec3> {
    let surf = body.faces.get(face)?.surface.as_ref()?;
    let (u, v) = surf.geometry.closest_point(p)?;
    let n = surf.geometry.normal_at(u, v);
    Some(if surf.orientation { n } else { n * -1.0 })
}

#[cfg(test)]
#[path = "imprint_classify_tests.rs"]
mod tests;
