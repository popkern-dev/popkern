//! Where the operands of a boolean touch, reported in its diagnostics
//! (the tangency foundations' observability, design §7.3).
//!
//! The boolean acts on intersection curves; a touch, where two faces
//! are tangent without crossing, leaves no curve, and whether it was
//! seen is otherwise visible only in the result's topology. Each pair
//! of faces, one from each operand, whose boxes come within the
//! tolerance is searched for the critical points of its signed gap
//! (`find_surface_contacts`) over each face's parameter box, and a
//! contact on both trims is reported with its classification (a touch
//! at a point with its side, branches crossing, a tangent curve, a
//! coincidence, or undecided with the reason).

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::contact_classify::{ContactClassification, ContactScale};
use popkern_core::geometry::intersection::contact_discovery::find_surface_contacts;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

/// Boundary samples per edge for a face's parameter box.
const LOOP_SAMPLES: usize = 8;

/// A place where a face of the first operand touches a face of the
/// second without crossing it.
#[derive(Debug, Clone, PartialEq)]
pub struct OperandContact {
    /// The first operand's face.
    pub face_a: FaceHandle,
    /// The second operand's face.
    pub face_b: FaceHandle,
    /// Where they touch (on the first operand's face).
    pub location: Point3,
    /// The signed gap there, within the tolerance.
    pub gap: f64,
    /// What the surfaces do there, with the evidence.
    pub classification: ContactClassification,
}

/// Every contact of `a`'s faces with `b`'s, in face handle order.
pub(crate) fn operand_contacts(a: &BRep, b: &BRep, tolerance: &Tolerance) -> Vec<OperandContact> {
    let linear = tolerance.linear();
    let boxes_a = crate::face_bounding_boxes(a);
    let boxes_b = crate::face_bounding_boxes(b);
    let mut out = Vec::new();
    for (fa, box_a) in &boxes_a {
        for (fb, box_b) in &boxes_b {
            if !box_a.overlaps(box_b, linear) {
                continue;
            }
            let (Some(sa), Some(sb)) = (surface(a, *fa), surface(b, *fb)) else {
                continue;
            };
            let extent = box_a.diagonal().max(box_b.diagonal());
            let found = find_surface_contacts(
                sa,
                window(a, *fa, sa),
                sb,
                window(b, *fb, sb),
                ContactScale { linear, extent },
            );
            for c in found.contacts {
                if on_trim(a, *fa, sa, c.uv1) && on_trim(b, *fb, sb, c.uv2) {
                    out.push(OperandContact {
                        face_a: *fa,
                        face_b: *fb,
                        location: c.point,
                        gap: c.gap,
                        classification: c.classification,
                    });
                }
            }
        }
    }
    out
}

fn surface(brep: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// The face's parameter box, from samples along its outer loop.
fn window(brep: &BRep, face: FaceHandle, s: &SurfaceGeometry) -> ((f64, f64), (f64, f64)) {
    tight_bounded_parameter_range(
        &brep.face_outer_loop_sample_positions(face, LOOP_SAMPLES),
        s,
    )
}

fn on_trim(brep: &BRep, face: FaceHandle, s: &SurfaceGeometry, uv: (f64, f64)) -> bool {
    brep.face_uv_trim_containment(face, s, uv.0, uv.1)
        .is_some_and(|c| c != FaceUvContainment::Outside)
}

#[cfg(test)]
#[path = "diagnostics_contacts_tests.rs"]
mod tests;
