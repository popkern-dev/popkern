//! Orientation classification for coincident face pairs.
//!
//! Determines whether two coincident faces have the same or opposite
//! normal orientation. This is critical for correct boolean results:
//! the spec (SS3.3) prescribes different rules for each case.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::face_vertex_positions;

use super::types::CoincidenceOrientation;

/// Classify the orientation relationship between two coincident faces.
///
/// Samples normals at the centroid of each face's vertex positions and
/// computes the dot product. If normals point in the same direction
/// (dot > 0), the faces are `Same`-oriented. If they point in opposite
/// directions (dot < 0), they are `Opposite`-oriented.
///
/// Accounts for the face's `SurfaceOnFace::orientation` flag, which
/// indicates whether the surface normal agrees with the face normal.
pub(crate) fn classify_face_orientation(
    a: &BRep,
    face_a: FaceHandle,
    b: &BRep,
    face_b: FaceHandle,
) -> CoincidenceOrientation {
    let normal_a = face_outward_normal(a, face_a);
    let normal_b = face_outward_normal(b, face_b);

    if normal_a.dot(&normal_b) >= 0.0 {
        CoincidenceOrientation::Same
    } else {
        CoincidenceOrientation::Opposite
    }
}

/// Compute the outward-facing normal of a face at its centroid.
///
/// Takes the surface orientation flag into account: if `orientation`
/// is false, the surface normal is flipped.
fn face_outward_normal(brep: &BRep, face: FaceHandle) -> popkern_core::geometry::Vec3 {
    let sof = match brep.faces[face].surface.as_ref() {
        Some(s) => s,
        None => return popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
    };

    let positions = face_vertex_positions(brep, face);
    if positions.is_empty() {
        return popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0);
    }

    // Find UV centroid
    let ((u_min, u_max), (v_min, v_max)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, &sof.geometry);
    let u_mid = (u_min + u_max) * 0.5;
    let v_mid = (v_min + v_max) * 0.5;

    let mut normal = sof.geometry.normal_at(u_mid, v_mid);

    // Flip if surface orientation disagrees with face orientation
    if !sof.orientation {
        normal = normal * (-1.0);
    }

    normal
}
