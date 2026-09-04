//! Orientation classification for coincident face pairs.
//!
//! Determines whether two coincident faces have the same or opposite
//! normal orientation. This is critical for correct boolean results:
//! the spec (SS3.3) prescribes different rules for each case.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::face_vertex_positions;
use crate::error::PointClassification;

use super::types::CoincidenceOrientation;

/// Classify the orientation relationship between two coincident faces.
///
/// Computes each face's material-aware outward normal (see
/// [`face_outward_normal`]) and compares them: normals agreeing (dot ≥ 0)
/// are `Same`-oriented, opposing (dot < 0) are `Opposite`-oriented. The
/// outward direction is resolved by which side the operand's solid
/// material lies on, so it is correct regardless of how the cap plane's
/// normal and `orientation` flag were stored at construction.
pub(crate) fn classify_face_orientation(
    a: &BRep,
    face_a: FaceHandle,
    aabbs_a: &[(FaceHandle, Aabb)],
    b: &BRep,
    face_b: FaceHandle,
    aabbs_b: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> CoincidenceOrientation {
    let normal_a = face_outward_normal(a, face_a, aabbs_a, tolerance);
    let normal_b = face_outward_normal(b, face_b, aabbs_b, tolerance);

    if normal_a.dot(&normal_b) >= 0.0 {
        CoincidenceOrientation::Same
    } else {
        CoincidenceOrientation::Opposite
    }
}

/// Compute the outward-facing normal of a face at its centroid.
///
/// The stored surface normal and `SurfaceOnFace::orientation` flag do
/// not reliably encode the OUTWARD (away-from-solid) direction for
/// every face: a circle-extrude stores both of a cylinder's cap planes
/// with the same `+n` normal and `orientation = true`, so the bottom
/// cap — whose true outward normal is `-n` — would read `+n` (PK-435:
/// this flipped a bottom disc/rectangle coincident pair to `Opposite`
/// and dropped the kept union skin). Instead, take the plane normal's
/// axis and orient it by MATERIAL SIDE: the outward normal is the one
/// under which a point just off the face centroid lies OUTSIDE the
/// solid. This is construction-independent and correct for any closed
/// operand. Falls back to the stored (flag-adjusted) normal only when
/// the material probe is inconclusive (e.g. an open-sheet operand).
fn face_outward_normal(
    brep: &BRep,
    face: FaceHandle,
    aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Vec3 {
    let sof = match brep.faces[face].surface.as_ref() {
        Some(s) => s,
        None => return Vec3::new(0.0, 0.0, 1.0),
    };

    let positions = face_vertex_positions(brep, face);
    if positions.is_empty() {
        return Vec3::new(0.0, 0.0, 1.0);
    }

    // Find UV centroid
    let ((u_min, u_max), (v_min, v_max)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, &sof.geometry);
    let u_mid = (u_min + u_max) * 0.5;
    let v_mid = (v_min + v_max) * 0.5;

    let mut normal = sof.geometry.normal_at(u_mid, v_mid);
    if !sof.orientation {
        normal = normal * (-1.0);
    }
    let axis_len = normal.length();
    if axis_len < 1.0e-12 {
        return normal;
    }
    let axis = normal * (1.0 / axis_len);

    // Material-side probe from the 3D centroid of the face boundary.
    let n = positions.len() as f64;
    let (mut cx, mut cy, mut cz) = (0.0, 0.0, 0.0);
    for p in &positions {
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    let centroid = Point3::new(cx / n, cy / n, cz / n);
    let eps = (tolerance.linear() * 100.0).max(1.0e-5);
    let probe_plus = Point3::new(
        centroid.x + axis.x * eps,
        centroid.y + axis.y * eps,
        centroid.z + axis.z * eps,
    );
    let probe_minus = Point3::new(
        centroid.x - axis.x * eps,
        centroid.y - axis.y * eps,
        centroid.z - axis.z * eps,
    );
    let plus_out = matches!(
        crate::classify::point_in_solid(&probe_plus, brep, aabbs, tolerance),
        Ok(PointClassification::Outside)
    );
    let minus_out = matches!(
        crate::classify::point_in_solid(&probe_minus, brep, aabbs, tolerance),
        Ok(PointClassification::Outside)
    );
    match (plus_out, minus_out) {
        (true, false) => axis,  // material on -axis side; outward is +axis
        (false, true) => -axis, // material on +axis side; outward is -axis
        _ => normal,            // inconclusive — keep the flag-adjusted normal
    }
}
