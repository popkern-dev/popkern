//! Coplanar residual support-family discovery.

use super::*;

pub(super) fn support_family_faces(brep: &BRep, support_face: FaceHandle) -> Vec<FaceHandle> {
    let mut faces = Vec::new();
    faces.push(support_face);
    faces.extend(coplanar_residual_support_faces(brep, support_face));
    faces
}

pub(super) fn coplanar_residual_support_faces(
    brep: &BRep,
    support_face: FaceHandle,
) -> Vec<FaceHandle> {
    let Some(source_face) = brep.faces.get(support_face) else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: source_origin,
        normal: source_normal,
        ..
    }) = source_face
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == support_face || face_data.shell != source_face.shell {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let normals_parallel = source_normal.cross(normal).length() <= tol;
            let same_plane = (*origin - *source_origin).dot(source_normal).abs() <= tol;
            (normals_parallel && same_plane).then_some(face)
        })
        .collect()
}
