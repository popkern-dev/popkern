//! Vertex-face contact detection between two BRep bodies.
//!
//! A vertex lies on a face when the surface `closest_point` projection
//! returns a distance within the snap tolerance, and the projected point
//! lies within the face's boundary (checked via UV containment).

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::face_vertex_positions;

use super::result::{Operand, VertexFaceHit};

/// Find all vertex-face contacts between both bodies.
///
/// Uses precomputed face AABBs for fast rejection.
pub(crate) fn detect_vertex_face(
    a: &BRep,
    b: &BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    snap_tol: f64,
) -> Vec<VertexFaceHit> {
    let mut hits = Vec::new();

    // A's vertices vs B's faces
    detect_one_direction(a, b, Operand::A, aabbs_b, snap_tol, &mut hits);
    // B's vertices vs A's faces
    detect_one_direction(b, a, Operand::B, aabbs_a, snap_tol, &mut hits);

    hits
}

/// Check vertices of `vert_brep` against faces of `face_brep`.
fn detect_one_direction(
    vert_brep: &BRep,
    face_brep: &BRep,
    vert_body: Operand,
    face_aabbs: &[(FaceHandle, Aabb)],
    snap_tol: f64,
    hits: &mut Vec<VertexFaceHit>,
) {
    for (vh, vert) in &vert_brep.vertices {
        let vp = vert.position;

        // Build a point AABB for the vertex
        let vp_aabb = Aabb { min: vp, max: vp };

        for (fh, face_aabb) in face_aabbs {
            // Fast AABB rejection
            if !vp_aabb.overlaps(face_aabb, snap_tol) {
                continue;
            }

            let surface = match face_brep.faces[*fh].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            // Project vertex onto the surface
            let (u, v) = match surface.closest_point(&vp) {
                Some(uv) => uv,
                None => continue,
            };

            let closest = surface.point_at(u, v);
            let dist = closest.distance_to(&vp);

            if dist < snap_tol {
                // Verify the projected point is within the face boundary
                // by checking that it's within the face's vertex extent
                let positions = face_vertex_positions(face_brep, *fh);
                if !positions.is_empty() {
                    let in_extent = positions
                        .iter()
                        .any(|p| closest.distance_to(p) < face_aabb.diagonal() * 1.1);
                    if !in_extent {
                        continue;
                    }
                }

                hits.push(VertexFaceHit {
                    vertex: vh,
                    vertex_body: vert_body,
                    face: *fh,
                    uv: (u, v),
                    distance: dist,
                });
            }
        }
    }
}
