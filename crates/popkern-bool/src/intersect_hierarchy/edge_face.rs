//! Edge-face intersection (pierce point) detection between two BRep bodies.
//!
//! An edge pierces a face when the edge's curve intersects the face's surface
//! at a point interior to both the edge (not at endpoints) and the face.
//! Uses `intersect_surface_curve` from popkern-core.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::intersection::{SurfaceCurveIntersection, intersect_surface_curve};
use popkern_core::geometry::{CurveEval, curve::CurveGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::clip::face_vertex_positions;

use super::result::{EdgeFaceHit, Operand};

/// Find all edge-face intersections between both bodies.
pub(crate) fn detect_edge_face(
    a: &BRep,
    b: &BRep,
    aabbs_a: &[(popkern_core::topology::FaceHandle, Aabb)],
    aabbs_b: &[(popkern_core::topology::FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Vec<EdgeFaceHit> {
    let mut hits = Vec::new();

    // A's edges vs B's faces
    detect_one_direction(a, b, Operand::A, aabbs_b, tolerance, &mut hits);
    // B's edges vs A's faces
    detect_one_direction(b, a, Operand::B, aabbs_a, tolerance, &mut hits);

    hits
}

/// Check edges of `edge_brep` against faces of `face_brep`.
fn detect_one_direction(
    edge_brep: &BRep,
    face_brep: &BRep,
    edge_body: Operand,
    face_aabbs: &[(popkern_core::topology::FaceHandle, Aabb)],
    tolerance: &Tolerance,
    hits: &mut Vec<EdgeFaceHit>,
) {
    for (edge_h, edge) in &edge_brep.edges {
        // Get endpoint positions
        let he = match edge_brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let v_start = match edge_brep.vertices.get(he.origin) {
            Some(v) => v.position,
            None => continue,
        };
        let twin = match edge_brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };
        let v_end = match edge_brep.vertices.get(twin.origin) {
            Some(v) => v.position,
            None => continue,
        };

        let edge_aabb = match Aabb::from_points(&[v_start, v_end]) {
            Some(a) => a,
            None => continue,
        };

        // Construct a line curve from vertex positions
        let direction = v_end - v_start;
        if direction.length_squared() < 1e-30 {
            continue;
        }
        let edge_curve = CurveGeometry::Line {
            origin: v_start,
            direction,
        };
        let t_range = (0.0, 1.0);

        for (face_h, face_aabb) in face_aabbs {
            if !edge_aabb.overlaps(face_aabb, tolerance.linear) {
                continue;
            }

            let surface = match face_brep.faces[*face_h].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let positions = face_vertex_positions(face_brep, *face_h);
            let ((u_min, u_max), (v_min, v_max)) =
                popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, surface);

            let result = intersect_surface_curve(
                surface,
                (u_min, u_max),
                (v_min, v_max),
                &edge_curve,
                t_range,
                tolerance,
            );

            if let Ok(result) = result {
                for isect in &result.intersections {
                    if let SurfaceCurveIntersection::Point { t, u, v, .. } = *isect {
                        // Filter endpoints (VV/VE hits, not EF piercings)
                        let eps = tolerance.linear;
                        if (t - t_range.0).abs() < eps || (t - t_range.1).abs() < eps {
                            continue;
                        }
                        hits.push(EdgeFaceHit {
                            edge: edge_h,
                            edge_body,
                            face: *face_h,
                            point: edge_curve.point_at(t),
                            t_edge: t,
                            uv: (u, v),
                        });
                    }
                }
            }
        }
    }
}
