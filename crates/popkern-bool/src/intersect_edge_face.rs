//! Edge-face intersection detection (diagnostic only).
//!
//! Split out of `intersect.rs` to keep the parent under the 500-line
//! structural cap. These hits are informational and not yet integrated
//! into the split pipeline.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use super::intersect::bounded_parameter_range;

/// A detected edge-face intersection point.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields used for diagnostic inspection and future integration.
pub(crate) struct EdgeFaceHit {
    /// The edge handle from the source solid.
    pub edge: popkern_core::topology::EdgeHandle,
    /// The face handle from the target solid.
    pub face: FaceHandle,
    /// The 3D intersection point.
    pub point: popkern_core::geometry::Point3,
    /// Parameter t on the edge curve at the intersection.
    pub t_edge: f64,
    /// UV parameters (u, v) on the face surface at the intersection.
    pub uv: (f64, f64),
}

/// Detect edge-face intersections between two solids (diagnostic only).
///
/// For each edge of solid A, checks if it intersects any face of solid B,
/// and vice versa. Uses `intersect_surface_curve` with the edge's curve
/// geometry and the face's surface geometry.
///
/// Returns a list of [`EdgeFaceHit`] records. These are informational and
/// are not yet integrated into the split pipeline.
pub(crate) fn detect_edge_face_intersections(
    a: &BRep,
    b: &BRep,
    aabbs_b: &[(FaceHandle, Aabb)],
    aabbs_a: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Vec<EdgeFaceHit> {
    let mut hits = Vec::new();

    // Edges of A vs faces of B.
    detect_edges_vs_faces(a, b, aabbs_b, tolerance, &mut hits);

    // Edges of B vs faces of A.
    detect_edges_vs_faces(b, a, aabbs_a, tolerance, &mut hits);

    if !hits.is_empty() {
        #[cfg(debug_assertions)]
        eprintln!(
            "[intersect] detected {} edge-face intersection(s)",
            hits.len()
        );
    }

    hits
}

/// Check edges of `edge_brep` against faces of `face_brep`.
fn detect_edges_vs_faces(
    edge_brep: &BRep,
    face_brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
    hits: &mut Vec<EdgeFaceHit>,
) {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::intersection::{SurfaceCurveIntersection, intersect_surface_curve};

    use popkern_core::geometry::curve::CurveGeometry;

    for (edge_h, edge) in &edge_brep.edges {
        // Get endpoint positions from vertices (always up-to-date, even if
        // the edge curve geometry is stale after vertex translation).
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

        // Construct a line curve from vertex positions (handles cases where
        // stored CurveOnEdge geometry is out of sync with vertex positions).
        let direction = v_end - v_start;
        if direction.length_squared() < tolerance.numerical_floor().powi(2) {
            continue;
        }
        let edge_curve = CurveGeometry::Line {
            origin: v_start,
            direction,
        };
        let t_range = (0.0, 1.0);

        for (face_h, face_aabb) in face_aabbs {
            if !edge_aabb.overlaps(face_aabb, tolerance.linear()) {
                continue;
            }

            let face = &face_brep.faces[*face_h];
            let surface_on_face = match face.surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let ((u_min, u_max), (v_min, v_max)) =
                bounded_parameter_range(face_brep, *face_h, surface_on_face);

            let result = intersect_surface_curve(
                surface_on_face,
                (u_min, u_max),
                (v_min, v_max),
                &edge_curve,
                t_range,
                tolerance,
            );

            if let Ok(result) = result {
                for isect in &result.intersections {
                    if let SurfaceCurveIntersection::Point { t, u, v, .. } = *isect {
                        // Filter out intersections at edge endpoints (these are
                        // vertex-face touches, not edge-face piercings).
                        let eps = tolerance.linear();
                        if (t - t_range.0).abs() < eps || (t - t_range.1).abs() < eps {
                            continue;
                        }
                        hits.push(EdgeFaceHit {
                            edge: edge_h,
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
