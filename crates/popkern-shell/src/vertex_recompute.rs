//! Vertex position recomputation for inner skin.

use std::collections::{HashMap, HashSet};

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Recompute all vertex positions in the inner skin.
///
/// For each vertex, finds the adjacent faces, maps them back to the outer
/// BRep, retrieves the corresponding offset surfaces, and computes the
/// new vertex position as the intersection of those surfaces.
///
/// For vertices adjacent to 3+ planes, uses Cramer's rule for an exact
/// solution. Otherwise falls back to a Newton-iteration approach from
/// an initial guess.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if adjacency lookups fail.
pub fn recompute_vertices(
    inner: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
) -> Result<(), ShellError> {
    let reverse_face = skin_map.reverse_face_map();
    let reverse_vertex = skin_map.reverse_vertex_map();

    let inner_vertex_handles: Vec<_> = inner.vertices.keys().collect();

    for inner_vh in inner_vertex_handles {
        // Find adjacent faces by walking outgoing half-edges.
        let adjacent_outer_faces = find_adjacent_outer_faces(inner, inner_vh, &reverse_face);

        // Collect offset surfaces for adjacent faces.
        let mut planes: Vec<(Vec3, f64)> = Vec::new();
        let mut surfaces: Vec<&SurfaceGeometry> = Vec::new();

        for outer_fh in &adjacent_outer_faces {
            if let Some(surf) = offset_surfaces.get(outer_fh) {
                surfaces.push(surf);
                if let SurfaceGeometry::Plane { origin, normal, .. } = surf {
                    let d = normal.dot(&Vec3::new(origin.x, origin.y, origin.z));
                    planes.push((*normal, d));
                }
            }
        }

        if surfaces.is_empty() {
            continue;
        }

        // Try 3-plane intersection if we have 3+ planes.
        let new_pos = if planes.len() >= 3 {
            solve_3_planes(&planes[0], &planes[1], &planes[2])
        } else {
            // Fall back to Newton iteration from initial guess.
            let outer_vh = reverse_vertex
                .get(&inner_vh)
                .ok_or_else(|| ShellError::Topology("vertex not in reverse map".into()))?;
            let outer_pos = outer.vertices[*outer_vh].position;
            let avg_normal = compute_average_normal(&surfaces);
            let guess = outer_pos + avg_normal;
            newton_vertex_position(&surfaces, guess, tolerance)
        };

        if let Some(pos) = new_pos {
            inner.vertices[inner_vh].position = pos;
        }
    }

    Ok(())
}

/// Find outer face handles adjacent to an inner vertex.
fn find_adjacent_outer_faces(
    inner: &BRep,
    vertex: VertexHandle,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
) -> Vec<FaceHandle> {
    let mut faces = HashSet::new();
    let start_he = inner.vertices[vertex].outgoing;

    if !inner.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    // Walk the half-edge fan around the vertex.
    let mut current = start_he;
    loop {
        let he = &inner.half_edges[current];
        let loop_h = he.loop_;
        if inner.loops.contains_key(loop_h) {
            let face_h = inner.loops[loop_h].face;
            if let Some(&outer_fh) = reverse_face.get(&face_h) {
                faces.insert(outer_fh);
            }
        }

        // Move to the next half-edge around the vertex:
        // twin's next gives the next outgoing half-edge from this vertex.
        let twin = he.twin;
        if !inner.half_edges.contains_key(twin) {
            break;
        }
        let next_outgoing = inner.half_edges[twin].next;
        if next_outgoing == start_he || !inner.half_edges.contains_key(next_outgoing) {
            break;
        }
        current = next_outgoing;
    }

    faces.into_iter().collect()
}

/// Solve for the intersection of 3 planes using Cramer's rule.
///
/// Given planes `n_i . x = d_i`, builds the 3x3 system and solves.
/// Returns `None` if the system is singular.
fn solve_3_planes(
    plane0: &(Vec3, f64),
    plane1: &(Vec3, f64),
    plane2: &(Vec3, f64),
) -> Option<Point3> {
    let (n0, d0) = plane0;
    let (n1, d1) = plane1;
    let (n2, d2) = plane2;

    // Compute determinant using triple product: n0 . (n1 x n2).
    let cross12 = n1.cross(n2);
    let det = n0.dot(&cross12);

    // Use a tight numerical threshold for the 3x3 determinant singularity
    // check.  This is a pure linear-algebra guard — not a geometric tolerance —
    // so 1e-12 is appropriate regardless of the tolerance object.
    if det.abs() < 1e-12 {
        return None;
    }

    let cross20 = n2.cross(n0);
    let cross01 = n0.cross(n1);

    let inv_det = 1.0 / det;
    let x = (*d0 * cross12.x + *d1 * cross20.x + *d2 * cross01.x) * inv_det;
    let y = (*d0 * cross12.y + *d1 * cross20.y + *d2 * cross01.y) * inv_det;
    let z = (*d0 * cross12.z + *d1 * cross20.z + *d2 * cross01.z) * inv_det;

    Some(Point3::new(x, y, z))
}

/// Compute the average outward normal direction from offset surfaces.
fn compute_average_normal(surfaces: &[&SurfaceGeometry]) -> Vec3 {
    let mut sum = Vec3::zero();
    for surf in surfaces {
        if let SurfaceGeometry::Plane { normal, .. } = surf {
            sum = sum + *normal;
        }
    }
    let len = sum.length();
    if len > 1e-12 {
        sum * (1.0 / len)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

/// Newton iteration to find the point closest to all surfaces.
///
/// Uses a simple averaging approach: starting from `guess`, project onto
/// each surface and average the projections. This works well for near-planar
/// configurations.
fn newton_vertex_position(
    surfaces: &[&SurfaceGeometry],
    guess: Point3,
    tolerance: &Tolerance,
) -> Option<Point3> {
    let mut pos = guess;
    let max_iter = 20;

    for _ in 0..max_iter {
        let mut sum = Vec3::zero();
        let mut count = 0;

        for surf in surfaces {
            if let SurfaceGeometry::Plane { origin, normal, .. } = surf {
                // Project pos onto the plane.
                let d = pos - *origin;
                let dist = d.dot(normal);
                let projected = pos - *normal * dist;
                let delta = projected - pos;
                sum = sum + delta;
                count += 1;
            }
        }

        if count == 0 {
            break;
        }

        let correction = sum * (1.0 / count as f64);
        pos = pos + correction;

        if correction.length() < tolerance.linear() {
            break;
        }
    }

    Some(pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::classify_faces;
    use crate::inner_skin::build_inner_skin;
    use crate::surface_offset::compute_offset_surfaces;
    use crate::types::ShellDirection;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn cube_vertices_shift_by_thickness() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let thickness = 0.1;
        let classification =
            classify_faces(&brep, thickness, &[], ShellDirection::Inward, &[]).unwrap();
        let offset_result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(&brep, &[]).unwrap();

        recompute_vertices(
            &mut inner,
            &offset_result.surfaces,
            &skin_map,
            &brep,
            &tol(),
        )
        .unwrap();

        // Each vertex of the inner cube should have moved from its original
        // position. Verify all vertices actually changed.
        let reverse_vertex = skin_map.reverse_vertex_map();
        let mut moved_count = 0;
        for (inner_vh, v) in &inner.vertices {
            if let Some(&outer_vh) = reverse_vertex.get(&inner_vh) {
                let outer_pos = brep.vertices[outer_vh].position;
                let dist = outer_pos.distance_to(&v.position);
                if dist > 1e-8 {
                    moved_count += 1;
                }
                // The offset distance along the vertex diagonal of a cube corner
                // (where 3 orthogonal planes meet) should be thickness * sqrt(3).
                // This holds regardless of normal direction.
                let expected_diag = thickness * 3.0_f64.sqrt();
                assert!(
                    (dist - expected_diag).abs() < 0.05,
                    "vertex moved by {dist}, expected ~{expected_diag}"
                );
            }
        }
        assert_eq!(moved_count, 8, "all 8 vertices should have moved");
    }

    #[test]
    fn solve_3_planes_basic() {
        // Three coordinate planes offset by 1.
        let p0 = (Vec3::new(1.0, 0.0, 0.0), 1.0); // x = 1
        let p1 = (Vec3::new(0.0, 1.0, 0.0), 2.0); // y = 2
        let p2 = (Vec3::new(0.0, 0.0, 1.0), 3.0); // z = 3
        let result = solve_3_planes(&p0, &p1, &p2).unwrap();
        assert!((result.x - 1.0).abs() < 1e-12);
        assert!((result.y - 2.0).abs() < 1e-12);
        assert!((result.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn solve_3_planes_singular() {
        // Two parallel planes + one perpendicular → singular.
        let p0 = (Vec3::new(0.0, 0.0, 1.0), 1.0);
        let p1 = (Vec3::new(0.0, 0.0, 1.0), 2.0);
        let p2 = (Vec3::new(1.0, 0.0, 0.0), 0.0);
        let result = solve_3_planes(&p0, &p1, &p2);
        assert!(result.is_none());
    }

    #[test]
    fn newton_vertex_converges_for_planes() {
        let plane_x = SurfaceGeometry::Plane {
            origin: Point3::new(1.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let plane_y = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 2.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let surfaces: Vec<&SurfaceGeometry> = vec![&plane_x, &plane_y];
        let guess = Point3::new(0.0, 0.0, 0.0);
        let result = newton_vertex_position(&surfaces, guess, &tol()).unwrap();
        assert!((result.x - 1.0).abs() < 1e-4);
        assert!((result.y - 2.0).abs() < 1e-4);
    }
}
