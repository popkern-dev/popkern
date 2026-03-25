//! Vertex position recomputation for drafted faces.
//!
//! For each vertex adjacent to both drafted and undrafted faces, computes
//! the new position from the intersection of the relevant surfaces.

use std::collections::HashMap;

use crate::error::DraftError;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Recompute vertex positions affected by drafted face surfaces.
///
/// For each vertex adjacent to drafted faces, finds the adjacent surfaces
/// (using the new drafted surface where available, the original otherwise)
/// and computes the intersection point.
///
/// # Errors
///
/// Returns [`DraftError::Topology`] if adjacency lookups fail.
pub fn recompute_draft_vertices(
    brep: &mut BRep,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
    tolerance: &Tolerance,
) -> Result<(), DraftError> {
    if drafted.is_empty() {
        return Ok(());
    }

    let vertex_handles: Vec<VertexHandle> = brep.vertices.keys().collect();

    for vh in vertex_handles {
        let adjacent_faces = find_adjacent_faces(brep, vh);

        // Check if this vertex is adjacent to any drafted face.
        let touches_drafted = adjacent_faces.iter().any(|fh| drafted.contains_key(fh));
        if !touches_drafted {
            continue;
        }

        // Collect surfaces (drafted version if available, otherwise original).
        let mut planes: Vec<(Vec3, f64)> = Vec::new();
        let mut surfaces: Vec<SurfaceGeometry> = Vec::new();

        for fh in &adjacent_faces {
            let surf = if let Some(drafted_surf) = drafted.get(fh) {
                drafted_surf.clone()
            } else if let Some(face) = brep.faces.get(*fh) {
                if let Some(ref sof) = face.surface {
                    sof.geometry.clone()
                } else {
                    continue;
                }
            } else {
                continue;
            };

            if let SurfaceGeometry::Plane { origin, normal, .. } = &surf {
                let d = normal.dot(&Vec3::new(origin.x, origin.y, origin.z));
                planes.push((*normal, d));
            }
            surfaces.push(surf);
        }

        if surfaces.is_empty() {
            continue;
        }

        // Try 3-plane intersection if we have 3+ planes.
        let new_pos = if planes.len() >= 3 {
            solve_3_planes(&planes[0], &planes[1], &planes[2])
        } else {
            // Newton iteration from current position.
            let current_pos = brep.vertices[vh].position;
            let surf_refs: Vec<&SurfaceGeometry> = surfaces.iter().collect();
            newton_vertex_position(&surf_refs, current_pos, tolerance)
        };

        if let Some(pos) = new_pos {
            // Guard against NaN / infinity from degenerate 3-plane intersections
            // or Newton convergence failures. Keep the original position if the
            // computed point is non-finite.
            if pos.x.is_finite() && pos.y.is_finite() && pos.z.is_finite() {
                brep.vertices[vh].position = pos;
            }
        }
    }

    Ok(())
}

/// Find all faces adjacent to a vertex.
fn find_adjacent_faces(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut faces = Vec::new();
    let start_he = brep.vertices[vertex].outgoing;

    if !brep.half_edges.contains_key(start_he) {
        return faces;
    }

    let mut current = start_he;
    let mut visited = std::collections::HashSet::new();
    loop {
        if !visited.insert(current) {
            break;
        }
        let he = &brep.half_edges[current];
        let loop_h = he.loop_;
        if brep.loops.contains_key(loop_h) {
            let face_h = brep.loops[loop_h].face;
            if !faces.contains(&face_h) {
                faces.push(face_h);
            }
        }

        let twin = he.twin;
        if !brep.half_edges.contains_key(twin) {
            break;
        }
        let next_outgoing = brep.half_edges[twin].next;
        if next_outgoing == start_he || !brep.half_edges.contains_key(next_outgoing) {
            break;
        }
        current = next_outgoing;
    }

    faces
}

/// Solve for the intersection of 3 planes using Cramer's rule.
fn solve_3_planes(
    plane0: &(Vec3, f64),
    plane1: &(Vec3, f64),
    plane2: &(Vec3, f64),
) -> Option<Point3> {
    let (n0, d0) = plane0;
    let (n1, d1) = plane1;
    let (n2, d2) = plane2;

    let cross12 = n1.cross(n2);
    let det = n0.dot(&cross12);

    if det.abs() < 1e-12 {
        return None;
    }

    let cross20 = n2.cross(n0);
    let cross01 = n0.cross(n1);

    let inv_det = 1.0 / det;
    let x = (*d0 * cross12.x + *d1 * cross20.x + *d2 * cross01.x) * inv_det;
    let y = (*d0 * cross12.y + *d1 * cross20.y + *d2 * cross01.y) * inv_det;
    let z = (*d0 * cross12.z + *d1 * cross20.z + *d2 * cross01.z) * inv_det;

    // Guard against NaN / infinity from near-singular systems.
    if !x.is_finite() || !y.is_finite() || !z.is_finite() {
        return None;
    }

    Some(Point3::new(x, y, z))
}

/// Newton iteration to find the point closest to all surfaces.
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

        if correction.length() < tolerance.linear {
            break;
        }
    }

    Some(pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn solve_3_planes_basic() {
        let p0 = (Vec3::new(1.0, 0.0, 0.0), 1.0);
        let p1 = (Vec3::new(0.0, 1.0, 0.0), 2.0);
        let p2 = (Vec3::new(0.0, 0.0, 1.0), 3.0);
        let result = solve_3_planes(&p0, &p1, &p2).unwrap();
        assert!((result.x - 1.0).abs() < 1e-12);
        assert!((result.y - 2.0).abs() < 1e-12);
        assert!((result.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn solve_3_planes_singular() {
        let p0 = (Vec3::new(0.0, 0.0, 1.0), 1.0);
        let p1 = (Vec3::new(0.0, 0.0, 1.0), 2.0);
        let p2 = (Vec3::new(1.0, 0.0, 0.0), 0.0);
        assert!(solve_3_planes(&p0, &p1, &p2).is_none());
    }

    #[test]
    fn newton_converges_for_planes() {
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

    #[test]
    fn recompute_empty_drafted_is_noop() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let drafted = HashMap::new();
        recompute_draft_vertices(&mut brep, &drafted, &tol()).unwrap();
        // Vertices should be unchanged.
        assert_eq!(brep.vertices.len(), 8);
    }

    #[test]
    fn recompute_drafted_cube_face() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        // Create a slightly tilted plane as the drafted surface.
        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.087, 0.996),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let result = recompute_draft_vertices(&mut brep, &drafted, &tol());
        assert!(result.is_ok());
    }
}
