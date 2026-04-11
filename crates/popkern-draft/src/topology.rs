//! Topology modification orchestrator for draft operations.
//!
//! Coordinates the sequence of operations needed to apply draft changes
//! to a BRep: surface replacement, orientation fixup, edge/vertex
//! recomputation, and adjacent face handling.

use std::collections::HashMap;

use crate::adjacent;
use crate::error::{DraftError, DraftWarning};
use crate::recompute;
use crate::types::DraftOptions;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Apply draft modifications to a BRep clone.
///
/// Clones the input BRep, replaces each drafted face's surface geometry,
/// fixes orientation if the normal flipped, recomputes edges and vertices,
/// and handles adjacent faces.
///
/// # Errors
///
/// Returns [`DraftError`] if topology modification fails (e.g., degenerate
/// geometry or inconsistent adjacency).
pub fn apply_draft_topology(
    brep: &BRep,
    drafted_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    options: &DraftOptions,
    tolerance: &Tolerance,
) -> Result<(BRep, Vec<DraftWarning>), DraftError> {
    let mut result = brep.clone();
    let mut warnings = Vec::new();

    // Step 1: Replace surface geometry on drafted faces and fix orientation.
    for (face_h, new_surf) in drafted_surfaces {
        if let Some(face) = result.faces.get_mut(*face_h)
            && let Some(ref mut sof) = face.surface
        {
            // Check if the draft reverses the normal direction.
            let old_normal = sof.geometry.normal_at(0.0, 0.0);
            let new_normal = new_surf.normal_at(0.0, 0.0);
            let dot = old_normal.dot(&new_normal);
            if dot < 0.0 {
                sof.orientation = !sof.orientation;
            }
            sof.geometry = new_surf.clone();
        }
    }

    // Capture original vertex positions before recomputation so that
    // transition face insertion can compare old vs new positions.
    let original_positions: HashMap<VertexHandle, Point3> = result
        .vertices
        .iter()
        .map(|(vh, v)| (vh, v.position))
        .collect();

    // Step 2: Recompute vertex positions.
    recompute::recompute_draft_vertices(&mut result, drafted_surfaces, tolerance)?;

    // Step 3: Recompute edge curves.
    recompute::recompute_draft_edges(&mut result, drafted_surfaces, tolerance)?;

    // Step 4: Vertex consistency — when a vertex is shared by multiple
    // drafted faces, each face may produce a slightly different position.
    // Snap vertices that are very close to each other to reduce geometric gaps.
    snap_shared_vertices(&mut result, drafted_surfaces, tolerance);

    // Step 5: Handle adjacent faces (transition face insertion if needed).
    let adj_warnings = adjacent::handle_adjacent_faces(
        &mut result,
        drafted_surfaces,
        &original_positions,
        options.adjacent_mode,
        tolerance,
    )?;
    warnings.extend(adj_warnings);

    // Step 6: Recompute pcurves for all faces.
    result
        .recompute_all_pcurves()
        .map_err(|e| DraftError::Topology(e.to_string()))?;

    Ok((result, warnings))
}

/// Snap vertices shared by multiple drafted faces to reduce geometric gaps.
///
/// After vertex recomputation, vertices adjacent to 2+ drafted faces may
/// have been solved independently and ended up at slightly different
/// positions depending on which 3-plane combination was used. This function
/// re-solves each such vertex as the average of projections onto all
/// adjacent surfaces, ensuring consistency.
fn snap_shared_vertices(
    brep: &mut BRep,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
    tolerance: &Tolerance,
) {
    use popkern_core::geometry::Vec3;

    let vertex_handles: Vec<VertexHandle> = brep.vertices.keys().collect();

    for vh in vertex_handles {
        // Find how many drafted faces touch this vertex.
        let adjacent_faces = collect_vertex_faces(brep, vh);
        let drafted_count = adjacent_faces
            .iter()
            .filter(|fh| drafted.contains_key(fh))
            .count();

        if drafted_count < 2 {
            continue;
        }

        // Re-project: average the nearest-surface projections.
        let pos = brep.vertices[vh].position;
        let mut projected_sum = Vec3::zero();
        let mut count = 0;

        for fh in &adjacent_faces {
            let surf = if let Some(s) = drafted.get(fh) {
                s.clone()
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
                let d = (pos - *origin).dot(normal);
                let projected = pos - *normal * d;
                projected_sum = projected_sum + (projected - Point3::new(0.0, 0.0, 0.0));
                count += 1;
            }
        }

        if count >= 2 {
            let avg = Point3::new(
                projected_sum.x / count as f64,
                projected_sum.y / count as f64,
                projected_sum.z / count as f64,
            );
            let delta = pos.distance_to(&avg);
            if delta < tolerance.linear() * 100.0 && delta > tolerance.linear() * 0.01 {
                brep.vertices[vh].position = avg;
            }
        }
    }
}

/// Collect all faces adjacent to a vertex by walking the half-edge fan.
fn collect_vertex_faces(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn draft_one_cube_face_valid_brep() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.087, 0.996),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let (result, warnings) =
            apply_draft_topology(&brep, &drafted, &DraftOptions::default(), &tol()).unwrap();

        // Euler: V - E + F = 2 for a closed solid.
        let v = result.vertices.len() as i64;
        let e = result.edges.len() as i64;
        let f = result.faces.len() as i64;
        assert_eq!(v - e + f, 2, "Euler: V={v} E={e} F={f}");
        assert!(warnings.is_empty());
    }

    #[test]
    fn drafted_face_has_new_surface() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        let new_normal = Vec3::new(0.0, 0.087, 0.996);
        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: new_normal,
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let (result, _) =
            apply_draft_topology(&brep, &drafted, &DraftOptions::default(), &tol()).unwrap();

        // Verify the face now has the drafted surface.
        let result_face = &result.faces[face];
        if let Some(ref sof) = result_face.surface {
            if let SurfaceGeometry::Plane { normal, .. } = &sof.geometry {
                let dot = normal.dot(&new_normal);
                assert!(dot > 0.99, "Drafted face should have new normal");
            } else {
                panic!("Expected plane surface after draft");
            }
        }
    }

    #[test]
    fn empty_drafted_is_clone() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let drafted = HashMap::new();

        let (result, warnings) =
            apply_draft_topology(&brep, &drafted, &DraftOptions::default(), &tol()).unwrap();

        assert_eq!(result.vertices.len(), brep.vertices.len());
        assert_eq!(result.edges.len(), brep.edges.len());
        assert_eq!(result.faces.len(), brep.faces.len());
        assert!(warnings.is_empty());
    }
}
