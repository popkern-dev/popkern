//! Inner skin construction: clone, map, and geometry replacement.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle};

/// Clone the outer BRep to create an inner skin, building a handle map.
///
/// Faces listed in `removed_faces` are deleted from the inner skin
/// (along with their exclusive loops). Shared edges and vertices are kept.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if the clone or removal encounters
/// inconsistent topology.
pub fn build_inner_skin(
    outer: &BRep,
    removed_faces: &[FaceHandle],
) -> Result<(BRep, SkinMap), ShellError> {
    // 1. Clone outer BRep -> inner.
    let inner = outer.clone();

    // 2. Build SkinMap by walking outer and inner slotmaps in parallel.
    //    Since clone preserves insertion order, zip gives matching pairs.
    let vertices: HashMap<_, _> = outer.vertices.keys().zip(inner.vertices.keys()).collect();
    let edges: HashMap<_, _> = outer.edges.keys().zip(inner.edges.keys()).collect();
    let half_edges: HashMap<_, _> = outer
        .half_edges
        .keys()
        .zip(inner.half_edges.keys())
        .collect();
    let faces: HashMap<_, _> = outer.faces.keys().zip(inner.faces.keys()).collect();

    let mut skin_map = SkinMap {
        vertices,
        edges,
        half_edges,
        faces,
    };

    // 3. Remove faces from inner skin.
    let mut inner = inner;
    for &outer_face in removed_faces {
        let inner_face = skin_map
            .faces
            .get(&outer_face)
            .copied()
            .ok_or_else(|| ShellError::Topology("removed face not in skin map".into()))?;

        // Remove the face's loops from the inner BRep.
        if inner.faces.contains_key(inner_face) {
            let face = inner.faces[inner_face].clone();
            // Remove outer loop.
            inner.loops.remove(face.outer_loop);
            // Remove inner loops.
            for loop_handle in &face.inner_loops {
                inner.loops.remove(*loop_handle);
            }
            // Remove the face itself.
            inner.faces.remove(inner_face);
        }

        // Remove from the skin map.
        skin_map.faces.remove(&outer_face);
    }

    Ok((inner, skin_map))
}

/// Replace inner skin face geometry with offset surfaces.
///
/// For each inner face, looks up the corresponding outer face via the
/// reverse skin map, retrieves the offset surface, and replaces the
/// face's surface geometry. Orientation flags are preserved: skins
/// that must mirror their material side do so through
/// [`BRep::reverse_body_orientation`], which flips winding and flag
/// together — a flag-only flip ships the winding-dishonest mirror
/// class (PK-244/PK-245).
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if a face mapping is missing.
pub fn apply_offset_geometry(
    inner: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
) -> Result<(), ShellError> {
    // Build reverse map: inner face -> outer face.
    let reverse = skin_map.reverse_face_map();

    let inner_face_handles: Vec<_> = inner.faces.keys().collect();
    for inner_fh in inner_face_handles {
        let outer_fh = match reverse.get(&inner_fh) {
            Some(&fh) => fh,
            None => continue, // Face may have been removed.
        };

        if let Some(offset_surface) = offset_surfaces.get(&outer_fh) {
            let face = &mut inner.faces[inner_fh];
            let orientation = face
                .surface
                .as_ref()
                .map(|sof| sof.orientation)
                .unwrap_or(true);

            face.surface = Some(SurfaceOnFace {
                geometry: offset_surface.clone(),
                orientation,
            });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn clone_preserves_topology_counts() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        assert_eq!(inner.vertices.len(), 8);
        assert_eq!(inner.edges.len(), 12);
        assert_eq!(inner.faces.len(), 6);
        assert_eq!(skin_map.vertices.len(), 8);
        assert_eq!(skin_map.edges.len(), 12);
        assert_eq!(skin_map.faces.len(), 6);
    }

    #[test]
    fn vertex_pairs_match() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        for (&outer_vh, &inner_vh) in &skin_map.vertices {
            let outer_pos = outer.vertices[outer_vh].position;
            let inner_pos = inner.vertices[inner_vh].position;
            assert!(
                outer_pos.distance_to(&inner_pos) < 1e-12,
                "vertex positions should match after clone"
            );
        }
    }

    #[test]
    fn remove_one_face_reduces_count() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let first_face = outer.faces.keys().next().unwrap();
        let (inner, skin_map) = build_inner_skin(&outer, &[first_face]).unwrap();
        assert_eq!(inner.faces.len(), 5);
        assert_eq!(skin_map.faces.len(), 5);
        // Vertices and edges are kept (shared with other faces).
        assert_eq!(inner.vertices.len(), 8);
        assert_eq!(inner.edges.len(), 12);
    }

    #[test]
    fn half_edge_map_has_correct_count() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (_, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        assert_eq!(skin_map.half_edges.len(), outer.half_edges.len());
    }

    #[test]
    fn apply_offset_geometry_replaces_surfaces() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();

        // Create fake offset surfaces (same as original for test).
        let mut offset_surfaces = HashMap::new();
        for &outer_fh in skin_map.faces.keys() {
            if let Some(sof) = &outer.faces[outer_fh].surface {
                offset_surfaces.insert(outer_fh, sof.geometry.clone());
            }
        }

        apply_offset_geometry(&mut inner, &offset_surfaces, &skin_map).unwrap();

        // Geometry replacement preserves the authored orientation;
        // mirroring is reverse_body_orientation's job.
        for (_, face) in &inner.faces {
            let sof = face.surface.as_ref().unwrap();
            assert!(sof.orientation, "orientation must be preserved");
        }
    }

    #[test]
    fn apply_offset_geometry_uses_offset_surface() {
        use popkern_core::geometry::{Point3, Vec3};
        let outer = BRep::unit_cube(tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();

        // Create offset surfaces that are distinct from the originals
        // (shifted planes).
        let mut offset_surfaces = HashMap::new();
        for &outer_fh in skin_map.faces.keys() {
            let shifted = SurfaceGeometry::Plane {
                origin: Point3::new(99.0, 99.0, 99.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            };
            offset_surfaces.insert(outer_fh, shifted);
        }

        apply_offset_geometry(&mut inner, &offset_surfaces, &skin_map).unwrap();

        // After applying, inner faces should have the offset geometry.
        for (_, face) in &inner.faces {
            let sof = face.surface.as_ref().unwrap();
            match &sof.geometry {
                SurfaceGeometry::Plane { origin, .. } => {
                    assert!(
                        (origin.x - 99.0).abs() < 1e-12,
                        "inner face should use offset surface"
                    );
                }
                _ => panic!("expected Plane"),
            }
        }
    }

    #[test]
    fn apply_with_removed_faces_skips_missing() {
        let outer = BRep::unit_cube(tol()).unwrap();
        let first_face = outer.faces.keys().next().unwrap();
        let (mut inner, skin_map) = build_inner_skin(&outer, &[first_face]).unwrap();

        // Only provide offset surfaces for faces that weren't removed.
        let mut offset_surfaces = HashMap::new();
        for &outer_fh in skin_map.faces.keys() {
            if let Some(sof) = &outer.faces[outer_fh].surface {
                offset_surfaces.insert(outer_fh, sof.geometry.clone());
            }
        }

        // Should not error even though one face was removed.
        apply_offset_geometry(&mut inner, &offset_surfaces, &skin_map).unwrap();
        assert_eq!(inner.faces.len(), 5);
    }

    #[test]
    fn full_pipeline_geometry_replacement() {
        use crate::classify::classify_faces;
        use crate::surface_offset::compute_offset_surfaces;
        use crate::types::ShellDirection;
        use crate::vertex_recompute::recompute_vertices;

        let brep = BRep::unit_cube(tol()).unwrap();
        let classification = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]).unwrap();
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

        apply_offset_geometry(&mut inner, &offset_result.surfaces, &skin_map).unwrap();

        // All inner faces should have offset geometry with preserved
        // orientation (the mirror happens in reverse_body_orientation).
        for (_, face) in &inner.faces {
            let sof = face.surface.as_ref().unwrap();
            assert!(sof.orientation, "orientation must be preserved");
            // Geometry should be a plane (offset cube faces are still planes).
            assert!(
                matches!(sof.geometry, SurfaceGeometry::Plane { .. }),
                "offset cube face should still be a plane"
            );
        }
    }
}
