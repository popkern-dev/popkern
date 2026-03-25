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

/// Replace inner skin face geometry with offset surfaces and flip normals.
///
/// For each inner face, looks up the corresponding outer face via the
/// reverse skin map, retrieves the offset surface, and replaces the
/// face's surface geometry. The orientation is flipped so the inner
/// skin normal points inward.
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
            // Get the current orientation, then flip it.
            let current_orientation = face
                .surface
                .as_ref()
                .map(|sof| sof.orientation)
                .unwrap_or(true);

            face.surface = Some(SurfaceOnFace {
                geometry: offset_surface.clone(),
                orientation: !current_orientation,
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
        for (&outer_fh, _) in &skin_map.faces {
            if let Some(sof) = &outer.faces[outer_fh].surface {
                offset_surfaces.insert(outer_fh, sof.geometry.clone());
            }
        }

        apply_offset_geometry(&mut inner, &offset_surfaces, &skin_map).unwrap();

        // After applying, all inner faces should have flipped orientation.
        for (_, face) in &inner.faces {
            let sof = face.surface.as_ref().unwrap();
            // Original cube has orientation=true, after flip it should be false.
            assert!(!sof.orientation, "orientation should be flipped");
        }
    }
}
