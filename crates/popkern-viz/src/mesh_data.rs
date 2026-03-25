//! Shared mesh data types and extraction for GPU and CPU rendering paths.
//!
//! [`MeshVertex`] is the canonical per-vertex representation used by both
//! the GPU renderer and the CPU software renderer.  [`build_mesh_vertices`]
//! converts a tessellated [`popkern_mesh::TriMesh`] (with face-map) into
//! vertex + index arrays, and [`extract_edge_pairs`] pulls position pairs
//! from BRep edges for wireframe drawing.

use popkern_core::topology::BRep;
use popkern_mesh::TriMesh;

use crate::colors::surface_color;

// ---------------------------------------------------------------------------
// Vertex format
// ---------------------------------------------------------------------------

/// Per-vertex data for mesh rendering (GPU or CPU).
///
/// Layout: position (vec3), normal (vec3), color (vec4) = 40 bytes.
/// This struct is GPU-layout-compatible but does **not** depend on wgpu.
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct MeshVertex {
    /// World-space position.
    pub position: [f32; 3],
    /// World-space normal (should be unit-length).
    pub normal: [f32; 3],
    /// RGBA color.
    pub color: [f32; 4],
}

// ---------------------------------------------------------------------------
// Mesh vertex construction
// ---------------------------------------------------------------------------

/// Default face color used when a face has no surface geometry.
const DEFAULT_COLOR: [f32; 4] = [0.6, 0.6, 0.6, 1.0];

/// Build vertex and index arrays from a BRep and its tessellation.
///
/// For each triangle in `mesh`, the originating BRep face is looked up
/// via [`TriMesh::face_map`] and the face's surface type is mapped to a
/// color using [`crate::colors::surface_color`].
///
/// Returns `(vertices, indices)` with three vertices per triangle.
pub fn build_mesh_vertices(brep: &BRep, mesh: &TriMesh) -> (Vec<MeshVertex>, Vec<u32>) {
    let face_color = |face_idx: usize| -> [f32; 4] {
        let face_handle = mesh.face_map.get(face_idx).copied();
        let color = face_handle.and_then(|fh| {
            let face = brep.faces.get(fh)?;
            let surf = face.surface.as_ref()?;
            Some(surface_color(&surf.geometry))
        });
        color.unwrap_or(DEFAULT_COLOR)
    };

    let mut vertices: Vec<MeshVertex> = Vec::with_capacity(mesh.triangles.len() * 3);
    let mut indices: Vec<u32> = Vec::with_capacity(mesh.triangles.len() * 3);

    for (tri_idx, tri) in mesh.triangles.iter().enumerate() {
        let color = face_color(tri_idx);
        for &vi in tri {
            let mv = &mesh.vertices[vi];
            let idx = vertices.len() as u32;
            vertices.push(MeshVertex {
                position: [
                    mv.position.x as f32,
                    mv.position.y as f32,
                    mv.position.z as f32,
                ],
                normal: [mv.normal.x as f32, mv.normal.y as f32, mv.normal.z as f32],
                color,
            });
            indices.push(idx);
        }
    }

    (vertices, indices)
}

// ---------------------------------------------------------------------------
// Edge extraction
// ---------------------------------------------------------------------------

/// Extract edge position pairs from a BRep for wireframe rendering.
///
/// Walks all edges, resolves both endpoint vertex positions, and returns
/// one `([f32; 3], [f32; 3])` pair per edge.  Edges whose half-edge
/// topology is incomplete are silently skipped.
pub fn extract_edge_pairs(brep: &BRep) -> Vec<([f32; 3], [f32; 3])> {
    let mut pairs: Vec<([f32; 3], [f32; 3])> = Vec::with_capacity(brep.edges.len());

    for (_edge_handle, edge) in &brep.edges {
        let he_a_handle = edge.half_edge;
        let Some(he_a) = brep.half_edges.get(he_a_handle) else {
            continue;
        };
        let he_b_handle = he_a.twin;
        let Some(he_b) = brep.half_edges.get(he_b_handle) else {
            continue;
        };

        let Some(va) = brep.vertices.get(he_a.origin) else {
            continue;
        };
        let Some(vb) = brep.vertices.get(he_b.origin) else {
            continue;
        };

        let pa = [
            va.position.x as f32,
            va.position.y as f32,
            va.position.z as f32,
        ];
        let pb = [
            vb.position.x as f32,
            vb.position.y as f32,
            vb.position.z as f32,
        ];

        pairs.push((pa, pb));
    }

    pairs
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_mesh::{TessellationOptions, tessellate};

    #[test]
    fn build_vertices_from_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
        let (vertices, indices) = build_mesh_vertices(&brep, &mesh);
        assert!(!vertices.is_empty());
        assert!(!indices.is_empty());
        assert_eq!(vertices.len(), mesh.triangles.len() * 3);
        assert_eq!(indices.len(), mesh.triangles.len() * 3);
        for v in &vertices {
            assert_eq!(v.color[3], 1.0);
        }
    }

    #[test]
    fn extract_edges_from_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let edges = extract_edge_pairs(&brep);
        assert_eq!(edges.len(), 12);
        for (p0, p1) in &edges {
            let dist =
                ((p0[0] - p1[0]).powi(2) + (p0[1] - p1[1]).powi(2) + (p0[2] - p1[2]).powi(2))
                    .sqrt();
            assert!(dist > 0.001);
        }
    }

    #[test]
    fn empty_brep_produces_no_vertices() {
        let brep = BRep::new(Tolerance::default());
        let mesh = popkern_mesh::TriMesh::default();
        let (vertices, indices) = build_mesh_vertices(&brep, &mesh);
        assert!(vertices.is_empty());
        assert!(indices.is_empty());
    }

    #[test]
    fn empty_brep_produces_no_edges() {
        let brep = BRep::new(Tolerance::default());
        let edges = extract_edge_pairs(&brep);
        assert!(edges.is_empty());
    }

    #[test]
    fn mesh_vertex_size_is_40_bytes() {
        assert_eq!(std::mem::size_of::<MeshVertex>(), 40);
    }

    #[test]
    fn mesh_vertex_is_pod() {
        let v = MeshVertex {
            position: [1.0, 2.0, 3.0],
            normal: [0.0, 1.0, 0.0],
            color: [1.0, 0.0, 0.0, 1.0],
        };
        let bytes = bytemuck::bytes_of(&v);
        assert_eq!(bytes.len(), 40);
        let v2: &MeshVertex = bytemuck::from_bytes(bytes);
        assert_eq!(v.position, v2.position);
        assert_eq!(v.normal, v2.normal);
        assert_eq!(v.color, v2.color);
    }

    #[test]
    fn default_color_is_opaque_gray() {
        assert!((DEFAULT_COLOR[0] - 0.6).abs() < f32::EPSILON);
        assert!((DEFAULT_COLOR[3] - 1.0).abs() < f32::EPSILON);
    }
}
