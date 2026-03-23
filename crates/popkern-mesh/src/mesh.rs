//! Mesh data types for tessellated BRep geometry.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::FaceHandle;

/// A single vertex in a tessellated mesh, carrying position, normal, and UV coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshVertex {
    /// World-space position of the vertex.
    pub position: Point3,
    /// Surface normal at the vertex.
    pub normal: Vec3,
    /// UV (texture) coordinates at the vertex.
    pub uv: (f64, f64),
}

/// The tessellated mesh for a single BRep face.
///
/// Triangles index into the local `vertices` buffer using zero-based indices.
#[derive(Debug, Clone)]
pub struct FaceMesh {
    /// The BRep face this mesh was generated from.
    pub face: FaceHandle,
    /// Vertices local to this face mesh.
    pub vertices: Vec<MeshVertex>,
    /// Triangles as triples of indices into `vertices`.
    pub triangles: Vec<[usize; 3]>,
}

/// A combined triangle mesh for an entire BRep solid or shell.
///
/// All face meshes are merged into a single flat vertex buffer with reindexed triangles.
/// Vertices are **not** deduplicated across faces.
#[derive(Debug, Clone, Default)]
pub struct TriMesh {
    /// All vertices from every contributing face mesh, in face order.
    pub vertices: Vec<MeshVertex>,
    /// Triangles as triples of indices into `vertices`.
    pub triangles: Vec<[usize; 3]>,
    /// For each entry in `triangles`, the originating face handle.
    ///
    /// `face_map[i]` gives the [`FaceHandle`] for `triangles[i]`.
    pub face_map: Vec<FaceHandle>,
    /// The individual per-face meshes that were merged to form this combined mesh.
    pub face_meshes: Vec<FaceMesh>,
}

impl TriMesh {
    /// Merge a collection of [`FaceMesh`] values into a single combined [`TriMesh`].
    ///
    /// Vertices are concatenated in order and triangle indices are offset accordingly.
    /// Vertices are **not** deduplicated across face boundaries.
    pub fn from_face_meshes(face_meshes: Vec<FaceMesh>) -> Self {
        let mut vertices: Vec<MeshVertex> = Vec::new();
        let mut triangles: Vec<[usize; 3]> = Vec::new();
        let mut face_map: Vec<FaceHandle> = Vec::new();

        for fm in &face_meshes {
            let offset = vertices.len();
            vertices.extend(fm.vertices.iter().cloned());
            for tri in &fm.triangles {
                triangles.push([tri[0] + offset, tri[1] + offset, tri[2] + offset]);
                face_map.push(fm.face);
            }
        }

        TriMesh {
            vertices,
            triangles,
            face_map,
            face_meshes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_vertex(x: f64, y: f64, z: f64) -> MeshVertex {
        MeshVertex {
            position: Point3::new(x, y, z),
            normal: Vec3::new(0.0, 0.0, 1.0),
            uv: (0.0, 0.0),
        }
    }

    #[test]
    fn from_face_meshes_empty_input_gives_empty_mesh() {
        let mesh = TriMesh::from_face_meshes(vec![]);
        assert!(mesh.vertices.is_empty());
        assert!(mesh.triangles.is_empty());
        assert!(mesh.face_map.is_empty());
        assert!(mesh.face_meshes.is_empty());
    }

    #[test]
    fn from_face_meshes_single_face_correct_counts() {
        let fm = FaceMesh {
            face: FaceHandle::default(),
            vertices: vec![
                make_vertex(0.0, 0.0, 0.0),
                make_vertex(1.0, 0.0, 0.0),
                make_vertex(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
        };
        let mesh = TriMesh::from_face_meshes(vec![fm]);
        assert_eq!(mesh.vertices.len(), 3);
        assert_eq!(mesh.triangles.len(), 1);
        assert_eq!(mesh.triangles[0], [0, 1, 2]);
        assert_eq!(mesh.face_map.len(), 1);
    }

    #[test]
    fn from_face_meshes_two_faces_triangles_reindexed() {
        let fm1 = FaceMesh {
            face: FaceHandle::default(),
            vertices: vec![
                make_vertex(0.0, 0.0, 0.0),
                make_vertex(1.0, 0.0, 0.0),
                make_vertex(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
        };
        let fm2 = FaceMesh {
            face: FaceHandle::default(),
            vertices: vec![
                make_vertex(2.0, 0.0, 0.0),
                make_vertex(3.0, 0.0, 0.0),
                make_vertex(2.0, 1.0, 0.0),
                make_vertex(3.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
        };
        let mesh = TriMesh::from_face_meshes(vec![fm1, fm2]);
        // 3 vertices from fm1 + 4 from fm2
        assert_eq!(mesh.vertices.len(), 7);
        // 1 triangle from fm1 + 2 from fm2
        assert_eq!(mesh.triangles.len(), 3);
        // fm1 triangle unchanged (offset 0)
        assert_eq!(mesh.triangles[0], [0, 1, 2]);
        // fm2 triangles offset by 3
        assert_eq!(mesh.triangles[1], [3, 4, 5]);
        assert_eq!(mesh.triangles[2], [4, 6, 5]);
        assert_eq!(mesh.face_map.len(), 3);
    }
}
