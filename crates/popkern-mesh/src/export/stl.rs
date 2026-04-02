//! STL (STereoLithography) format writer — ASCII and binary variants.

use crate::{MeshError, TriMesh};
use std::io::Write;

/// Write `mesh` in STL ASCII format to `writer`.
///
/// The solid name is `popkern`.  Every triangle gets a face normal computed
/// as the average of its three vertex normals.
pub fn write_stl_ascii<W: Write>(mesh: &TriMesh, writer: &mut W) -> Result<(), MeshError> {
    writeln!(writer, "solid popkern")?;

    for tri in &mesh.triangles {
        let [i0, i1, i2] = *tri;
        let n0 = &mesh.vertices[i0].normal;
        let n1 = &mesh.vertices[i1].normal;
        let n2 = &mesh.vertices[i2].normal;

        let nx = (n0.x + n1.x + n2.x) / 3.0;
        let ny = (n0.y + n1.y + n2.y) / 3.0;
        let nz = (n0.z + n1.z + n2.z) / 3.0;

        writeln!(writer, "  facet normal {nx} {ny} {nz}")?;
        writeln!(writer, "    outer loop")?;

        for &idx in &[i0, i1, i2] {
            let p = &mesh.vertices[idx].position;
            writeln!(writer, "      vertex {} {} {}", p.x, p.y, p.z)?;
        }

        writeln!(writer, "    endloop")?;
        writeln!(writer, "  endfacet")?;
    }

    writeln!(writer, "endsolid popkern")?;
    Ok(())
}

/// Write `mesh` in STL binary format to `writer`.
///
/// Layout: 80-byte zero header, 4-byte LE triangle count, then per triangle:
/// 12 × f32 LE (normal xyz + 3 × vertex xyz) followed by 2 zero attribute bytes
/// — 50 bytes per triangle total.
pub fn write_stl_binary<W: Write>(mesh: &TriMesh, writer: &mut W) -> Result<(), MeshError> {
    // 80-byte header (all zeros)
    writer.write_all(&[0u8; 80])?;

    // Triangle count
    let count = mesh.triangles.len() as u32;
    writer.write_all(&count.to_le_bytes())?;

    for tri in &mesh.triangles {
        let [i0, i1, i2] = *tri;
        let n0 = &mesh.vertices[i0].normal;
        let n1 = &mesh.vertices[i1].normal;
        let n2 = &mesh.vertices[i2].normal;

        let nx = ((n0.x + n1.x + n2.x) / 3.0) as f32;
        let ny = ((n0.y + n1.y + n2.y) / 3.0) as f32;
        let nz = ((n0.z + n1.z + n2.z) / 3.0) as f32;

        // Normal
        writer.write_all(&nx.to_le_bytes())?;
        writer.write_all(&ny.to_le_bytes())?;
        writer.write_all(&nz.to_le_bytes())?;

        // Three vertices
        for &idx in &[i0, i1, i2] {
            let p = &mesh.vertices[idx].position;
            writer.write_all(&(p.x as f32).to_le_bytes())?;
            writer.write_all(&(p.y as f32).to_le_bytes())?;
            writer.write_all(&(p.z as f32).to_le_bytes())?;
        }

        // Attribute byte count (zeros)
        writer.write_all(&[0u8; 2])?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{FaceMesh, MeshVertex};
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::topology::FaceHandle;

    fn single_triangle_mesh() -> TriMesh {
        let vertices = vec![
            MeshVertex {
                position: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (1.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 1.0),
            },
        ];
        let fm = FaceMesh {
            face: FaceHandle::default(),
            vertices: vertices.clone(),
            triangles: vec![[0, 1, 2]],
        };
        TriMesh::from_face_meshes(vec![fm])
    }

    #[test]
    fn ascii_starts_with_solid_popkern() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.starts_with("solid popkern\n"), "got: {text:?}");
        assert!(text.ends_with("endsolid popkern\n"), "got: {text:?}");
    }

    #[test]
    fn ascii_correct_facet_count() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let count = text.matches("facet normal").count();
        assert_eq!(count, mesh.triangles.len());
    }

    #[test]
    fn binary_correct_byte_length() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        let expected = 80 + 4 + mesh.triangles.len() * 50;
        assert_eq!(buf.len(), expected);
    }

    #[test]
    fn binary_triangle_count_in_header() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count as usize, mesh.triangles.len());
    }

    #[test]
    fn ascii_empty_mesh_valid_structure() {
        let mesh = TriMesh::default();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text, "solid popkern\nendsolid popkern\n");
    }

    #[test]
    fn binary_empty_mesh_correct_length() {
        let mesh = TriMesh::default();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        assert_eq!(buf.len(), 80 + 4);
    }

    #[test]
    fn binary_header_is_all_zeros() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        assert!(
            buf[..80].iter().all(|&b| b == 0),
            "header should be 80 zero bytes"
        );
    }

    #[test]
    fn binary_vertex_data_matches_input() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        // After 80 header + 4 count + 12 normal bytes, vertex data starts at offset 96.
        // First vertex: (0,0,0)
        let v0x = f32::from_le_bytes(buf[96..100].try_into().unwrap());
        let v0y = f32::from_le_bytes(buf[100..104].try_into().unwrap());
        let v0z = f32::from_le_bytes(buf[104..108].try_into().unwrap());
        assert!((v0x - 0.0).abs() < 1e-6);
        assert!((v0y - 0.0).abs() < 1e-6);
        assert!((v0z - 0.0).abs() < 1e-6);
        // Second vertex: (1,0,0)
        let v1x = f32::from_le_bytes(buf[108..112].try_into().unwrap());
        assert!((v1x - 1.0).abs() < 1e-6);
        // Third vertex: (0,1,0)
        let v2y = f32::from_le_bytes(buf[124..128].try_into().unwrap());
        assert!((v2y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn binary_normal_is_average_of_vertex_normals() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        // Normal at offset 84 (after 80 header + 4 count).
        let nx = f32::from_le_bytes(buf[84..88].try_into().unwrap());
        let ny = f32::from_le_bytes(buf[88..92].try_into().unwrap());
        let nz = f32::from_le_bytes(buf[92..96].try_into().unwrap());
        // All vertex normals are (0,0,1), so average is (0,0,1).
        assert!((nx - 0.0).abs() < 1e-6);
        assert!((ny - 0.0).abs() < 1e-6);
        assert!((nz - 1.0).abs() < 1e-6);
    }

    #[test]
    fn binary_attribute_bytes_are_zero() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        // Last 2 bytes of the single triangle record (at offset 80+4+48 = 132).
        assert_eq!(buf[132], 0);
        assert_eq!(buf[133], 0);
    }

    fn multi_triangle_mesh() -> TriMesh {
        let vertices = vec![
            MeshVertex {
                position: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (1.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(1.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (1.0, 1.0),
            },
            MeshVertex {
                position: Point3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 1.0),
            },
        ];
        let fm = FaceMesh {
            face: FaceHandle::default(),
            vertices: vertices.clone(),
            triangles: vec![[0, 1, 2], [0, 2, 3]],
        };
        TriMesh::from_face_meshes(vec![fm])
    }

    #[test]
    fn binary_multi_triangle_correct_count_and_size() {
        let mesh = multi_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count, 2);
        assert_eq!(buf.len(), 80 + 4 + 2 * 50);
    }

    #[test]
    fn ascii_multi_triangle_has_two_facets() {
        let mesh = multi_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text.matches("facet normal").count(), 2);
        assert_eq!(text.matches("endloop").count(), 2);
        assert_eq!(text.matches("endfacet").count(), 2);
    }

    #[test]
    fn ascii_vertex_values_match_input() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Check that the three vertex positions appear in the output.
        assert!(text.contains("vertex 0 0 0"));
        assert!(text.contains("vertex 1 0 0"));
        assert!(text.contains("vertex 0 1 0"));
    }

    #[test]
    fn ascii_normal_direction_correct() {
        let mesh = single_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // All vertex normals are (0,0,1), average is (0,0,1).
        assert!(text.contains("facet normal 0 0 1"));
    }

    #[test]
    fn binary_empty_mesh_triangle_count_is_zero() {
        let mesh = TriMesh::default();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        let count = u32::from_le_bytes(buf[80..84].try_into().unwrap());
        assert_eq!(count, 0);
    }

    fn degenerate_triangle_mesh() -> TriMesh {
        // Triangle where two vertices are the same point.
        let vertices = vec![
            MeshVertex {
                position: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (1.0, 0.0),
            },
        ];
        let fm = FaceMesh {
            face: FaceHandle::default(),
            vertices: vertices.clone(),
            triangles: vec![[0, 1, 2]],
        };
        TriMesh::from_face_meshes(vec![fm])
    }

    #[test]
    fn binary_degenerate_triangle_still_exports() {
        let mesh = degenerate_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        assert_eq!(buf.len(), 80 + 4 + 50);
    }

    #[test]
    fn ascii_degenerate_triangle_still_exports() {
        let mesh = degenerate_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_ascii(&mesh, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text.matches("facet normal").count(), 1);
    }

    #[test]
    fn binary_second_triangle_record_correct_offset() {
        let mesh = multi_triangle_mesh();
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        // Second triangle starts at 80 + 4 + 50 = 134.
        // Its first vertex (index 0) is at offset 134 + 12 (normal) = 146.
        let v0x = f32::from_le_bytes(buf[146..150].try_into().unwrap());
        let v0y = f32::from_le_bytes(buf[150..154].try_into().unwrap());
        let v0z = f32::from_le_bytes(buf[154..158].try_into().unwrap());
        // Second triangle is [0,2,3] -> first vertex is (0,0,0).
        assert!((v0x - 0.0).abs() < 1e-6);
        assert!((v0y - 0.0).abs() < 1e-6);
        assert!((v0z - 0.0).abs() < 1e-6);
    }

    #[test]
    fn binary_mixed_normals_averaged_correctly() {
        let vertices = vec![
            MeshVertex {
                position: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(1.0, 0.0, 0.0),
                uv: (0.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 1.0, 0.0),
                uv: (1.0, 0.0),
            },
            MeshVertex {
                position: Point3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                uv: (0.0, 1.0),
            },
        ];
        let fm = FaceMesh {
            face: FaceHandle::default(),
            vertices: vertices.clone(),
            triangles: vec![[0, 1, 2]],
        };
        let mesh = TriMesh::from_face_meshes(vec![fm]);
        let mut buf = Vec::new();
        write_stl_binary(&mesh, &mut buf).unwrap();
        let nx = f32::from_le_bytes(buf[84..88].try_into().unwrap());
        let ny = f32::from_le_bytes(buf[88..92].try_into().unwrap());
        let nz = f32::from_le_bytes(buf[92..96].try_into().unwrap());
        // Average of (1,0,0), (0,1,0), (0,0,1) is (1/3, 1/3, 1/3).
        let expected = 1.0_f32 / 3.0;
        assert!((nx - expected).abs() < 1e-6);
        assert!((ny - expected).abs() < 1e-6);
        assert!((nz - expected).abs() < 1e-6);
    }
}
