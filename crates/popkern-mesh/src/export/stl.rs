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
}
