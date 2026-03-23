//! glTF 2.0 format writer — JSON manifest + binary buffer.
//!
//! Produces a `.gltf` JSON file referencing a separate `.bin` binary buffer.
//! The binary buffer contains interleaved vertex data (position, normal, UV)
//! preceded by triangle index data.

use crate::{MeshError, TriMesh};
use std::io::Write;

/// Write `mesh` in glTF 2.0 format as a JSON manifest and a binary buffer.
///
/// `gltf_writer` receives the JSON manifest, `bin_writer` receives the raw
/// binary buffer, and `bin_filename` is embedded in the JSON as the buffer URI.
///
/// # Binary buffer layout
///
/// ```text
/// [indices: u32 LE × 3 × num_triangles] [vertices: f32 LE × 8 × num_vertices]
/// ```
///
/// Each vertex is 8 floats: 3 position + 3 normal + 2 UV (32 bytes stride).
pub fn write_gltf<W1: Write, W2: Write>(
    mesh: &TriMesh,
    gltf_writer: &mut W1,
    bin_writer: &mut W2,
    bin_filename: &str,
) -> Result<(), MeshError> {
    let num_tris = mesh.triangles.len();
    let num_verts = mesh.vertices.len();
    let num_indices = num_tris * 3;

    // --- Write binary buffer ---

    // Indices (u32 LE)
    let mut max_index: u32 = 0;
    for tri in &mesh.triangles {
        for &idx in tri {
            let idx32 = idx as u32;
            bin_writer.write_all(&idx32.to_le_bytes())?;
            if idx32 > max_index {
                max_index = idx32;
            }
        }
    }

    // Vertices (interleaved f32 LE: pos.x, pos.y, pos.z, n.x, n.y, n.z, u, v)
    let mut min_pos = [f64::MAX; 3];
    let mut max_pos = [f64::MIN; 3];

    for v in &mesh.vertices {
        let px = v.position.x;
        let py = v.position.y;
        let pz = v.position.z;

        if px < min_pos[0] {
            min_pos[0] = px;
        }
        if py < min_pos[1] {
            min_pos[1] = py;
        }
        if pz < min_pos[2] {
            min_pos[2] = pz;
        }
        if px > max_pos[0] {
            max_pos[0] = px;
        }
        if py > max_pos[1] {
            max_pos[1] = py;
        }
        if pz > max_pos[2] {
            max_pos[2] = pz;
        }

        bin_writer.write_all(&(px as f32).to_le_bytes())?;
        bin_writer.write_all(&(py as f32).to_le_bytes())?;
        bin_writer.write_all(&(pz as f32).to_le_bytes())?;
        bin_writer.write_all(&(v.normal.x as f32).to_le_bytes())?;
        bin_writer.write_all(&(v.normal.y as f32).to_le_bytes())?;
        bin_writer.write_all(&(v.normal.z as f32).to_le_bytes())?;
        bin_writer.write_all(&(v.uv.0 as f32).to_le_bytes())?;
        bin_writer.write_all(&(v.uv.1 as f32).to_le_bytes())?;
    }

    // --- Compute sizes ---
    let index_bytes = num_indices * 4; // u32 = 4 bytes
    let vertex_bytes = num_verts * 32; // 8 floats × 4 bytes
    let total_bytes = index_bytes + vertex_bytes;

    // Handle empty mesh — clamp min/max to zero
    if num_verts == 0 {
        min_pos = [0.0; 3];
        max_pos = [0.0; 3];
    }

    // --- Build JSON ---
    let json = format!(
        r#"{{
  "asset": {{ "version": "2.0", "generator": "popkern" }},
  "buffers": [{{ "uri": "{bin_filename}", "byteLength": {total_bytes} }}],
  "bufferViews": [
    {{ "buffer": 0, "byteOffset": 0, "byteLength": {index_bytes}, "target": 34963 }},
    {{ "buffer": 0, "byteOffset": {index_bytes}, "byteLength": {vertex_bytes}, "target": 34962, "byteStride": 32 }}
  ],
  "accessors": [
    {{ "bufferView": 0, "componentType": 5125, "count": {num_indices}, "type": "SCALAR", "max": [{max_index}], "min": [0] }},
    {{ "bufferView": 1, "byteOffset": 0, "componentType": 5126, "count": {num_verts}, "type": "VEC3", "max": [{mx}, {my}, {mz}], "min": [{nx}, {ny}, {nz}] }},
    {{ "bufferView": 1, "byteOffset": 12, "componentType": 5126, "count": {num_verts}, "type": "VEC3" }},
    {{ "bufferView": 1, "byteOffset": 24, "componentType": 5126, "count": {num_verts}, "type": "VEC2" }}
  ],
  "meshes": [{{ "primitives": [{{ "attributes": {{ "POSITION": 1, "NORMAL": 2, "TEXCOORD_0": 3 }}, "indices": 0 }}] }}],
  "nodes": [{{ "mesh": 0 }}],
  "scenes": [{{ "nodes": [0] }}],
  "scene": 0
}}
"#,
        mx = max_pos[0],
        my = max_pos[1],
        mz = max_pos[2],
        nx = min_pos[0],
        ny = min_pos[1],
        nz = min_pos[2],
    );

    gltf_writer.write_all(json.as_bytes())?;
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
    fn gltf_json_has_required_fields() {
        let mesh = single_triangle_mesh();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
        let json = String::from_utf8(json_buf).unwrap();
        assert!(json.contains("\"asset\""));
        assert!(json.contains("\"version\""));
        assert!(json.contains("\"buffers\""));
        assert!(json.contains("\"accessors\""));
        assert!(json.contains("\"meshes\""));
    }

    #[test]
    fn gltf_binary_buffer_correct_size() {
        let mesh = single_triangle_mesh();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
        let num_tris = mesh.triangles.len();
        let num_verts = mesh.vertices.len();
        let expected = num_tris * 3 * 4 + num_verts * 8 * 4;
        assert_eq!(bin_buf.len(), expected);
    }

    #[test]
    fn gltf_json_contains_bin_filename() {
        let mesh = single_triangle_mesh();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "my_model.bin").unwrap();
        let json = String::from_utf8(json_buf).unwrap();
        assert!(json.contains("my_model.bin"));
    }

    #[test]
    fn gltf_position_bounds_correct() {
        let mesh = single_triangle_mesh();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
        let json = String::from_utf8(json_buf).unwrap();
        // max position is (1, 1, 0), min is (0, 0, 0)
        assert!(json.contains("\"max\": [1, 1, 0]"));
        assert!(json.contains("\"min\": [0, 0, 0]"));
    }

    #[test]
    fn gltf_empty_mesh_produces_valid_output() {
        let mesh = TriMesh::default();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "empty.bin").unwrap();
        assert_eq!(bin_buf.len(), 0);
        let json = String::from_utf8(json_buf).unwrap();
        assert!(json.contains("\"asset\""));
    }

    #[test]
    fn gltf_byte_stride_present() {
        let mesh = single_triangle_mesh();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
        let json = String::from_utf8(json_buf).unwrap();
        assert!(json.contains("\"byteStride\": 32"));
    }

    #[test]
    fn gltf_indices_correctly_written() {
        let mesh = single_triangle_mesh();
        let mut json_buf = Vec::new();
        let mut bin_buf = Vec::new();
        write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
        // First 12 bytes are 3 u32 indices: 0, 1, 2
        let i0 = u32::from_le_bytes(bin_buf[0..4].try_into().unwrap());
        let i1 = u32::from_le_bytes(bin_buf[4..8].try_into().unwrap());
        let i2 = u32::from_le_bytes(bin_buf[8..12].try_into().unwrap());
        assert_eq!((i0, i1, i2), (0, 1, 2));
    }
}
