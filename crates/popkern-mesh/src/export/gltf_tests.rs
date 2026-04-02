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

#[test]
fn gltf_json_is_valid_structure() {
    // Verify the JSON can be parsed and has key structural fields.
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).expect("JSON should parse");
    assert_eq!(json["asset"]["version"], "2.0");
    assert_eq!(json["asset"]["generator"], "popkern");
    assert_eq!(json["scene"], 0);
    assert!(json["scenes"].is_array());
    assert!(json["nodes"].is_array());
    assert!(json["meshes"].is_array());
    assert!(json["accessors"].is_array());
    assert!(json["bufferViews"].is_array());
    assert!(json["buffers"].is_array());
}

#[test]
fn gltf_buffer_byte_length_matches_binary() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    let declared_len = json["buffers"][0]["byteLength"].as_u64().unwrap();
    assert_eq!(declared_len as usize, bin_buf.len());
}

#[test]
fn gltf_buffer_views_cover_full_buffer() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    let bv0_offset = json["bufferViews"][0]["byteOffset"].as_u64().unwrap();
    let bv0_len = json["bufferViews"][0]["byteLength"].as_u64().unwrap();
    let bv1_offset = json["bufferViews"][1]["byteOffset"].as_u64().unwrap();
    let bv1_len = json["bufferViews"][1]["byteLength"].as_u64().unwrap();
    // Index view starts at 0
    assert_eq!(bv0_offset, 0);
    // Vertex view starts where index view ends
    assert_eq!(bv1_offset, bv0_len);
    // Total coverage
    assert_eq!(bv0_len + bv1_len, bin_buf.len() as u64);
}

#[test]
fn gltf_accessor_counts_correct() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    // Accessor 0 = indices: count should be 3 (1 triangle * 3)
    assert_eq!(json["accessors"][0]["count"], 3);
    assert_eq!(json["accessors"][0]["type"], "SCALAR");
    assert_eq!(json["accessors"][0]["componentType"], 5125); // UNSIGNED_INT
    // Accessor 1 = position: count should be 3 vertices
    assert_eq!(json["accessors"][1]["count"], 3);
    assert_eq!(json["accessors"][1]["type"], "VEC3");
    assert_eq!(json["accessors"][1]["componentType"], 5126); // FLOAT
    // Accessor 2 = normal
    assert_eq!(json["accessors"][2]["count"], 3);
    assert_eq!(json["accessors"][2]["type"], "VEC3");
    // Accessor 3 = texcoord
    assert_eq!(json["accessors"][3]["count"], 3);
    assert_eq!(json["accessors"][3]["type"], "VEC2");
}

#[test]
fn gltf_mesh_primitive_attributes_correct() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    let attrs = &json["meshes"][0]["primitives"][0]["attributes"];
    assert_eq!(attrs["POSITION"], 1);
    assert_eq!(attrs["NORMAL"], 2);
    assert_eq!(attrs["TEXCOORD_0"], 3);
    assert_eq!(json["meshes"][0]["primitives"][0]["indices"], 0);
}

#[test]
fn gltf_vertex_data_decodes_correctly() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    // Vertex data starts after indices: 3 indices * 4 bytes = 12 bytes offset.
    let vdata_start = 12;
    // First vertex: pos=(0,0,0), normal=(0,0,1), uv=(0,0) -> 8 floats
    let px = f32::from_le_bytes(bin_buf[vdata_start..vdata_start + 4].try_into().unwrap());
    let py = f32::from_le_bytes(
        bin_buf[vdata_start + 4..vdata_start + 8]
            .try_into()
            .unwrap(),
    );
    let pz = f32::from_le_bytes(
        bin_buf[vdata_start + 8..vdata_start + 12]
            .try_into()
            .unwrap(),
    );
    assert!((px - 0.0).abs() < 1e-6);
    assert!((py - 0.0).abs() < 1e-6);
    assert!((pz - 0.0).abs() < 1e-6);
    let nx = f32::from_le_bytes(
        bin_buf[vdata_start + 12..vdata_start + 16]
            .try_into()
            .unwrap(),
    );
    let ny = f32::from_le_bytes(
        bin_buf[vdata_start + 16..vdata_start + 20]
            .try_into()
            .unwrap(),
    );
    let nz = f32::from_le_bytes(
        bin_buf[vdata_start + 20..vdata_start + 24]
            .try_into()
            .unwrap(),
    );
    assert!((nx - 0.0).abs() < 1e-6);
    assert!((ny - 0.0).abs() < 1e-6);
    assert!((nz - 1.0).abs() < 1e-6);
    let u = f32::from_le_bytes(
        bin_buf[vdata_start + 24..vdata_start + 28]
            .try_into()
            .unwrap(),
    );
    let v = f32::from_le_bytes(
        bin_buf[vdata_start + 28..vdata_start + 32]
            .try_into()
            .unwrap(),
    );
    assert!((u - 0.0).abs() < 1e-6);
    assert!((v - 0.0).abs() < 1e-6);
}

#[test]
fn gltf_second_vertex_data_correct() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    // Second vertex starts at 12 + 32 = 44.
    let off = 44;
    let px = f32::from_le_bytes(bin_buf[off..off + 4].try_into().unwrap());
    assert!((px - 1.0).abs() < 1e-6); // position.x = 1.0
}

#[test]
fn gltf_multi_triangle_buffer_size() {
    let mesh = multi_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    // 2 triangles * 3 indices * 4 bytes = 24 index bytes
    // 4 vertices * 32 bytes = 128 vertex bytes
    assert_eq!(bin_buf.len(), 24 + 128);
}

#[test]
fn gltf_multi_triangle_index_accessor() {
    let mesh = multi_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    // 6 indices total
    assert_eq!(json["accessors"][0]["count"], 6);
    // max index should be 3
    assert_eq!(json["accessors"][0]["max"][0], 3);
    assert_eq!(json["accessors"][0]["min"][0], 0);
}

#[test]
fn gltf_multi_triangle_indices_decode() {
    let mesh = multi_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let indices: Vec<u32> = (0..6)
        .map(|i| {
            let off = i * 4;
            u32::from_le_bytes(bin_buf[off..off + 4].try_into().unwrap())
        })
        .collect();
    assert_eq!(indices, vec![0, 1, 2, 0, 2, 3]);
}

#[test]
fn gltf_empty_mesh_json_is_valid() {
    let mesh = TriMesh::default();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "empty.bin").unwrap();
    let json: serde_json::Value =
        serde_json::from_slice(&json_buf).expect("empty mesh JSON should parse");
    assert_eq!(json["asset"]["version"], "2.0");
    assert_eq!(json["buffers"][0]["byteLength"], 0);
    assert_eq!(json["accessors"][0]["count"], 0);
    assert_eq!(json["accessors"][1]["count"], 0);
}

#[test]
fn gltf_position_bounds_multi_triangle() {
    let mesh = multi_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    // max position is (1, 1, 0)
    let max_arr = &json["accessors"][1]["max"];
    assert_eq!(max_arr[0].as_f64().unwrap(), 1.0);
    assert_eq!(max_arr[1].as_f64().unwrap(), 1.0);
    assert_eq!(max_arr[2].as_f64().unwrap(), 0.0);
    // min position is (0, 0, 0)
    let min_arr = &json["accessors"][1]["min"];
    assert_eq!(min_arr[0].as_f64().unwrap(), 0.0);
    assert_eq!(min_arr[1].as_f64().unwrap(), 0.0);
    assert_eq!(min_arr[2].as_f64().unwrap(), 0.0);
}

#[test]
fn gltf_buffer_view_targets_correct() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    // 34963 = ELEMENT_ARRAY_BUFFER (indices)
    assert_eq!(json["bufferViews"][0]["target"], 34963);
    // 34962 = ARRAY_BUFFER (vertices)
    assert_eq!(json["bufferViews"][1]["target"], 34962);
}

#[test]
fn gltf_node_references_mesh() {
    let mesh = single_triangle_mesh();
    let mut json_buf = Vec::new();
    let mut bin_buf = Vec::new();
    write_gltf(&mesh, &mut json_buf, &mut bin_buf, "test.bin").unwrap();
    let json: serde_json::Value = serde_json::from_slice(&json_buf).unwrap();
    assert_eq!(json["nodes"][0]["mesh"], 0);
    assert_eq!(json["scenes"][0]["nodes"][0], 0);
}
