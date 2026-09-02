//! BRep → GPU-shaped buffers: the default-quality tessellation every
//! gate uses, flattened to `f32` positions/normals, `u32` indices and a
//! per-triangle face id, plus the welded open-edge census (the same
//! measure `Part::open_edge_count` reports).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::{MeshError, TessellationOptions, TriMesh};

use crate::envelope::BodyMesh;

/// Tessellate `brep` into buffers. `label_face` spells a face's stable
/// name (or answers `None` for an unnamed face).
pub fn body_mesh(
    brep: &BRep,
    id: u64,
    label: Option<String>,
    label_face: impl Fn(FaceHandle) -> Option<String>,
) -> Result<BodyMesh, MeshError> {
    let mesh = popkern_mesh::tessellate(brep, &TessellationOptions::default())?;
    let mut face_index: HashMap<FaceHandle, u32> = HashMap::new();
    let mut face_labels = Vec::with_capacity(brep.faces.len());
    for (handle, _) in brep.faces.iter() {
        face_index.insert(handle, face_labels.len() as u32);
        face_labels.push(label_face(handle).unwrap_or_else(|| "unnamed".to_string()));
    }
    let mut positions = Vec::with_capacity(mesh.vertices.len() * 3);
    let mut normals = Vec::with_capacity(mesh.vertices.len() * 3);
    let mut bounds = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
    for v in &mesh.vertices {
        let p = [v.position.x, v.position.y, v.position.z];
        for axis in 0..3 {
            bounds[0][axis] = bounds[0][axis].min(p[axis]);
            bounds[1][axis] = bounds[1][axis].max(p[axis]);
        }
        positions.extend(p.iter().map(|&c| c as f32));
        normals.extend([v.normal.x as f32, v.normal.y as f32, v.normal.z as f32]);
    }
    if mesh.vertices.is_empty() {
        bounds = [[0.0; 3]; 2];
    }
    let mut indices = Vec::with_capacity(mesh.triangles.len() * 3);
    let mut face_ids = Vec::with_capacity(mesh.triangles.len());
    for (t, tri) in mesh.triangles.iter().enumerate() {
        indices.extend(tri.iter().map(|&i| i as u32));
        let face = mesh
            .face_map
            .get(t)
            .and_then(|h| face_index.get(h))
            .copied()
            .unwrap_or(u32::MAX);
        face_ids.push(face);
    }
    Ok(BodyMesh {
        id,
        label,
        vertex_count: mesh.vertices.len(),
        triangle_count: mesh.triangles.len(),
        face_count: brep.faces.len(),
        face_labels,
        bounds,
        open_edges: welded_open_edges(&mesh),
        positions,
        normals,
        indices,
        face_ids,
    })
}

/// Welded open-edge census: weld vertices at 1e-7, count edges used by
/// exactly one triangle (0 = watertight). Mirrors the part facade's
/// gate so a body the page calls watertight is one the CLI would too.
pub fn welded_open_edges(mesh: &TriMesh) -> usize {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: Point3| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}

/// FNV-1a over the little-endian bytes of `f32` values: the page-side
/// smoke test hashes the same bytes so wasm and native tessellations
/// can be compared bit for bit without shipping the buffers back.
pub fn fnv1a_f32(values: &[f32]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for v in values {
        for byte in v.to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

/// FNV-1a over the little-endian bytes of `u32` values.
pub fn fnv1a_u32(values: &[u32]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for v in values {
        for byte in v.to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}
