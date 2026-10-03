//! BRep → GPU-shaped buffers: the display-quality tessellation
//! ([`TessellationOptions::display`] — the gate chord height plus an
//! angular facet bound, so small blends read smoothly in shaded and
//! zebra views), flattened to `f32` positions/normals, `u32` indices
//! and a per-triangle face id, plus the welded open-edge census (the
//! same measure `Part::open_edge_count` reports).

use std::collections::HashMap;

use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::{MeshError, TessellationOptions, TriMesh};

use crate::envelope::BodyMesh;

/// Angular facet bound of the ANALYSIS mesh the zebra view renders
/// from (PK-577). At the 12° display bound the zebra flips the stripe
/// band on 3–10% of pixels on curved NURBS parts. At 6° the worst
/// gallery part is 5.6% (pyramid_apex 0.8%). Shaded views stay on the
/// display mesh.
pub const ANALYSIS_MAX_FACET_ANGLE_DEG: f64 = 6.0;

/// Largest estimated analysis mesh the zebra view builds: 400k
/// triangles, about 11 MB of position, normal and index buffers. It
/// covers every gallery part; coil_spring, the largest, meshes 361k
/// analysis triangles in about 7 s native since PK-589.
pub const ANALYSIS_TRIANGLE_BUDGET: usize = 400_000;

/// The analysis tessellation: the display preset at
/// [`ANALYSIS_MAX_FACET_ANGLE_DEG`].
pub fn analysis_options() -> TessellationOptions {
    TessellationOptions {
        max_facet_angle: ANALYSIS_MAX_FACET_ANGLE_DEG.to_radians(),
        ..TessellationOptions::display()
    }
}

/// What a body mesh is for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshPurpose {
    /// Shaded and curvature views: the display preset (12°) with the
    /// analytic per-vertex curvature fields.
    Display,
    /// The zebra view (PK-577): [`analysis_options`] (6°), with no
    /// curvature fields. Zebra shading reads normals only, and the
    /// analytic fields over a 6° mesh dominated its cost (coil_spring
    /// 691 s against about 60 s for the tessellation alone).
    ZebraAnalysis,
}

/// Tessellate `brep` into buffers for `purpose`. `label_face` spells a
/// face's stable name (or answers `None` for an unnamed face).
pub fn body_mesh(
    brep: &BRep,
    id: u64,
    label: Option<String>,
    label_face: impl Fn(FaceHandle) -> Option<String>,
    purpose: MeshPurpose,
) -> Result<BodyMesh, MeshError> {
    let options = match purpose {
        MeshPurpose::Display => TessellationOptions::display(),
        MeshPurpose::ZebraAnalysis => analysis_options(),
    };
    let mesh = popkern_mesh::tessellate(brep, &options)?;
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
    // Analysis curvature fields (PK-520): analytic per-vertex values,
    // computed once at mesh time so the viewer's mode switch is a
    // buffer swap, not a recompute round-trip.
    let field = |kind| {
        popkern_measure::curvature_field(brep, &mesh, kind)
            .map(|f| f.per_vertex.iter().map(|&v| v as f32).collect())
            .unwrap_or_default()
    };
    let (gaussian, mean): (Vec<f32>, Vec<f32>) = match purpose {
        MeshPurpose::Display => (
            field(popkern_measure::CurvatureKind::Gaussian),
            field(popkern_measure::CurvatureKind::Mean),
        ),
        MeshPurpose::ZebraAnalysis => (Vec::new(), Vec::new()),
    };
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
        gaussian,
        mean,
    })
}

/// Welded open-edge census: the part facade's gate
/// ([`popkern_mesh::welded_open_edges`]), so a body the page calls
/// watertight is one the CLI would too.
pub fn welded_open_edges(mesh: &TriMesh) -> usize {
    popkern_mesh::welded_open_edges(mesh)
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
