//! Which way each face's surface chart faces relative to the mesh
//! (split from the zebra rasterizer under the 500-line budget).

use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TriMesh;

/// Per-face facing sign of the surface chart relative to the mesh:
/// the length-weighted vote of `chart_normal · stored_normal` over the
/// face mesh's vertices. The mesh's facing is the authority (the
/// tessellator orients faces by global winding coherence, because
/// `SurfaceOnFace::orientation` is not trustworthy on its own).
pub(super) fn facing_signs(
    brep: &BRep,
    mesh: &TriMesh,
) -> (Vec<(FaceHandle, f64)>, Vec<FaceHandle>) {
    let mut signs = Vec::with_capacity(mesh.face_meshes.len());
    let mut mixed = Vec::new();
    for fm in &mesh.face_meshes {
        let Some(sof) = brep.faces.get(fm.face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let (mut agree, mut oppose) = (0usize, 0usize);
        for v in &fm.vertices {
            if let Some(n) = popkern_mesh::surface_limit_normal(&sof.geometry, v.uv.0, v.uv.1) {
                let d = n.dot(&v.normal);
                if d > 0.0 {
                    agree += 1;
                } else if d < 0.0 {
                    oppose += 1;
                }
            }
        }
        if agree > 0 && oppose > 0 {
            mixed.push(fm.face);
        }
        signs.push((fm.face, if oppose > agree { -1.0 } else { 1.0 }));
    }
    (signs, mixed)
}
