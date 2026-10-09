//! Tessellations kept across interference checks (PK-722).
//!
//! A drive-train sweep runs the interference check at every step of a
//! slider. Only placements change between steps, but each check
//! re-tessellated every definition body, and re-meshed near-contact
//! pairs finer, from scratch. A cycloid disc spent seconds per step on
//! meshes identical to the last step's.
//!
//! The cache is keyed by CONTENT: a hash of the body's topology and
//! geometry, with the chord height and subdivision cap. An edited or
//! rebuilt definition hashes differently and misses; nothing has to
//! invalidate it. The hash covers vertices, half-edges with their
//! pcurves, edges with their curves, loops, and faces with their
//! surfaces, in full precision. It leaves out names, which the
//! tessellation never reads.

use std::collections::HashMap;
use std::sync::Arc;

use popkern_core::topology::BRep;

use crate::mesh_bvh::LocalMesh;

/// Entries kept before the cache is cleared (a long session of
/// edits must not grow it without bound).
const CAPACITY: usize = 512;

type Key = (u64, u64, usize, usize, u64);

/// Tessellations by body content and mesh settings.
#[derive(Debug, Clone, Default)]
pub(crate) struct MeshCache {
    meshes: HashMap<Key, Arc<LocalMesh>>,
}

/// The content key of `brep` (see the module docs).
pub(crate) fn fingerprint(brep: &BRep) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::hash::DefaultHasher::new();
    format!(
        "{:?}{:?}{:?}{:?}{:?}",
        brep.vertices, brep.half_edges, brep.edges, brep.loops, brep.faces
    )
    .hash(&mut hasher);
    hasher.finish()
}

impl MeshCache {
    /// The tessellation of `brep` (content key `content`) under
    /// `options`, built on a miss.
    pub(crate) fn get_or_build(
        &mut self,
        content: u64,
        brep: &BRep,
        options: &popkern_mesh::TessellationOptions,
    ) -> Result<Arc<LocalMesh>, popkern_mesh::MeshError> {
        let key = (
            content,
            options.chord_height.to_bits(),
            options.max_subdivisions,
            options.min_subdivisions,
            options.max_facet_angle.to_bits(),
        );
        if let Some(mesh) = self.meshes.get(&key) {
            return Ok(mesh.clone());
        }
        let mesh = popkern_mesh::tessellate(brep, options)?;
        let local = Arc::new(LocalMesh::new(brep, &mesh));
        if self.meshes.len() >= CAPACITY {
            self.meshes.clear();
        }
        self.meshes.insert(key, local.clone());
        Ok(local)
    }
}
