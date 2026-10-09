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
use std::sync::{Arc, Mutex};

use popkern_core::topology::BRep;

use crate::mesh_bvh::{LocalMesh, WorldMesh};
use crate::placement::Placement;

/// Entries kept before the cache is cleared (a long session of
/// edits must not grow it without bound).
const CAPACITY: usize = 512;

type Key = (u64, u64, usize, usize, u64);

/// Tessellations by body content and mesh settings. Shared by the
/// interference check's worker threads.
#[derive(Debug, Default)]
pub(crate) struct MeshCache {
    meshes: Mutex<HashMap<Key, Arc<LocalMesh>>>,
    /// Placed meshes by occurrence path and body, with the local mesh
    /// and placement they were built from: an occurrence that has not
    /// moved since the last check (half a swing's) reuses its BVH.
    placed: Mutex<HashMap<(String, usize), Placed>>,
}

type Placed = (Arc<LocalMesh>, [u64; 12], Arc<WorldMesh>);

impl Clone for MeshCache {
    fn clone(&self) -> Self {
        MeshCache {
            meshes: Mutex::new(self.meshes.lock().map(|m| m.clone()).unwrap_or_default()),
            placed: Mutex::new(self.placed.lock().map(|m| m.clone()).unwrap_or_default()),
        }
    }
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
        &self,
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
        if let Some(mesh) = self.meshes.lock().ok().and_then(|m| m.get(&key).cloned()) {
            return Ok(mesh);
        }
        // Built outside the lock: two workers may build the same mesh
        // once each; the second insert is a no-op in effect.
        let mesh = popkern_mesh::tessellate(brep, options)?;
        let local = Arc::new(LocalMesh::new(brep, &mesh));
        if let Ok(mut meshes) = self.meshes.lock() {
            if meshes.len() >= CAPACITY {
                meshes.clear();
            }
            meshes.insert(key, local.clone());
        }
        Ok(local)
    }
}

impl MeshCache {
    /// `mesh` placed at `world` for occurrence `path` (body `body`),
    /// reused when neither changed since the last call.
    pub(crate) fn placed(
        &self,
        path: &str,
        body: usize,
        mesh: &Arc<LocalMesh>,
        world: &Placement,
    ) -> Arc<WorldMesh> {
        let m = &world.transform().matrix;
        let mut bits = [0u64; 12];
        for i in 0..3 {
            for j in 0..4 {
                bits[i * 4 + j] = m[i][j].to_bits();
            }
        }
        let key = (path.to_string(), body);
        if let Some((local, b, placed)) = self.placed.lock().ok().and_then(|p| p.get(&key).cloned())
            && Arc::ptr_eq(&local, mesh)
            && b == bits
        {
            return placed;
        }
        let placed = Arc::new(WorldMesh::new(mesh, world));
        if let Ok(mut all) = self.placed.lock() {
            if all.len() >= CAPACITY * 8 {
                all.clear();
            }
            all.insert(key, (mesh.clone(), bits, placed.clone()));
        }
        placed
    }
}
