//! The occurrences an interference check measures, placed lazily
//! (PK-723).
//!
//! A check used to flatten the whole assembly first: every part body
//! moved to world coordinates (a BRep transform) and every mesh placed
//! with its own hierarchy, for every occurrence, at every call. On the
//! six-axis arm (741 occurrences, 124 of them threaded screws) that was
//! seconds per call, and a swing calls it at every step while only a
//! few dozen pairs come within its clearance.
//!
//! The broad phase now runs on world bounds alone: each body's local
//! mesh box, its eight corners placed. Only the occurrences in a
//! candidate pair get their world BRep and placed mesh.

use std::collections::HashMap;
use std::sync::Arc;

use crate::assembly::{Assembly, DefinitionId};
use crate::interference::Occurrence;
use popkern_core::topology::BRep;

use crate::mesh_bvh::{Aabb, LocalMesh};
use crate::placement::Placement;

/// A body moved to its world placement on first use: only pairs that
/// refine onto the true faces or build the overlap need it.
#[derive(Debug)]
pub(crate) struct WorldBrep {
    parts: Arc<crate::parts_cache::CachedParts>,
    index: usize,
    world: Placement,
    cell: std::sync::OnceLock<Option<BRep>>,
}

impl WorldBrep {
    /// The world BRep (`None` when the placement is refused).
    pub(crate) fn get(&self) -> Option<&BRep> {
        self.cell
            .get_or_init(|| {
                let local = self.parts.parts.get(self.index)?.brep();
                if self.world.is_identity() {
                    Some(local.clone())
                } else {
                    local.transform(self.world.transform()).ok()
                }
            })
            .as_ref()
    }

    /// One vertex of the body, in world coordinates (a containment
    /// probe), without placing the whole body.
    pub(crate) fn first_vertex(&self) -> Option<popkern_core::geometry::Point3> {
        let local = self.parts.parts.get(self.index)?.brep();
        let v = local.vertices.values().next()?;
        Some(self.world.apply_point(&v.position))
    }
}

/// One part occurrence before placement: its bodies' local meshes and
/// world bounds (with the meshes' largest deviation).
pub(crate) struct Candidate {
    pub(crate) path: String,
    pub(crate) definition: DefinitionId,
    pub(crate) world: Placement,
    pub(crate) bodies: Vec<(usize, Aabb, f64)>,
}

/// `local`'s box placed at `world`: the box of its eight corners.
pub(crate) fn world_bounds(local: &LocalMesh, world: &Placement) -> Aabb {
    let (lo, hi) = (local.bounds.lo, local.bounds.hi);
    let mut out = Aabb {
        lo: [f64::INFINITY; 3],
        hi: [f64::NEG_INFINITY; 3],
    };
    for k in 0..8 {
        let corner = popkern_core::geometry::Point3::new(
            if k & 1 == 0 { lo[0] } else { hi[0] },
            if k & 2 == 0 { lo[1] } else { hi[1] },
            if k & 4 == 0 { lo[2] } else { hi[2] },
        );
        let q = world.apply_point(&corner);
        for (i, c) in [q.x, q.y, q.z].into_iter().enumerate() {
            out.lo[i] = out.lo[i].min(c);
            out.hi[i] = out.hi[i].max(c);
        }
    }
    out
}

impl Assembly {
    /// The selected part occurrences with their bodies' world bounds,
    /// in depth-first instance order.
    pub(crate) fn interference_candidates(
        &self,
        selected: &dyn Fn(&str) -> bool,
        local: &HashMap<(DefinitionId, usize), Arc<LocalMesh>>,
    ) -> Vec<Candidate> {
        let mut out = Vec::new();
        for row in self.leaves() {
            let path = row.path.to_string();
            if !row.is_part || !selected(&path) {
                continue;
            }
            let mut bodies = Vec::new();
            for body in 0.. {
                let Some(mesh) = local.get(&(row.definition, body)) else {
                    break;
                };
                bodies.push((body, world_bounds(mesh, &row.world), mesh.max_eps));
            }
            if !bodies.is_empty() {
                out.push(Candidate {
                    path,
                    definition: row.definition,
                    world: row.world,
                    bodies,
                });
            }
        }
        out
    }

    /// Place candidate `c`: its bodies' meshes (reused while the
    /// occurrence has not moved) and lazily placed BReps.
    pub(crate) fn place_candidate(
        &self,
        c: &Candidate,
        parts: &HashMap<DefinitionId, Arc<crate::parts_cache::CachedParts>>,
        local: &HashMap<(DefinitionId, usize), Arc<LocalMesh>>,
    ) -> Occurrence {
        let mut bodies = Vec::with_capacity(c.bodies.len());
        for &(body, _, _) in &c.bodies {
            let (Some(cached), Some(mesh)) =
                (parts.get(&c.definition), local.get(&(c.definition, body)))
            else {
                continue;
            };
            let world = self.meshes.placed(&c.path, body, mesh, &c.world);
            bodies.push((
                Arc::new(WorldBrep {
                    parts: cached.clone(),
                    index: body,
                    world: c.world,
                    cell: std::sync::OnceLock::new(),
                }),
                world,
                body,
            ));
        }
        Occurrence {
            path: c.path.clone(),
            definition: c.definition,
            world: c.world,
            bodies,
        }
    }
}
