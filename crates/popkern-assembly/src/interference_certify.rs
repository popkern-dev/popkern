//! Adaptive certification of near-contact pairs (PK-747).
//!
//! The interference check meshes every definition at one chord height,
//! a fraction of the part's size: about 0.08 mm on an 80 mm part. A
//! fit whose clearance is below that cannot be decided at that
//! resolution. The meshes cross, so the exact boolean runs and finds no
//! overlap, and the certified lower bound (mesh distance minus both
//! meshes' deviations) goes negative. A cycloid disc 0.02 mm off its
//! ring pins read as touching every pin.
//!
//! An unsettled pair is re-meshed finer and checked again, only that
//! pair and only as fine as it needs. A pair is unsettled when:
//!
//! - its contact rests on mesh distances, not on refined points of the
//!   true faces within the tolerance;
//! - its refined gap is above the tolerance but the lower bound is not;
//! - or the boolean declined on crossing meshes, which a near-contact
//!   fit can provoke.
//!
//! Each level divides the chord height by 8, up to [`MAX_LEVELS`]; the
//! meshes join the assembly's content-keyed cache, so a sweep re-meshes
//! a body at most once per level. The first decided report wins. A pair still undecided at the finest
//! level keeps its last report, so a decline is never hidden.

use std::collections::HashMap;
use std::sync::Arc;

use popkern_part::Part;

use crate::assembly::{Assembly, DefinitionId};
use crate::interference::{Occurrence, PairReport};
use crate::mesh_bvh::{LocalMesh, WorldMesh};

/// Finer world meshes per occurrence path, body and level, shared by
/// every pair (and worker) certifying against that occurrence.
pub(crate) type WorldMeshes = std::sync::Mutex<HashMap<(String, usize, u32), Arc<WorldMesh>>>;

/// Refinement levels past the base chord (each 1/8 of the last).
pub(crate) const MAX_LEVELS: u32 = 3;

impl Assembly {
    /// Re-check pair (`a`, `b`) on finer meshes until it is decided;
    /// the overlap boolean runs only at the finest level. `None` when
    /// a finer mesh could not be built (the caller decides at the base
    /// level, boolean included); `Some(report)` otherwise — the first
    /// decided one, or the finest level's.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn certify_pair(
        &self,
        a: &Occurrence,
        b: &Occurrence,
        parts: &HashMap<DefinitionId, Vec<Part>>,
        contents: &HashMap<(DefinitionId, usize), u64>,
        worlds: &WorldMeshes,
        chord_fraction: f64,
        tolerance: f64,
        clearance: f64,
    ) -> Option<Option<PairReport>> {
        let mut last: Option<Option<PairReport>> = None;
        for level in 1..=MAX_LEVELS {
            let fraction = chord_fraction / 8f64.powi(level as i32);
            let mut local: HashMap<(DefinitionId, usize), Arc<LocalMesh>> = HashMap::new();
            let meshes = &self.meshes;
            let mut finer_occurrence = |o: &Occurrence| -> Option<Occurrence> {
                let mut bodies = Vec::with_capacity(o.bodies.len());
                for (brep, _, body) in &o.bodies {
                    let part = parts.get(&o.definition)?.get(*body)?;
                    let local_brep = part.brep();
                    let content = *contents.get(&(o.definition, *body))?;
                    let size = crate::interference_refine::bounds_diagonal(local_brep).max(1e-12);
                    // The subdivision cap rises with the chord: at the
                    // default 64 an r = 10 bore keeps a 0.012 sagitta
                    // whatever chord is asked for.
                    let base = popkern_mesh::TessellationOptions::default();
                    let options = popkern_mesh::TessellationOptions {
                        chord_height: fraction * size,
                        max_subdivisions: base.max_subdivisions << (3 * level),
                        ..base
                    };
                    let mesh = meshes.get_or_build(content, local_brep, &options).ok()?;
                    let key = (o.path.clone(), *body, level);
                    let cached = worlds.lock().ok().and_then(|w| w.get(&key).cloned());
                    let world = match cached {
                        Some(w) => w,
                        None => {
                            let w = Arc::new(WorldMesh::new(&mesh, &o.world));
                            if let Ok(mut all) = worlds.lock() {
                                all.insert(key, w.clone());
                            }
                            w
                        }
                    };
                    bodies.push((brep.clone(), world, *body));
                    local.insert((o.definition, *body), mesh);
                }
                Some(Occurrence {
                    path: o.path.clone(),
                    definition: o.definition,
                    world: o.world,
                    bodies,
                })
            };
            // No finer mesh: the caller decides at the base level.
            let (Some(fa), Some(fb)) = (finer_occurrence(a), finer_occurrence(b)) else {
                return None;
            };
            let check =
                self.check_pair(&fa, &fb, &local, tolerance, clearance, level == MAX_LEVELS);
            if check.settled {
                return Some(check.report);
            }
            last = Some(check.report);
        }
        last
    }
}
