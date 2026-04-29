//! Phase 6c-C.2d — apex-region rebuild orchestration.
//!
//! Wires Phase 6c-C.1 (detection), 6c-C.2 (plan), 6c-C.2b
//! (mef execute), and 6c-C.2c (lens collapse) into a single
//! pass that runs from the fillet pipeline.
//!
//! For each cap face produced by the fillet pipeline:
//! 1. Build the deletion plan via `plan_apex_region_deletion`.
//! 2. If actionable (uniform-radius case with the canonical
//!    6-edge cap loop), run `execute_deletion_plan` to insert
//!    the cap-tangent arc edges and create the lens sub-faces.
//! 3. Collapse each lens via `collapse_apex_lens`.
//! 4. If any step fails, the snapshot-rollback in each
//!    primitive ensures the BRep is restored to its pre-call
//!    state for that cap.  The pass continues with the next
//!    cap face.
//!
//! Activation: env-gated by `POPKERN_APEX_TOPOLOGY_REBUILD=1`.
//! Default off.  Phase 6c-F flips the default once the
//! validation matrix passes on the full fixture set.

use popkern_core::topology::{BRep, FaceHandle};

use super::collapse::collapse_apex_lens;
use super::execute::execute_deletion_plan;
use super::plan::plan_apex_region_deletion;

/// Run apex topology rebuild on every cap face in `cap_faces`.
/// Returns the count of successfully-rebuilt apexes.  Failures
/// are silent — the BRep is unchanged for caps that don't fit
/// the canonical pattern.
#[allow(dead_code)]
pub(crate) fn rebuild_all_apex_caps(brep: &mut BRep, cap_faces: &[FaceHandle]) -> usize {
    let mut succeeded = 0usize;
    for &cap_face in cap_faces {
        if rebuild_one_apex_cap(brep, cap_face) {
            succeeded += 1;
        }
    }
    succeeded
}

/// Rebuild one apex cap.  Returns `true` on full success
/// (plan actionable + execute OK + all lenses collapsed),
/// `false` otherwise.  Snapshot-rollback inside each primitive
/// keeps the BRep clean on partial failure.
#[allow(dead_code)]
pub(crate) fn rebuild_one_apex_cap(brep: &mut BRep, cap_face: FaceHandle) -> bool {
    let trace = std::env::var_os("POPKERN_APEX_REBUILD_TRACE").is_some();
    let plan = plan_apex_region_deletion(brep, cap_face);
    if trace {
        eprintln!(
            "[rebuild] cap={cap_face:?} keep={} delete={} gens={} arcs={} actionable={}",
            plan.keep_vertices.len(),
            plan.delete_vertices.len(),
            plan.delete_cap_cylinder_generators.len(),
            plan.insert_cap_tangent_arcs.len(),
            plan.is_actionable(),
        );
    }
    if !plan.is_actionable() {
        return false;
    }

    // Snapshot before mef so we can fully roll back on failure
    // (the per-step primitives already rollback on validate
    // failure, but a partial collapse — e.g., 2 of 3 lenses
    // collapsed — is still in an intermediate state we don't
    // want to ship by default).
    let snapshot = brep.clone();

    let exec = match execute_deletion_plan(brep, &plan) {
        Ok(r) => r,
        Err(e) => {
            if trace {
                eprintln!("[rebuild] execute_deletion_plan FAILED: {e:?}");
            }
            *brep = snapshot;
            return false;
        }
    };
    if trace {
        eprintln!(
            "[rebuild] execute OK: lens_faces={} inserted_arcs={}",
            exec.lens_faces.len(),
            exec.inserted_arcs.len()
        );
    }

    // Each lens face's outer loop has 3 edges.  Two of them are
    // the obsolete generators; one is the cap-tangent arc (just
    // inserted by mef).  Find the gen pair for each lens and
    // pass to `collapse_apex_lens`.
    for &lens_face in &exec.lens_faces {
        let lens_loop = match brep.faces.get(lens_face) {
            Some(f) => f.outer_loop,
            None => {
                *brep = snapshot;
                return false;
            }
        };
        let start_he = brep.loops[lens_loop].half_edge;
        let mut gen_pair: Vec<popkern_core::topology::HalfEdgeHandle> = Vec::new();
        let mut cursor = start_he;
        for _ in 0..3 {
            let he_edge = brep.half_edges[cursor].edge;
            let he_next = brep.half_edges[cursor].next;
            if !exec.inserted_arcs.contains(&he_edge) {
                gen_pair.push(cursor);
            }
            cursor = he_next;
            if cursor == start_he {
                break;
            }
        }
        if gen_pair.len() != 2 {
            *brep = snapshot;
            return false;
        }
        if let Err(e) = collapse_apex_lens(brep, lens_face, gen_pair[0], gen_pair[1]) {
            if trace {
                eprintln!("[rebuild] collapse_apex_lens FAILED on {lens_face:?}: {e:?}");
            }
            *brep = snapshot;
            return false;
        }
    }

    if trace {
        eprintln!("[rebuild] success");
    }
    true
}
