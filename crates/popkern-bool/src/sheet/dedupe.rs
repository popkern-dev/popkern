//! Face-subset removal with topology pruning.
//!
//! [`drop_faces_and_prune`] removes a specified set of faces from a
//! BRep in place, along with the loops they own, and then prunes any
//! orphaned half-edges, edges, and vertices. Used by the sheet
//! dispatcher to deduplicate coincident faces before concatenating
//! two open shells for `SheetSheetUnion`.

use popkern_core::topology::{BRep, FaceHandle};
use std::collections::HashSet;

/// Remove every face in `faces_to_remove` from `brep`, then prune
/// orphaned half-edges / edges / vertices left behind.
///
/// The face's owned outer and inner loops are removed first, then
/// half-edges whose loop is gone are dropped, then edges whose
/// canonical half-edge is gone, then vertices no longer referenced
/// by any half-edge. The result is a topologically consistent
/// sub-BRep containing only the surviving faces.
///
/// This is structurally identical to the orphan-cleanup tail of
/// [`super::trim::trim_sheet`] — extracted into its own helper so
/// [`super::dispatch::boolean_with_sheets`] can use it to drop
/// coincident duplicate faces during [`crate::sheet::SheetOperationKind::SheetSheetUnion`].
pub(super) fn drop_faces_and_prune(brep: &mut BRep, faces_to_remove: &HashSet<FaceHandle>) {
    // Step 1: remove the targeted faces + their loops, in canonical
    // handle order — HashSet iteration order would permute the face
    // and loop slotmap free lists of a served result BRep, making
    // handle allocation in any later insertion run-varying (PK-300
    // hash-iteration audit).
    let mut face_order: Vec<FaceHandle> = faces_to_remove.iter().copied().collect();
    face_order.sort_unstable();
    for fh in &face_order {
        if let Some(face) = brep.faces.get(*fh) {
            let outer_loop = face.outer_loop;
            let inner_loops = face.inner_loops.clone();
            brep.loops.remove(outer_loop);
            for il in inner_loops {
                brep.loops.remove(il);
            }
            brep.faces.remove(*fh);
        }
    }

    // Step 2: prune orphaned half-edges.
    let orphan_he: Vec<_> = brep
        .half_edges
        .iter()
        .filter_map(|(heh, he)| {
            if brep.loops.get(he.loop_).is_none() {
                Some(heh)
            } else {
                None
            }
        })
        .collect();
    for heh in orphan_he {
        brep.remove_half_edge_reanchored(heh);
    }

    // Step 3: prune orphaned edges.
    let orphan_edges: Vec<_> = brep
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            if brep.half_edges.get(edge.half_edge).is_none() {
                Some(eh)
            } else {
                None
            }
        })
        .collect();
    for eh in orphan_edges {
        brep.edges.remove(eh);
    }

    // Step 4: prune orphaned vertices.
    let mut vertex_in_use = HashSet::new();
    for (_, he) in &brep.half_edges {
        vertex_in_use.insert(he.origin);
    }
    let orphan_vertices: Vec<_> = brep
        .vertices
        .keys()
        .filter(|vh| !vertex_in_use.contains(vh))
        .collect();
    for vh in orphan_vertices {
        brep.vertices.remove(vh);
    }
}
