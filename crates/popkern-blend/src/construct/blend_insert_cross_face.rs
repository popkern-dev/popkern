//! Merging a corner's cross-face V-edges into one blend face.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::blend_topo::{collect_edges_through_vertex, repair_vertex_outgoing};

/// Kef every V-edge that separates two distinct faces.
///
/// Iterates until no cross-face V-edges remain. For trihedral
/// corners this merges (N-1) face-pairs out of N; the last V-edge
/// becomes a same-face bridge which `kef` cannot remove (its
/// precondition requires distinct loops). The bridge is handled
/// by `kev_v_bridge` in the kev pass, after the surface has been
/// assigned to the merged face.
///
/// Returns the surviving merged-face handle.
///
/// `pub(crate)` because the §6c-DirectApex.2.e apex builder
/// (`crate::apex_construction::construct`) calls this after its
/// own mef sequence to merge the corner sub-faces.
pub(crate) fn kef_cross_face_v_edges(
    brep: &mut BRep,
    vertex: VertexHandle,
) -> Result<FaceHandle, String> {
    let mut blend_face: Option<FaceHandle> = None;
    let mut safety = brep.edges.len() + 16;

    loop {
        if safety == 0 {
            return Err("kef_cross_face_v_edges: iteration cap hit".to_string());
        }
        safety -= 1;

        if !brep.vertices.contains_key(vertex) {
            break;
        }

        // `kef` invalidates V's `outgoing` field when the removed
        // half-edge *was* V's outgoing. Re-pick a valid outgoing
        // before every incidence walk.
        repair_vertex_outgoing(brep, vertex)?;

        let edges_through_v = collect_edges_through_vertex(brep, vertex)?;
        if edges_through_v.is_empty() {
            break;
        }

        let mut kef_target: Option<HalfEdgeHandle> = None;
        for edge in &edges_through_v {
            if !brep.edges.contains_key(*edge) {
                continue;
            }
            let he = brep.edges[*edge].half_edge;
            let he_loop = brep.half_edges[he].loop_;
            let twin = brep.require_mate(he).map_err(|e| e.to_string())?;
            let twin_loop = brep.half_edges[twin].loop_;
            if he_loop != twin_loop {
                kef_target = Some(he);
                break;
            }
        }

        let Some(he) = kef_target else {
            // Only same-face bridges remain — the cross-face kef pass is done.
            // Record the pinched face as the surviving blend face
            // (it's the loop that currently holds V).
            if blend_face.is_none() {
                let outgoing = brep.vertices[vertex].outgoing;
                if brep.half_edges.contains_key(outgoing) {
                    let loop_h = brep.half_edges[outgoing].loop_;
                    blend_face = Some(brep.loops[loop_h].face);
                }
            }
            break;
        };
        let he_face = brep.loops[brep.half_edges[he].loop_].face;
        if popkern_core::trace::trace_enabled("corner-carve") {
            let twin = brep.require_mate(he).map_err(|e| e.to_string())?;
            let twin_face = brep.loops[brep.half_edges[twin].loop_].face;
            let count = |fh: FaceHandle| {
                let start = brep.loops[brep.faces[fh].outer_loop].half_edge;
                let mut he2 = start;
                let mut n = 0;
                for _ in 0..256 {
                    n += 1;
                    he2 = brep.half_edges[he2].next;
                    if he2 == start {
                        break;
                    }
                }
                n
            };
            eprintln!(
                "  [vblend] kef at {vertex:?}: merging {he_face:?} (loop {}) with {twin_face:?} (loop {})",
                count(he_face),
                count(twin_face)
            );
        }
        brep.kef(he)
            .map_err(|e| format!("kef failed through vertex {vertex:?}: {e}"))?;
        blend_face = Some(he_face);
    }

    blend_face.ok_or_else(|| {
        "kef_cross_face_v_edges: no merge performed — blend face not created".to_string()
    })
}
