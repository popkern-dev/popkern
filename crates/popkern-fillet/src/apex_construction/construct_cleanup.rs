//! Cleanup helpers for the legacy A2 mutating builder.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

pub(super) fn dissolve_cap_face_slits(
    brep: &mut BRep,
    cap_face: FaceHandle,
    protected_edges: &[EdgeHandle],
) {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let max_iters = brep.half_edges.len() + 16;
    for _ in 0..max_iters {
        let Some(face_data) = brep.faces.get(cap_face) else {
            return;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            return;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_walk = brep.half_edges.len() + 1;
        let mut slit: Option<HalfEdgeHandle> = None;
        while let Some(he) = brep.half_edges.get(cursor) {
            let twin = he.twin;
            let twin_loop = brep.half_edges.get(twin).map(|t| t.loop_);
            if twin_loop == Some(outer_loop) && !protected_edges.contains(&he.edge) {
                slit = Some(cursor);
                break;
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_walk {
                break;
            }
        }
        let Some(slit_he) = slit else {
            return;
        };
        if trace {
            eprintln!("[a2-trace §28.H/3.8] dissolving cap-face slit {slit_he:?} on {cap_face:?}");
        }
        if brep.dissolve_same_loop_bridge(slit_he).is_err() {
            return;
        }
    }
}

pub(super) fn dissolve_v_bridges_until_valence_le_1(brep: &mut BRep, vertex: VertexHandle) {
    if !brep.vertices.contains_key(vertex) {
        return;
    }
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let max_iters = brep.half_edges.len() + 16;
    for _ in 0..max_iters {
        if !brep.vertices.contains_key(vertex) {
            return;
        }
        let outgoings: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, h)| h.origin == vertex)
            .map(|(h, _)| h)
            .collect();
        if outgoings.len() <= 1 {
            return;
        }
        let mut to_dissolve: Option<HalfEdgeHandle> = None;
        for &he in &outgoings {
            let Some(he_data) = brep.half_edges.get(he) else {
                continue;
            };
            let twin = he_data.twin;
            let Some(twin_data) = brep.half_edges.get(twin) else {
                continue;
            };
            if he_data.loop_ == twin_data.loop_ {
                to_dissolve = Some(he);
                break;
            }
        }
        let Some(he) = to_dissolve else {
            return;
        };
        if trace {
            eprintln!(
                "[a2-trace §28.H/3.5] dissolving same-loop bridge {he:?} at V={vertex:?} (valence={})",
                outgoings.len()
            );
        }
        if brep.dissolve_same_loop_bridge(he).is_err() {
            return;
        }
    }
}
