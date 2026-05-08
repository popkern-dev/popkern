//! Cleanup helpers for the legacy A2 mutating builder.

use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

pub(super) fn dissolve_v_bridges_until_valence_le_1(brep: &mut BRep, vertex: VertexHandle) {
    if !brep.vertices.contains_key(vertex) {
        return;
    }
    let trace = crate::a2_trace::enabled();
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
