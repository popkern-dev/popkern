//! Face-use traversal for STEP pcurve export.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct EdgeFaceUse {
    pub face: FaceHandle,
    pub half_edge: HalfEdgeHandle,
}

pub(super) fn edge_face_uses(brep: &BRep, edge: EdgeHandle) -> Vec<EdgeFaceUse> {
    let mut uses = Vec::new();
    for (face_handle, face) in &brep.faces {
        collect_loop_edge_face_uses(brep, face_handle, face.outer_loop, edge, &mut uses);
        for &inner_loop in &face.inner_loops {
            collect_loop_edge_face_uses(brep, face_handle, inner_loop, edge, &mut uses);
        }
    }
    uses.sort_by_key(|face_use| format!("{:?}:{:?}", face_use.face, face_use.half_edge));
    uses.dedup();
    uses
}

pub(super) fn unique_edge_half_edges(brep: &BRep, edge: EdgeHandle) -> Vec<HalfEdgeHandle> {
    brep.half_edges
        .iter()
        .filter_map(|(half_edge, data)| (data.edge == edge).then_some(half_edge))
        .collect()
}

fn collect_loop_edge_face_uses(
    brep: &BRep,
    face: FaceHandle,
    loop_handle: LoopHandle,
    edge: EdgeHandle,
    uses: &mut Vec<EdgeFaceUse>,
) {
    let Some(loop_data) = brep.loops.get(loop_handle) else {
        return;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return;
    }
    let mut cursor = start;
    let bound = brep.half_edges.len().max(1);
    for _ in 0..=bound {
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            break;
        };
        if half_edge.edge == edge {
            uses.push(EdgeFaceUse {
                face,
                half_edge: cursor,
            });
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
}
