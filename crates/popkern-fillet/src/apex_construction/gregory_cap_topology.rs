//! Topology helpers for A2 Gregory cap assembly.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

pub(super) fn merge_cap_side_accumulators(
    brep: &mut BRep,
    apex_edge_handles: &[EdgeHandle],
) -> Result<(), String> {
    use std::collections::HashSet;
    if apex_edge_handles.is_empty() {
        return Ok(());
    }
    let perimeter: HashSet<EdgeHandle> = apex_edge_handles.iter().copied().collect();
    let cap_iters = apex_edge_handles.len() * 4 + 32;
    for _ in 0..cap_iters {
        let face_apex_count = count_apex_edges_per_face(brep, &perimeter);
        let mut victim: Option<HalfEdgeHandle> = None;
        for (eh, edge) in brep.edges.iter() {
            if perimeter.contains(&eh) {
                continue;
            }
            let primary = edge.half_edge;
            let Some(he) = brep.half_edges.get(primary) else {
                continue;
            };
            let Some(twin_he) = brep.half_edges.get(he.twin) else {
                continue;
            };
            if he.loop_ == twin_he.loop_ {
                continue;
            }
            let face_a = match brep.loops.get(he.loop_) {
                Some(l) => l.face,
                None => continue,
            };
            let face_b = match brep.loops.get(twin_he.loop_) {
                Some(l) => l.face,
                None => continue,
            };
            let count_a = face_apex_count.get(&face_a).copied().unwrap_or(0);
            let count_b = face_apex_count.get(&face_b).copied().unwrap_or(0);
            if count_a == 0 || count_b == 0 {
                continue;
            }
            if count_a < 2 && count_b < 2 {
                continue;
            }
            victim = Some(primary);
            break;
        }
        let Some(he) = victim else {
            return Ok(());
        };
        brep.kef(he)
            .map_err(|e| format!("kef on cap-side accumulator merge: {e:?}"))?;
    }
    Err("merge_cap_side_accumulators: iteration cap hit".to_string())
}

fn count_apex_edges_per_face(
    brep: &BRep,
    perimeter: &std::collections::HashSet<EdgeHandle>,
) -> std::collections::HashMap<FaceHandle, usize> {
    let mut out: std::collections::HashMap<FaceHandle, usize> = std::collections::HashMap::new();
    for (fh, face) in &brep.faces {
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut count = 0usize;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if perimeter.contains(&he.edge) {
                count += 1;
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        if count > 0 {
            out.insert(fh, count);
        }
    }
    out
}

#[allow(dead_code)]
pub(super) fn merge_unmerged_apex_subfaces(
    brep: &mut BRep,
    apex_edge_handles: &[EdgeHandle],
) -> Result<(), String> {
    let cap_iters = apex_edge_handles.len() * 2 + 8;
    let apex_set: std::collections::HashSet<EdgeHandle> =
        apex_edge_handles.iter().copied().collect();
    let apex_vertex_set = apex_vertex_set(brep, &apex_set);
    for _ in 0..cap_iters {
        let cap_side_faces = faces_with_multi_apex_edges(brep, &apex_set);
        let mut victim: Option<HalfEdgeHandle> = None;
        if cap_side_faces.len() >= 2 {
            victim = find_kef_candidate_between(brep, &cap_side_faces, &apex_set);
        }
        if victim.is_none() {
            victim = find_bisected_spoke_kef_candidate(brep, &apex_set, &apex_vertex_set);
        }
        let Some(he) = victim else {
            return Ok(());
        };
        brep.kef(he)
            .map_err(|e| format!("kef on cap-merging edge: {e:?}"))?;
    }
    Err("merge_unmerged_apex_subfaces: iteration cap hit".to_string())
}

fn apex_vertex_set(
    brep: &BRep,
    apex_set: &std::collections::HashSet<EdgeHandle>,
) -> std::collections::HashSet<VertexHandle> {
    let mut out = std::collections::HashSet::new();
    for &eh in apex_set {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let primary = edge.half_edge;
        if let Some(he) = brep.half_edges.get(primary) {
            out.insert(he.origin);
            if let Some(twin) = brep.half_edges.get(he.twin) {
                out.insert(twin.origin);
            }
        }
    }
    out
}

fn find_kef_candidate_between(
    brep: &BRep,
    cap_side_faces: &[FaceHandle],
    apex_set: &std::collections::HashSet<EdgeHandle>,
) -> Option<HalfEdgeHandle> {
    for fh in cap_side_faces {
        let face = brep.faces.get(*fh)?;
        let loop_data = brep.loops.get(face.outer_loop)?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if !apex_set.contains(&he.edge) {
                let twin_loop = brep.half_edges.get(he.twin).map(|t| t.loop_);
                let twin_face = twin_loop.and_then(|tl| brep.loops.get(tl).map(|l| l.face));
                if let Some(tf) = twin_face
                    && tf != *fh
                    && cap_side_faces.contains(&tf)
                {
                    return Some(cursor);
                }
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    None
}

fn find_bisected_spoke_kef_candidate(
    brep: &BRep,
    apex_set: &std::collections::HashSet<EdgeHandle>,
    apex_vertex_set: &std::collections::HashSet<VertexHandle>,
) -> Option<HalfEdgeHandle> {
    let containing_faces: Vec<FaceHandle> = brep
        .faces
        .iter()
        .filter_map(|(fh, face)| {
            let loop_data = brep.loops.get(face.outer_loop)?;
            let start = loop_data.half_edge;
            let mut cursor = start;
            let mut iters = 0;
            let max_iters = brep.half_edges.len() + 1;
            while let Some(he) = brep.half_edges.get(cursor) {
                if apex_set.contains(&he.edge) {
                    return Some(fh);
                }
                cursor = he.next;
                iters += 1;
                if cursor == start || iters >= max_iters {
                    return None;
                }
            }
            None
        })
        .collect();
    for fh in &containing_faces {
        let face = brep.faces.get(*fh)?;
        let loop_data = brep.loops.get(face.outer_loop)?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if !apex_set.contains(&he.edge) && apex_vertex_set.contains(&he.origin) {
                let twin = he.twin;
                let twin_origin = brep.half_edges.get(twin).map(|t| t.origin);
                let twin_loop = brep.half_edges.get(twin).map(|t| t.loop_);
                let twin_face = twin_loop.and_then(|tl| brep.loops.get(tl).map(|l| l.face));
                if let (Some(to), Some(tf)) = (twin_origin, twin_face)
                    && apex_vertex_set.contains(&to)
                    && tf != *fh
                    && containing_faces.contains(&tf)
                {
                    return Some(cursor);
                }
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    None
}

fn faces_with_multi_apex_edges(
    brep: &BRep,
    apex_set: &std::collections::HashSet<EdgeHandle>,
) -> Vec<FaceHandle> {
    let mut out = Vec::new();
    for (fh, face) in &brep.faces {
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        let mut count = 0usize;
        while let Some(he) = brep.half_edges.get(cursor) {
            if apex_set.contains(&he.edge) {
                count += 1;
                if count >= 2 {
                    out.push(fh);
                    break;
                }
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    out
}

pub(super) fn find_face_bounded_by_all_apex_edges(
    brep: &BRep,
    apex_edge_handles: &[EdgeHandle],
) -> Option<FaceHandle> {
    if apex_edge_handles.is_empty() {
        return None;
    }
    for (fh, face) in &brep.faces {
        let outer_loop = face.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut loop_edges: Vec<EdgeHandle> = Vec::new();
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            loop_edges.push(he.edge);
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        if apex_edge_handles.iter().all(|eh| loop_edges.contains(eh)) {
            return Some(fh);
        }
    }
    None
}
