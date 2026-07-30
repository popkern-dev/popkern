//! Trace-only cap-neighborhood dumps for the repeated-apex partition's
//! cap-to-cylinder seam census failure diagnosis.

use super::*;

/// Trace-only ground-truth dump of the cap outer loop: for every half-edge,
/// report the twin's owning face/surface and whether the edge record's
/// canonical half-edge still references a live side of this edge. Diagnoses
/// seams that exist in the loop but are invisible to edge-map iteration.
pub(super) fn trace_cap_loop_walk(brep: &BRep, cap_face: FaceHandle) {
    let mut visited = std::collections::BTreeSet::new();
    let mut frontier = vec![cap_face];
    for ring in 0..2 {
        let mut next_frontier = Vec::new();
        for face in frontier {
            if !visited.insert(face) {
                continue;
            }
            for neighbor in trace_face_loop_walk(brep, face, ring) {
                next_frontier.push(neighbor);
            }
        }
        frontier = next_frontier;
    }
}

/// Dump one face's outer loop (he/edge/twin-face/surface) and return its
/// distinct neighbor faces so the caller can walk outward ring by ring.
fn trace_face_loop_walk(brep: &BRep, face: FaceHandle, ring: usize) -> Vec<FaceHandle> {
    let Some(outer_loop) = brep.faces.get(face).map(|face_data| face_data.outer_loop) else {
        eprintln!("[direct-apex-trace census-walk] ring{ring} face {face:?} is stale");
        return Vec::new();
    };
    let surface = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| format!("{:?}", std::mem::discriminant(&surface.geometry)))
        .unwrap_or_else(|| "NO-SURFACE".to_owned());
    let Ok(half_edges) = collect_loop_half_edges(brep, outer_loop) else {
        eprintln!(
            "[direct-apex-trace census-walk] ring{ring} face {face:?} ({surface}) loop \
             {outer_loop:?} does not close"
        );
        return Vec::new();
    };
    let mut neighbors = Vec::new();
    let mut entries = Vec::new();
    for he in half_edges {
        let Some(side) = brep.half_edges.get(he) else {
            continue;
        };
        let twin_desc = match brep.half_edges.get(side.twin) {
            None => format!("twin {:?} DEAD", side.twin),
            Some(twin_side) => match brep.loops.get(twin_side.loop_) {
                None => format!("twin loop {:?} stale", twin_side.loop_),
                Some(loop_data) => {
                    neighbors.push(loop_data.face);
                    format!("{:?}", loop_data.face)
                }
            },
        };
        entries.push(format!("{he:?}/{:?}->{twin_desc}", side.edge));
    }
    eprintln!(
        "[direct-apex-trace census-walk] ring{ring} face {face:?} ({surface}): {}",
        entries.join(", ")
    );
    neighbors.sort();
    neighbors.dedup();
    neighbors
}

/// Trace-gated face census for partition-stage observability: total
/// faces and surfaceless faces after each late apply pass.
pub(in crate::apex_construction) fn trace_face_census(brep: &BRep, stage: &str) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    let surfaceless = brep
        .faces
        .iter()
        .filter(|(_, face)| face.surface.is_none())
        .count();
    eprintln!(
        "[direct-apex-trace §14.9/face-census] {stage}: faces={} surfaceless={surfaceless}",
        brep.faces.len()
    );
}
