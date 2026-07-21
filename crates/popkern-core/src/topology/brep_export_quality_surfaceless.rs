//! Surfaceless-face census for export-quality diagnostics (PK-156).

use super::BrepSurfacelessFaceRecord;
use crate::topology::BRep;

/// Census of faces with no surface geometry at the export boundary.
pub(super) fn surfaceless_face_records(brep: &BRep) -> Vec<BrepSurfacelessFaceRecord> {
    let mut records: Vec<BrepSurfacelessFaceRecord> = brep
        .faces
        .iter()
        .filter(|(_, face)| face.surface.is_none())
        .map(|(face, face_data)| {
            let outer_ring_edges = brep
                .loops
                .get(face_data.outer_loop)
                .map(|_| {
                    let mut count = 0usize;
                    if let Some(loop_data) = brep.loops.get(face_data.outer_loop) {
                        let start = loop_data.half_edge;
                        let mut current = start;
                        for _ in 0..brep.half_edges.len() {
                            count += 1;
                            let Some(half_edge) = brep.half_edges.get(current) else {
                                break;
                            };
                            current = half_edge.next;
                            if current == start {
                                break;
                            }
                        }
                    }
                    count
                })
                .unwrap_or(0);
            BrepSurfacelessFaceRecord {
                face: format!("{face:?}"),
                role: format!("{:?}", face_data.role),
                outer_ring_edges,
                inner_loops: face_data.inner_loops.len(),
            }
        })
        .collect();
    records.sort_by(|a, b| a.face.cmp(&b.face));
    records
}
