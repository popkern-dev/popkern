//! Phase B/C data structures and helpers for single-pass fillet topology.

use popkern_core::geometry::NurbsCurve;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::error::BlendFailure;

/// A trim to be applied to a face during Phase B.
pub(super) struct FaceTrim {
    /// Index into the edge_geometries array (determines sort order).
    pub edge_idx: usize,
    /// Which side of the fillet (1 or 2).
    pub side: u8,
    /// Anchor vertices (start, end) created during Phase A.
    pub anchors: (VertexHandle, VertexHandle),
}

/// Result of Phase B face splitting for one trim.
pub(super) struct TrimSplitResult {
    /// The fillet edge index this trim belongs to.
    pub edge_idx: usize,
    /// Which side (1 or 2).
    pub side: u8,
    /// The consumed sub-face (contains the fillet edge).
    pub consumed: FaceHandle,
    /// The trim edge(s) created by split_face.
    pub trim_edges: Vec<EdgeHandle>,
    /// Pcurve on the kept face's surface, computed during Phase B before kef
    /// changes face identity.
    pub pcurve: NurbsCurve,
    /// 3D trim curve between actual anchor vertex positions.
    pub trim_3d: NurbsCurve,
}

/// Find which face in `candidates` has both vertices on its outer loop.
pub(super) fn find_face_containing_vertices(
    brep: &BRep,
    candidates: &[FaceHandle],
    v1: VertexHandle,
    v2: VertexHandle,
) -> Option<FaceHandle> {
    for &fh in candidates {
        let face = match brep.faces.get(fh) {
            Some(f) => f,
            None => continue,
        };
        let start_he = match brep.loops.get(face.outer_loop) {
            Some(l) => l.half_edge,
            None => continue,
        };
        if !brep.half_edges.contains_key(start_he) {
            continue;
        }

        let mut found_v1 = false;
        let mut found_v2 = false;
        let mut cursor = start_he;
        let limit = brep.half_edges.len() + 1;
        for _ in 0..limit {
            let origin = brep.half_edges[cursor].origin;
            if origin == v1 {
                found_v1 = true;
            }
            if origin == v2 {
                found_v2 = true;
            }
            if found_v1 && found_v2 {
                return Some(fh);
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }
    }
    None
}

/// Split a face along a trim curve, handling errors gracefully.
pub(super) fn split_face_with_trim(
    brep: &mut BRep,
    face: FaceHandle,
    pcurve: &NurbsCurve,
    trim_3d: &NurbsCurve,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<popkern_core::topology::SplitFaceResult, BlendFailure> {
    brep.split_face(face, pcurve, trim_3d, tolerance)
        .map_err(|e| BlendFailure::InvalidTopology {
            reason: format!("face split failed: {e}"),
        })
}
