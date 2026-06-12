//! Face-use-local trim authority storage.

use crate::geometry::pcurve::Pcurve;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};
use serde::{Deserialize, Serialize};

/// Stable key for trim authority owned by one face's use of a half-edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct FaceUseTrimKey {
    /// Face that owns the trim use.
    pub face: FaceHandle,
    /// Half-edge used by that face loop.
    pub half_edge: HalfEdgeHandle,
}

impl BRep {
    /// Store a face-use-local pcurve. Returns false for stale handles.
    pub fn set_face_use_pcurve(
        &mut self,
        face: FaceHandle,
        half_edge: HalfEdgeHandle,
        pcurve: Pcurve,
    ) -> bool {
        if !self.faces.contains_key(face) || !self.half_edges.contains_key(half_edge) {
            return false;
        }
        self.face_use_pcurves
            .insert(FaceUseTrimKey { face, half_edge }, pcurve);
        true
    }

    /// Return explicit face-use-local pcurve authority for a live face use.
    pub fn face_use_pcurve(&self, face: FaceHandle, half_edge: HalfEdgeHandle) -> Option<&Pcurve> {
        if !self.faces.contains_key(face) || !self.half_edges.contains_key(half_edge) {
            return None;
        }
        self.face_use_pcurves
            .get(&FaceUseTrimKey { face, half_edge })
    }
}
