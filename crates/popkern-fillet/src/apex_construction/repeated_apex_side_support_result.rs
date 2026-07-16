//! Result payloads from executable side-support continuation.

use popkern_core::topology::{FaceHandle, LoopHandle};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct SideSupportRegionContinuationResult {
    primitive: SideSupportContinuationPrimitive,
    bridge_faces: Vec<FaceHandle>,
    bridge_loops: Vec<LoopHandle>,
    /// For the cross-owner paired primitive: the fillet index each
    /// bridge face serves, aligned with `bridge_faces`. The bridge IS
    /// that fillet's new strip and must receive its cylinder surface
    /// (a naked bridge face makes the partition's cap-to-cylinder
    /// seam census read zero — measured on the strict V-free gate).
    bridge_fillet_indices: Vec<usize>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) enum SideSupportContinuationPrimitive {
    #[default]
    None,
    FaceFamilyBridge,
    OutsideOutsideCylinderMediated,
}

impl SideSupportRegionContinuationResult {
    pub(super) fn not_applied() -> Self {
        Self::default()
    }

    pub(super) fn from_face_family_bridge_payload(
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self::from_bridge_payload(
            SideSupportContinuationPrimitive::FaceFamilyBridge,
            bridge_faces,
            bridge_loops,
        )
    }

    pub(super) fn from_outside_outside_cylinder_payload(
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self::from_bridge_payload(
            SideSupportContinuationPrimitive::OutsideOutsideCylinderMediated,
            bridge_faces,
            bridge_loops,
        )
    }

    pub(super) fn from_cross_owner_paired_payload(
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
        bridge_fillet_indices: Vec<usize>,
    ) -> Self {
        let mut result = Self::from_bridge_payload(
            SideSupportContinuationPrimitive::OutsideOutsideCylinderMediated,
            bridge_faces,
            bridge_loops,
        );
        result.bridge_fillet_indices = bridge_fillet_indices;
        result
    }

    /// `(bridge_face, fillet_index)` pairs for the cross-owner paired
    /// primitive; empty for every other outcome.
    pub(super) fn bridge_fillet_assignments(&self) -> Vec<(FaceHandle, usize)> {
        self.bridge_faces
            .iter()
            .copied()
            .zip(self.bridge_fillet_indices.iter().copied())
            .collect()
    }

    fn from_bridge_payload(
        primitive: SideSupportContinuationPrimitive,
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self {
            primitive,
            bridge_faces,
            bridge_loops,
            bridge_fillet_indices: Vec::new(),
        }
    }

    pub(super) fn was_applied(&self) -> bool {
        !self.bridge_faces.is_empty() || !self.bridge_loops.is_empty()
    }

    pub(super) fn bridge_loops(&self) -> &[LoopHandle] {
        &self.bridge_loops
    }

    pub(super) fn bridge_face_count(&self) -> usize {
        self.bridge_faces.len()
    }

    pub(super) fn bridge_loop_count(&self) -> usize {
        self.bridge_loops.len()
    }

    pub(super) fn primitive_code(&self) -> &'static str {
        self.primitive.code()
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "side_support_result=primitive={} applied={} bridge_faces={:?} bridge_loops={:?}",
            self.primitive_code(),
            self.was_applied(),
            self.bridge_faces,
            self.bridge_loops
        )
    }
}

impl SideSupportContinuationPrimitive {
    fn code(&self) -> &'static str {
        match self {
            SideSupportContinuationPrimitive::None => "none",
            SideSupportContinuationPrimitive::FaceFamilyBridge => "face_family_bridge",
            SideSupportContinuationPrimitive::OutsideOutsideCylinderMediated => {
                "outside_outside_cylinder_mediated"
            }
        }
    }
}
